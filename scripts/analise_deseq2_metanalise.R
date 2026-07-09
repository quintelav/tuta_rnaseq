#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline de Metanálise Transcriptómica (Nativo em R - Padrão-Ouro)
# ------------------------------------------------------------------------------
# 1. Importação de dados do Salmon via tximport.
# 2. Pré-filtragem de genes com baixa contagem para evitar inflação de zeros.
# 3. Modelagem Binomial Negativa (DESeq2) para cada tratamento independentemente.
# 4. Extração de p-valores ajustados (FDR) e Log2FoldChange.
# 5. Geração de matriz consolidada de Genes Diferencialmente Expressos (DEGs).
# ==============================================================================

cat("==================================================\n")
cat("   METANÁLISE DE RNA-SEQ VIA DESEQ2 (R)           \n")
cat("==================================================\n")

# Carregar bibliotecas silenciosamente
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

# Configuração de caminhos
dir_quant <- "quant"
arquivo_meta <- "metadata_geral.csv"
# Novo nome do arquivo de saída (livre da nomenclatura antiga)
arquivo_saida <- "quant/deseq2_resultados_metanalise.csv"

# 1. Ler Metadados
if (!file.exists(arquivo_meta)) {
  stop("❌ [ERRO] Arquivo 'metadata_geral.csv' não encontrado!")
}

meta <- read.csv(arquivo_meta, stringsAsFactors = FALSE)
meta$condition <- tolower(meta$condition)
meta$study_id <- tolower(meta$study_id)

# Mapear ficheiros do Salmon
arquivos_salmon <- file.path(dir_quant, meta$sample_id, "quant.sf")
names(arquivos_salmon) <- meta$sample_id

# Filtrar apenas as amostras que realmente existem no disco
existe <- file.exists(arquivos_salmon)
if (sum(existe) == 0) {
  stop("❌ [ERRO] Nenhuma pasta do Salmon encontrada!")
}
meta <- meta[existe, ]
arquivos_salmon <- arquivos_salmon[existe]

cat(sprintf("📂 Foram encontradas %d amostras processadas pelo Salmon.\n", length(arquivos_salmon)))
cat("⏳ Importando contagens com tximport...\n")

# txOut=TRUE indica que queremos análise ao nível de transcrito (não gene)
txi <- tximport(arquivos_salmon, type = "salmon", txOut = TRUE)

estudos <- unique(meta$study_id)
genes_totais <- rownames(txi$counts)

# Dataframe matriz para guardar os resultados puros
df_resultados <- data.frame(row.names = genes_totais)
matriz_up <- numeric(length(genes_totais))
matriz_down <- numeric(length(genes_totais))

colunas_log2fc <- c()
colunas_padj <- c()

# 2. Executar DESeq2 iterativamente por Inseticida
for (estudo in estudos) {
  cat("\n==================================================\n")
  cat(sprintf("📊 Ajustando Modelo DESeq2: '%s'\n", toupper(estudo)))
  cat("==================================================\n")
  
  meta_sub <- meta[meta$study_id == estudo, ]
  
  if (nrow(meta_sub) < 4) {
    cat("⚠️ Amostras insuficientes para estatística. Pulando...\n")
    next
  }
  
  # Subconjunto dos dados do tximport
  txi_sub <- list(
    abundance = txi$abundance[, meta_sub$sample_id],
    counts = txi$counts[, meta_sub$sample_id],
    length = txi$length[, meta_sub$sample_id],
    countsFromAbundance = txi$countsFromAbundance
  )
  
  # Preparar colData (As linhas devem bater com as colunas do txi)
  rownames(meta_sub) <- meta_sub$sample_id
  meta_sub$condition <- factor(meta_sub$condition, levels = c("control", "treatment"))
  
  # Criar objeto DESeq
  dds <- DESeqDataSetFromTximport(txi_sub, colData = meta_sub, design = ~ condition)
  
  # ---------------------------------------------------------------------------
  # PRÉ-FILTRAGEM RIGOROSA
  # Manter transcritos com pelo menos 10 reads em N réplicas do menor grupo
  # ---------------------------------------------------------------------------
  min_rep <- floor(nrow(meta_sub) / 2)
  keep <- rowSums(counts(dds) >= 10) >= min_rep
  dds <- dds[keep, ]
  cat(sprintf("   🧹 Pré-filtragem: Restaram %d transcritos expressos.\n", sum(keep)))
  
  # Executar o pipeline core do DESeq2
  cat("   ⏳ Estimando tamanho de bibliotecas e dispersão (Binomial Negativa)...\n")
  dds <- DESeq(dds, quiet = TRUE)
  
  cat("   ⏳ Extraindo Contrastes e P-valores ajustados (FDR < 0.05)...\n")
  res <- results(dds, contrast = c("condition", "treatment", "control"))
  
  # Preencher vetores de resultados com zeros/uns para genes filtrados
  lfc_vec <- rep(0, length(genes_totais))
  padj_vec <- rep(1.0, length(genes_totais))
  names(lfc_vec) <- genes_totais
  names(padj_vec) <- genes_totais
  
  # Mapear resultados estatísticos reais
  valid_genes <- rownames(res)
  lfc_vec[valid_genes] <- ifelse(is.na(res$log2FoldChange), 0, res$log2FoldChange)
  padj_vec[valid_genes] <- ifelse(is.na(res$padj), 1.0, res$padj)
  
  # Gravar no dataframe mestre com nomenclatura padrão
  col_lfc <- paste0("Log2FC_", estudo)
  col_p <- paste0("padj_", estudo)
  
  df_resultados[[col_lfc]] <- lfc_vec
  df_resultados[[col_p]] <- padj_vec
  
  colunas_log2fc <- c(colunas_log2fc, col_lfc)
  colunas_padj <- c(colunas_padj, col_p)
  
  # Identificação de DEGs estritos (FDR < 0.05 e |Log2FC| >= 1.0)
  is_up <- (lfc_vec >= 1.0) & (padj_vec < 0.05)
  is_down <- (lfc_vec <= -1.0) & (padj_vec < 0.05)
  
  cat(sprintf("   ✅ DEGs Identificados: %d Upregulated | %d Downregulated\n", 
              sum(is_up), sum(is_down)))
              
  matriz_up <- matriz_up + as.numeric(is_up)
  matriz_down <- matriz_down + as.numeric(is_down)
}

# 3. CONSOLIDAÇÃO DA METANÁLISE (Sem dependência de metodologias obscuras)
cat("\n📈 A consolidar a intersecção de DEGs entre os inseticidas...\n")

df_resultados$Total_Studies_Upregulated <- matriz_up
df_resultados$Total_Studies_Downregulated <- matriz_down

if (length(colunas_log2fc) > 0) {
  df_resultados$Mean_Log2FC <- rowMeans(df_resultados[, colunas_log2fc, drop = FALSE])
}

if (length(colunas_padj) > 0) {
  # Extrai o p-valor mais significativo entre os estudos como indicador global de confiança
  df_resultados$Min_FDR <- apply(df_resultados[, colunas_padj, drop = FALSE], 1, min)
}

# 4. SALVAR RESULTADOS
if (!dir.exists(dir_quant)) {
  dir.create(dir_quant)
}

write.csv(df_resultados, file = arquivo_saida, quote = FALSE, row.names = TRUE)

cat(sprintf("\n💾 Tabela Master de Expressão Diferencial salva em '%s'!\n", arquivo_saida))
cat("\n==================================================\n")
cat("📊 RESUMO DE DEGs CONSERVADOS (|Log2FC| >= 1.0 & FDR < 0.05):\n")
cat("==================================================\n")
cat(sprintf("   🌟 Upregulated em todos os tratamentos (Core): %d transcritos\n", sum(matriz_up == 3)))
cat(sprintf("   ❄️ Downregulated em todos os tratamentos (Core): %d transcritos\n", sum(matriz_down == 3)))
cat("==================================================\n")
cat("Pipeline Estatístico Concluído. Prossiga para a anotação funcional!\n")
