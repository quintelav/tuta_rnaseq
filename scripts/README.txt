Projeto RNA-seq Tuta absoluta

Status atual:
- WSL instalado
- Conda funcionando
- Ambiente rnaseq criado
- SRA Toolkit instalado
- Salmon instalado
- FastQC instalado
- MultiQC instalado
- Estrutura do projeto criada
- Metadados obtidos do BioProject PRJNA1359366
- Testes iniciais de download realizados

Próximo passo:
- recriar scripts
- baixar poucas amostras
- rodar FastQC

📑 Protocolo de Bioinformática: Meta-análise de RNA-Seq (Tuta absoluta)
Este documento registra todas as etapas executadas no terminal para a preparação do ambiente, reconstrução do transcriptoma de referência unificado e execução do pipeline de quantificação por blocos.

👥 Dados do Projeto
Organismo-alvo: Tuta absoluta (Traça-do-tomateiro)

Objetivo: Meta-análise de expressão gênica diferencial (Abamectina vs. Spinosad vs. Emamectina)

Ambiente Virtual: Conda (rnaseq)

🚀 Passo a Passo dos Comandos
Fase 1: Ativação do Ambiente e Instalação de Dependências
Bash
# 1. Ativar o ambiente com as ferramentas de RNA-Seq
conda activate rnaseq

# 2. Instalar as ferramentas oficiais do NCBI e o gerador de transcritos
conda install -c conda-forge ncbi-datasets-cli unzip -y
conda install -c bioconda gffread -y
Fase 2: Obtenção do Genoma Estrutural (GenBank GCA)
Como a montagem oficial do artigo (GCA_027580185.1) não possui um arquivo de transcritos prontos (rna.fna), baixamos o DNA genômico bruto e as coordenadas de anotação para reconstruí-lo localmente.

Bash
# 1. Garantir que a pasta de referência existe
mkdir -p reference

# 2. Baixar o pacote contendo o Genoma (DNA) e a Anotação estrutural (GFF3)
datasets download genome accession GCA_027580185.1 --include genome,gff3 --filename reference/tuta_genomic.zip

# 3. Extrair o pacote compactado do NCBI
unzip -o reference/tuta_genomic.zip -d reference/tuta_genomic_extracted
Fase 3: Extração Cirúrgica do Transcriptoma (gffread)
Para evitar erros de formatação com linhas não-codificantes, entramos na pasta profunda e usamos filtros estruturais para costurar apenas as regiões de exons que codificam mRNAs válidos.

Bash
# 1. Entrar na pasta onde os arquivos de 615MB e 74MB foram extraídos
cd reference/tuta_genomic_extracted/ncbi_dataset/data/GCA_027580185.1/

# 2. Executar o gffread para gerar o arquivo de transcritos com o ID correto (Tabs)
gffread genomic.gff \
        -g GCA_027580185.1_ASM2758018v1_genomic.fna \
        -w tuta_transcripts.fna \
        -y tuta_proteins.faa \
        -x tuta_cds.fna

# 3. Mover o transcriptoma final para a pasta de referência principal do projeto
mv tuta_transcripts.fna ~/projects/tuta_rnaseq/reference/

# 4. Voltar para a pasta raiz do projeto
cd ~/projects/tuta_rnaseq

# 5. Faxina de segurança: apagar os gigabytes de DNA bruto e zips temporários
rm -rf reference/tuta_genomic.zip reference/tuta_genomic_extracted
Fase 4: Validação e Indexação no Salmon
Bash
# 1. Testar o cabeçalho para garantir que o intruso do bicudo foi removido e a Tuta (Tabs) assumiu
head -n 5 reference/tuta_transcripts.fna

# 2. Deletar qualquer índice misturado antigo
rm -rf reference/salmon_index

# 3. Criar o novo índice molecular definitivo para o alinhamento seletivo
salmon index -t reference/tuta_transcripts.fna -i reference/salmon_index
🛠️ Fase 5: Estrutura do Script de Execução por Bloco (nano)
Para editar ou visualizar o script que automatiza o processamento das amostras sem estourar o armazenamento do seu computador, usamos o editor nano:

Bash
nano rodar_saqe_por_bloco.sh
📝 Conteúdo Lógico do Script Interno (rodar_saqe_por_bloco.sh)
O seu script executa uma lógica estruturada por blocos (ex: spinosad), contendo a seguinte automação protetora de hardware:

Bash
#!/bin/bash
# Script de automação do pipeline SAQE por Bloco de Inseticida

BLOCO=$1
echo "Iniciando o processamento do bloco: $BLOCO"

# Exemplo de loop interno processando os IDs do SRA do bloco selecionado
for ACESSATION in SRR33779431 SRR33779430 SRR33779429 # ... suas outras amostras
do
    echo "--------------------------------------------------"
    echo "Processando amostra: $ACESSATION"
    echo "--------------------------------------------------"
    
    # 1. Download do arquivo binário SRA
    prefetch $ACESSATION
    
    # 2. Extração paralela para FastQ pareado (R1 e R2)
    fasterq-dump --split-files --threads 6 $ACESSATION
    
    # 3. Compactação imediata para economizar espaço em disco
    pigz ${ACESSATION}_1.fastq ${ACESSATION}_2.fastq
    
    # 4. Controle de qualidade e remoção de adaptadores/Poly-G
    trim_galore --paired --quality 20 --fastqc ${ACESSATION}_1.fastq.gz ${ACESSATION}_2.fastq.gz
    
    # 5. Pseudo-alinhamento e Quantificação de Abundância via Salmon
    salmon quant -i reference/salmon_index -l A \
                 -1 ${ACESSATION}_1_val_1.fq.gz \
                 -2 ${ACESSATION}_2_val_2.fq.gz \
                 -p 6 --validateMappings \
                 -o quant_results/${ACESSATION}_quant
                 
    # 6. FAXINA AUTOMÁTICA: Remove os FastQs pesados deixando apenas as matrizes de contagem
    echo "Limpando arquivos intermediários para proteger os 16GB de RAM e o HD..."
    rm -rf $ACESSATION ${ACESSATION}_1.fastq.gz ${ACESSATION}_2.fastq.gz ${ACESSATION}_1_val_1.fq.gz ${ACESSATION}_2_val_2.fq.gz
done

echo "Bloco $BLOCO finalizado com sucesso!"
Para sair do nano salvando as alterações, lembre-se dos comandos: Ctrl + O (para gravar), Enter (para confirmar o nome) e Ctrl + X (para sair).

🚀 Execução do Bloco
Bash
# Dar permissão de execução ao script (só precisa fazer uma vez)
chmod +x rodar_saqe_por_bloco.sh

# Disparar o pipeline para o grupo do Spinosad
./rodar_saqe_por_bloco.sh spinosad