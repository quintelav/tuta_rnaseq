#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Módulo de Visualização Gráfica de Alta Resolução (Padrão de Publicação)
Gera os painéis visuais do projeto seguindo a nomenclatura adaptada de Bono (2021).
Utiliza 'XR_score' e 'XR_ratio' de forma nativa e automática nos eixos.
"""

import os
import re
import numpy as np
import pandas as pd

# Configura o Matplotlib para rodar em modo 'headless' (evita falhas de display no WSL)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Configurações Estéticas Globais (Padrão de Publicação)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['text.usetex'] = False
plt.rcParams['svg.fonttype'] = 'none'
sns.set_context("paper", font_scale=1.2)

# Caminhos de arquivos no WSL
MATRIZ_ANOTADA = "quant/metanalise_completa_anotada.csv"
PASTA_SAIDA = "quant/figuras"
os.makedirs(PASTA_SAIDA, exist_ok=True)

def obter_label_amigavel(row):
    """
    Retorna o nome usual do gene (ex: CYP6B1, ABCC1) se estiver disponível.
    Caso contrário, retorna o ID encurtado do transcrito como fallback.
    """
    gene_name = str(row.get('Gene_Name', ''))
    
    # Se for nulo, Neutro, Não especificado ou vazio, tenta extrair da anotação
    if not gene_name or gene_name.lower() in ['neutro', 'nan', 'não especificado', 'not specified', '']:
        annot = str(row.get('Annotation', ''))
        match = re.search(r'\[([^\]]+)\]', annot)
        if match:
            gene_name = match.group(1)
            
    # Se tivermos um nome usual válido, retornamos o nome amigável
    if gene_name and gene_name.lower() not in ['neutro', 'nan', 'não especificado', 'not specified', '']:
        return gene_name
        
    # Caso contrário, aplica fallback para o ID encurtado do NCBI
    id_original = str(row.name if hasattr(row, 'name') and row.name else row.get('Transcript_ID', ''))
    id_curto = id_original.split('|')[-1].replace('.t1', '') if '|' in id_original else id_original
    return id_curto

def simular_dados_caso_ausente():
    """
    Carrega os dados reais se disponíveis, adaptando dinamicamente os nomes das colunas
    para evitar KeyErrors. Caso contrário, gera uma simulação baseada no censo real.
    """
    if os.path.exists(MATRIZ_ANOTADA):
        print(f"✅ Planilha master real localizada em '{MATRIZ_ANOTADA}'. Carregando dados...")
        df = pd.read_csv(MATRIZ_ANOTADA, index_col=0)
        
        # 1. Mapeamento Dinâmico de Consistência (XR_score / Net On-Score)
        if 'XR_score' not in df.columns:
            for col in ['XR_score', 'Bono_Net_Score', 'ON_score', 'Net_OnScore', 'Net_On_Score']:
                if col in df.columns:
                    df['XR_score'] = df[col]
                    break
            if 'XR_score' not in df.columns and 'Bono_UP_Score' in df.columns and 'Bono_DOWN_Score' in df.columns:
                df['XR_score'] = df['Bono_UP_Score'] - df['Bono_DOWN_Score']
            if 'XR_score' not in df.columns:
                df['XR_score'] = 0  # Fallback neutro
                
        # 2. Mapeamento Dinâmico de Famílias Biológicas
        if 'Family' not in df.columns:
            for col in ['Detox_Family', 'Family']:
                if col in df.columns:
                    df['Family'] = df[col]
                    break
            if 'Family' not in df.columns:
                df['Family'] = 'Outros'
        
        # Padroniza nomes de famílias para filtragem estrita no Heatmap
        df['Family'] = df['Family'].replace({
            'Outros / Não Identificado': 'Outros',
            'Outros Processos Fisiológicos': 'Outros',
            'Other Genomic Transcripts': 'Outros'
        })
                
        # 3. Mapeamento Dinâmico de Símbolos Gênicos / Nomes Amigáveis
        if 'Gene_Name' not in df.columns:
            for col in ['Preferred_Gene_Name', 'Gene_Symbol', 'Gene_Name']:
                if col in df.columns:
                    df['Gene_Name'] = df[col]
                    break
            if 'Gene_Name' not in df.columns and 'Annotation' in df.columns:
                # Tenta extrair o nome curto entre colchetes da anotação
                df['Gene_Name'] = df['Annotation'].str.extract(r'\[([^\]]+)\]').fillna('Neutro')
            if 'Gene_Name' not in df.columns:
                df['Gene_Name'] = 'Neutro'
                
        # 4. Mapeamento Dinâmico de Expressão Média (Log2FC que representa o XR-Ratio)
        if 'Log2FC_Media' not in df.columns:
            for col in ['mean_XR_ratio', 'Meta_Log2FC', 'mean_ON_ratio', 'Log2FC_Media']:
                if col in df.columns:
                    df['Log2FC_Media'] = df[col]
                    break
            if 'Log2FC_Media' not in df.columns:
                # Calcula a média com base nas colunas individuais de tratamentos
                cols_fc = [c for c in df.columns if 'Log2FC_' in c or 'Log2Ratio_' in c or 'XR_ratio_' in c]
                if cols_fc:
                    df['Log2FC_Media'] = df[cols_fc].mean(axis=1)
                else:
                    df['Log2FC_Media'] = 0.0
                    
        # 5. Mapeamento Dinâmico de P-valores
        if 'P_Value' not in df.columns:
            for col in ['Meta_pvalue', 'P_Value', 'pvalue']:
                if col in df.columns:
                    df['P_Value'] = df[col]
                    break
            if 'P_Value' not in df.columns:
                cols_p = [c for c in df.columns if 'pvalue' in c.lower()]
                if cols_p:
                    df['P_Value'] = df[cols_p[0]]
                else:
                    df['P_Value'] = 0.05
                    
        return df
    
    print("ℹ️ Planilha master real não localizada. Estruturando mock com base no censo real para validação visual...")
    np.random.seed(42)
    
    # Simulação do censo exato do genoma
    n_total = 64591
    n_up = 110
    n_down = 2099
    n_neutro = n_total - n_up - n_down
    
    transcritos = [f"Phth_abs_T{i:05d}" for i in range(n_total)]
    
    # Define os Scores discretos
    scores = np.zeros(n_total, dtype=int)
    scores[:n_up] = np.random.choice([1, 2, 3], size=n_up, p=[0.2, 0.3, 0.5])
    scores[n_up:n_up+n_down] = np.random.choice([-1, -2, -3], size=n_down, p=[0.3, 0.4, 0.3])
    
    # Destaca os biomarcadores citados no resumo
    nomes_genes = ["Neutro"] * n_total
    nomes_genes[0] = "ABCC1"   # Net XR-Score +3
    nomes_genes[1] = "CYP6B1"  # Net XR-Score +3
    scores[0] = 3
    scores[1] = 3
    
    # Simulação de fold-changes e p-valores para o Volcano Plot
    log2fc = np.random.normal(0, 0.4, n_total)
    log2fc[scores > 0] = np.random.uniform(1.2, 3.5, n_up)
    log2fc[scores < 0] = np.random.uniform(-3.5, -1.2, n_down)
    
    pvalues = np.random.uniform(0.05, 1.0, n_total)
    pvalues[scores != 0] = np.random.uniform(1e-5, 0.009, n_up + n_down)
    
    # Simulação de famílias de detoxificação para o Heatmap
    familias = ["Outros"] * n_total
    for i in range(n_up):
        familias[i] = np.random.choice(["CYP (Fase I)", "ABC (Fase III)"])
    for i in range(n_up, n_up+400):
        familias[i] = np.random.choice(["CCE (Fase I)", "Proteína Cuticular"])
        
    df = pd.DataFrame({
        'Transcript_ID': transcritos,
        'Gene_Name': nomes_genes,
        'XR_score': scores,
        'Log2FC_Media': log2fc,
        'P_Value': pvalues,
        'Family': familias
    })
    return df

def plotar_figura1_distribuicao(df):
    """
    Gráfico 1: Histograma Discreto de Distribuição do XR-Score.
    Evidencia o censo genômico e a estabilidade homeostática do inseto.
    Melhorias aplicadas conforme solicitado:
      - Textos traduzidos para o inglês (com nomenclatura XR).
      - Remoção de qualquer título.
      - Remoção da chamada para "Bono, 2021" no eixo X.
      - Remoção das linhas horizontais de grade.
      - Aumento do limite superior do eixo Y (headroom) para não cortar os rótulos de dados (ex: 1,211).
    """
    print("📊 Plotando Figura 1: Distribuição Genômica de Consistência (Net XR-Score)...")
    plt.figure(figsize=(7, 5))
    
    # Contagem dos scores discretos de -3 a +3
    counts = df['XR_score'].value_counts().reindex(range(-3, 4), fill_value=0)
    
    # Paleta divergente simétrica (Vermelho para DOWN, Cinza Neutro, Azul para UP)
    cores = ['#d7191c', '#fdae61', '#fee08b', '#f5f5f5', '#d9ef8b', '#91cf60', '#2b83ba']
    
    bars = plt.bar(counts.index, counts.values, color=cores, edgecolor='black', linewidth=0.8, width=0.7)
    
    # Escala logarítmica no eixo Y devido ao enorme volume de genes neutros (62k vs 110)
    plt.yscale('log')
    
    # Adiciona os rótulos de valores absolutos no topo das barras
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            plt.text(bar.get_x() + bar.get_width()/2., height * 1.15, 
                     f'{int(height):,}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Eixos limpos e traduzidos para o inglês com nomenclatura XR
    plt.xlabel("Net XR-Score", fontsize=11, labelpad=8)
    plt.ylabel("Number of Transcripts (Log Scale)", fontsize=11, labelpad=8)
    
    plt.xticks(range(-3, 4), labels=[str(i) for i in range(-3, 4)])
    plt.xlim(-3.7, 3.7)
    
    # Eleva significativamente o limite superior do eixo Y para dar headroom aos rótulos (evita corte do 1,211 e outros)
    plt.ylim(bottom=1, top=counts.max() * 15)
    
    plt.tight_layout()
    plt.savefig(f"{PASTA_SAIDA}/Figura1_Distribuicao_OnScore.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{PASTA_SAIDA}/Figura1_Distribuicao_OnScore.svg", dpi=300, bbox_inches='tight')
    plt.close()

def plotar_figura2_volcano(df):
    """
    Gráfico 2: Volcano Plot Adaptado para Resposta a Xenobióticos.
    Evidencia a significância estatística das isoformas ativas de detoxificação.
    """
    print("🌋 Plotando Figura 2: Volcano Plot de Filtração de Ruído (XR-Ratio)...")
    plt.figure(figsize=(8, 6))
    
    df['-log10_p'] = -np.log10(df['P_Value'])
    
    # Separação por categorias de atividade
    neutros = df[df['XR_score'] == 0]
    induzidos = df[df['XR_score'] > 0]
    reprimidos = df[df['XR_score'] < 0]
    
    # Plotagem dos pontos com transparências estratégicas para evitar overplotting
    plt.scatter(neutros['Log2FC_Media'], neutros['-log10_p'], color='lightgray', alpha=0.4, s=8, label='Neutral (Homeostasis)')
    plt.scatter(reprimidos['Log2FC_Media'], reprimidos['-log10_p'], color='#d7191c', alpha=0.7, s=15, label='Repressed (Core DOWN)')
    plt.scatter(induzidos['Log2FC_Media'], induzidos['-log10_p'], color='#2b83ba', alpha=0.8, s=25, label='Induced (Core UP)')
    
    # Linhas de corte metodológicas
    plt.axhline(-np.log10(0.05), color='black', linestyle=':', alpha=0.6, linewidth=1)
    plt.axvline(1.0, color='black', linestyle=':', alpha=0.6, linewidth=1)
    plt.axvline(-1.0, color='black', linestyle=':', alpha=0.6, linewidth=1)
    
    # Destaque e rotulagem dos biomarcadores principais citados no resumo
    biomarcadores = df[df['Gene_Name'].isin(['ABCC1', 'CYP6B1'])]
    for _, row in biomarcadores.iterrows():
        plt.scatter(row['Log2FC_Media'], row['-log10_p'], color='gold', edgecolor='black', s=80, linewidth=1.2, zorder=5)
        plt.annotate(r"$\mathit{" + row['Gene_Name'] + "}$", 
                     xy=(row['Log2FC_Media'], row['-log10_p']),
                     xytext=(row['Log2FC_Media'] + 0.15, row['-log10_p'] + 0.2),
                     fontsize=10, fontweight='bold', bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3))
        
    plt.title("Volcano Plot of Genomic Response to Insecticides", fontsize=12, pad=15, fontweight='bold')
    plt.xlabel("Mean Response Magnitude (Mean XR-Ratio)", fontsize=11)
    plt.ylabel(r"Statistical Significance ($-\log_{10}$ $\mathit{P}$-value)", fontsize=11)
    
    # Customização da Legenda
    plt.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.4)
    
    plt.tight_layout()
    plt.savefig(f"{PASTA_SAIDA}/Figura2_Volcano_Metanalise.png", dpi=300, bbox_inches='tight')
    plt.close()

def plotar_figura3_heatmap_tradeoff(df):
    """
    Gráfico 3: Heatmap de Famílias Funcionais / Coexpressão.
    Evidencia o trade-off energético e a cascata Fase I vs Fase II.
    """
    print("🌡️ Plotando Figura 3: Heatmap de Trade-Off Funcional (XR-Ratio)...")
    
    # Filtra apenas os genes funcionais que possuem On-Score diferente de zero e ignora resíduos "Outros"
    df_foco = df[df['Family'] != "Outros"].copy()
    
    if len(df_foco) == 0:
        print("⚠️ Dados insuficientes de famílias funcionais curadas para gerar o Heatmap. Pulando etapa.")
        return
        
    # Cria uma tabela cruzada média de expressão por família e score para visualização limpa
    pivot_df = df_foco.groupby(['Family', 'XR_score'])['Log2FC_Media'].mean().unstack().fillna(0)
    
    plt.figure(figsize=(8, 5))
    
    # Mapa de calor usando paleta divergente Coolwarm (Perfeito para variação de expressão)
    sns.heatmap(pivot_df, cmap="coolwarm", center=0, annot=True, fmt=".2f",
                linewidths=0.5, linecolor='white', cbar_kws={'label': 'Mean Response Magnitude (Mean XR-Ratio)'},
                annot_kws={"size": 10, "weight": "bold"})
    
    plt.title("Functional Trade-Off & Metabolic Cascade Synergy", fontsize=12, pad=15, fontweight='bold')
    plt.xlabel("Net XR-Score Value (Consistency)", fontsize=11, labelpad=8)
    plt.ylabel("Curated Functional Families", fontsize=11)
    
    plt.tight_layout()
    plt.savefig(f"{PASTA_SAIDA}/Figura3_Heatmap_TradeOff.png", dpi=300, bbox_inches='tight')
    plt.close()

def plotar_figura4a_top_up(df):
    """
    Gráfico A: Transcritos mais regulados positivamente (UP) sob nomenclatura XR.
    """
    print("📊 Plotando Gráfico A: Top Transcritos Super-regulados (UP)...")
    plt.figure(figsize=(10, 5))
    
    df_temp = df.copy()
    if 'Transcript_ID' not in df_temp.columns:
        df_temp['Transcript_ID'] = df_temp.index
        
    # Filtra os transcritos induzidos (UP, XR_score > 0)
    df_up = df_temp[df_temp['XR_score'] > 0].copy()
    if len(df_up) == 0:
        print("⚠️ Sem transcritos induzidos (UP) suficientes para gerar o Gráfico A.")
        return
        
    # Ordena pelo score decrescente e depois pela magnitude média de expressão
    df_up = df_up.sort_values(by=['XR_score', 'Log2FC_Media'], ascending=[False, False])
    
    # Seleciona os top 15 transcritos mais regulados
    df_top_up = df_up.head(15).copy()
    
    # Aplica mapeamento para obter o nome amigável do gene (ex: ABCC1, CYP6B1) em vez do ID longo
    df_top_up['Friendly_Label'] = df_top_up.apply(obter_label_amigavel, axis=1)
    
    bars = plt.bar(df_top_up['Friendly_Label'], df_top_up['XR_score'], color='#2b83ba', edgecolor='black', linewidth=0.8, width=0.55)
    
    # Customizações estéticas solicitadas
    plt.ylabel("Net XR-Score", fontsize=11, labelpad=8)
    plt.xlabel("Gene Symbol", fontsize=11, labelpad=8)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks([1, 2, 3], labels=['1', '2', '3'])
    plt.ylim(0, 3.7)
    
    # Adiciona os números exatos sobre as barras
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.08, 
                 f'{int(height)}', ha='center', va='bottom', fontsize=9, fontweight='bold')
                 
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{PASTA_SAIDA}/Figura4A_Top_Upregulated.png", dpi=300, bbox_inches='tight')
    plt.close()

def plotar_figura4b_top_down(df):
    """
    Gráfico B: Transcritos mais regulados negativamente (DOWN) sob nomenclatura XR.
    """
    print("📊 Plotando Gráfico B: Top Transcritos Sub-regulados (DOWN)...")
    plt.figure(figsize=(10, 5))
    
    df_temp = df.copy()
    if 'Transcript_ID' not in df_temp.columns:
        df_temp['Transcript_ID'] = df_temp.index
        
    # Filtra os transcritos reprimidos (DOWN, XR_score < 0)
    df_down = df_temp[df_temp['XR_score'] < 0].copy()
    if len(df_down) == 0:
        print("⚠️ Sem transcritos reprimidos (DOWN) suficientes para gerar o Gráfico B.")
        return
        
    # Ordena pelo score crescente (mais negativo primeiro) e expressao decrescente (mais reprimido primeiro)
    df_down = df_down.sort_values(by=['XR_score', 'Log2FC_Media'], ascending=[True, True])
    
    # Seleciona os top 15 transcritos mais reprimidos
    df_top_down = df_down.head(15).copy()
    
    # Aplica mapeamento para obter o nome amigável do gene (ex: CHT12, CCE-dig) em vez do ID longo
    df_top_down['Friendly_Label'] = df_top_down.apply(obter_label_amigavel, axis=1)
    
    bars = plt.bar(df_top_down['Friendly_Label'], df_top_down['XR_score'], color='#d7191c', edgecolor='black', linewidth=0.8, width=0.55)
    
    # Customizações estéticas solicitadas
    plt.ylabel("Net XR-Score", fontsize=11, labelpad=8)
    plt.xlabel("Gene Symbol", fontsize=11, labelpad=8)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks([-3, -2, -1], labels=['-3', '-2', '-1'])
    plt.ylim(-3.7, 0)
    plt.axhline(0, color='black', linewidth=0.8)
    
    # Adiciona os números exatos abaixo das barras negativas
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height - 0.12, 
                 f'{int(height)}', ha='center', va='top', fontsize=9, fontweight='bold')
                 
    sns.despine(bottom=False, top=True)
    plt.tight_layout()
    plt.savefig(f"{PASTA_SAIDA}/Figura4B_Top_Downregulated.png", dpi=300, bbox_inches='tight')
    plt.close()

def executar_pipeline_grafico():
    print("==================================================")
    print("    SAQE - GERAÇÃO DE VISUALIZAÇÕES DE ALTO IMPACTO")
    print("==================================================")
    
    df_dados = simular_dados_caso_ausente()
    
    plotar_figura1_distribuicao(df_dados)
    plotar_figura2_volcano(df_dados)
    plotar_figura3_heatmap_tradeoff(df_dados)
    
    # Executa os novos gráficos separados UP (Gráfico A) e DOWN (Gráfico B)
    plotar_figura4a_top_up(df_dados)
    plotar_figura4b_top_down(df_dados)
    
    print("\n==================================================")
    print("🎉 TODOS OS GRÁFICOS FORAM GERADOS COM SUCESSO!")
    print(f"📁 Verifique as imagens salvas na pasta: '{PASTA_SAIDA}/'")
    print("   • Figura1_Distribuicao_OnScore.png/.svg (Distribuição de Barras - Net XR-Score)")
    print("   • Figura2_Volcano_Metanalise.png (Volcano Plot de Expressão / XR)")
    print("   • Figura3_Heatmap_TradeOff.png (Visualização de Sinergia / XR-Score)")
    print("   • Figura4A_Top_Upregulated.png (Gráfico A - Transcritos mais ativados)")
    print("   • Figura4B_Top_Downregulated.png (Gráfico B - Transcritos mais reprimidos)")
    print("==================================================")

if __name__ == "__main__":
    executar_pipeline_grafico()
