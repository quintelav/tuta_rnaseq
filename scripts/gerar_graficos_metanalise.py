#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Gerador de Gráficos Científicos de Alta Resolução (Meta-Volcano Plot & Clustered Heatmaps)
Executa uma metanálise estatística combinando p-valores via Método de Fisher e calcula o Meta-Log2FC.
Gera figuras profissionais (300 DPI) em inglês e sem sobreposição de rótulos.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

# Configura o Matplotlib para rodar em modo 'headless' (sem necessidade de servidor X11/interface gráfica)
# Essencial para WSL, servidores SSH e ambientes de nuvem.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Dicionário de tradução para garantir que todos os elementos gráficos finais estejam em inglês técnico
MAPA_TRADUCAO_FAMILIAS = {
    "Fase I: Cytochrome P450 (CYP)": "Phase I: Cytochrome P450 (CYP)",
    "Fase I: Carboxylesterase (CCE)": "Phase I: Carboxylesterase (CCE)",
    "Fase II: Glutathione S-transferase (GST)": "Phase II: Glutathione S-transferase (GST)",
    "Fase II: UDP-glucuronosyltransferase (UGT)": "Phase II: UDP-glucuronosyltransferase (UGT)",
    "Fase III: Transportadores ABC": "Phase III: ABC Transporters",
    "Antioxidante: Superoxide Dismutase (SOD)": "Antioxidant: Superoxide Dismutase (SOD)",
    "Antioxidante: Peroxidase": "Antioxidant: Peroxidase",
    "Antioxidante: Catalase (CAT)": "Antioxidant: Catalase (CAT)",
    "Outros / Não Identificado": "Others / Unidentified",
    "Outros Processos Fisiológicos": "Other Physiological Processes"
}

# Paleta de cores sofisticada em inglês, adaptada para periódicos internacionais de alto impacto
CORES_FAMILIAS = {
    "Phase I: Cytochrome P450 (CYP)": "#8e44ad",          # Roxo
    "Phase I: Carboxylesterase (CCE)": "#2980b9",         # Azul
    "Phase II: Glutathione S-transferase (GST)": "#27ae60", # Verde
    "Phase II: UDP-glucuronosyltransferase (UGT)": "#d35400",# Laranja
    "Phase III: ABC Transporters": "#c0392b",           # Vermelho
    "Antioxidant: Superoxide Dismutase (SOD)": "#16a085",# Teal
    "Antioxidant: Peroxidase": "#1abc9c",                # Turquesa claro
    "Antioxidant: Catalase (CAT)": "#2ecc71",            # Verde claro
    "Others / Unidentified": "#7f8c8d",               # Cinza escuro
    "Other Physiological Processes": "#bdc3c7"            # Cinza claro
}

def carregar_e_calcular_pvalores():
    """
    Lê a matriz de TPM unificada com TODOS os transcritos e calcula p-valores de forma
    vetorizada extremamente rápida usando Welch's t-test. Combina-os via método de Fisher.
    """
    print("🧪 A ler dados e a calcular p-valores estatísticos vetorizados para todo o genoma...")
    
    matriz_tpm_path = "quant/matriz_tpm_geral.csv"
    metadata_path = "metadata_geral.csv"
    anotado_path = "quant/metanalise_completa_anotada.csv"
    bono_path = "quant/metanalise_completa_bono.csv"
    
    if not all(os.path.exists(f) for f in [matriz_tpm_path, metadata_path, anotado_path, bono_path]):
        print("❌ [ERRO] Arquivos de entrada necessários não encontrados na pasta 'quant/' ou 'reference/'!")
        return None
        
    df_tpm = pd.read_csv(matriz_tpm_path, index_col=0)
    df_meta = pd.read_csv(metadata_path)
    df_anotado = pd.read_csv(anotado_path, index_col=0)
    df_bono_completo = pd.read_csv(bono_path, index_col=0)
    
    # Alinhamento e interseção de índices extremamente segura para evitar KeyError de desalinhamento
    comum_index = df_bono_completo.index.intersection(df_tpm.index)
    if len(comum_index) == 0:
        print("❌ [ERRO] Sem correspondência de índices entre a matriz de TPM e os resultados de Bono!")
        return None
        
    # Iniciamos com a matriz completa de todos os genes correspondentes para termos o background real
    df_volcano = pd.DataFrame(index=comum_index)
    
    # Mapeamos as anotações das famílias de detox para todos os genes (padrão: Outros transcritos)
    df_volcano['Detox_Family'] = df_volcano.index.map(df_anotado['Detox_Family']).fillna("Other Genomic Transcripts")
    df_volcano['Detox_Family'] = df_volcano['Detox_Family'].map(MAPA_TRADUCAO_FAMILIAS).fillna(df_volcano['Detox_Family'])
    
    # Mapeamos também o Bono_UP_Score e Bono_DOWN_Score para filtragens posteriores
    df_volcano['Bono_UP_Score'] = df_bono_completo.loc[comum_index, 'Bono_UP_Score']
    df_volcano['Bono_DOWN_Score'] = df_bono_completo.loc[comum_index, 'Bono_DOWN_Score']
    
    estudos = df_meta['study_id'].unique()
    colunas_p = []
    colunas_l2fc = []
    
    for estudo in estudos:
        print(f"   🔹 A processar estatísticas de todo o genoma para o estudo: '{estudo}'...")
        meta_sub = df_meta[df_meta['study_id'] == estudo]
        
        amostras_control = meta_sub[meta_sub['condition'] == 'control']['sample_id'].tolist()
        amostras_treat = meta_sub[meta_sub['condition'] == 'treatment']['sample_id'].tolist()
        
        amostras_control = [x for x in amostras_control if x in df_tpm.columns]
        amostras_treat = [x for x in amostras_treat if x in df_tpm.columns]
        
        if len(amostras_control) < 2 or len(amostras_treat) < 2:
            continue
            
        # Filtra de forma totalmente segura a matriz de TPM para as amostras e genes correspondentes
        tpm_ctrl_sub = df_tpm.loc[comum_index, amostras_control].values
        tpm_treat_sub = df_tpm.loc[comum_index, amostras_treat].values
        
        # Teste t de Welch vetorizado (calcula os 64.000 genes instantaneamente)
        _, p_vals = stats.ttest_ind(tpm_treat_sub, tpm_ctrl_sub, axis=1, equal_var=False)
        p_vals = np.nan_to_num(p_vals, nan=1.0)
        
        # Log2FC vetorizado com pseudocount de +1.0
        mean_ctrl = tpm_ctrl_sub.mean(axis=1)
        mean_treat = tpm_treat_sub.mean(axis=1)
        l2fc = np.log2((mean_treat + 1.0) / (mean_ctrl + 1.0))
        
        p_col = f'pvalue_{estudo}'
        l2fc_col = f'Log2FC_{estudo}'
        
        df_volcano[p_col] = p_vals
        df_volcano[l2fc_col] = l2fc
        
        colunas_p.append(p_col)
        colunas_l2fc.append(l2fc_col)

    # --- METANÁLISE CONSOLIDADA (Método de Fisher) ---
    print("\n🧬 A consolidar metanálise via Método de Fisher para combinação de p-valores...")
    
    # Prepara matrizes para cálculo vetorizado de Fisher
    p_matrix = df_volcano[colunas_p].values
    p_matrix = np.clip(p_matrix, 1e-300, 1.0)
    
    # Qui-quadrado de Fisher: -2 * sum(ln(p))
    chi2_stats = -2.0 * np.sum(np.log(p_matrix), axis=1)
    combined_p = stats.chi2.sf(chi2_stats, df=2 * len(colunas_p))
    combined_p = np.nan_to_num(combined_p, nan=1.0)
    
    mean_l2fc = df_volcano[colunas_l2fc].mean(axis=1).values
    
    df_volcano['Meta_pvalue'] = combined_p
    df_volcano['Meta_Log2FC'] = mean_l2fc
    df_volcano['Meta_pvalue'] = df_volcano['Meta_pvalue'].replace(0, 1e-300)
    df_volcano['neg_log10_meta_p'] = -np.log10(df_volcano['Meta_pvalue'])
    
    return df_volcano

def evitar_colisao_labels(x, y, placed_coords, min_dx=0.8, min_dy=2.5):
    """
    Algoritmo de Prevenção de Colisão de etiquetas de texto.
    """
    for px, py in placed_coords:
        if abs(x - px) < min_dx and abs(y - py) < min_dy:
            return True
    return False

def plotar_meta_volcano_plot(df_volcano):
    """
    Gera um único e deslumbrante Volcano Plot Consolidado (Meta-Volcano Plot)
    mostrando a totalidade dos 64k+ genes de fundo e destacando os mais regulados.
    """
    print("\n🚀 A desenhar o Meta-Volcano Plot consolidado de alta resolução...")
    
    plt.figure(figsize=(10, 8), dpi=300)
    sns.set_theme(style="ticks")
    
    # 1. Desenha o fundo de todos os transcritos gerais sem atividade marcante (Genoma Completo)
    plt.scatter(df_volcano['Meta_Log2FC'], df_volcano['neg_log10_meta_p'], 
                c='#bdc3c7', alpha=0.25, s=6, label='Other genomic transcripts')
    
    etiquetas_desenhadas = []
    
    # 2. Desenha e destaca as famílias de detoxificação de interesse científico
    for familia, cor in CORES_FAMILIAS.items():
        if "Others" in familia or "Other" in familia:
            continue
        sub_fam = df_volcano[df_volcano['Detox_Family'] == familia]
        if len(sub_fam) > 0:
            plt.scatter(sub_fam['Meta_Log2FC'], sub_fam['neg_log10_meta_p'], 
                        c=cor, alpha=0.9, s=40, edgecolors='black', linewidths=0.3, label=familia)
            
            # Rotula de forma inteligente os IDs das maiores ativações e genes mais regulados de forma extrema
            sub_fam_sorted = sub_fam.sort_values(by='neg_log10_meta_p', ascending=False).head(3)
            for idx, row in sub_fam_sorted.iterrows():
                x_coord = row['Meta_Log2FC']
                y_coord = row['neg_log10_meta_p']
                
                if y_coord > 1.3: # p-value < 0.05
                    if not evitar_colisao_labels(x_coord, y_coord, etiquetas_desenhadas, min_dx=1.0, min_dy=4.0):
                        etiquetas_desenhadas.append((x_coord, y_coord))
                        
                        gene_name = idx.split('|')[-1] if '|' in idx else idx
                        gene_name = gene_name.replace(".t1", "")
                        
                        plt.annotate(
                            gene_name, 
                            (x_coord, y_coord),
                            textcoords="offset points", 
                            xytext=(0, 6), 
                            ha='center', 
                            fontsize=6.5, 
                            fontweight='bold',
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="slategray", alpha=0.85, lw=0.3)
                        )
                        
    # 3. Adiciona rotulagem extra para os genes com os maiores Log2FC absolutos do genoma de detoxificação para destacar extremos
    familias_detox_validas = [
        "Phase I: Cytochrome P450 (CYP)",
        "Phase I: Carboxylesterase (CCE)",
        "Phase II: Glutathione S-transferase (GST)",
        "Phase II: UDP-glucuronosyltransferase (UGT)",
        "Phase III: ABC Transporters",
        "Antioxidant: Superoxide Dismutase (SOD)",
        "Antioxidant: Peroxidase",
        "Antioxidant: Catalase (CAT)"
    ]
    df_extremos = df_volcano[df_volcano['Detox_Family'].isin(familias_detox_validas)].copy()
    top_up = df_extremos.sort_values(by='Meta_Log2FC', ascending=False).head(2)
    top_down = df_extremos.sort_values(by='Meta_Log2FC', ascending=True).head(2)
    
    for df_ext in [top_up, top_down]:
        for idx, row in df_ext.iterrows():
            x_coord = row['Meta_Log2FC']
            y_coord = row['neg_log10_meta_p']
            if not evitar_colisao_labels(x_coord, y_coord, etiquetas_desenhadas, min_dx=1.0, min_dy=4.0):
                etiquetas_desenhadas.append((x_coord, y_coord))
                gene_name = idx.split('|')[-1] if '|' in idx else idx
                gene_name = gene_name.replace(".t1", "")
                plt.annotate(
                    gene_name, 
                    (x_coord, y_coord),
                    textcoords="offset points", 
                    xytext=(0, 6), 
                    ha='center', 
                    fontsize=6.5, 
                    fontweight='bold',
                    color='red' if x_coord > 0 else 'blue',
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="slategray", alpha=0.85, lw=0.3)
                )
    
    # Linhas de corte biológicas padrões
    plt.axhline(y=-np.log10(0.05), color='#7f8c8d', linestyle='--', linewidth=0.8, alpha=0.7) # p-value = 0.05
    plt.axvline(x=1.0, color='#7f8c8d', linestyle='--', linewidth=0.8, alpha=0.7)  # Fold Change = 2 (super-expresso)
    plt.axvline(x=-1.0, color='#7f8c8d', linestyle='--', linewidth=0.8, alpha=0.7) # Fold Change = 2 (sub-expresso)
    
    # Títulos e labels sofisticados em Inglês Técnico
    plt.title("Meta-Volcano Plot - Complete Genomic Landscape under Insecticide Stress\n(Fisher's Combined P-Value Method across Abamectin, Spinosad & Emamectin)", 
              fontsize=11, fontweight='bold', pad=15)
    plt.xlabel(r"Consolidated Magnitude of Expression ($\text{Mean }\log_2 \text{ Fold Change}$)", fontsize=10)
    plt.ylabel(r"Combined Statistical Significance ($-\log_{10} \text{ Combined p-value}$)", fontsize=10)
    
    plt.xlim(-8, 8) # Limita eixo X de forma simétrica para realçar a dispersão dos pontos
    plt.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, edgecolor='#e2e8f0', fontsize=8)
    sns.despine(trim=True)
    plt.tight_layout()
    
    fig_path = "quant/meta_volcano_plot.png"
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Meta-Volcano Plot completo guardado em '{fig_path}'!")

def plotar_heatmaps_detox(df_dados):
    """
    Gera duas abordagens de Heatmaps profissionais livres de qualquer sobreposição:
    1. Family-Level Heatmap: Agrupado por classe funcional para resumos visuais.
    2. Dot/Bubble Plot de Co-expressão (Top 25): A melhor alternativa moderna a heatmaps embolados!
    """
    print("\n🚀 A desenhar os Heatmaps de alta resolução livres de sobreposição...")
    
    # Filtra os candidatos principais (On-Score >= 2)
    df_candidatos = df_dados[df_dados['Detox_Family'] != 'Other Genomic Transcripts'].copy()
    df_candidatos = df_candidatos[(df_candidatos['Bono_UP_Score'] >= 2) | (df_candidatos['Bono_DOWN_Score'] >= 2)].copy()
    
    if len(df_candidatos) == 0:
        print("⚠️  Nenhum candidato com On-Score >= 2 encontrado para compor os Heatmaps.")
        return
        
    colunas_ratios = ['Log2FC_abamectina', 'Log2FC_spinosad', 'Log2FC_emamectina']
    colunas_pvalues = ['pvalue_abamectina', 'pvalue_spinosad', 'pvalue_emamectina']
    nomes_colunas_ingles = ['Abamectin', 'Spinosad', 'Emamectin']
    
    # -------------------------------------------------------------------------
    # ABORDAGEM 1: HEATMAP AGRUPADO POR FAMÍLIA (Zero sobreposição, ideal para o texto)
    # -------------------------------------------------------------------------
    print("   🔹 A compilar Heatmap agrupado por famílias metabólicas...")
    df_grouped = df_candidatos.groupby('Detox_Family')[colunas_ratios].mean()
    df_grouped.columns = nomes_colunas_ingles
    
    df_grouped = df_grouped.drop(index=["Other Physiological Processes", "Others / Unidentified"], errors='ignore')
    
    plt.figure(figsize=(8, 4.5), dpi=300)
    sns.heatmap(
        df_grouped,
        cmap='RdYlBu_r',
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        center=0.0,
        vmin=-3.0, vmax=3.0,
        cbar_kws={'label': 'Mean ' + r'$\log_2$ Fold Change'}
    )
    plt.title("Family-Level Consolidated Detoxification Heatmap\n(Mean Expression Pattern of Candidate Genes per Family)", 
              fontsize=10, fontweight='bold', pad=15)
    plt.xlabel("Insecticide Stress Treatments", fontsize=9)
    plt.ylabel("Cellular Detoxification Families", fontsize=9)
    plt.xticks(fontsize=9, fontweight='bold')
    plt.yticks(fontsize=9)
    plt.tight_layout()
    
    fig_path_group = "quant/heatmap_detox_families.png"
    plt.savefig(fig_path_group, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Heatmap agrupado por famílias guardado em '{fig_path_group}'!")
    
    # -------------------------------------------------------------------------
    # ABORDAGEM 2: BUBBLE PLOT (DOT PLOT) DOS TOP 25 CANDIDATOS (Sem qualquer embolado!)
    # Representa os dados em uma grade cartesiana perfeita eliminando sobreposições.
    # -------------------------------------------------------------------------
    print("   🔹 A compilar Dot Plot (Bubble Heatmap) de alta clareza para os Top 25 candidatos...")
    
    df_top = df_candidatos.sort_values(by='Meta_pvalue', ascending=True).head(25).copy()
    
    # Cria os nomes de exibição explicativos combinando a família e o ID encurtado
    nomes_rotulados = []
    for idx, row in df_top.iterrows():
        id_curto = idx.split('|')[-1] if '|' in idx else idx
        id_curto = id_curto.replace(".t1", "")
        
        fam = row['Detox_Family']
        abrev = "CYP" if "P450" in fam else "CCE" if "Carboxylesterase" in fam else "GST" if "Glutathione" in fam else "UGT" if "UDP" in fam else "ABC" if "ABC" in fam else "SOD" if "Dismutase" in fam else "PRX" if "Peroxidase" in fam else "DET"
        nomes_rotulados.append(f"[{abrev}] {id_curto}")
        
    df_top['Gene_Label'] = nomes_rotulados
    
    # Estruturamos os dados em formato longo (Melt) para plotar com Scatter de forma perfeita
    dados_melted = []
    for i, row in df_top.iterrows():
        for col_idx, (col_ratio, col_p) in enumerate(zip(colunas_ratios, colunas_pvalues)):
            treatment_name = nomes_colunas_ingles[col_idx]
            fc_val = row[col_ratio]
            pval = row[col_p]
            neg_log_p = -np.log10(pval) if pval > 0 else 0
            # Adiciona pequeno limite visual de tamanho mínimo para bolhas visíveis
            bubble_size = np.clip(neg_log_p * 45, 10, 450) 
            
            dados_melted.append({
                'Gene': row['Gene_Label'],
                'Treatment': treatment_name,
                'Log2FC': fc_val,
                'BubbleSize': bubble_size,
                'pvalue': pval
            })
            
    df_melt = pd.DataFrame(dados_melted)
    
    # Configura a estrutura da figura do Dot Plot
    plt.figure(figsize=(9, 11), dpi=300)
    sns.set_theme(style="whitegrid")
    
    # Plota os círculos representativos com mapa divergente de cores
    scatter = plt.scatter(
        x=df_melt['Treatment'],
        y=df_melt['Gene'],
        s=df_melt['BubbleSize'],
        c=df_melt['Log2FC'],
        cmap='RdYlBu_r',
        alpha=0.9,
        edgecolors='black',
        linewidths=0.5,
        vmin=-4.0, vmax=4.0
    )
    
    # Formata a barra lateral de cores (Legenda de intensidade de expressão)
    cbar = plt.colorbar(scatter, shrink=0.5, aspect=12, pad=0.03)
    cbar.set_label(r"$\log_2$ Fold Change", fontsize=9, fontweight='bold')
    
    # Cria uma legenda manual para ilustrar os tamanhos das bolhas (Significância p-value) - Corrigido para 'slategray'
    for size_label, size_val in [("p < 0.05", 1.3), ("p < 0.01", 2.0), ("p < 0.001", 3.0)]:
        plt.scatter([], [], c='slategray', alpha=0.7, s=size_val * 45, label=size_label, edgecolors='black', linewidths=0.5)
        
    plt.legend(
        title="Significance Threshold", 
        title_fontsize=9, 
        loc="upper right", 
        bbox_to_anchor=(1.35, 1.0),
        frameon=True, 
        facecolor='white', 
        edgecolor='#e2e8f0',
        fontsize=8
    )
    
    plt.title("Co-expression Profile of Top 25 Cellular Detoxification Candidate Transcripts\n(Consolidated Metanalysis - Dot Plot Visualization)", 
              fontsize=10, fontweight='bold', pad=20)
    plt.xlabel("Insecticide Stress Treatments", fontsize=9, fontweight='bold', labelpad=10)
    plt.ylabel("Candidate Genes & Isoforms", fontsize=9, fontweight='bold')
    plt.xticks(fontsize=9, fontweight='bold')
    plt.yticks(fontsize=8, fontname='monospace')
    
    plt.tight_layout()
    fig_path_top = "quant/heatmap_top_candidates.png" # Salvamos com o mesmo nome para manter compatibilidade com o pipeline
    plt.savefig(fig_path_top, bbox_inches='tight')
    plt.close()
    print(f"   ✅ Dot Plot (Bubble Heatmap) de alta legibilidade guardado em '{fig_path_top}'!")

def executar_pipeline_visualizacao():
    df_volcano = carregar_e_calcular_pvalores()
    if df_volcano is not None:
        plotar_meta_volcano_plot(df_volcano)
        plotar_heatmaps_detox(df_volcano)
        print("\n==================================================")
        print("🎉 TODOS OS GRÁFICOS CIENTÍFICOS FORAM GERADOS!")
        print("==================================================")
        print("📂 Guardados com sucesso na sua pasta 'quant/':")
        print("   👉 'meta_volcano_plot.png'      -> Volcano Plot Unificado do Genoma Inteiro (English)")
        print("   👉 'heatmap_detox_families.png' -> Heatmap Simplificado por Família (Sem sobreposição)")
        print("   👉 'heatmap_top_candidates.png' -> Dot Plot de Alta Resolução das Top 25 Isoformas (Sem embolar)")

if __name__ == "__main__":
    executar_pipeline_visualizacao()
