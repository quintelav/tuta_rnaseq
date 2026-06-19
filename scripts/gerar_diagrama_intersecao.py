#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Gerador de Diagrama de Interseção de Venn de Detoxificação
Calcula a interseção de genes ativos (Log2Ratio/FC >= 1.0) entre os três 
inseticidas utilizando diretamente os dados anotados reais e plota o gráfico de alta resolução.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Dicionário de cores sofisticado para publicação
CORES_ESTUDOS = {
    'Abamectin': '#8e44ad',  # Roxo
    'Spinosad': '#2980b9',   # Azul
    'Emamectin': '#27ae60',  # Verde
    'Mixed': '#f1c40f',      # Amarelo/Dourado para interseções
}

def extrair_nome_gene(annotation_text):
    """
    Extrai o símbolo ou nome curto do gene dos colchetes (ex: [CYP6B1] -> CYP6B1).
    """
    if not isinstance(annotation_text, str):
        return "Uncharacterized"
    match = re.search(r'\[([^\]]+)\]', annotation_text)
    if match:
        return match.group(1)
    
    # Fallback
    text_lower = annotation_text.lower()
    if "cytochrome p450" in text_lower or "cyp" in text_lower:
        return "CYP-like"
    if "carboxylesterase" in text_lower or "esterase" in text_lower:
        return "CCE-like"
    if "glutathione s-transferase" in text_lower or "gst" in text_lower:
        return "GST-like"
    return "Detox-gene"

def carregar_dados_e_definir_sets():
    """
    Lê o ficheiro da metanálise anotada real e filtra os genes ativos sob cada inseticida.
    Garante mapeamento robusto tolerante a diferenças de nomes de colunas (Log2Ratio vs Log2FC).
    """
    caminho_graficos = "quant/metanalise_completa_anotada.csv"
    if not os.path.exists(caminho_graficos):
        # Fallback caso tenha sido rodado apenas o pipeline de graficos
        caminho_graficos = "quant/metanalise_completa_bono.csv"
        
    if not os.path.exists(caminho_graficos):
        print("❌ [ERRO] Não foi encontrada nenhuma tabela de resultados na pasta 'quant/'.")
        print("Execute primeiro a anotação funcional: 'python3 anotar_resultados.py'")
        return None

    df = pd.read_csv(caminho_graficos, index_col=0)
    
    # Filtra apenas famílias de detoxificação válidas
    familias_validas = [
        "Fase I: Cytochrome P450 (CYP)", "Phase I: Cytochrome P450 (CYP)",
        "Fase I: Carboxylesterase (CCE)", "Phase I: Carboxylesterase (CCE)",
        "Fase II: Glutathione S-transferase (GST)", "Phase II: Glutathione S-transferase (GST)",
        "Fase II: UDP-glucuronosyltransferase (UGT)", "Phase II: UDP-glucuronosyltransferase (UGT)",
        "Fase III: Transportadores ABC", "Phase III: ABC Transporters",
        "Antioxidante: Superoxide Dismutase (SOD)", "Antioxidant: Superoxide Dismutase (SOD)",
        "Antioxidante: Peroxidase", "Antioxidant: Peroxidase",
        "Antioxidante: Catalase (CAT)", "Antioxidant: Catalase (CAT)"
    ]
    df_detox = df[df['Detox_Family'].isin(familias_validas)].copy()
    
    # Define as colunas de fold change de forma robusta e dinâmica (Log2Ratio ou Log2FC)
    col_abam = 'Log2Ratio_abamectina' if 'Log2Ratio_abamectina' in df_detox.columns else 'Log2FC_abamectina'
    col_spin = 'Log2Ratio_spinosad' if 'Log2Ratio_spinosad' in df_detox.columns else 'Log2FC_spinosad'
    col_emam = 'Log2Ratio_emamectina' if 'Log2Ratio_emamectina' in df_detox.columns else 'Log2FC_emamectina'
    
    # Se as colunas de p-valor calculadas existem, usamos o filtro estrito (p < 0.05), caso contrário usamos apenas o fold change (FC >= 1.0)
    p_abam = 'pvalue_abamectina' if 'pvalue_abamectina' in df_detox.columns else None
    p_spin = 'pvalue_spinosad' if 'pvalue_spinosad' in df_detox.columns else None
    p_emam = 'pvalue_emamectina' if 'pvalue_emamectina' in df_detox.columns else None
    
    # Define ativação de super-expressão para cada inseticida
    if p_abam and p_spin and p_emam:
        set_abam = set(df_detox[(df_detox[col_abam] >= 1.0) & (df_detox[p_abam] < 0.05)].index)
        set_spin = set(df_detox[(df_detox[col_spin] >= 1.0) & (df_detox[p_spin] < 0.05)].index)
        set_emam = set(df_detox[(df_detox[col_emam] >= 1.0) & (df_detox[p_emam] < 0.05)].index)
    else:
        set_abam = set(df_detox[df_detox[col_abam] >= 1.0].index)
        set_spin = set(df_detox[df_detox[col_spin] >= 1.0].index)
        set_emam = set(df_detox[df_detox[col_emam] >= 1.0].index)
    
    return set_abam, set_spin, set_emam, df_detox, col_abam, col_spin, col_emam

def plotar_diagrama_venn_custom(set_abam, set_spin, set_emam, df_detox, col_abam, col_spin, col_emam):
    """
    Desenha um diagrama de Venn de 3 círculos totalmente customizado e sem 
    sobreposição de rótulos de texto, ideal para publicações.
    """
    print("🚀 A calcular interseções genómicas de detoxificação...")
    
    # Cálculo das interseções exclusivas
    abam_only = set_abam - set_spin - set_emam
    spin_only = set_spin - set_abam - set_emam
    emam_only = set_emam - set_abam - set_spin
    
    abam_spin = (set_abam & set_spin) - set_emam
    abam_emam = (set_abam & set_emam) - set_spin
    spin_emam = (set_spin & set_emam) - set_abam
    
    shared_all = set_abam & set_spin & set_emam
    
    print(f"   🔹 Abamectin only: {len(abam_only)}")
    print(f"   🔹 Spinosad only: {len(spin_only)}")
    print(f"   🔹 Emamectin only: {len(emam_only)}")
    print(f"   🔹 Shared by All Three (Core): {len(shared_all)}")
    
    # Configura a estrutura da figura (Layout de duas colunas: Venn à esquerda, Tabela à direita)
    fig, (ax_venn, ax_table) = plt.subplots(1, 2, figsize=(14, 8), dpi=300, gridspec_kw={'width_ratios': [1.2, 1]})
    
    # --- COLUNA 1: VENN DIAGRAM ---
    ax_venn.set_xlim(-3.5, 3.5)
    ax_venn.set_ylim(-3.5, 3.5)
    ax_venn.set_aspect('equal')
    ax_venn.axis('off')
    
    # Desenha os círculos transparentes com posições perfeitas
    circ_abam = patches.Circle((0, 0.8), 1.8, facecolor=CORES_ESTUDOS['Abamectin'], alpha=0.35, edgecolor='black', linewidth=1)
    circ_spin = patches.Circle((-1.0, -0.6), 1.8, facecolor=CORES_ESTUDOS['Spinosad'], alpha=0.35, edgecolor='black', linewidth=1)
    circ_emam = patches.Circle((1.0, -0.6), 1.8, facecolor=CORES_ESTUDOS['Emamectin'], alpha=0.35, edgecolor='black', linewidth=1)
    
    ax_venn.add_patch(circ_abam)
    ax_venn.add_patch(circ_spin)
    ax_venn.add_patch(circ_emam)
    
    # Rótulos externos dos círculos (Títulos dos tratamentos)
    ax_venn.text(0, 2.8, "Abamectin\n(Stress A)", fontsize=11, fontweight='bold', ha='center', va='center', color='#5e2a74')
    ax_venn.text(-2.4, -1.8, "Spinosad\n(Stress B)", fontsize=11, fontweight='bold', ha='center', va='center', color='#1b4f72')
    ax_venn.text(2.4, -1.8, "Emamectin\n(Stress C)", fontsize=11, fontweight='bold', ha='center', va='center', color='#145a32')
    
    # Injeta os valores numéricos em posições espaciais calculadas
    # 1. Valores Únicos
    ax_venn.text(0, 1.7, f"{len(abam_only)}", fontsize=14, fontweight='bold', ha='center', va='center')
    ax_venn.text(-1.8, -0.9, f"{len(spin_only)}", fontsize=14, fontweight='bold', ha='center', va='center')
    ax_venn.text(1.8, -0.9, f"{len(emam_only)}", fontsize=14, fontweight='bold', ha='center', va='center')
    
    # 2. Interseções Duplas
    ax_venn.text(-0.9, 0.5, f"{len(abam_spin)}", fontsize=12, fontweight='bold', ha='center', va='center', color='#2c3e50')
    ax_venn.text(0.9, 0.5, f"{len(abam_emam)}", fontsize=12, fontweight='bold', ha='center', va='center', color='#2c3e50')
    ax_venn.text(0, -1.3, f"{len(spin_emam)}", fontsize=12, fontweight='bold', ha='center', va='center', color='#2c3e50')
    
    # 3. Interseção Tripla (O Núcleo "Core")
    ax_venn.text(0, -0.1, f"{len(shared_all)}", fontsize=16, fontweight='bold', ha='center', va='center', color='#7b241c',
                 bbox=dict(boxstyle="circle,pad=0.2", fc="#f9ebd2", ec="#7b241c", lw=1.2))
    
    ax_venn.set_title("Venn Diagram of Active Detoxification Genes\n(Log2FC >= 1.0)", 
                     fontsize=12, fontweight='bold', pad=10)
    
    # --- COLUNA 2: TABELA DE CANDIDATOS NÚCLEO (Core Shared com Fallback Inteligente) ---
    ax_table.axis('off')
    
    # Define a lista de genes para preencher a tabela aplicando o critério de hierarquia de fallback
    genes_para_tabela = list(shared_all)
    titulo_tabela = f"Core Shared Candidates (All 3 Insecticides) ({len(shared_all)} total)"
    
    if len(genes_para_tabela) == 0:
        # Fallback 1: Genes compartilhados por pelo menos 2 inseticidas (interseções duplas)
        compartilhados_dois = (set_abam & set_spin) | (set_abam & set_emam) | (set_spin & set_emam)
        if len(compartilhados_dois) > 0:
            genes_para_tabela = list(compartilhados_dois)
            titulo_tabela = f"Shared in at least 2 Insecticides (Top 12 of {len(compartilhados_dois)})"
        else:
            # Fallback 2: Todos os genes ativos gerais de detox com maior expressão média
            ativos_geral = set_abam | set_spin | set_emam
            if len(ativos_geral) > 0:
                genes_para_tabela = list(ativos_geral)
                titulo_tabela = f"Top Active Detoxification Candidates (Top 12 of {len(ativos_geral)})"

    if len(genes_para_tabela) > 0:
        # Recupera os dados detalhados para os genes selecionados
        df_shared = df_detox.loc[genes_para_tabela].copy()
        df_shared['Gene'] = df_shared['Annotation'].apply(extrair_nome_gene)
        df_shared['Mean_FC'] = df_shared[[col_abam, col_spin, col_emam]].mean(axis=1)
        
        # Seleciona e ordena os top 12 mais ativos para caber perfeitamente na figura
        df_shared_sorted = df_shared.sort_values(by='Mean_FC', ascending=False).head(12)
        
        col_labels = ['Short ID', 'Gene Symbol', 'Detox Family', 'Mean Log2FC']
        cell_text = []
        for idx, row in df_shared_sorted.iterrows():
            id_curto = idx.split('|')[-1] if '|' in idx else idx
            id_curto = id_curto.replace(".t1", "")
            id_curto = id_curto[:18] if len(id_curto) > 18 else id_curto
            fam_curta = row['Detox_Family'].split(': ')[-1] if ': ' in row['Detox_Family'] else row['Detox_Family']
            cell_text.append([id_curto, row['Gene'], fam_curta, f"{row['Mean_FC']:.2f}"])
            
        tabela = ax_table.table(
            cellText=cell_text,
            colLabels=col_labels,
            loc='center',
            cellLoc='center',
            colColours=['#f2f4f4'] * 4
        )
        tabela.auto_set_font_size(False)
        tabela.set_fontsize(7.5)
        tabela.scale(1.05, 1.8)
        
        ax_table.set_title(titulo_tabela, fontsize=10, fontweight='bold', pad=15)
    else:
        ax_table.text(0.5, 0.5, "No active detox genes found\nwith Log2FC >= 1.0", 
                      ha='center', va='center', fontsize=11, color='gray')

    plt.tight_layout()
    output_path = "quant/venn_intersection_detox.png"
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"📊 [SUCESSO] Diagrama de Venn académico guardado em '{output_path}'!")
    
    # --- GERAÇÃO DE RELATÓRIO TSV DETALHADO DA INTERSEÇÃO ---
    gerar_relatorio_intersecao(abam_only, spin_only, emam_only, abam_spin, abam_emam, spin_emam, shared_all, df_detox, col_abam, col_spin, col_emam)

def gerar_relatorio_intersecao(abam_only, spin_only, emam_only, abam_spin, abam_emam, spin_emam, shared_all, df_detox, col_abam, col_spin, col_emam):
    """
    Grava um relatório em formato Markdown detalhado com a identidade biológica 
    de todos os transcritos em cada compartimento da interseção de Venn.
    """
    path_report = "quant/venn_detailed_report.md"
    print(f"📝 A gerar relatório analítico em '{path_report}'...")
    
    linhas = []
    linhas.append("# Detailed Intersection Analysis of Detoxification Genomes")
    linhas.append("")
    linhas.append("This document catalogs the identity of *Phthorimaea absoluta* transcripts inside each Venn compartment.")
    linhas.append("")
    
    # Secção 1: Core Shared (All 3)
    linhas.append("## 1. Core Shared Detoxification Genes (Active in All 3 Treatments)")
    if len(shared_all) > 0:
        linhas.append("| Transcript ID | Gene | Family | Abamectin FC | Spinosad FC | Emamectin FC |")
        linhas.append("|---|---|---|---|---|---|")
        for idx in shared_all:
            row = df_detox.loc[idx]
            gene = extrair_nome_gene(row['Annotation'])
            linhas.append(f"| `{idx}` | **{gene}** | {row['Detox_Family']} | {row[col_abam]:.2f} | {row[col_spin]:.2f} | {row[col_emam]:.2f} |")
    else:
        linhas.append("*No genes were strictly shared across all three treatments under the current threshold (Log2FC >= 1.0).*")
    linhas.append("")
    
    # Secção 2: Abamectin Exclusive
    linhas.append("## 2. Abamectin-Specific Detoxification Responders")
    if len(abam_only) > 0:
        linhas.append("| Transcript ID | Gene | Family | Abamectin FC |")
        linhas.append("|---|---|---|---|")
        for idx in abam_only:
            row = df_detox.loc[idx]
            gene = extrair_nome_gene(row['Annotation'])
            linhas.append(f"| `{idx}` | {gene} | {row['Detox_Family']} | {row[col_abam]:.2f} |")
    else:
        linhas.append("*No genes found exclusively active under Abamectin stress.*")
    linhas.append("")

    # Secção 3: Spinosad Exclusive
    linhas.append("## 3. Spinosad-Specific Detoxification Responders")
    if len(spin_only) > 0:
        linhas.append("| Transcript ID | Gene | Family | Spinosad FC |")
        linhas.append("|---|---|---|---|")
        for idx in spin_only:
            row = df_detox.loc[idx]
            gene = extrair_nome_gene(row['Annotation'])
            linhas.append(f"| `{idx}` | {gene} | {row['Detox_Family']} | {row[col_spin]:.2f} |")
    else:
        linhas.append("*No genes found exclusively active under Spinosad stress.*")
    linhas.append("")

    # Secção 4: Emamectin Exclusive
    linhas.append("## 4. Emamectin-Specific Detoxification Responders")
    if len(emam_only) > 0:
        linhas.append("| Transcript ID | Gene | Family | Emamectin FC |")
        linhas.append("|---|---|---|---|")
        for idx in emam_only:
            row = df_detox.loc[idx]
            gene = extrair_nome_gene(row['Annotation'])
            linhas.append(f"| `{idx}` | {gene} | {row['Detox_Family']} | {row[col_emam]:.2f} |")
    else:
        linhas.append("*No genes found exclusively active under Emamectin stress.*")
    linhas.append("")

    # Secção 5: Interseções Duplas (Pairwise)
    linhas.append("## 5. Shared Responders (Pairwise Intersections)")
    linhas.append("")
    
    linhas.append("### Abamectin & Spinosad")
    if len(abam_spin) > 0:
        linhas.append("| Transcript ID | Gene | Family | Abamectin FC | Spinosad FC |")
        linhas.append("|---|---|---|---|---|")
        for idx in abam_spin:
            row = df_detox.loc[idx]
            gene = extrair_nome_gene(row['Annotation'])
            linhas.append(f"| `{idx}` | {gene} | {row['Detox_Family']} | {row[col_abam]:.2f} | {row[col_spin]:.2f} |")
    else:
        linhas.append("*No shared genes found between Abamectin & Spinosad.*")
    linhas.append("")

    linhas.append("### Abamectin & Emamectin")
    if len(abam_emam) > 0:
        linhas.append("| Transcript ID | Gene | Family | Abamectin FC | Emamectin FC |")
        linhas.append("|---|---|---|---|---|")
        for idx in abam_emam:
            row = df_detox.loc[idx]
            gene = extrair_nome_gene(row['Annotation'])
            linhas.append(f"| `{idx}` | {gene} | {row['Detox_Family']} | {row[col_abam]:.2f} | {row[col_emam]:.2f} |")
    else:
        linhas.append("*No shared genes found between Abamectin & Emamectin.*")
    linhas.append("")

    linhas.append("### Spinosad & Emamectin")
    if len(spin_emam) > 0:
        linhas.append("| Transcript ID | Gene | Family | Spinosad FC | Emamectin FC |")
        linhas.append("|---|---|---|---|---|")
        for idx in spin_emam:
            row = df_detox.loc[idx]
            gene = extrair_nome_gene(row['Annotation'])
            linhas.append(f"| `{idx}` | {gene} | {row['Detox_Family']} | {row[col_spin]:.2f} | {row[col_emam]:.2f} |")
    else:
        linhas.append("*No shared genes found between Spinosad & Emamectin.*")
    linhas.append("")

    with open(path_report, 'w', encoding='utf-8') as f:
        f.write("\n".join(linhas))
    print("   ✅ Relatório detalhado guardado com sucesso!")

def executar_venn():
    sets = carregar_dados_e_definir_sets()
    if sets is not None:
        set_abam, set_spin, set_emam, df_detox, col_abam, col_spin, col_emam = sets
        plotar_diagrama_venn_custom(set_abam, set_spin, set_emam, df_detox, col_abam, col_spin, col_emam)

if __name__ == "__main__":
    executar_venn()
