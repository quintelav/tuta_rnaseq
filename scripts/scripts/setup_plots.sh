#!/bin/bash
# ==============================================================================
# SAQE - Script de Configuração Automatizada das Figuras Globais
# Cria de forma robusta os arquivos independentes de visualização do transcriptoma.
# ==============================================================================

echo "===================================================================="
echo "⚙️  SAQE - Configurando Ficheiros de Visualização Separados..."
echo "===================================================================="

# --- 1. CRIAÇÃO DO SCRIPT: gerar_donut_global.py ---
echo "📝 Gravando 'gerar_donut_global.py'..."
cat << 'EOF' > gerar_donut_global.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Analisador e Visualizador de Partição Global do Transcriptoma (Publication Ready)
Focado exclusivamente na geração do gráfico de Donut da partição do transcriptoma
com base no On-Score original de Bono (2021): Net On-Score = count(UP) - count(DOWN).
Gera uma figura de alta resolução (300 DPI) em inglês, sem sobreposições.
"""

import os
import pandas as pd
import numpy as np

# Configura o Matplotlib para rodar em modo 'headless' (sem interface gráfica/X11)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def gerar_grafico_donut():
    print("==================================================")
    print("     SAQE - GERADOR DO GRÁFICO DE DONUT GLOBAL     ")
    print("==================================================")
    
    bono_path = "quant/metanalise_completa_bono.csv"
    
    if not os.path.exists(bono_path):
        print(f"❌ [ERRO] O ficheiro bruto da metanálise '{bono_path}' não foi encontrado!")
        return
        
    # Lê todos os transcritos do genoma completo
    df = pd.read_csv(bono_path, index_col=0)
    total_transcritos = len(df)
    
    print(f"🧬 Genoma analisado com sucesso! Total de transcritos: {total_transcritos:,}")
    
    # --- 1. METODOLOGIA ESTREITA DE BONO (2021) ---
    if 'Bono_Net_Score' not in df.columns:
        df['Bono_Net_Score'] = df['Bono_UP_Score'] - df['Bono_DOWN_Score']
        
    count_up = len(df[df['Bono_Net_Score'] > 0])
    count_down = len(df[df['Bono_Net_Score'] < 0])
    count_neutral = len(df[df['Bono_Net_Score'] == 0])
    
    pct_up = (count_up / total_transcritos) * 100
    pct_down = (count_down / total_transcritos) * 100
    pct_neutral = (count_neutral / total_transcritos) * 100

    print(f"\n📊 PARTIÇÃO GLOBAL DO GENOMA (NET ON-SCORE):")
    print(f"   🔹 Induced (On-Score > 0): {count_up:,} ({pct_up:.2f}%)")
    print(f"   🔹 Repressed (On-Score < 0): {count_down:,} ({pct_down:.2f}%)")
    print(f"   🔹 Neutral (On-Score == 0): {count_neutral:,} ({pct_neutral:.2f}%)")

    # --- 2. DESENHO DO GRÁFICO CIRCULAR (DONUT CHART) ---
    fig, ax = plt.subplots(figsize=(8, 7.5), dpi=300)
    sns.set_theme(style="ticks")
    
    sizes_global = [count_up, count_down, count_neutral]
    colors_global = ['#e74c3c', '#3498db', '#ecf0f1']
    
    # Desenho do donut limpo (labels=None evita amontoados de texto ao redor do círculo)
    wedges, texts, autotexts = ax.pie(
        sizes_global, 
        labels=None, 
        autopct=lambda pct: f'{pct:.2f}%' if pct > 0.5 else '',
        startangle=140, 
        colors=colors_global,
        textprops=dict(color="black", fontsize=9, fontweight='bold'),
        wedgeprops=dict(width=0.4, edgecolor='black', linewidth=0.5)
    )
    
    for autotext in autotexts:
        autotext.set_fontsize(9)
        
    ax.set_title("Global Transcriptome Partition", fontsize=11, fontweight='bold', pad=15)
    
    # Legenda estruturada na base da figura para evitar qualquer sobreposição
    labels_legenda = [
        f"Induced (On-Score > 0): {count_up:,} ({pct_up:.2f}%)",
        f"Repressed (On-Score < 0): {count_down:,} ({pct_down:.2f}%)",
        f"Neutral (On-Score = 0): {count_neutral:,} ({pct_neutral:.2f}%)"
    ]
    ax.legend(
        wedges, 
        labels_legenda,
        title="Transcriptome Cohorts",
        title_fontsize=9.5,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=1,
        fontsize=8.5,
        frameon=True,
        facecolor='#f8f9fa',
        edgecolor='#e2e8f0'
    )
    
    # --- NOTA DE RODAPÉ METODOLÓGICA SÓBRIA ---
    legend_note = (
        "Net On-Score = count(upregulated studies) - count(downregulated studies).\n"
        "UP (Induction) threshold: log2 Fold Change >= 1.0 and adjusted p < 0.05.\n"
        "DOWN (Repression) threshold: log2 Fold Change <= -1.0 and adjusted p < 0.05."
    )
    fig.text(0.5, 0.02, legend_note, ha='center', va='bottom', fontsize=7.5, style='italic',
             bbox=dict(boxstyle='round,pad=0.5', fc='#f8f9fa', ec='#e2e8f0', lw=0.8))
    
    # Reserva de espaço para a legenda e nota de rodapé
    plt.tight_layout(rect=[0, 0.18, 1, 0.98])
    
    output_png = "quant/global_transcriptome_partition_donut.png"
    plt.savefig(output_png, bbox_inches='tight')
    plt.close()
    
    print(f"\n📊 [SUCESSO] Gráfico Donut guardado em '{output_png}'!")

if __name__ == "__main__":
    gerar_grafico_donut()
EOF

# --- 2. CRIAÇÃO DO SCRIPT: gerar_distribuicao_net_onscore.py ---
echo "📝 Gravando 'gerar_distribuicao_net_onscore.py'..."
cat << 'EOF' > gerar_distribuicao_net_onscore.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Analisador e Visualizador de Frequência do Net On-Score (Publication Ready)
Focado exclusivamente na geração do gráfico de barras do Net On-Score de Bono (2021).
Gera uma figura de publicação de alta resolução (300 DPI) em inglês, com eixos tangentes.
"""

import os
import pandas as pd
import numpy as np

# Configura o Matplotlib para rodar em modo 'headless' (sem interface gráfica/X11)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator

def gerar_grafico_barras_onscore():
    print("==================================================")
    print("   SAQE - GERADOR DO GRÁFICO DE BARRAS ON-SCORE   ")
    print("==================================================")
    
    bono_path = "quant/metanalise_completa_bono.csv"
    
    if not os.path.exists(bono_path):
        print(f"❌ [ERRO] O ficheiro bruto da metanálise '{bono_path}' não foi encontrado!")
        return
        
    # Lê todos os transcritos do genoma completo
    df = pd.read_csv(bono_path, index_col=0)
    
    if 'Bono_Net_Score' not in df.columns:
        df['Bono_Net_Score'] = df['Bono_UP_Score'] - df['Bono_DOWN_Score']
        
    # Contabilização de cada nível de consistência líquida do On-Score (Exclui os neutros 0)
    score_neg3 = len(df[df['Bono_Net_Score'] == -3])
    score_neg2 = len(df[df['Bono_Net_Score'] == -2])
    score_neg1 = len(df[df['Bono_Net_Score'] == -1])
    score_pos1 = len(df[df['Bono_Net_Score'] == 1])
    score_pos2 = len(df[df['Bono_Net_Score'] == 2])
    score_pos3 = len(df[df['Bono_Net_Score'] == 3])
    
    print(f"\n📋 DISTRIBUIÇÃO DETALHADA DO NET ON-SCORE:")
    print(f"   🔹 Score -3 (Universal Down): {score_neg3:,} transcritos")
    print(f"   🔹 Score -2 (Consistent Down): {score_neg2:,} transcritos")
    print(f"   🔹 Score -1 (Specific Down): {score_neg1:,} transcritos")
    print(f"   🔹 Score +1 (Specific Up): {score_pos1:,} transcritos")
    print(f"   🔹 Score +2 (Consistent Up): {score_pos2:,} transcritos")
    print(f"   🔹 Score +3 (Universal Up): {score_pos3:,} transcritos")

    # --- DESENHO DO GRÁFICO (BAR CHART) ---
    fig, ax = plt.subplots(figsize=(8.5, 6.5), dpi=300)
    sns.set_theme(style="ticks")
    
    categories = ['-3', '-2', '-1', '+1', '+2', '+3']
    scores_values = [score_neg3, score_neg2, score_neg1, score_pos1, score_pos2, score_pos3]
    
    # Paleta simétrica clássica: tons de azul para repressão (-) e tons de laranja/vermelho para indução (+)
    colors_bars = ['#2e86c1', '#5dade2', '#85c1e9', '#f5b041', '#eb984e', '#cb4335']
    
    rects = ax.bar(categories, scores_values, color=colors_bars, edgecolor='black', linewidth=0.5)
    
    ax.set_title("On-Score Distribution (Active Transcripts)", fontsize=11, fontweight='bold', pad=12)
    ax.set_ylabel("Transcript Count", fontsize=9.5)
    ax.set_xlabel("Net On-Score Value", fontsize=9.5)
    ax.tick_params(axis='both', which='both', bottom=True, left=True, labelsize=9)
    
    # FORÇA INTERSEÇÃO PERFEITA NO ZERO ABSOLUTO (Remove recuos automáticos das bordas)
    ax.set_xlim(left=-0.5, right=len(categories)-0.5)
    ax.set_ylim(bottom=0, top=max(scores_values) * 1.15)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Função auxiliar de anotação de valores sobre as barras
    def autolabel(rect_set):
        for rect in rect_set:
            height = rect.get_height()
            if height > 0:
                ax.annotate(f'{height:,}',
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),
                            textcoords="offset points",
                            ha='center', va='bottom', fontsize=8, fontweight='bold')
                            
    autolabel(rects)
    
    # --- NOTA DE RODAPÉ METODOLÓGICA SÓBRIA ---
    legend_note = (
        "Net On-Score = count(upregulated studies) - count(downregulated studies).\n"
        "UP (Induction) threshold: log2 Fold Change >= 1.0 and adjusted p < 0.05.\n"
        "DOWN (Repression) threshold: log2 Fold Change <= -1.0 and adjusted p < 0.05."
    )
    fig.text(0.5, 0.02, legend_note, ha='center', va='bottom', fontsize=7.5, style='italic',
             bbox=dict(boxstyle='round,pad=0.5', fc='#f8f9fa', ec='#e2e8f0', lw=0.8))
    
    sns.despine(trim=True)
    # Reserva de espaço para a nota explicativa
    plt.tight_layout(rect=[0, 0.12, 1, 0.98])
    
    output_png = "quant/global_onscore_distribution_bar.png"
    plt.savefig(output_png, bbox_inches='tight')
    plt.close()
    
    print(f"\n📊 [SUCESSO] Gráfico de barras guardado em '{output_png}'!")

if __name__ == "__main__":
    gerar_grafico_barras_onscore()
