#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Analisador e Visualizador de Partição Global do Transcriptoma
Este script analisa a tabela de metanálise contendo todo o genoma para quantificar
quantos genes foram ativados (UP), reprimidos (DOWN) ou mantiveram-se neutros.
Gera uma figura de publicação de dois painéis (Donut + Bar Chart) em inglês.
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

def analisar_transcriptoma_global():
    print("==================================================")
    print("     SAQE - ANALISADOR GLOBAL DO TRANSCRIPTOMA     ")
    print("==================================================")
    
    bono_path = "quant/metanalise_completa_bono.csv"
    
    if not os.path.exists(bono_path):
        print(f"❌ [ERRO] O ficheiro bruto da metanálise '{bono_path}' não foi encontrado!")
        return
        
    # Lê todos os transcritos do genoma completo
    df = pd.read_csv(bono_path, index_col=0)
    total_transcritos = len(df)
    
    print(f"🧬 Genoma analisado com sucesso! Total de transcritos: {total_transcritos:,}")
    
    # --- 1. QUANTIFICAÇÃO DOS GRUPOS GERAIS ---
    # UP: Ativo em pelo menos 1 inseticida (Bono_UP_Score >= 1)
    # DOWN: Ativo em pelo menos 1 inseticida (Bono_DOWN_Score >= 1)
    # Neutro: Não cruzou os limiares em nenhum dos 3 estudos
    
    up_total = df[df['Bono_UP_Score'] >= 1]
    down_total = df[df['Bono_DOWN_Score'] >= 1]
    
    # Transcritos que são estritamente neutros (Score = 0 para UP e DOWN)
    neutros = df[(df['Bono_UP_Score'] == 0) & (df['Bono_DOWN_Score'] == 0)]
    
    # Transcritos ambíguos (raros casos onde o gene é UP em um estudo e DOWN em outro)
    ambiguos = df[(df['Bono_UP_Score'] >= 1) & (df['Bono_DOWN_Score'] >= 1)]
    
    count_up = len(up_total) - len(ambiguos)
    count_down = len(down_total) - len(ambiguos)
    count_neutral = len(neutros)
    count_ambiguous = len(ambiguos)
    
    print(f"\n📊 PARTIÇÃO GLOBAL DO GENOMA:")
    print(f"   🔹 Up-regulated (Induzidos): {count_up:,} ({count_up/total_transcritos*100:.2f}%)")
    print(f"   🔹 Down-regulated (Reprimidos): {count_down:,} ({count_down/total_transcritos*100:.2f}%)")
    print(f"   🔹 Neutral (Sem alteração): {count_neutral:,} ({count_neutral/total_transcritos*100:.2f}%)")
    if count_ambiguous > 0:
        print(f"   🔹 Dual/Ambiguous (UP e DOWN): {count_ambiguous:,} ({count_ambiguous/total_transcritos*100:.2f}%)")

    # --- 2. DETALHAMENTO DA CONSISTÊNCIA DE REGULAÇÃO (SCORES 1, 2, 3) ---
    up_s1 = len(df[df['Bono_UP_Score'] == 1])
    up_s2 = len(df[df['Bono_UP_Score'] == 2])
    up_s3 = len(df[df['Bono_UP_Score'] == 3])
    
    down_s1 = len(df[df['Bono_DOWN_Score'] == 1])
    down_s2 = len(df[df['Bono_DOWN_Score'] == 2])
    down_s3 = len(df[df['Bono_DOWN_Score'] == 3])
    
    print(f"\n📋 CONSISTÊNCIA DE ATIVAÇÃO (ON-SCORE):")
    print(f"   🔹 UP-Score 1 (Específico): {up_s1:,} transcritos")
    print(f"   🔹 UP-Score 2 (Consistente): {up_s2:,} transcritos")
    print(f"   🔹 UP-Score 3 (Core Shared): {up_s3:,} transcritos")
    print(f"   🔹 DOWN-Score 1 (Específico): {down_s1:,} transcritos")
    print(f"   🔹 DOWN-Score 2 (Consistente): {down_s2:,} transcritos")
    print(f"   🔹 DOWN-Score 3 (Core Repressed): {down_s3:,} transcritos")

    # --- 3. DESENHO DO GRÁFICO GLOBAL (DONUT + BAR CHART) ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.5), dpi=300)
    sns.set_theme(style="ticks")
    
    # Painel Esquerdo: Donut Chart da Partição Global do Genoma
    labels_global = ['Induced (UP)\n' + f'{count_up:,}', 
                     'Repressed (DOWN)\n' + f'{count_down:,}', 
                     'Neutral\n' + f'{count_neutral:,}']
    sizes_global = [count_up, count_down, count_neutral]
    colors_global = ['#e74c3c', '#3498db', '#ecf0f1']
    
    # Desenha a pizza
    wedges, texts, autotexts = ax1.pie(
        sizes_global, 
        labels=labels_global, 
        autopct='%1.2f%%',
        startangle=140, 
        colors=colors_global,
        textprops=dict(color="black", fontsize=8.5, fontweight='bold'),
        wedgeprops=dict(width=0.4, edgecolor='black', linewidth=0.5) # Cria o furo central (donut)
    )
    
    # Ajusta os textos de percentual dentro do donut
    for autotext in autotexts:
        autotext.set_fontsize(8)
        autotext.set_color('black')
        
    ax1.set_title("Global Transcriptome Partition", fontsize=11, fontweight='bold', pad=15)
    
    # Painel Direito: Barras Agrupadas dos Níveis de Consistência (Scores 1, 2, 3)
    categories = ['Score 1\n(Specific)', 'Score 2\n(Consistent)', 'Score 3\n(Core Shared)']
    up_scores = [up_s1, up_s2, up_s3]
    down_scores = [down_s1, down_s2, down_s3]
    
    x = np.arange(len(categories))
    width = 0.35
    
    rects1 = ax2.bar(x - width/2, up_scores, width, label='Induced (UP)', color='#e74c3c', edgecolor='black', linewidth=0.5)
    rects2 = ax2.bar(x + width/2, down_scores, width, label='Repressed (DOWN)', color='#3498db', edgecolor='black', linewidth=0.5)
    
    ax2.set_title("Consistency Score Distribution", fontsize=11, fontweight='bold', pad=15)
    ax2.set_ylabel("Transcript Count", fontsize=9)
    ax2.set_xticks(x)
    ax2.set_xticklabels(categories, fontsize=8.5, fontweight='bold')
    ax2.legend(fontsize=9, loc='upper right')
    
    # Força a interseção exata no zero absoluto e adiciona ticks secundários
    ax2.set_ylim(bottom=0)
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.tick_params(axis='both', which='both', bottom=True, left=True, labelsize=8.5)
    
    # Adiciona os valores numéricos exatos no topo das barras
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            if height > 0:
                ax2.annotate(f'{height:,}',
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),  # Deslocamento vertical de 3 pontos
                            textcoords="offset points",
                            ha='center', va='bottom', fontsize=7.5, fontweight='bold')
                            
    autolabel(rects1)
    autolabel(rects2)
    
    # Título geral sóbrio da figura com nome científico em itálico
    plt.suptitle(r"Global Landscape of Differential Transcript Regulation in $\mathit{Phthorimaea\ absoluta}$", 
                 fontsize=12, fontweight='bold', y=0.98)
    
    # Nota de rodapé técnica
    legend_note = (
        "UP (Induction) filter: log2 Fold Change >= 1.0 and adjusted p < 0.05. "
        "DOWN (Repression) filter: log2 Fold Change <= -1.0 and adjusted p < 0.05.\n"
        "Consistency Score represents the number of insecticide treatments (Abamectin, Spinosad, Emamectin) meeting the significance threshold."
    )
    fig.text(0.5, 0.02, legend_note, ha='center', va='bottom', fontsize=7.5, style='italic',
             bbox=dict(boxstyle='round,pad=0.5', fc='#f8f9fa', ec='#e2e8f0', lw=0.8))
    
    sns.despine(trim=True)
    plt.tight_layout(rect=[0, 0.08, 1, 0.94])
    
    output_png = "quant/global_transcriptome_partition.png"
    plt.savefig(output_png, bbox_inches='tight')
    plt.close()
    
    print(f"\n📊 [SUCESSO] Gráfico de partição global gravado em '{output_png}'!")

if __name__ == "__main__":
    analisar_transcriptoma_global()
