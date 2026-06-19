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
