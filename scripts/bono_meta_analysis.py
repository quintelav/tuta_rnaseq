#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
📖 TUTORIAL DE USO DO NANO (Como Criar, Editar e Executar este Script)
==============================================================================
1. Como criar ou abrir este arquivo no terminal usando o editor de texto nano:
   $ nano bono_meta_analysis.py

2. Como colar o código copiado do Canvas no terminal (WSL ou Linux Nativo):
   - Use o atalho de teclado: Ctrl + Shift + V
   - Ou simplesmente clique com o botão direito do mouse/rato sobre a janela.

3. Como Gravar/Salvar as alterações dentro do editor nano:
   - Pressione o atalho: Ctrl + O
   - Confirme o nome do ficheiro pressionando: Enter

4. Como Sair/Fechar o editor de texto nano de volta ao terminal Bash:
   - Pressione o atalho: Ctrl + X

5. Como dar permissão de execução para o script no terminal:
   - Rode o comando:
     $ chmod +x bono_meta_analysis.py

6. Como executar a Metanálise de Bono (2021) com todos os dados integrados:
   - Rode o comando:
     $ python3 bono_meta_analysis.py
==============================================================================

SAQE - Implementação da Metanálise de Transcriptomas de Insetos (Bono, 2021)
Calcula as métricas de ON-ratio (escala log10 com pseudocount de 0.01) e 
ON-score (frequência líquida de regulação: contagem_UP - contagem_DOWN).
Gera uma figura de publicação (300 DPI) com a distribuição do ON-score (Figura 2B).
Realiza varredura física completa e mapeia dinamicamente os 26 SRAs dos 3 inseticidas.
"""

import os
import glob
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def escanear_e_ajustar_todos_metadados(quant_dir, metadata_path):
    """
    Varre a pasta 'quant' procurando por todas as pastas de quantificação físicas (SRR*).
    Mapeia os SRRs detectados para os seus respectivos tratamentos e condições baseando-se
    no mapa mestre de 26 amostras (6 Spinosad, 8 Emamectina, 12 Abamectina).
    Atualiza o arquivo de metadados de forma robusta e dinâmica.
    """
    print("\n🔍 A iniciar varredura e mapeamento físico de todas as amostras SRA no disco...")
    
    # Mapa mestre de todas as 26 amostras possíveis dos 3 estudos
    master_sra_map = {
        # Spinosad (6 amostras)
        'SRR33779431': {'condition': 'control', 'study_id': 'spinosad'},
        'SRR33779430': {'condition': 'control', 'study_id': 'spinosad'},
        'SRR33779429': {'condition': 'control', 'study_id': 'spinosad'},
        'SRR33779428': {'condition': 'treatment', 'study_id': 'spinosad'},
        'SRR33779427': {'condition': 'treatment', 'study_id': 'spinosad'},
        'SRR33779426': {'condition': 'treatment', 'study_id': 'spinosad'},
        
        # Emamectina (8 amostras fisicamente presentes)
        'SRR15248447': {'condition': 'control', 'study_id': 'emamectina'},
        'SRR15248446': {'condition': 'control', 'study_id': 'emamectina'},
        'SRR15248441': {'condition': 'control', 'study_id': 'emamectina'},
        'SRR15248440': {'condition': 'control', 'study_id': 'emamectina'},
        'SRR15248443': {'condition': 'treatment', 'study_id': 'emamectina'},
        'SRR15248442': {'condition': 'treatment', 'study_id': 'emamectina'},
        'SRR15248439': {'condition': 'treatment', 'study_id': 'emamectina'},
        'SRR15248438': {'condition': 'treatment', 'study_id': 'emamectina'},
        
        # Abamectina (12 amostras)
        'SRR35985991': {'condition': 'control', 'study_id': 'abamectina'},
        'SRR35985990': {'condition': 'control', 'study_id': 'abamectina'},
        'SRR35985989': {'condition': 'control', 'study_id': 'abamectina'},
        'SRR35985995': {'condition': 'control', 'study_id': 'abamectina'},
        'SRR35985988': {'condition': 'treatment', 'study_id': 'abamectina'},
        'SRR35985987': {'condition': 'treatment', 'study_id': 'abamectina'},
        'SRR35985986': {'condition': 'treatment', 'study_id': 'abamectina'},
        'SRR35985994': {'condition': 'treatment', 'study_id': 'abamectina'},
        'SRR35985985': {'condition': 'treatment', 'study_id': 'abamectina'},
        'SRR35985984': {'condition': 'treatment', 'study_id': 'abamectina'},
        'SRR35985993': {'condition': 'treatment', 'study_id': 'abamectina'},
        'SRR35985992': {'condition': 'treatment', 'study_id': 'abamectina'}
    }
    
    # Varre a pasta de quantificação para achar pastas que iniciam com SRR
    padrao = os.path.join(quant_dir, "SRR*")
    pastas_reais = [os.path.basename(p) for p in glob.glob(padrao) if os.path.isdir(p)]
    
    print(f"   📂 Total de pastas físicas SRA encontradas no disco: {len(pastas_reais)}")
    
    # Filtra e mapeia apenas as amostras fisicamente presentes que estão no nosso mapa mestre
    amostras_mapeadas = []
    for sra in sorted(pastas_reais):
        if sra in master_sra_map:
            info = master_sra_map[sra]
            amostras_mapeadas.append({
                'sample_id': sra,
                'condition': info['condition'],
                'study_id': info['study_id']
            })
    
    # Estatísticas do auto-descobrimento e gravação
    df_novas = pd.DataFrame(amostras_mapeadas)
    if not df_novas.empty:
        print("\n📊 Resumo das Amostras Detectadas e Mapeadas para o Pipeline:")
        contagem_estudos = df_novas['study_id'].value_counts()
        for estudo in sorted(contagem_estudos.index):
            total_sra = contagem_estudos[estudo]
            ctrls = len(df_novas[(df_novas['study_id'] == estudo) & (df_novas['condition'] == 'control')])
            trats = len(df_novas[(df_novas['study_id'] == estudo) & (df_novas['condition'] == 'treatment')])
            print(f"   🔹 Inseticida: {estudo:<12} | {total_sra:>2} amostras totais ({ctrls} Controles vs {trats} Tratados)")
        
        # Sobrescreve o arquivo de metadados para garantir total sincronia com o disco físico
        df_novas.to_csv(metadata_path, index=False)
        print(f"\n   💾 Tabela de metadados '{metadata_path}' atualizada e sincronizada com sucesso!")
    else:
        print("   ⚠️  Aviso: Nenhuma pasta de amostra válida correspondente ao mapa mestre foi encontrada no disco.")

def carregar_dados_salmon_fallback(quant_dir):
    """
    Função de segurança: Caso o arquivo matriz_tpm_geral.csv não exista,
    varre a pasta 'quant' e reconstrói a matriz de TPM real consolidada.
    """
    print("📂 Procurando pastas de quantificação do Salmon para reconstruir a matriz de TPM real...")
    dados_tpm = {}
    padrao_busca = os.path.join(quant_dir, "*", "quant.sf")
    arquivos = glob.glob(padrao_busca)
    
    if not arquivos:
        return None
        
    print(f"   🔹 Encontrados {len(arquivos)} ficheiros quant.sf. Unificando dados...")
    for caminho in arquivos:
        sra_id = os.path.basename(os.path.dirname(caminho))
        try:
            df = pd.read_csv(caminho, sep='\t')
            if 'Name' in df.columns and 'TPM' in df.columns:
                dados_tpm[sra_id] = df.set_index('Name')['TPM']
        except Exception as e:
            print(f"   ⚠️  Erro ao ler o ficheiro {caminho}: {e}")
            
    if not dados_tpm:
        return None
        
    df_expressao = pd.DataFrame(dados_tpm)
    df_expressao.to_csv(f"{quant_dir}/matriz_tpm_geral.csv")
    print(f"   ✅ Matriz de TPM geral unificada guardada em '{quant_dir}/matriz_tpm_geral.csv'.")
    return df_expressao

def gerar_dados_simulados():
    """
    Gera dados sintéticos apenas se não houver absolutamente nenhum dado real
    no ambiente, evitando que o script falhe.
    """
    print("✨ Nenhum dado real de RNA-Seq detectado. Gerando dados de simulação para teste...")
    np.random.seed(42)
    genes = [f"Gene_{i:03d}" for i in range(1, 501)]
    
    amostras = []
    dados = {}
    for est in [1, 2, 3]:
        for cond in ['Ctrl', 'Str']:
            for rep in [1, 2, 3]:
                nome_amostra = f"E{est}_{cond}_R{rep}"
                amostras.append(nome_amostra)
                if cond == 'Str':
                    tpm = np.random.exponential(scale=10, size=500)
                    tpm[:50] *= np.random.uniform(5, 50, size=50)
                    tpm[50:100] /= np.random.uniform(5, 50, size=50)
                else:
                    tpm = np.random.exponential(scale=10, size=500)
                
                zero_mask = np.random.rand(500) < 0.1
                tpm[zero_mask] = 0.0
                dados[nome_amostra] = tpm
                
    df_tpm = pd.DataFrame(dados, index=genes)
    df_tpm.index.name = 'Gene_ID'
    
    os.makedirs("quant", exist_ok=True)
    df_tpm.to_csv("quant/matriz_tpm_geral.csv")
    
    meta_rows = []
    for est in [1, 2, 3]:
        for cond in ['control', 'treatment']:
            for rep in [1, 2, 3]:
                cond_prefix = 'Ctrl' if cond == 'control' else 'Str'
                meta_rows.append({
                    'sample_id': f"E{est}_{cond_prefix}_R{rep}",
                    'condition': cond,
                    'study_id': f"Estudo_{est}"
                })
    df_meta = pd.DataFrame(meta_rows)
    df_meta.to_csv("metadata_geral.csv", index=False)
    print("✅ Ficheiros de teste 'quant/matriz_tpm_geral.csv' e 'metadata_geral.csv' criados!")

def executar_metanalise_bono(tpm_path, metadata_path, output_path="quant/resultado_metanalise_bono.csv"):
    """
    Executa a metanálise quantitativa de Bono (2021) com dados reais:
    1. Calcula o log10 fold change (ON-ratio) com pseudocount de 0.01 por estudo.
    2. Aplica os limiares de corte (UP se FC > 10, DOWN se FC < 0.1).
    3. Consolida as métricas: média do ON-ratio e o ON-score líquido (UP - DOWN).
    """
    print("\n🚀 A iniciar processamento da Metanálise de Bono (2021)...")
    
    # Carrega dados
    df_tpm = pd.read_csv(tpm_path, index_col=0)
    df_meta = pd.read_csv(metadata_path)
    
    # Validações básicas das colunas exigidas nos metadados
    for col in ['sample_id', 'condition', 'study_id']:
        if col not in df_meta.columns:
            raise KeyError(f"A coluna '{col}' não foi encontrada no ficheiro de metadados!")
            
    estudos = df_meta['study_id'].unique()
    print(f"📊 Detetados {len(estudos)} estudos independentes para integração: {list(estudos)}")
    
    # Inicializa DataFrame de resultados
    df_resultados = pd.DataFrame(index=df_tpm.index)
    colunas_ratios = []
    
    total_up_matrix = np.zeros(len(df_tpm))
    total_down_matrix = np.zeros(len(df_tpm))
    
    # Processa cada estudo individualmente
    for estudo in estudos:
        meta_sub = df_meta[df_meta['study_id'] == estudo]
        
        amostras_ctrl = meta_sub[meta_sub['condition'].str.lower() == 'control']['sample_id'].tolist()
        amostras_stress = meta_sub[meta_sub['condition'].str.lower().isin(['treatment', 'stress'])]['sample_id'].tolist()
        
        # Filtra amostras que realmente existem na matriz de TPM
        amostras_ctrl = [x for x in amostras_ctrl if x in df_tpm.columns]
        amostras_stress = [x for x in amostras_stress if x in df_tpm.columns]
        
        if len(amostras_ctrl) == 0 or len(amostras_stress) == 0:
            print(f"⚠️  Estudo '{estudo}' ignorado devido a amostras em falta na matriz de TPM.")
            continue
            
        print(f"   🔹 A processar '{estudo}': {len(amostras_ctrl)} controles vs {len(amostras_stress)} estresses")
        
        # Média das TPMs para controle e estresse dentro deste estudo
        mean_ctrl = df_tpm[amostras_ctrl].mean(axis=1)
        mean_stress = df_tpm[amostras_stress].mean(axis=1)
        
        # FÓRMULA DE BONO (2021) - ON-ratio individual com pseudocount de 0.01 em base log10
        on_ratio_estudo = np.log10(mean_stress + 0.01) - np.log10(mean_ctrl + 0.01)
        
        col_name = f"ON_ratio_{estudo}"
        df_resultados[col_name] = on_ratio_estudo
        colunas_ratios.append(col_name)
        
        up_mask = on_ratio_estudo > 1.0
        down_mask = on_ratio_estudo < -1.0
        
        total_up_matrix += up_mask.astype(int)
        total_down_matrix += down_mask.astype(int)
        
    # --- CONSOLIDAÇÃO DA METANÁLISE ---
    df_resultados['ON_score'] = (total_up_matrix - total_down_matrix).astype(int)
    df_resultados['mean_ON_ratio'] = df_resultados[colunas_ratios].mean(axis=1)
    df_resultados['count_UP'] = total_up_matrix.astype(int)
    df_resultados['count_DOWN'] = total_down_matrix.astype(int)
    
    # Salva a tabela final
    colunas_finais = ['mean_ON_ratio', 'ON_score', 'count_UP', 'count_DOWN'] + colunas_ratios
    df_saida = df_resultados[colunas_finais]
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_saida.to_csv(output_path)
    
    print(f"\n💾 Tabela final da metanálise real guardada com sucesso em: '{output_path}'!")
    return df_saida

def plotar_distribuicao_onscore(df_resultados, N_estudos, fig_path="quant/onscore_distribution.png"):
    """
    Gera o histograma de distribuição do ON-score (Semelhante à Figura 2B de Bono, 2021).
    Utiliza eixos tangentes limpos e sóbrios.
    """
    print("📊 A gerar figura científica da distribuição do ON-score...")
    
    plt.figure(figsize=(7.5, 5), dpi=300)
    
    score_range = np.arange(-N_estudos, N_estudos + 1)
    contagens = [len(df_resultados[df_resultados['ON_score'] == s]) for s in score_range]
    
    cores = []
    for s in score_range:
        if s < 0:
            cores.append('#3498db')
        elif s > 0:
            cores.append('#e74c3c')
        else:
            cores.append('#bdc3c7')
            
    bars = plt.bar(score_range, contagens, color=cores, edgecolor='black', linewidth=0.6, width=0.6)
    
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
    plt.xlim(-N_estudos - 0.6, N_estudos + 0.6)
    plt.ylim(bottom=0, top=max(contagens) * 1.15)
    
    plt.xticks(score_range, [f"+{s}" if s > 0 else str(s) for s in score_range], fontsize=9, fontweight='bold')
    plt.yticks(fontsize=9)
    
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            plt.annotate(f'{height:,}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),  
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8, fontweight='bold')
            
    plt.title("Distribution of Net ON-score Values across Insecticide Stress Trials", 
              fontsize=10, fontweight='bold', pad=15)
    plt.xlabel("Net ON-score Value (count UP - count DOWN)", fontsize=9, fontweight='bold')
    plt.ylabel("Number of Transcripts", fontsize=9, fontweight='bold')
    
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(fig_path, bbox_inches='tight')
    plt.close()
    print(f"🖼️  Gráfico de publicação guardado em: '{fig_path}'!")

if __name__ == "__main__":
    tpm_file = "quant/matriz_tpm_geral.csv"
    metadata_file = "metadata_geral.csv"
    
    # 1. Garante que as pastas básicas existam
    os.makedirs("quant", exist_ok=True)
    
    # 2. Executa a varredura inteligente de TODOS os 3 blocos no disco físico ANTES de qualquer processamento
    escanear_e_ajustar_todos_metadados("quant", metadata_file)
    
    # 3. Força a reconstrução da matriz de TPM caso o metadata_geral.csv tenha sido atualizado ou não exista
    if os.path.exists(tpm_file):
        os.remove(tpm_file)
        
    if not os.path.exists(tpm_file):
        carregar_dados_salmon_fallback("quant")
        
    # 4. Se os arquivos continuarem ausentes, executa a simulação automática para segurança
    if not os.path.exists(tpm_file) or not os.path.exists(metadata_file):
        gerar_dados_simulados()
        
    # 5. Executa a metanálise com os dados reais/ajustados
    df_meta_res = executar_metanalise_bono(tpm_file, metadata_file)
    
    # 6. Obtém o número real de estudos dinamicamente para calibrar o gráfico
    df_meta = pd.read_csv(metadata_file)
    N_estudos = len(df_meta['study_id'].unique())
    
    # 7. Gera e grava o gráfico de publicação
    plotar_distribuicao_onscore(df_meta_res, N_estudos)
    
    print("\n====================================================================")
    print("🎉 PIPELINE DE METANÁLISE DE BONO EXECUTADO COM SUCESSO!")
    print("====================================================================")
    print("   👉 Resultados salvos em: 'quant/resultado_metanalise_bono.csv'")
    print("   👉 Figura gerada em:     'quant/onscore_distribution.png'")
    print("====================================================================")
