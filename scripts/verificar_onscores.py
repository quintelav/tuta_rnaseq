#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
📖 TUTORIAL DE USO DO NANO (Como Criar, Editar e Executar este Script)
==============================================================================
1. Como criar ou abrir este arquivo no terminal usando o editor de texto nano:
   $ nano verificar_onscores.py

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
     $ chmod +x verificar_onscores.py

6. Como executar a ferramenta interativa de busca e verificação de On-Scores:
   - Rode o comando:
     $ python3 verificar_onscores.py   (ou simplesmente  ./verificar_onscores.py)
==============================================================================

SAQE - Navegador e Inspecionador Interativo de On-Scores de Bono (2021)
Este utilitário permite pesquisar, classificar e auditar a tabela de expressão 
geral do transcriptoma diretamente pela linha de comando de forma sóbria e rápida.
"""

import os
import sys
import pandas as pd

def limpar_terminal():
    """
    Limpa a tela do terminal de forma compatível com Windows (cmd) e Linux/WSL.
    """
    os.system('cls' if os.name == 'nt' else 'clear')

def obter_dados():
    """
    Localiza os arquivos de resultados de metanálise e calcula de forma dinâmica
    o Net On-Score de Bono (2021): Net On-Score = count(UP) - count(DOWN).
    """
    caminhos = [
        "quant/metanalise_completa_anotada.csv",
        "quant/metanalise_completa_bono.csv",
        "metanalise_completa_anotada.csv",
        "metanalise_completa_bono.csv"
    ]
    
    df = None
    caminho_final = ""
    for c in caminhos:
        if os.path.exists(c):
            df = pd.read_csv(c, index_col=0)
            caminho_final = c
            break
            
    if df is None:
        print("❌ [ERRO] Nenhum arquivo de resultados de metanálise encontrado nos diretórios padrão!")
        print("Caminhos verificados:")
        for c in caminhos:
            print(f"   • {c}")
        print("\nDica: Certifique-se de executar o pipeline diferencial primeiro:")
        print("   $ python3 analise_meta_bono.py")
        sys.exit(1)
        
    # Calcula dinamicamente o Net On-Score real alinhado com Bono (2021)
    if 'Bono_Net_Score' not in df.columns:
        df['Bono_Net_Score'] = df['Bono_UP_Score'] - df['Bono_DOWN_Score']
        
    return df, caminho_final

def menu_principal():
    """
    Exibe a interface gráfica baseada em texto para navegação simplificada.
    """
    print("\n" + "="*65)
    print("       SAQE - INSPEÇÃO DE ON-SCORES DE BONO (TUTA ABSOLUTA)     ")
    print("="*65)
    print("  [1] Exibir Resumo Estatístico do Transcriptoma Completo")
    print("  [2] Listar Transcritos Super-expressos Core (Net On-Score = +3)")
    print("  [3] Listar Transcritos Reprimidos Core (Net On-Score = -3)")
    print("  [4] Filtrar por Net On-Score Líquido Específico (-3 a +3)")
    print("  [5] Buscar Transcrito por ID ou Gene Symbol")
    print("  [6] Sair")
    print("="*65)

def main():
    df, caminho = obter_dados()
    limpar_terminal()
    
    while True:
        menu_principal()
        try:
            opcao = input("👉 Escolha uma opção [1-6]: ").strip()
        except (KeyboardInterrupt, EOFError):
            print("\n\nSaindo... Bons estudos e bom trabalho com o manuscrito!")
            break
            
        if opcao == '1':
            limpar_terminal()
            total = len(df)
            up = len(df[df['Bono_Net_Score'] > 0])
            down = len(df[df['Bono_Net_Score'] < 0])
            neutral = len(df[df['Bono_Net_Score'] == 0])
            
            print("\n" + "-"*55)
            print("📊 RESUMO GERAL DO TRANSCRIPTOMA (METODOLOGIA BONO)")
            print("-"*55)
            print(f"Ficheiro de origem:   {caminho}")
            print(f"Total de Transcritos: {total:,}")
            print(f"Induzidos (UP):       {up:,} ({up/total*100:.2f}%)")
            print(f"Reprimidos (DOWN):    {down:,} ({down/total*100:.2f}%)")
            print(f"Neutros (On-Score=0): {neutral:,} ({neutral/total*100:.2f}%)")
            print("-"*55)
            
            print("\nContagem Detalhada por Net On-Score Líquido:")
            for val in [-3, -2, -1, 0, 1, 2, 3]:
                cont = len(df[df['Bono_Net_Score'] == val])
                print(f"   • Net On-Score {val: >2}: {cont: >6,} transcritos ({cont/total*100: >5.2f}%)")
            print("-"*55)
            input("\nPressione [Enter] para continuar...")
            limpar_terminal()
            
        elif opcao == '2':
            limpar_terminal()
            df_core = df[df['Bono_Net_Score'] == 3]
            print(f"\n🔥 TRANSCRITOS COM INDUÇÃO UNIVERSAL (Net On-Score = +3) | Total: {len(df_core):,}")
            print("-"*110)
            if len(df_core) > 0:
                cols = ['Bono_UP_Score', 'Bono_DOWN_Score', 'Bono_Net_Score']
                if 'Annotation' in df_core.columns:
                    cols.append('Annotation')
                if 'Detox_Family' in df_core.columns:
                    cols.append('Detox_Family')
                
                print(df_core[cols].head(30))
                if len(df_core) > 30:
                    print(f"\n... e mais {len(df_core)-30} transcritos no ficheiro.")
            else:
                print("Nenhum transcrito com Net On-Score de +3 encontrado nos thresholds estipulados.")
            input("\nPressione [Enter] para continuar...")
            limpar_terminal()
            
        elif opcao == '3':
            limpar_terminal()
            df_core_down = df[df['Bono_Net_Score'] == -3]
            print(f"\n❄️  TRANSCRITOS COM REPRESSÃO UNIVERSAL (Net On-Score = -3) | Total: {len(df_core_down):,}")
            print("-"*110)
            if len(df_core_down) > 0:
                cols = ['Bono_UP_Score', 'Bono_DOWN_Score', 'Bono_Net_Score']
                if 'Annotation' in df_core_down.columns:
                    cols.append('Annotation')
                if 'Detox_Family' in df_core_down.columns:
                    cols.append('Detox_Family')
                    
                print(df_core_down[cols].head(30))
                if len(df_core_down) > 30:
                    print(f"\n... e mais {len(df_core_down)-30} transcritos no ficheiro.")
            else:
                print("Nenhum transcrito com Net On-Score de -3 encontrado nos thresholds estipulados.")
            input("\nPressione [Enter] para continuar...")
            limpar_terminal()
            
        elif opcao == '4':
            limpar_terminal()
            try:
                score_alvo = int(input("Digite o Net On-Score desejado [-3 a 3]: ").strip())
                if score_alvo < -3 or score_alvo > 3:
                    print("❌ Score inválido. Digite um valor entre -3 e 3.")
                    continue
            except ValueError:
                print("❌ Entrada inválida. Digite um número inteiro.")
                continue
                
            df_sub = df[df['Bono_Net_Score'] == score_alvo]
            print(f"\n🔍 Transcritos com Net On-Score = {score_alvo} | Total: {len(df_sub):,}")
            print("-"*110)
            if len(df_sub) > 0:
                cols = ['Bono_UP_Score', 'Bono_DOWN_Score', 'Bono_Net_Score']
                if 'Annotation' in df_sub.columns:
                    cols.append('Annotation')
                if 'Detox_Family' in df_sub.columns:
                    cols.append('Detox_Family')
                    
                print(df_sub[cols].head(25))
                if len(df_sub) > 25:
                    print(f"\n... e mais {len(df_sub)-25} transcritos.")
            else:
                print("Nenhum transcrito encontrado para este score de consistência líquida.")
            input("\nPressione [Enter] para continuar...")
            limpar_terminal()
            
        elif opcao == '5':
            limpar_terminal()
            termo = input("Digite o ID do Transcrito ou Gene de busca (ex: Tabs_g, CYP6, GST): ").strip()
            
            # Filtro robusto tolerante a maiúsculas e buscas parciais
            df_match = df[df.index.str.contains(termo, case=False, na=False)]
            if 'Annotation' in df.columns:
                df_match_annot = df[df['Annotation'].str.contains(termo, case=False, na=False)]
                df_match = pd.concat([df_match, df_match_annot]).drop_duplicates()
                
            print(f"\n🔍 Correspondências encontradas para '{termo}' | Total: {len(df_match):,}")
            print("-"*110)
            if len(df_match) > 0:
                cols = ['Bono_UP_Score', 'Bono_DOWN_Score', 'Bono_Net_Score']
                if 'Annotation' in df_match.columns:
                    cols.append('Annotation')
                if 'Detox_Family' in df_match.columns:
                    cols.append('Detox_Family')
                    
                print(df_match[cols].head(40))
                if len(df_match) > 40:
                    print(f"\n... e mais {len(df_match)-40} correspondências listadas no ficheiro.")
            else:
                print("Nenhuma correspondência molecular ou descrição coincidente encontrada.")
            input("\nPressione [Enter] para continuar...")
            limpar_terminal()
            
        elif opcao == '6':
            print("\nSaindo da ferramenta... Bons estudos e bom trabalho com o seu manuscrito!")
            break
        else:
            print("❌ Opção inválida. Escolha um número de 1 a 6.")

if __name__ == "__main__":
    main()
