#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Script de Anotação Funcional e Classificação de Famílias de Detoxificação
Versão Metodologia Bono (2021) - Suporta mapeamento por GFF3, Homologia BLAST (Swiss-Prot),
anotações ricas do eggNOG-mapper (.emapper.annotations) ou mapeamento online via NCBI.
Inclui descarregador automático de GFF3 oficial do NCBI e gerador de gráficos.
"""

import os
import glob
import re
import sys
import json
import time
import gzip
import shutil
import urllib.request
import urllib.parse
import pandas as pd
import numpy as np

def normalizar_id(gene_id):
    """
    Remove prefixos estruturais comuns do NCBI para unificar chaves de busca.
    """
    if not isinstance(gene_id, str):
        return ""
    id_limpo = gene_id.strip()
    for prefixo in ['rna-', 'cds-', 'gene-', 'id-', 'transcript-']:
        if id_limpo.lower().startswith(prefixo):
            id_limpo = id_limpo[len(prefixo):]
    return id_limpo

def descarregar_gff3_ncbi():
    """
    Descarrega e descompacta automaticamente o GFF3 oficial de referência da Tuta absoluta 
    (Phthorimaea absoluta - GCF_022605935.1_ASM2260593v1) direto do FTP do NCBI.
    """
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/605/935/GCF_022605935.1_ASM2260593v1/GCF_022605935.1_ASM2260593v1_genomic.gff.gz"
    out_gz = "reference/GCF_022605935.1_genomic.gff.gz"
    out_gff = "reference/GCF_022605935.1_genomic.gff"
    
    print("\n🌐 A descarregar o ficheiro GFF3 oficial do NCBI para Phthorimaea absoluta (GCF_022605935.1)...")
    try:
        os.makedirs("reference", exist_ok=True)
        
        # Download do arquivo compactado (.gz)
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=40) as response, open(out_gz, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
            
        # Descompactação nativa
        print("📦 A descompactar o ficheiro genomic.gff.gz...")
        with gzip.open(out_gz, 'rb') as f_in, open(out_gff, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
        os.remove(out_gz) # Remove o arquivo .gz temporário para poupar espaço
        print("✅ Ficheiro GFF3 descompactado com sucesso na pasta 'reference/'!")
        return out_gff
    except Exception as e:
        print(f"❌ Erro ao descarregar o GFF3 automaticamente: {e}")
        print("Dica: Descarregue manualmente o GFF3 de GCF_022605935.1 e coloque em 'reference/'")
        return None

def parser_avancado_gff3(gff_path):
    """
    Lê o ficheiro GFF3/GTF hierarquicamente e mapeia as descrições funcionais.
    """
    print(f"📖 A abrir e a processar o ficheiro de anotação local: '{gff_path}'...")
    mapeamento_final = {}
    relacao_parentesco = {}
    funcoes_diretas = {}
    
    with open(gff_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attributes_str = parts[8]
            
            attrs = {}
            for attr in attributes_str.split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attrs[k.strip()] = urllib.parse.unquote(v.strip())
            
            feature_id = attrs.get('ID')
            parent_id = attrs.get('Parent')
            product = attrs.get('product') or attrs.get('Note') or attrs.get('description')
            transcript_id = attrs.get('transcript_id') or attrs.get('Name')
            
            if product:
                if feature_id:
                    funcoes_diretas[feature_id] = product
                if transcript_id:
                    funcoes_diretas[transcript_id] = product
                if parent_id:
                    for parent in parent_id.split(','):
                        funcoes_diretas[parent.strip()] = product
            
            if feature_type == 'CDS' and feature_id and parent_id:
                relacao_parentesco[feature_id] = parent_id

    for cds_id, parent_ids in relacao_parentesco.items():
        produto_cds = funcoes_diretas.get(cds_id)
        if produto_cds:
            for parent in parent_ids.split(','):
                parent_clean = parent.strip()
                if parent_clean not in funcoes_diretas:
                    funcoes_diretas[parent_clean] = produto_cds

    for k, v in funcoes_diretas.items():
        mapeamento_final[k] = v
        k_norm = normalizar_id(k)
        if k_norm:
            mapeamento_final[k_norm] = v
            mapeamento_final[f"rna-{k_norm}"] = v

    return mapeamento_final

def parse_eggnog_mapper(emapper_path):
    """
    NOVO PARSER REVOLUCIONÁRIO: Realiza a leitura estrita do arquivo .emapper.annotations
    gerado pelo eggNOG-mapper (localmente ou online).
    Mapeia de forma robusta o ID do transcrito para obter nome de gene, descrição, termos GO e KEGG.
    """
    print(f"📖 A abrir e a processar o relatório do eggNOG-mapper: '{emapper_path}'...")
    mapeamento_eggnog = {}
    
    try:
        colunas_relevantes = {}
        with open(emapper_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                # Localiza a linha de cabeçalho que define as posições das colunas
                if line.startswith('#query') or line.startswith('# query'):
                    cabecalho = line.strip().lstrip('#').strip().split('\t')
                    for idx, col in enumerate(cabecalho):
                        colunas_relevantes[col.strip().lower()] = idx
                    break
                    
        if not colunas_relevantes:
            # Layout padrão v2 do eggNOG-mapper caso o cabeçalho não seja encontrado
            colunas_relevantes = {
                'query': 0, 'preferred_name': 8, 'description': 7, 'gos': 9, 'kegg_ko': 11
            }
            
        idx_query = colunas_relevantes.get('query', 0)
        idx_name = colunas_relevantes.get('preferred_name', 8)
        idx_desc = colunas_relevantes.get('description', 7)
        idx_go = colunas_relevantes.get('gos', 9)
        idx_kegg = colunas_relevantes.get('kegg_ko', 11)
        
        contagem = 0
        with open(emapper_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) <= max(idx_query, idx_name, idx_desc):
                    continue
                    
                q_id = parts[idx_query]
                pref_name = parts[idx_name].strip()
                description = parts[idx_desc].strip()
                gos = parts[idx_go].strip() if len(parts) > idx_go else ""
                kegg = parts[idx_kegg].strip() if len(parts) > idx_kegg else ""
                
                # Formata uma descrição robusta combinando nome do gene e função molecular
                anotacao_detalhada = ""
                if pref_name and pref_name != "-":
                    anotacao_detalhada += f"[{pref_name}] "
                if description and description != "-":
                    anotacao_detalhada += description
                else:
                    anotacao_detalhada += "Proteína identificada via ortologia"
                    
                dados_completos = {
                    'annotation': anotacao_detalhada,
                    'gene_name': pref_name if pref_name != "-" else "Não especificado",
                    'go_terms': gos if gos != "-" else "Sem termos GO",
                    'kegg_pathways': kegg if kegg != "-" else "Sem KEGG"
                }
                
                mapeamento_eggnog[q_id] = dados_completos
                id_norm = normalizar_id(q_id)
                if id_norm:
                    mapeamento_eggnog[id_norm] = dados_completos
                    mapeamento_eggnog[f"rna-{id_norm}"] = dados_completos
                contagem += 1
                
        print(f"   ✅ Processadas {contagem} anotações de ortologia do eggNOG com sucesso!")
        return mapeamento_eggnog
    except Exception as e:
        print(f"⚠️  Erro crítico ao processar o arquivo eggNOG: {e}")
        return {}

def extrair_anotacoes_fasta_local(fasta_path):
    """
    Fallback Local: Lê o FASTA de referência e extrai anotações dos cabeçalhos.
    """
    print(f"🧬 A vasculhar o FASTA de referência local em busca de anotações: '{fasta_path}'...")
    mapeamento_fasta = {}
    if not os.path.exists(fasta_path):
        return mapeamento_fasta
        
    try:
        with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
            for linha in f:
                if linha.startswith('>'):
                    cabecalho = linha[1:].strip()
                    partes = cabecalho.split(maxsplit=1)
                    if partes:
                        id_seq = partes[0]
                        if len(partes) > 1:
                            desc = partes[1]
                            
                            lower_desc = desc.lower()
                            if any(x in lower_desc for x in ['cds=', 'location=', 'gbkey=', 'parent=']) and not any(x in lower_desc for x in ['protein=', 'product=']):
                                continue
                                
                            protein_match = re.search(r'\[protein=([^\]]+)\]', desc)
                            product_match = re.search(r'\[product=([^\]]+)\]', desc)
                            
                            if protein_match:
                                desc_limpa = protein_match.group(1)
                            elif product_match:
                                desc_limpa = product_match.group(1)
                            else:
                                desc_limpa = re.sub(r'^[A-Z0-9_]+\.[0-9]+\s+', '', desc)
                                
                            if len(desc_limpa) > 5 and not desc_limpa.lower().startswith("predicted protein"):
                                mapeamento_fasta[id_seq] = desc_limpa
                                id_norm = normalizar_id(id_seq)
                                mapeamento_fasta[id_norm] = desc_limpa
                                mapeamento_fasta[f"rna-{id_norm}"] = desc_limpa
    except Exception as e:
        print(f"⚠️  Não foi possível extrair anotações do FASTA de referência: {e}")
        
    print(f"   ✅ Extraídas {len(mapeamento_fasta)} descrições funcionais reais do FASTA.")
    return mapeamento_fasta

def buscar_uniprot_api(protein_id):
    """
    Consulta a API pública do UniProtKB para obter a anotação funcional com Retry automático.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    retries = 3
    delay = 1.0
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=5) as response:
                data = json.loads(response.read().decode('utf-8'))
                desc = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value')
                if desc:
                    return desc
        except Exception:
            time.sleep(delay)
            delay *= 2
    return "Proteína homóloga (UniProt ID: " + protein_id + ")"

def parse_blast_swiss_prot(blast_path):
    """
    Processa os resultados de uma pesquisa BLAST/DIAMOND (outfmt 6) contra o Swiss-Prot/UniProt.
    """
    print(f"📖 A processar os resultados de homologia BLAST/DIAMOND: '{blast_path}'...")
    mapeamento_blast = {}
    hits_temp = {}
    
    with open(blast_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            qseqid = parts[0]
            sseqid = parts[1]
            
            if qseqid not in hits_temp:
                uniprot_id = sseqid
                if '|' in sseqid:
                    uniprot_id = sseqid.split('|')[1]
                hits_temp[qseqid] = uniprot_id
                
    print(f"🧬 A procurar descrições funcionais na API do UniProt para os {len(hits_temp)} melhores hits...")
    total = len(hits_temp)
    for idx, (qseqid, uniprot_id) in enumerate(hits_temp.items(), 1):
        if idx % 10 == 0 or idx == total:
            print(f"   🔹 Progresso: {idx}/{total} proteínas consultadas...")
        desc = buscar_uniprot_api(uniprot_id)
        mapeamento_blast[qseqid] = desc
        mapeamento_blast[normalizar_id(qseqid)] = desc
        time.sleep(0.1)
        
    return mapeamento_blast

def buscar_ncbi_online_api_batch(transcript_ids):
    """
    FUNÇÃO EM LOTE (Batching) CORRIGIDA: Consulta o NCBI esummary para múltiplos accessions.
    """
    mapeamento_batch = {}
    ids_limpos = []
    mapa_original = {}
    
    for tid in transcript_ids:
        limpo = normalizar_id(tid)
        if '.' in limpo:
            limpo = limpo.split('.')[0]
        if '|' in limpo:
            limpo = limpo.split('|')[-1]
        
        if limpo.startswith("Tabs_g") or "gnl" in tid:
            continue
            
        ids_limpos.append(limpo)
        mapa_original[limpo] = tid

    if not ids_limpos:
        return mapeamento_batch

    tamanho_lote = 50
    for i in range(0, len(ids_limpos), tamanho_lote):
        lote = ids_limpos[i:i+tamanho_lote]
        
        term = " OR ".join([f"{x}[Accession]" for x in lote])
        query_params = urllib.parse.urlencode({'db': 'nuccore', 'term': term, 'retmode': 'json'})
        url_search = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{query_params}"
        
        retries = 3
        delay = 1.0
        uids = []
        
        for r in range(retries):
            try:
                req = urllib.request.Request(url_search, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=10) as response:
                    res_json = json.loads(response.read().decode('utf-8'))
                    uids = res_json.get('esearchresult', {}).get('idlist', [])
                    break
            except Exception:
                time.sleep(delay)
                delay *= 2
                
        if not uids:
            continue
            
        uids_str = ",".join(uids)
        query_params_sum = urllib.parse.urlencode({'db': 'nuccore', 'id': uids_str, 'retmode': 'json'})
        url_summary = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{query_params_sum}"
        
        summary_data = {}
        for r in range(retries):
            try:
                time.sleep(0.35)
                req_sum = urllib.request.Request(url_summary, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req_sum, timeout=10) as summary_response:
                    summary_data = json.loads(summary_response.read().decode('utf-8')).get('result', {})
                    break
            except Exception:
                time.sleep(delay)
                delay *= 2
                
        for uid in summary_data:
            if uid == "uids":
                continue
            registro = summary_data[uid]
            ac_ver = registro.get('accessionversion', '')
            ac_sem_ver = ac_ver.split('.')[0] if ac_ver else ''
            title = registro.get('title', '')
            
            if title:
                if ", mRNA" in title:
                    title = title.split(", mRNA")[0]
                if "Tuta absoluta" in title:
                    title = title.replace("Tuta absoluta", "").strip()
                title = re.sub(r'^[A-Z0-9_]+\.[0-9]+\s+', '', title)
                
                for chave_ac in [ac_ver, ac_sem_ver]:
                    if chave_ac in mapa_original:
                        orig_id = mapa_original[chave_ac]
                        mapeamento_batch[orig_id] = title

    return mapeamento_batch

def classificar_familia_detox(annotation_text):
    """
    Classifica o texto da anotação funcional nas famílias metabólicas clássicas de detoxificação.
    """
    if not isinstance(annotation_text, str) or annotation_text in ["Sem anotação disponível", "Sem descrição molecular disponível", "Pendente", "Sem anotação funcional encontrada", "Sem anotação funcional encontrada"]:
        return "Outros / Não Identificado"
        
    text = annotation_text.lower()
    
    # Fase I: Oxidação e Hidrólise
    if any(p in text for p in ["cytochrome p450", "cyp", "monooxygenase"]):
        return "Fase I: Cytochrome P450 (CYP)"
    if any(p in text for p in ["carboxylesterase", "esterase", "acetylcholinesterase"]):
        return "Fase I: Carboxylesterase (CCE)"
        
    # Fase II: Conjugação
    if any(p in text for p in ["glutathione s-transferase", "gst"]):
        return "Fase II: Glutathione S-transferase (GST)"
    if any(p in text for p in ["udp-glucuronosyltransferase", "ugt", "glucuronyltransferase"]):
        return "Fase II: UDP-glucuronosyltransferase (UGT)"
        
    # Fase III: Excreção/Transporte
    if any(p in text for p in ["abc transporter", "p-glycoprotein", "multidrug resistance"]):
        return "Fase III: Transportadores ABC"
        
    # Antioxidantes e Estresse Oxidativo (Bono, 2021)
    if any(p in text for p in ["superoxide dismutase", "sod"]):
        return "Antioxidante: Superoxide Dismutase (SOD)"
    if "catalase" in text:
        return "Antioxidante: Catalase (CAT)"
    if "peroxidase" in text:
        return "Antioxidante: Peroxidase"
        
    return "Outros Processos Fisiológicos"

def gerar_grafico_distribuicao(df_resultados):
    """
    Gera automaticamente um gráfico de barras com alta resolução técnica (300 DPI) para SSH/WSL.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        df_plot = df_resultados[df_resultados['Detox_Family'] != "Outros Processos Fisiológicos"]
        df_plot = df_plot[df_plot['Detox_Family'] != "Outros / Não Identificado"]
        
        if len(df_plot) == 0:
            return
            
        plt.figure(figsize=(10, 6), dpi=300)
        sns.set_theme(style="ticks")
        
        paleta = sns.color_palette("viridis", len(df_plot['Detox_Family'].unique()))
        order = df_plot['Detox_Family'].value_counts().index
        ax = sns.countplot(y='Detox_Family', data=df_plot, order=order, palette=paleta)
        
        plt.title("Distribuição das Famílias Metabólicas de Detoxificação\n(Genes Candidatos da Metanálise - Bono, 2021)", 
                  fontsize=12, fontweight='bold', pad=15)
        plt.xlabel("Contagem de Transcritos Diferencialmente Ativos", fontsize=10)
        plt.ylabel("Família Funcional de Detoxificação", fontsize=10)
        
        for p in ax.patches:
            width = p.get_width()
            ax.text(width + 0.3, p.get_y() + p.get_height()/2 + 0.07, 
                    f'{int(width)}', ha="left", va="center", fontsize=9, fontweight='bold')
            
        sns.despine(trim=True)
        plt.tight_layout()
        
        fig_path = "quant/distribuicao_familias_detox.png"
        plt.savefig(fig_path, bbox_inches='tight')
        plt.close()
        print(f"📊 [SUCESSO] Gráfico de publicação guardado em '{fig_path}'!")
    except ImportError:
        print("\n💡 Nota: Matplotlib ou Seaborn não estão instalados no seu ambiente Conda.")

def executar_anotacao():
    print("==================================================")
    print("      SAQE - ANOTAÇÃO FUNCIONAL E CURADORIA       ")
    print("           (MÉTODO DE BONO / GFF3-FREE)           ")
    print("==================================================")
    
    results_file = "quant/metanalise_completa_bono.csv"
    fasta_file = "reference/tuta_transcripts.fna"
    
    if not os.path.exists(results_file):
        print(f"❌ [ERRO] O ficheiro de resultados '{results_file}' não foi encontrado!")
        print("Dica: Execute primeiro o pipeline de expressão diferencial: 'python3 analise_meta_bono.py'")
        return
        
    df_resultados = pd.read_csv(results_file, index_col=0)
    print(f"📊 Carregados {len(df_resultados)} transcritos totais da metanálise.")
    
    # --- FILTRAGEM INTERATIVA ---
    print("\n💡 Opções de filtragem para acelerar a anotação:")
    print("   [1] Apenas genes consistentes (On-Score >= 2) - Recomendado e Rápido (~394 genes)")
    print("   [2] Apenas genes core universais (On-Score = 3) - Ultra rápido (~16 genes)")
    print("   [3] Todos os transcritos ativos gerais (On-Score >= 1) - Completo (~2424 genes)")
    
    escolha = input("\n👉 Escolha uma opção (1, 2 ou 3) [Padrão: 1]: ").strip()
    if escolha == "2":
        df_filtrado = df_resultados[(df_resultados['Bono_UP_Score'] == 3) | (df_resultados['Bono_DOWN_Score'] == 3)].copy()
    elif escolha == "3":
        df_filtrado = df_resultados[(df_resultados['Bono_UP_Score'] >= 1) | (df_resultados['Bono_DOWN_Score'] >= 1)].copy()
    else:
        df_filtrado = df_resultados[(df_resultados['Bono_UP_Score'] >= 2) | (df_resultados['Bono_DOWN_Score'] >= 2)].copy()
        
    print(f"\n🎯 Filtrados {len(df_filtrado)} transcritos para receberem anotação funcional.")
    
    # Dicionário mestre
    anotacoes_acumuladas = {}
    
    # --- PASSO A: DETECTAR RELATÓRIO DO EGGNOG-MAPPER (PRIORIDADE MÁXIMA) ---
    arquivos_eggnog = glob.glob(os.path.join("quant", "*.emapper.annotations")) + \
                      glob.glob(os.path.join("reference", "*.emapper.annotations")) + \
                      glob.glob("*.emapper.annotations")
                      
    usando_eggnog = False
    if arquivos_eggnog:
        print(f"\n🎯 [EGGNOG DETETADO] Encontrado ficheiro de anotações do eggNOG-mapper: '{arquivos_eggnog[0]}'")
        mapa_eggnog = parse_eggnog_mapper(arquivos_eggnog[0])
        if mapa_eggnog:
            anotacoes_acumuladas.update(mapa_eggnog)
            usando_eggnog = True
            
    # --- PASSO B: VERIFICAÇÃO AUTOMÁTICA DE IDS 'Tabs_g' E DOWNLOAD DO GFF3 DO NCBI ---
    arquivos_gff = glob.glob(os.path.join("reference", "*"))
    gff_selecionado = [f for f in arquivos_gff if f.endswith(('.gff', '.gff3', '.gtf'))]
    
    possui_ids_tuta = any("Tabs_g" in str(x) or "JAQAHW" in str(x) for x in df_filtrado.index)
    
    if possui_ids_tuta and not gff_selecionado and not usando_eggnog:
        print("\n⚠️  [DETETADO] Os seus identificadores de transcritos ('Tabs_g...') pertencem à montagem de referência de Phthorimaea absoluta (Tuta absoluta).")
        print("Para anotar estes IDs customizados locais, é necessário o ficheiro GFF3 de referência oficial.")
        resposta = input("👉 Deseja descarregar o ficheiro GFF3 oficial (GCF_022605935.1) automaticamente do NCBI? (S/N) [Padrão: S]: ").strip().lower()
        if resposta != 'n':
            gff_baixado = descarregar_gff3_ncbi()
            if gff_baixado:
                gff_selecionado = [gff_baixado]
                
    # --- PASSO C: BUSCA LOCAL NO FASTA DE REFERÊNCIA (Se não usou eggNOG) ---
    if not usando_eggnog:
        anotacoes_fasta_local = extrair_anotacoes_fasta_local(fasta_file)
        anotacoes_acumuladas.update(anotacoes_fasta_local)
    
    # --- PASSO D: DECISION TREE PARA ANOTAÇÕES COMPLEMENTARES ---
    arquivos_blast = glob.glob(os.path.join("quant", "*blast*")) + glob.glob("*.txt") + glob.glob("*.tsv")
    blast_selecionado = [f for f in arquivos_blast if "blast" in f.lower() or "swissprot" in f.lower() or "diamond" in f.lower()]
    
    metodo_complementar = "eggNOG-mapper (Ortologia)" if usando_eggnog else "Nenhum (FASTA apenas)"
    
    if not usando_eggnog:
        if gff_selecionado:
            print("\n📁 [MODO LOCAL] Ficheiro GFF3 detetado. Mapeando coordenadas locais...")
            mapa_gff = parser_avancado_gff3(gff_selecionado[0])
            anotacoes_acumuladas.update(mapa_gff)
            metodo_complementar = "GFF3 Local"
        elif blast_selecionado:
            print("\n🎯 [MODO LOCAL] Tabela BLAST/DIAMOND de homologia detetada. Consultando Swiss-Prot...")
            mapa_blast = parse_blast_swiss_prot(blast_selecionado[0])
            anotacoes_acumuladas.update(mapa_blast)
            metodo_complementar = "Homologia BLAST/DIAMOND"
        else:
            print("\n🌐 [MODO REMOTO] Procurando dados em lote via API do NCBI...")
            ids_para_consultar = df_filtrado.index.tolist()
            ids_em_falta = [x for x in ids_para_consultar if x not in anotacoes_acumuladas and normalizar_id(x) not in anotacoes_acumuladas]
            
            if ids_em_falta:
                print(f"   🔹 Consultando {len(ids_em_falta)} accessions pendentes no NCBI...")
                mapa_ncbi_batch = buscar_ncbi_online_api_batch(ids_em_falta)
                anotacoes_acumuladas.update(mapa_ncbi_batch)
                metodo_complementar = "NCBI Entrez Batch API"
            else:
                print("   ✅ Todas as descrições necessárias já foram recuperadas do FASTA/GFF3 local!")
                metodo_complementar = "FASTA/GFF3 Local Direto"
            
    # --- CRUZAMENTO DE DADOS FINAL ---
    print(f"\n🧬 Cruzando identificadores com a tabela de resultados...")
    anotacoes_finais = []
    nomes_genes = []
    termos_go = []
    vias_kegg = []
    hits = 0
    
    for id_transcrito in df_filtrado.index:
        func = None
        gene_name = "Não especificado"
        go = "Sem termos GO"
        kegg = "Sem KEGG"
        
        for chave in [id_transcrito, normalizar_id(id_transcrito), f"rna-{normalizar_id(id_transcrito)}"]:
            if chave in anotacoes_acumuladas:
                dados = anotacoes_acumuladas[chave]
                if isinstance(dados, dict):
                    # Se veio do eggNOG, extrai o dicionário completo de atributos
                    func = dados.get('annotation')
                    gene_name = dados.get('gene_name', 'Não especificado')
                    go = dados.get('go_terms', 'Sem termos GO')
                    kegg = dados.get('kegg_pathways', 'Sem KEGG')
                else:
                    func = dados
                break
                
        if func:
            anotacoes_finais.append(func)
            nomes_genes.append(gene_name)
            termos_go.append(go)
            vias_kegg.append(kegg)
            hits += 1
        else:
            anotacoes_finais.append("Sem anotação funcional encontrada")
            nomes_genes.append("Não especificado")
            termos_go.append("Sem termos GO")
            vias_kegg.append("Sem KEGG")
            
    df_filtrado['Annotation'] = anotacoes_finais
    
    # Se o eggNOG foi utilizado, injeta colunas ricas de bioinformática na tabela final
    if usando_eggnog:
        df_filtrado['Preferred_Gene_Name'] = nomes_genes
        df_filtrado['GO_Terms'] = termos_go
        df_filtrado['KEGG_Pathways'] = vias_kegg
        
    print(f"\n📊 Resultados do cruzamento de dados:")
    print(f"   🔹 Método principal de anotação: {metodo_complementar}")
    print(f"   🔹 Transcritos analisados: {len(df_filtrado)}")
    print(f"   🔹 Transcritos com anotação recuperada: {hits} ({hits/len(df_filtrado)*100:.1f}%)")
    
    # Classificar nas famílias metabólicas de detoxificação
    df_filtrado['Detox_Family'] = df_filtrado['Annotation'].apply(classificar_familia_detox)
    
    # Salvar a planilha master unificada anotada
    output_master = "quant/metanalise_completa_anotada.csv"
    df_filtrado.to_csv(output_master)
    print(f"💾 Tabela consolidada salva em: '{output_master}'")
    
    # Distribuição estatística das famílias
    print("\n==================================================")
    print("📊 DISTRIBUIÇÃO DAS FAMÍLIAS DE DETOXIFICAÇÃO:")
    print("==================================================")
    contagem_familias = df_filtrado['Detox_Family'].value_counts()
    for familia, total_cont in contagem_familias.items():
        print(f"   🔹 {familia}: {total_cont} transcritos")
        
    # Exportar datasets individuais por família de detoxificação para o Excel/R
    pasta_detox = os.path.join("quant", "detox_families")
    os.makedirs(pasta_detox, exist_ok=True)
    
    familias_alvo = [
        ("Fase I: Cytochrome P450 (CYP)", "cyp_p450.csv"),
        ("Fase I: Carboxylesterase (CCE)", "carboxylesterases.csv"),
        ("Fase II: Glutathione S-transferase (GST)", "gst_transferases.csv"),
        ("Fase II: UDP-glucuronosyltransferase (UGT)", "ugt_transferases.csv"),
        ("Fase III: Transportadores ABC", "abc_transporters.csv"),
        ("Antioxidante: Superoxide Sodumase (SOD)", "sod_antiox.csv"),
        ("Antioxidante: Catalase (CAT)", "catalase_antiox.csv")
    ]
    
    for nome_familia, nome_arquivo in familias_alvo:
        sub_df = df_filtrado[df_filtrado['Detox_Family'] == nome_familia]
        if len(sub_df) > 0:
            sub_df.to_csv(os.path.join(pasta_detox, nome_arquivo))
            
    # Geração de gráficos para publicação
    gerar_grafico_distribuicao(df_filtrado)
    
    print("\n🎉 Processo de anotação funcional concluído com sucesso!")

if __name__ == "__main__":
    executar_anotacao()