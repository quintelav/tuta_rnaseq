#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Alinhamento DIAMOND Local Otimizado contra Swiss-Prot
Descarrega a base de dados curada Swiss-Prot (~250MB), extrai apenas as sequências
dos transcritos candidatos da metanálise e corre o DIAMOND de forma ultra-rápida.
"""

import os
import sys
import shutil
import urllib.request
import subprocess
import pandas as pd

def verificar_diamond():
    """
    Valida se o utilitário DIAMOND está instalado e acessível no ambiente Conda/WSL.
    """
    if shutil.which("diamond") is None:
        print("❌ [ERRO] O utilitário 'diamond' não foi encontrado no seu terminal!")
        print("Dica: Instale-o no seu ambiente Conda ativo executando o comando:")
        print("      conda install -c bioconda diamond")
        return False
    return True

def extrair_sequencias_candidatas(csv_path, fasta_completo, fasta_saida):
    """
    Otimização de Disco: Lê a lista de IDs candidatos e extrai apenas as suas
    sequências correspondentes do FASTA de referência geral.
    """
    print("🧬 A ler a lista de transcritos candidatos da metanálise...")
    try:
        df = pd.read_csv(csv_path, index_col=0)
        # Considera candidatos ativos com On-Score >= 1
        candidatos = set(df.index.tolist())
    except Exception as e:
        print(f"❌ Erro ao ler planilha de candidatos: {e}")
        return False

    print(f"📊 Total de transcritos candidatos identificados: {len(candidatos)}")
    print(f"📖 A extrair sequências correspondentes de '{fasta_completo}'...")

    try:
        total_extraidos = 0
        escrevendo = False
        with open(fasta_completo, 'r', encoding='utf-8', errors='ignore') as f_in, \
             open(fasta_saida, 'w', encoding='utf-8') as f_out:
            for linha in f_in:
                if linha.startswith('>'):
                    cabecalho = linha[1:].strip()
                    id_seq = cabecalho.split()[0]
                    # Verifica correspondência direta ou com variações de prefixos rna-
                    id_seq_sem_rna = id_seq[4:] if id_seq.startswith("rna-") else id_seq
                    
                    if id_seq in candidatos or id_seq_sem_rna in candidatos or f"rna-{id_seq}" in candidatos:
                        escrevendo = True
                        f_out.write(linha)
                        total_extraidos += 1
                    else:
                        escrevendo = False
                elif escrevendo:
                    f_out.write(linha)
                    
        print(f"   ✅ Extraídas {total_extraidos} sequências candidatas em '{fasta_saida}'.")
        return total_extraidos > 0
    except Exception as e:
        print(f"❌ Falha crítica ao extrair sequências: {e}")
        return False

def descarregar_swissprot(db_fasta):
    """
    Descarrega de forma segura a base de dados de proteínas curadas Swiss-Prot.
    """
    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    out_gz = db_fasta + ".gz"
    
    if os.path.exists(db_fasta):
        print("   ✅ Base de dados Swiss-Prot em formato FASTA já existe localmente.")
        return True
        
    print("🌐 A descarregar a base curada Swiss-Prot do FTP oficial do UniProt (~250 MB)...")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=60) as response, open(out_gz, 'wb') as f_out:
            shutil.copyfileobj(response, f_out)
            
        print("📦 Decompressing uniprot_sprot.fasta.gz...")
        import gzip
        with gzip.open(out_gz, 'rb') as f_in, open(db_fasta, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
        os.remove(out_gz)
        print("   ✅ Download e descompactação concluídos!")
        return True
    except Exception as e:
        print(f"❌ Erro ao descarregar base de dados Swiss-Prot: {e}")
        return False

def executar_pipeline_diamond():
    print("==================================================")
    print("      SAQE - PIPELINE DIAMOND SWISS-PROT LOCAL     ")
    print("==================================================")
    
    if not verificar_diamond():
        return
        
    os.makedirs("reference", exist_ok=True)
    os.makedirs("quant", exist_ok=True)
    
    csv_candidatos = "quant/metanalise_completa_bono.csv"
    fasta_referencia = "reference/tuta_transcripts.fna"
    fasta_candidatos = "quant/candidatos_ativos.fasta"
    swissprot_fasta = "reference/uniprot_sprot.fasta"
    swissprot_dmnd = "reference/swissprot.dmnd"
    blast_out = "quant/blast_swissprot.txt"
    
    if not os.path.exists(csv_candidatos):
        print(f"❌ Planilha de candidatos '{csv_candidatos}' não encontrada!")
        print("Dica: Execute primeiro o pipeline diferencial: 'python3 analise_meta_bono.py'")
        return
        
    if not os.path.exists(fasta_referencia):
        print(f"❌ Ficheiro de sequências de referência '{fasta_referencia}' não encontrado!")
        return

    # Passo 1: Extrair sequências de interesse para poupar processamento
    if not extrair_sequencias_candidatas(csv_candidatos, fasta_referencia, fasta_candidatos):
        return
        
    # Passo 2: Descarregar a base curada
    if not descarregar_swissprot(swissprot_fasta):
        return
        
    # Passo 3: Compilar a base de dados no formato do DIAMOND
    if not os.path.exists(swissprot_dmnd):
        print("\n⚙️  A compilar e a indexar a base de dados Swiss-Prot para o DIAMOND...")
        cmd_makedb = ["diamond", "makedb", "--in", swissprot_fasta, "-d", "reference/swissprot"]
        res = subprocess.run(cmd_makedb)
        if res.returncode != 0:
            print("❌ Erro ao indexar base de dados com 'diamond makedb'!")
            return
    else:
        print("   ✅ Base de dados Swiss-Prot já indexada para o DIAMOND.")
        
    # Passo 4: Executar o alinhamento blastx local em lote
    print("\n🚀 A executar o alinhamento blastx ultra-rápido via DIAMOND...")
    # Usa modo de sensibilidade balanceada, outfmt 6 (tabular padrão), 6 threads para alta performance
    cmd_blastx = [
        "diamond", "blastx",
        "-q", fasta_candidatos,
        "-d", swissprot_dmnd,
        "-o", blast_out,
        "--outfmt", "6",
        "--threads", "6",
        "--max-target-seqs", "1",
        "--evalue", "1e-5"
    ]
    
    res_blast = subprocess.run(cmd_blastx)
    if res_blast.returncode == 0:
        print("\n==================================================")
        print("🎉 ALINHAMENTO CONCLUÍDO COM SUCESSO!")
        print("==================================================")
        print(f"📂 Ficheiro de homologia gerado em: '{blast_out}'")
        print("\n🧹 A limpar ficheiros temporários para libertar espaço...")
        
        # Elimina os fastas brutos para preservar o seu armazenamento em disco
        if os.path.exists(fasta_candidatos):
            os.remove(fasta_candidatos)
        if os.path.exists(swissprot_fasta):
            os.remove(swissprot_fasta)
            
        print("✅ Concluído! Agora execute o seu script de anotação funcional:")
        print("   👉 python3 anotar_resultados.py")
    else:
        print("❌ Erro durante a execução do DIAMOND blastx.")

if __name__ == "__main__":
    executar_pipeline_diamond()
