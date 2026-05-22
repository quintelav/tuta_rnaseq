import urllib.request
import json
import time

projetos = {
    "abamectina": "PRJNA1359366",
    "spinosad": "PRJNA1270399",
    "emamectina": "PRJNA749726"
}

for nome, bioproject in projetos.items():
    print(f"📥 Consultando IDs internos para o projeto {nome} ({bioproject})...")
    
    # ETAPA 1: Buscar os identificadores numéricos internos (UIDs) no banco SRA
    url_search = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={bioproject}&retmode=json&retmax=500"
    
    try:
        requisicao = urllib.request.urlopen(url_search)
        dados_busca = json.loads(requisicao.read().decode('utf-8'))
        lista_ids = dados_busca.get('esearchresult', {}).get('idlist', [])
        
        if not lista_ids:
            print(f"⚠️ Alerta: Nenhuma amostra encontrada para o projeto {bioproject}.\n")
            continue
            
        print(f"   -> {len(lista_ids)} runs localizadas. Extraindo tabela RunInfo...")
        
        # ETAPA 2: Fazer o efetch coletando o formato tabular runinfo real
        ids_formatados = ",".join(lista_ids)
        url_fetch = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={ids_formatados}&rettype=runinfo&retmode=text"
        
        urllib.request.urlretrieve(url_fetch, f"info_{nome}.csv")
        print(f"✅ Arquivo 'info_{nome}.csv' salvo com sucesso!\n")
        
        # Pausa estratégica para evitar que o NCBI bloqueie seu IP por excesso de requisições
        time.sleep(1)
        
    except Exception as e:
        print(f"❌ Falha crítica no projeto {nome}: {e}\n")
