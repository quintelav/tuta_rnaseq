#!/bin/bash

# Ativar o ambiente do conda de forma segura dentro do script
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# LISTA DE ACESSOS (Substitua ou adicione os SRRs reais do projeto aqui)
AMOSTRAS=("SRR35985985" "SRR35985986" "SRR35985987")

for SRR in "${AMOSTRAS[@]}"
do
    echo "=========================================================="
    echo "Iniciando processamento da amostra: $SRR"
    echo "=========================================================="

    # 1. Download do arquivo bruto (.sra)
    echo "[1/4] Baixando arquivo via prefetch..."
    prefetch $SRR

    # 2. Extração direta para FASTQ compactado (.gz) para poupar o HD
    echo "[2/4] Extraindo sequências para FASTQ.GZ..."
    fastq-dump --split-files --gzip -O raw_data $SRR

    # 3. Quantificação de alta performance no Salmon usando o RefSeq de 95MB
    echo "[3/4] Rodando a quantificação no Salmon..."
    salmon quant \
      -i reference/salmon_index \
      -l A \
      -1 raw_data/${SRR}_1.fastq.gz \
      -2 raw_data/${SRR}_2.fastq.gz \
      -p 6 \
      --validateMappings \
      -o quant/$SRR

    # 4. Limpeza cirúrgica do armazenamento
    echo "[4/4] Limpando arquivos brutos para liberar espaço..."
    rm -rf raw_data/${SRR}*
    rm -rf ~/ncbi/public/sra/${SRR}.sra
    rm -rf ${SRR}/

    echo "Concluído com sucesso para $SRR!"
    echo "----------------------------------------------------------"
done

echo "Todas as amostras selecionadas foram processadas!"
