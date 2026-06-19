#!/bin/bash
set -e

# Verificar se o usuário passou o argumento do bloco
if [ -z "$1" ]; then
    echo "❌ Erro: Você precisa especificar o bloco que deseja rodar!"
    echo "Uso correto: ./rodar_saqe_por_bloco.sh [abamectina | spinosad | emamectina]"
    exit 1
fi

BLOCO=$1
mkdir -p sra_bruto raw_data/trimmed quant

# Definir as amostras com base no bloco escolhido
if [ "$BLOCO" == "abamectina" ]; then
    AMOSTRAS=(
        "SRR35985995" "SRR35985994" "SRR35985991" "SRR35985990"
        "SRR35985989" "SRR35985988" "SRR35985987" "SRR35985986"
        "SRR35985985" "SRR35985984" "SRR35985993" "SRR35985992"
    )
elif [ "$BLOCO" == "spinosad" ]; then
    AMOSTRAS=("SRR33779431" "SRR33779430" "SRR33779429" "SRR33779428" "SRR33779427" "SRR33779426")
elif [ "$BLOCO" == "emamectina" ]; then
    AMOSTRAS=(
        "SRR15248447" "SRR15248446" "SRR15248443" "SRR15248442"
        "SRR15248441" "SRR15248440" "SRR15248439" "SRR15248438"
    )
else
    echo "❌ Bloco inválido! Escolha entre: abamectina, spinosad ou emamectina."
    exit 1
fi

TOTAL=${#AMOSTRAS[@]}
CONTADOR=1

echo "🧬 Iniciando SAQE - Bloco: $BLOCO ($TOTAL amostras)..."

for SRA in "${AMOSTRAS[@]}"; do
    echo "📊 Amostra [$CONTADOR/$TOTAL]: $SRA"

    # SAQE 01: Download
    prefetch $SRA --output-directory sra_bruto/

    # SAQE 02: Conversão e Compressão
    PATH_SRA=$(find sra_bruto/ -name "$SRA.sra")
    fasterq-dump --split-files --outdir raw_data/ "$PATH_SRA"
    pigz -f raw_data/${SRA}_1.fastq raw_data/${SRA}_2.fastq || gzip -f raw_data/${SRA}_1.fastq raw_data/${SRA}_2.fastq

    # SAQE 03: Trimming
    trim_galore --paired --fastqc --output_dir raw_data/trimmed/ raw_data/${SRA}_1.fastq.gz raw_data/${SRA}_2.fastq.gz

    # SAQE 05: Salmon Quant
    salmon quant -i reference/salmon_index -l A \
        -1 raw_data/trimmed/${SRA}_1_val_1.fq.gz \
        -2 raw_data/trimmed/${SRA}_2_val_2.fq.gz \
        -p 6 --validateMappings -o quant/${SRA}

    # Faxina
    rm -rf sra_bruto/${SRA}
    rm -f raw_data/${SRA}_1.fastq.gz raw_data/${SRA}_2.fastq.gz
    rm -f raw_data/trimmed/${SRA}_1_val_1.fq.gz raw_data/trimmed/${SRA}_2_val_2.fq.gz

    echo "✅ Concluído: $SRA"
    CONTADOR=$((CONTADOR + 1))
done

echo "🏆 Bloco $BLOCO finalizado com sucesso!"
