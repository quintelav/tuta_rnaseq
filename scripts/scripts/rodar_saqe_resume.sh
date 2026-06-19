#!/bin/bash
set -e

# Verificar argumento
if [ -z "$1" ]; then
    echo "❌ Erro: especifique o bloco!"
    echo "Uso: ./rodar_saqe_por_bloco.sh [abamectina | spinosad | emamectina]"
    exit 1
fi

BLOCO=$1
mkdir -p sra_bruto raw_data raw_data/trimmed quant

# Definir amostras
if [ "$BLOCO" == "abamectina" ]; then
    AMOSTRAS=(
        "SRR35985995" "SRR35985994" "SRR35985991" "SRR35985990"
        "SRR35985989" "SRR35985988" "SRR35985987" "SRR35985986"
        "SRR35985985" "SRR35985984" "SRR35985993" "SRR35985992"
    )
elif [ "$BLOCO" == "spinosad" ]; then
    AMOSTRAS=(
        "SRR33779431" "SRR33779430" "SRR33779429"
        "SRR33779428" "SRR33779427" "SRR33779426"
    )
elif [ "$BLOCO" == "emamectina" ]; then
    AMOSTRAS=(
        "SRR15248447" "SRR15248446" "SRR15248443" "SRR15248442"
        "SRR15248441" "SRR15248440" "SRR15248439" "SRR15248438"
    )
else
    echo "❌ Bloco inválido!"
    exit 1
fi

TOTAL=${#AMOSTRAS[@]}
CONTADOR=1

echo "🧬 Iniciando SAQE - Bloco: $BLOCO ($TOTAL amostras)"

for SRA in "${AMOSTRAS[@]}"; do

    echo "📊 Amostra [$CONTADOR/$TOTAL]: $SRA"

    # Pular se já quantificado
    if [ -f "quant/${SRA}/quant.sf" ]; then
        echo "⏩ $SRA já quantificado. Pulando..."
        CONTADOR=$((CONTADOR + 1))
        continue
    fi

    # Procurar .sra já existente
    PATH_SRA=$(find sra_bruto -name "${SRA}.sra" 2>/dev/null | head -1)

    # Download apenas se necessário
    if [ -z "$PATH_SRA" ]; then
        echo "⬇️ Baixando $SRA..."
        prefetch $SRA --output-directory sra_bruto --max-size 100G
        PATH_SRA=$(find sra_bruto -name "${SRA}.sra" | head -1)
    else
        echo "📦 $SRA já baixado. Pulando prefetch..."
    fi

    # Conversão
    if [ ! -f "raw_data/${SRA}_1.fastq.gz" ]; then
        echo "🔄 fasterq-dump..."
        fasterq-dump --split-files --outdir raw_data "$PATH_SRA"

        echo "🗜️ Compressão..."
        pigz -f raw_data/${SRA}_1.fastq raw_data/${SRA}_2.fastq || \
        gzip -f raw_data/${SRA}_1.fastq raw_data/${SRA}_2.fastq
    else
        echo "📦 FASTQ já existe. Pulando conversão..."
    fi

    # Trim
    if [ ! -f "raw_data/trimmed/${SRA}_1_val_1.fq.gz" ]; then
        echo "✂️ Trim Galore..."
        trim_galore --paired --fastqc \
            --output_dir raw_data/trimmed \
            raw_data/${SRA}_1.fastq.gz \
            raw_data/${SRA}_2.fastq.gz
    else
        echo "📦 Trim já existe. Pulando..."
    fi

    # Salmon
    echo "🐟 Salmon quant..."
    salmon quant -i reference/salmon_index -l A \
        -1 raw_data/trimmed/${SRA}_1_val_1.fq.gz \
        -2 raw_data/trimmed/${SRA}_2_val_2.fq.gz \
        -p 6 --validateMappings \
        -o quant/${SRA}

    # Faxina
    rm -rf sra_bruto/${SRA}
    rm -f raw_data/${SRA}_1.fastq.gz raw_data/${SRA}_2.fastq.gz
    rm -f raw_data/trimmed/${SRA}_1_val_1.fq.gz
    rm -f raw_data/trimmed/${SRA}_2_val_2.fq.gz

    echo "✅ Concluído: $SRA"

    CONTADOR=$((CONTADOR + 1))
done

echo "🏆 Bloco $BLOCO finalizado!"
