# 🚀 Main Processing Pipeline

The project uses a block-based RNA-seq processing strategy.

Each insecticide dataset is processed independently:

- abamectina
- spinosad
- emamectina

The main script automates:

1. SRA download
2. FASTQ conversion
3. gzip compression
4. adapter trimming
5. FastQC
6. transcript quantification with Salmon
7. temporary file cleanup

---

# 📂 Expected Project Structure

```text
tuta_rnaseq/
│
├── raw_data/
│   ├── trimmed/
│   └── fastqc/
│
├── quant/
│
├── reference/
│   └── salmon_index/
│
├── sra_bruto/
│
├── scripts/
│
├── rodar_saqe_por_bloco.sh
└── README.md
```

---

# ⚙️ Required Software

Install Miniconda first.

## Create environment

```bash
conda create -n rnaseq python=3.11
conda activate rnaseq
```

## Install dependencies

```bash
conda install -c bioconda \
sra-tools \
fastqc \
multiqc \
trim-galore \
salmon \
pigz
```

---

# 🧬 Reference Transcriptome

The project uses the annotated transcriptome from:

GCA_027580185.1_ASM2758018v1

Download:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/580/185/GCA_027580185.1_ASM2758018v1/GCA_027580185.1_ASM2758018v1_rna_from_genomic.fna.gz
```

Decompress:

```bash
gunzip GCA_027580185.1_ASM2758018v1_rna_from_genomic.fna.gz
```

Rename:

```bash
mv GCA_027580185.1_ASM2758018v1_rna_from_genomic.fna tuta_transcripts.fna
```

---

# 🧱 Build Salmon Index

```bash
salmon index \
-t reference/tuta_transcripts.fna \
-i reference/salmon_index \
-k 31
```
# 🧾 Main Processing Script

The entire RNA-seq preprocessing and quantification workflow is automated using the script below.

File:

```text
rodar_saqe_por_bloco.sh
```

```bash
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

    pigz -f raw_data/${SRA}_1.fastq raw_data/${SRA}_2.fastq || \
    gzip -f raw_data/${SRA}_1.fastq raw_data/${SRA}_2.fastq

    # SAQE 03: Trimming
    trim_galore \
        --paired \
        --fastqc \
        --output_dir raw_data/trimmed/ \
        raw_data/${SRA}_1.fastq.gz \
        raw_data/${SRA}_2.fastq.gz

    # SAQE 04: Salmon Quantification
    salmon quant \
        -i reference/salmon_index \
        -l A \
        -1 raw_data/trimmed/${SRA}_1_val_1.fq.gz \
        -2 raw_data/trimmed/${SRA}_2_val_2.fq.gz \
        -p 6 \
        --validateMappings \
        -o quant/${SRA}

    # Cleanup
    rm -rf sra_bruto/${SRA}

    rm -f raw_data/${SRA}_1.fastq.gz \
          raw_data/${SRA}_2.fastq.gz

    rm -f raw_data/trimmed/${SRA}_1_val_1.fq.gz \
          raw_data/trimmed/${SRA}_2_val_2.fq.gz

    echo "✅ Concluído: $SRA"

    CONTADOR=$((CONTADOR + 1))
done

echo "🏆 Bloco $BLOCO finalizado com sucesso!"
```
---

# ▶️ Running the Pipeline

Make the script executable:

```bash
chmod +x rodar_saqe_por_bloco.sh
```

Run one insecticide block at a time:

## Abamectin

```bash
./rodar_saqe_por_bloco.sh abamectina
```

## Spinosad

```bash
./rodar_saqe_por_bloco.sh spinosad
```

## Emamectin benzoate

```bash
./rodar_saqe_por_bloco.sh emamectina
```

---

# 🔄 Pipeline Workflow

For each SRA sample the script performs:

```text
SRA download
    ↓
FASTQ conversion
    ↓
gzip compression
    ↓
Trim Galore + FastQC
    ↓
Salmon quantification
    ↓
cleanup
```

---

# 📤 Main Outputs

## Quantification files

```text
quant/SRRXXXXXXX/
```

Main file:

```text
quant.sf
```

Contains:

- TPM
- estimated counts
- transcript length
- transcript identifiers

---

# 📊 Quality Control

FastQC reports are generated automatically during trimming.

MultiQC summary:

```bash
multiqc raw_data/trimmed -o multiqc
```

---

# 🧹 Automatic Cleanup

To reduce disk usage, the script automatically removes:

- downloaded `.sra`
- intermediate FASTQ files
- trimmed FASTQ files after quantification

Final quantification results are preserved in:

```text
quant/
```
