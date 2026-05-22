### 🧬 RNA-seq Reanalysis Pipeline in Phthorimaea absoluta

Reanalysis of publicly available RNA-seq datasets to investigate transcriptional responses associated with insecticide exposure and resistance in Phthorimaea absoluta using a lightweight reproducible workflow based on Salmon pseudoalignment.

### 📌 Project Overview

This repository contains a complete RNA-seq reanalysis workflow focused on publicly available sequencing datasets from Phthorimaea absoluta.

This bioinformatics pipeline for exploratory transcriptomic analyses aiming to identify patterns of gene expression potentially associated with insecticide response, detoxification mechanisms, and resistance-related pathways.

The workflow includes:

Download of public SRA datasets
Quality control and preprocessing
Adapter trimming and filtering
Transcriptome indexing
Transcript quantification using Salmon
Preparation for downstream differential expression analysis

The pipeline prioritizes:

Reproducibility
Lightweight computational requirements
Modular execution
Compatibility with Linux environments and Conda

### 🎯 Biological Questions

This reanalysis pipeline can support questions such as:

Which transcripts are highly expressed after insecticide exposure?
Are detoxification-related genes enriched among expressed transcripts?
Do RNA-seq datasets show evidence of strong transcriptional remodeling?
Are cytochrome P450s, GSTs, UGTs, ABC transporters, or esterases among the most abundant transcripts?
How reproducible are publicly available transcriptomic datasets?
Does sequencing quality influence downstream transcript quantification?

Although the original datasets were previously published, this project focuses on:

independent computational validation;
reproducible reanalysis;
exploratory biological interpretation;
methodological transparency.

### 🧪 Pipeline Structure

### 1️⃣ Download of Public RNA-seq Data

Raw sequencing data are retrieved from NCBI SRA using the SRA Toolkit.

Example
prefetch SRR33779430

fasterq-dump SRR33779430 \
  --split-files \
  -O raw_data \
  -e 6

Generated files:

raw_data/
├── SRR33779430_1.fastq.gz
└── SRR33779430_2.fastq.gz

### 2️⃣ Quality Control

Initial sequencing quality is evaluated using:

FastQC
MultiQC
Main Metrics Evaluated
Per-base quality
Adapter contamination
GC content distribution
Read length distribution
Overrepresented sequences
Poly-G artifacts
Duplication levels
Example
fastqc raw_data/*.fastq.gz

multiqc .

### 3️⃣ Read Trimming

Reads are processed using Trim Galore.

Example
trim_galore \
  --paired \
  raw_data/SRR33779430_1.fastq.gz \
  raw_data/SRR33779430_2.fastq.gz \
  -o raw_data/trimmed
Outputs
raw_data/trimmed/
├── SRR33779430_1_val_1.fq.gz
├── SRR33779430_2_val_2.fq.gz
├── trimming_report.txt
└── FastQC reports

### 4️⃣ Reference Transcriptome Preparation

The transcript reference was obtained from the NCBI genome assembly:

Assembly

GCA_027580185.1_ASM2758018v1

Downloaded Files
genomic FASTA
transcript FASTA
GTF annotation
Main Reference Used
GCA_027580185.1_ASM2758018v1_rna_from_genomic.fna.gz

### 5️⃣ Salmon Index Construction

Transcript quantification is performed using Salmon pseudoalignment.

Index Construction
salmon index \
-t tuta_transcripts.fna \
-i salmon_index \
-k 31

### 6️⃣ Transcript Quantification

Example
salmon quant \
-i reference/salmon_index \
-l A \
-1 raw_data/trimmed/SRR33779430_1_val_1.fq.gz \
-2 raw_data/trimmed/SRR33779430_2_val_2.fq.gz \
-p 6 \
--validateMappings \
-o quant/SRR33779430

### 📊 Outputs Generated

Each sample generates:

quant/
└── SRR33779430/
    ├── quant.sf
    ├── cmd_info.json
    ├── lib_format_counts.json
    └── logs/
📄 Main Output File
quant.sf

Contains transcript-level abundance estimates.

Columns
Column	Description
Name	Transcript ID
Length	Transcript length
EffectiveLength	Effective transcript length
TPM	Transcripts per million
NumReads	Estimated read counts

### 📈 Quality Metrics Interpreted During Analysis

FastQC / MultiQC
Adapter contamination

Moderate adapter presence (~35%) was detected in some datasets, justifying trimming.

Poly-G tails

Low-frequency poly-G artifacts consistent with two-color Illumina chemistry were detected.

Sequence quality

Overall Phred quality scores remained high across reads.

### 📌 Salmon Quantification Metrics

Typical mapping rates observed:

~45–50%

This is expected because:

datasets were not generated specifically for this transcriptome assembly;
transcriptome completeness may be limited;
pseudoalignment is highly stringent when using --validateMappings.
📂 Repository Structure
tuta_rnaseq/
├── raw_data/
├── quant/
├── reference/
├── scripts/
├── results/
├── environment.yml
├── rodar_amostras.sh
├── pegar_metadados.py
└── README.md
⚙️ Software Used
Software	Version
FastQC	latest
MultiQC	latest
Trim Galore	2.2.0
Salmon	1.10.3
SRA Toolkit	3.4.1

### 🧬 Planned Downstream Analyses

Future analyses may include:

Differential expression analysis
Functional annotation
GO enrichment
KEGG pathway analysis
Detoxification gene screening
Comparative transcriptomics between treatments

### 🚀 Environment Installation

Conda Environment
conda env create -f environment.yml
conda activate rnaseq

### 📚 How to Cite This Repository

Quintela, V. (2026). RNA-seq Reanalysis Pipeline in Phthorimaea absoluta: A lightweight workflow for exploratory transcriptomic analysis using public datasets. GitHub repository.

### 📚 Scientific References
