import pandas as pd
from Bio import SeqIO

# DEG table
deg = pd.read_csv("results/spinosad_sig_DEGs.csv")

# transcript IDs
ids = set(deg.iloc[:,0])

# transcriptome fasta
fasta_in = "reference/tuta_transcripts.fna"
fasta_out = "results/spinosad_sig_DEGs.fasta"

count = 0

with open(fasta_out, "w") as out_handle:
    for record in SeqIO.parse(fasta_in, "fasta"):
        if record.id in ids:
            SeqIO.write(record, out_handle, "fasta")
            count += 1

print("Sequences extracted:", count)
