import iv2py
import argparse
import time

parser = argparse.ArgumentParser(description="Naive search for pattern occurrences in reference sequences.")
parser.add_argument('--reference', type=str, required=True, help='Reference file (FASTA)')

args = parser.parse_args()

reference = ""
for record in iv2py.fasta.reader(file=args.reference): #"../data/hg38_partial.fasta.gz"
    reference = str(record.seq)

index = iv2py.fmindex(reference)
index.save("py_index.idx")