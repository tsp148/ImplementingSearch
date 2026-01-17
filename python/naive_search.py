import iv2py
import time
import argparse

def findOccurences(text, pattern):
	"""Finde alle Vorkommen von pattern in text und gebe die Startindizes zurück."""
	occurences = []
	n, m = len(text), len(pattern)
	for i in range(n - m + 1):
		if text[i:i+m] == pattern:
			occurences.append(i)
	return occurences

# def read_fasta_sequences(filepath):
	"""Liest alle Sequenzen aus einer FASTA-Datei als Liste von Strings ein."""



parser = argparse.ArgumentParser(description="Naive search for pattern occurrences in reference sequences.")
parser.add_argument('--reference', type=str, required=True, help='Reference file (FASTA)')
parser.add_argument('--query', type=str, required=True, help='Query file (FASTA)')
parser.add_argument('--query_ct', type=int, default=100, help='Number of queries (will be duplicated if needed)')

args = parser.parse_args()

start_time = time.time()

reference = ""
for record in iv2py.fasta.reader(file=args.reference): #"../data/hg38_partial.fasta.gz"
    reference = str(record.seq)

queries = []
for record in iv2py.fasta.reader(file=args.query):
    queries.append(str(record.seq))

# Extend queries if needed
while len(queries) < args.query_ct:
    old_count = len(queries)
    queries.extend(queries[:old_count])
queries = queries[:args.query_ct]


# Search for all queries the reference
# occurences = []
# for ref in reference:

for query in queries:
    occurences = []
    occurences = findOccurences(reference, query)
    # print(f"Query: \n", query, f"\nfound at positions: {occurences}")
    # print("ref:", reference[occurences[0]:occurences[0]+len(query)] if occurences else "No occurrence found.")

elapsed = time.time() - start_time
print(f"Elapsed time: {int(elapsed)} s")
