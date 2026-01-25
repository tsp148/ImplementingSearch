import iv2py
import argparse
import time

parser = argparse.ArgumentParser(description="Naive search for pattern occurrences in reference sequences.")
parser.add_argument('--index', type=str, required=True, help='Index file (created with fmindex_construct.py)')
parser.add_argument('--query', type=str, required=True, help='Query file (FASTA)')
parser.add_argument('--query_ct', type=int, default=100, help='Number of queries (will be duplicated if needed)')

args = parser.parse_args()

begin_time = time.time()

queries = []
for record in iv2py.fasta.reader(file=args.query):
    queries.append(str(record.seq))

# Extend queries if needed
while len(queries) < args.query_ct:
    old_count = len(queries)
    queries.extend(queries[:old_count])
queries = queries[:args.query_ct]

index = iv2py.fmindex(args.index)

for query in queries:
    results = index.search(query)
    # Uncomment to see the results and check correctness
    # for result in results:
    #     print(f"Found at: {result}")
elapsed = time.time() - begin_time
elapsed_ms = int(elapsed * 1000)
print(f"Elapsed time: {int(elapsed)} s")
print(f"Elapsed time: {elapsed_ms} ms")