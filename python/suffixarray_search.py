# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 15:54:03 2026

@author: camil
"""
import iv2py
import argparse
import time

parser = argparse.ArgumentParser(description="Naive search for pattern occurrences in reference sequences.")
parser.add_argument('--reference', type=str, required=True, help='Reference file (FASTA)')
parser.add_argument('--query', type=str, required=True, help='Query file (FASTA)')
parser.add_argument('--query_ct', type=int, default=100, help='Number of queries (will be duplicated if needed)')

args = parser.parse_args()


def suffix_array_find_range(s: str, sa: list[int], pattern: str) -> tuple[int, int]:
    n = len(sa)
    m = len(pattern)

    def lower_bound(p: str) -> int:
        low, high = 0, n
        while low < high:
            mid = (low + high) // 2
            start = sa[mid]
            pref = s[start:start + m]
            if pref < p:
                low = mid + 1
            else:
                high = mid
        return low

    L = lower_bound(pattern)


    low, high = L, n
    while low < high:
        mid = (low + high) // 2
        start = sa[mid]
        pref = s[start:start + m]
        if pref <= pattern:
            low = mid + 1
        else:
            high = mid
    R = low

 
    if L == n or s[sa[L]:sa[L] + m] != pattern:
        return (0, 0)

  
    while R > L and s[sa[R-1]:sa[R-1] + m] != pattern:
        R -= 1

    return (L, R)

def suffix_array_find_all(s: str, sa: list[int], pattern: str) -> list[int]:
    L, R = suffix_array_find_range(s, sa, pattern)
    return [sa[i] for i in range(L, R)]

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

sa = iv2py.create_suffixarray(reference)
elapsed = time.time() - start_time
print(f"Elapsed time for SA construction: {int(elapsed)} s")
start_time = time.time()
for query in queries:
    occurences = []
    occurences = suffix_array_find_all(reference, sa, query)
    # print(f"Query: \n", query, f"\nfound at positions: {occurences}")
    # print("ref:", reference[occurences[0]:occurences[0]+len(query)] if occurences else "No occurrence found.")

# print(suffix_array_find_all("banana$",  [6, 5, 3, 1, 0, 4, 2], "ana"))


elapsed = time.time() - start_time
print(f"Elapsed time: {elapsed} s")


# /usr/bin/time -v python suffixarray_search.py --reference ../data/hg38_partial.fasta.gz --query ../data/illumina_reads_100.fasta.gz --query_ct 1