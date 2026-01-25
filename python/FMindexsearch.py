# -*- coding: utf-8 -*-
"""
Created on Sat Jan 24 13:27:26 2026

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



# This funtion build the L column of the lexicographically index array  
def BWT(text, SA): 
    lenSA = len(SA)
    btw = [""]*lenSA
    for i in range(lenSA):
        if SA[i] > 0:
            btw[i] = text[SA[i]-1]
        if SA[i] == 0:
            btw[i] = "$"
    return btw


   # Build the cumlative counts funtion         
def C(text):
    frequency = {}
    for ch in text:
        frequency[ch] = frequency.get(ch, 0) + 1
                
    
    C = {}
    count = 0
    for ch in sorted(frequency):
        C[ch] = count
        count += frequency[ch]
        
    return C

# build the funtion Occ which indicate the indexes of a character in L
def BuildOcc(L):
    Occ = {}
    n = len(L)

    # initialize
    for c in set(L):
        Occ[c] = [0] * (n + 1)

    # build prefix sums
    for i in range(n):
        for c in Occ:
            Occ[c][i+1] = Occ[c][i]
        Occ[L[i]][i+1] += 1

    return Occ

   
            
def Occ(Occ_table, c, i):
    # number of c in L[0..i]
    return Occ_table[c][i+1]

    
def LF(L, C, charIndexes, i):
    # i is an index into L (0-based)
    c = L[i]
    return C[c] + Occ(charIndexes, c, i) - 1
    
def backward_search(pattern, L, C, charIndexes):
    n = len(L)
    l, r = 0, n - 1

    # process pattern from right to left
    for ch in reversed(pattern):
        if ch not in C:
            return (-1, -1)

        # Occ(c, l-1) should be 0 when l == 0
        occ_lm1 = 0 if l == 0 else Occ(charIndexes, ch, l - 1)
        occ_r   = Occ(charIndexes, ch, r)

        l = C[ch] + occ_lm1
        r = C[ch] + occ_r - 1

        if l > r:
            return (-1, -1)

    return (l, r)


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

if not reference.endswith("$"):
    reference += "$"
sa = iv2py.create_suffixarray(reference)
elapsed = time.time() - start_time
print(f"Elapsed time for SA construction: {int(elapsed)} s")
start_time = time.time()

L = BWT(reference, sa)
C_table = C(reference)
Occ_table = BuildOcc(L)


for query in queries:
    l, r = backward_search(query, L, C_table, Occ_table)
    if l == -1:
        occurrences = []
    else:
        occurrences = sa[l:r+1]
    # optional: print or count
    # print(query, len(occurrences))



elapsed = time.time() - start_time
print(f"Elapsed time: {elapsed} s")

# /usr/bin/time -v python suffixarray_search.py --reference ../data/hg38_partial.fasta.gz --query ../data/illumina_reads_100.fasta.gz --query_ct 1