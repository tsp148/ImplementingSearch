# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 15:54:03 2026

@author: camil
"""

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


#print(suffix_array_find_all("banana$",  [6, 5, 3, 1, 0, 4, 2], "ana"))


