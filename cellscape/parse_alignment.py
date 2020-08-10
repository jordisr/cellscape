from Bio import pairwise2
from Bio.PDB import *
from Bio.SubsMat.MatrixInfo import blosum62
import numpy as np
import os, sys, re

def identity_from_alignment(a):
    s1 = np.array(list(a[0]))
    s2 = np.array(list(a[1]))
    return np.sum(s1 == s2) / len(np.where( s1 != '-')[0])

def overlap_from_alignment(a):
    s1 = np.array(list(a[0]))
    s2 = np.array(list(a[1]))
    s1_nogap = np.where( s1 != '-')
    s2_nogap = np.where( s2 != '-')
    s1_start_align = np.min(s1_nogap)
    s1_end_align = np.max(s1_nogap)
    s2_start_align = np.min(s2_nogap)
    s2_end_align = np.max(s2_nogap)
    overlap_align = (max(s1_start_align, s2_start_align), min(s1_end_align, s2_end_align))
    return(
    np.where(s1_nogap == overlap_align[0])[1][0],
    np.where(s1_nogap == overlap_align[1])[1][0],
    np.where(s2_nogap == overlap_align[0])[1][0],
    np.where(s2_nogap == overlap_align[1])[1][0])

def align_pair(s1, s2):
    # wrapper for biopython pairwise alignment
    return pairwise2.align.localds(s1, s2, blosum62, -3, -3, one_alignment_only=True)[0]

def align_all_pairs(s):
    for i in range(len(s)):
        for j in range(i+1, len(s)):
            s1 = s[i][1]
            s2 = s[j][1]
            alignments = pairwise2.align.localds(s1, s2, blosum62, -3, -3, one_alignment_only=True)
            print(s[i][0], len(s1), s[j][0], len(s2), *overlap_from_alignment(alignments[0]), identity_from_alignment(alignments[0]))

if __name__ == '__main__':
    a1 = (
    '---------AAAAABBBBBBB',
    'BBBBBBBBBAAAAA-------'
    )

    print(overlap_from_alignment(a1))

    a2 = (
    'BBBBBBBBBAAAAABBBBBBB',
    '---------AAAAA-------'
    )

    print(overlap_from_alignment(a2))
