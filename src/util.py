LAH = 'CACAGTCGAAAGACTGTG'
MS2= 'GCATATGAGGATCACCCATATGC'

import numpy as np
LAH_idx = np.array([idx for idx, c in enumerate(list(LAH)) if c in ['A','C']])
MS2_idx = np.array([idx for idx, c in enumerate(list(MS2)) if c in ['A','C']])

def fasta_to_df(fasta_file):
    """Converts a fasta file to a pandas dataframe"""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines]
    names = [l[1:] for l in lines if l.startswith('>')]
    seqs = [l for l in lines if not l.startswith('>')]
    return pd.DataFrame({'construct': names, 'sequence': seqs})