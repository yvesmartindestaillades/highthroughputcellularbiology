import numpy as np

version = 'v0.0'
generate_plots = True

LAH = 'CACAGTCGAAAGACTGTG'
MS2= 'GCATATGAGGATCACCCATATGC'

LAH_idx = np.array([idx for idx, c in enumerate(list(LAH)) if c in ['A','C']])
MS2_idx = np.array([idx for idx, c in enumerate(list(MS2)) if c in ['A','C']])

# define boundaries for the sections
boundary = {
    'MS2': lambda s: [19, 42],
    'TC1': lambda s: [42, 44],
    'ROI': lambda s: [44, s.index(LAH)-2],
    'TC2': lambda s: [s.index(LAH)-2, s.index(LAH)],
    'LAH': lambda s: [s.index(LAH), s.index(LAH)+len(LAH)],
    'buffer': lambda s: [s.index(LAH)+len(LAH), 139],
    'barcode': lambda s: [139, 151],
    'full': lambda s: [0, len(s)],
}