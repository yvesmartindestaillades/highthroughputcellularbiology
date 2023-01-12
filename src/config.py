import numpy as np
import sys, os
import json
import pandas as pd
import matplotlib
sys.path.append('../../')

matplotlib.use('agg')

version = 'v0.1'
generate_plots = False
min_base_coverage = 1000

LAH = 'CACAGTCGAAAGACTGTG'
MS2= 'GCATATGAGGATCACCCATATGC'

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
sys.path.append('/Users/ymdt/src/dreem')
import dreem

path_data = '../../data/dreem_output/'
saved_feather = path_data+'df.feather'

if not os.path.exists(saved_feather):
    study = dreem.draw.study.Study(
        data = [json.load(open(path_data + f, 'r')) for f in os.listdir(path_data) if f.endswith('.json')]
    )
    study.df['deltaG'] = study.df['deltaG'].apply(lambda x: 0 if x == 'void' else x)
    study.df.to_feather(saved_feather)
else:
    study = dreem.draw.study.Study()
    study.df = pd.read_feather(saved_feather)

study.df = study.df[study.df['worst_cov_bases'] > min_base_coverage]

