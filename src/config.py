import numpy as np
import sys, os
import json
import pandas as pd
import matplotlib

sys.path.append('../../')

version = 'v0.2'
generate_plots = False
if generate_plots:
    matplotlib.use('agg')

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

bio_replicates_samples = [
['lauren470_S1','18_DMS_S7_L001'],
['lauren472_S3','19_DMS_S9_L001'],
['01_1_S22_reads','01_02_S23_reads'],
['05_1_S24_reads','05_2_S25_reads'],
['1_1_S26new_reads','1_2_S27new_reads'],
['2_1_s28new_reads','2_2_S29new_reads'],
['5_2_S31_reads','5_1_S30new_reads'],
['10_1_s32new_reads','10_2_s33new_reads']
 ] 

reaction_time_samples = [
    'Lauren_603_1min',
    'Lauren_603_3min',
    'Lauren_603_5min',
    'Lauren_603_10min',
]

dms_concentration_samples = [
 '01_1_S22_reads',
 '01_02_S23_reads',
 '05_1_S24_reads',
 '05_2_S25_reads'
 '1_1_S26new_reads',
 '1_2_S27new_reads',
 '2_1_s28new_reads',
 '2_2_S29new_reads',
 '5_1_S30new_reads'
 '5_2_S31_reads',
 '10_1_s32new_reads',
 '10_2_s33new_reads',
 ]
    
temperature_samples = [
    '65degrees_1_S20_L001',
    '5degrees_2_S9_L001', 
    '37degrees_01percent_2_S17_L001', 
    '25degrees_2_S13_L001', 
    '37degrees_1percent_2_S15_L001', 
    '37degrees_1percent_1_S14_L001', 
    '25degrees_1_S12_L001', 
    '37degrees_01percent_1_S16_L001', 
    '65degrees_2_S21_L001', 
    '5degrees_1_S8_L001', 
    '10degrees_2_S11_L001', 
    '45degrees_2_S19_L001', 
    '45degrees_1_S18_L001', 
    '10degrees_1_S10_L001'
    ]
# TODO

    
sys.path.append('/Users/ymdt/src/dreem')
import dreem

path_data = '../../data/'
saved_feather = path_data+'df.feather'

