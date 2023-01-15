import scipy.stats
import numpy as np
from util import findall, custom_pearsonr
import pandas as pd

def generate_barcode_replicates_pairs(study, sample):
    data = study.get_df(sample=sample, section = 'ROI')[['sequence', 'construct']]
    replicates = {}
    for _, g in data.groupby('sequence'):
        if len(g) > 1:
            for _, row in g.iterrows():
                replicates[row['construct']] = []
                for _, row2 in g.iterrows():
                    if row2['construct'] != row['construct']:
                        replicates[row['construct']].append(row2['construct'])
    return replicates


def compute_pearson_scores(study, sample, replicates_lists, sections):
    scores = []
    data = study.get_df(sample=sample, section = sections, base_type=['A','C'])[['sequence','construct','mut_rates']]

    # make sure that each construct has the same number of sections
    for c, g in data.groupby('construct'):
        if len(g) != len(sections) and c in replicates_lists:
            replicates_lists.pop(c)
    
    df = data.groupby('construct').agg({'mut_rates': lambda x: x}).reset_index()
    df['mut_rates'] = df['mut_rates'].apply(lambda x:np.concatenate(x))
    
    for construct, replicates in replicates_lists.items():
        scores_construct = []
        x = df[df['construct'] == construct]['mut_rates'].iloc[0]
        for replicate in replicates:
            if replicate not in replicates_lists.keys():
                continue
            y = df[df['construct'] == replicate]['mut_rates'].iloc[0]
            scores_construct.append(custom_pearsonr(x, y))
        scores.append(np.mean(scores_construct))

    sorted_idx = np.argsort(scores)[::-1]
    # return dict with construct names and scores as two sorted lists
    return {'constructs': [list(replicates_lists.keys())[i] for i in sorted_idx], 'scores': [scores[i] for i in sorted_idx]}


def find_frame_shift_ROI(study):
    
    if 'frame_shift_ROI' in study.df.columns:
        study.df.drop('frame_shift_ROI', axis=1, inplace=True)
            
    df = study.df[(study.df['section'] == 'ROI')][['sample','sequence', 'construct', 'family','section']].reset_index(drop=True)

    for _, g in df.groupby(['sample','family']):
        g.sort_values('sequence', key=lambda x: x.str.len(), inplace=True, ascending=False)
        reference = g['sequence'].iloc[0]
        for idx, row in g.iterrows():
            # assert sequence is unique in reference
            subsequence_matches = list(findall(row['sequence'], reference))
            most_likely_match = [-abs(len(reference)/2 - (m+len(row['sequence'])/2)) for m in subsequence_matches]
            assert len(most_likely_match) > 0, 'Sequence {} not found in reference {}'.format(row['sequence'], reference)
            df.loc[idx, 'frame_shift_ROI'] = subsequence_matches[most_likely_match.index(max(most_likely_match))]

    for _, g in df.groupby('family'):
        g.sort_values('sequence', key=lambda x: x.str.len(), inplace=True, ascending=False)
        assert g['frame_shift_ROI'].iloc[0] == 0, 'Frame shift is not 0 for reference sequence'

    df = study.df.merge(df[['sample','construct','section','sequence','frame_shift_ROI']], on=['sample','construct','section','sequence'], how='left')
    return df
