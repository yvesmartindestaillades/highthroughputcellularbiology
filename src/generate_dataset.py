import scipy.stats
import numpy as np
from util import findall, custom_pearsonr

def generate_barcode_replicates_pairs(study, sample):
    data = study.get_df(sample=sample, section = 'ROI')
    replicates = {}
    for construct in study.get_constructs(sample):
        sequence = study.get_df(sample=sample, construct=construct, section='ROI')['sequence'].values[0]
        replicates[construct] = []
        for _, row in data.iterrows():
            if row['sequence'] == sequence and row['construct'] != construct:
                replicates[construct].append(row['construct'])
        if len(replicates[construct]) == 0:
            del replicates[construct]
    return replicates


def compute_pearson_scores(study, sample, replicates_lists, sections):
    scores = []
    for construct, replicates in replicates_lists.items():
        scores_construct = []
        for replicate in replicates:
            x, y = [], []
            for section in sections:
                x += study.get_df(sample=sample, construct=construct, section=section, base_type=['A','C'])['mut_rates'].iloc[0].tolist()
                y += study.get_df(sample=sample, construct=replicate, section=section, base_type=['A','C'])['mut_rates'].iloc[0].tolist()
            scores_construct.append(custom_pearsonr(x, y))
        scores.append(np.mean(scores_construct))
    scores.sort(reverse=True)
    return scores


def find_frame_shift_ROI(study):
    df = study.get_df(section = 'ROI')
    for _, g in df.groupby('family'):
        g.sort_values('sequence', key=lambda x: x.str.len(), inplace=True, ascending=False)
        reference = g['sequence'].iloc[0]
        for idx, row in g.iterrows():
            # assert sequence is unique in reference
            subsequence_matches = list(findall(row['sequence'], reference))
            most_likely_match = [-abs(len(reference)/2 - (m+len(row['sequence'])/2)) for m in subsequence_matches]
            assert len(most_likely_match) > 0, 'Sequence {} not found in reference {}'.format(row['sequence'], reference)
            df.loc[idx, 'frame_shift_ROI'] = subsequence_matches[most_likely_match.index(max(most_likely_match))]
            
    return df['frame_shift_ROI'].astype(int)
