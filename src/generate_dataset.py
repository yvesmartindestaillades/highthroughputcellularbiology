import scipy.stats
import numpy as np
from util import findall, custom_pearsonr

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
    data = study.get_df(sample=sample, section = sections, base_type=['A','C'])[['construct','mut_rates']]
    df = data.groupby('construct').agg({'mut_rates': lambda x: x}).reset_index()
    df['mut_rates'] = df['mut_rates'].apply(lambda x:np.concatenate(x))
    
    for construct, replicates in replicates_lists.items():
        scores_construct = []
        x = df[df['construct'] == construct]['mut_rates'].iloc[0]
        for replicate in replicates:
            y = df[df['construct'] == replicate]['mut_rates'].iloc[0]
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
