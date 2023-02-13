import scipy.stats
import numpy as np
from util import findall, custom_pearsonr
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm

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


def compute_k_fold_fit(study, sample, family, stride='turner'):
    """Compute the K-fold fit for a given sample and family
    
    Parameters
    ----------
    study : Study
        Study object
    sample : str
        Sample name
    family : str
        Family name
    stride : str, optional
        Stride to use for the fit, by default 'turner'. Can be 'turner' or 'child#'.
        
    """
    
    data = study.get_df(sample=sample, family=family, section='ROI') #TODO remove unpaired bases

    data['deltaG'] = data['deltaG'].apply(lambda x: 0 if x == 'void' else float(x))

    assert len(data)>0, 'No data for sample {} and family {}'.format(sample, family)

    # turn it into a dataframe
    df = pd.DataFrame(
        columns= [base + str(idx+1) for base, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])],
        data = [int(offset)*[np.nan] + list(mr) for offset, mr in zip(data['frame_shift_ROI'], data['mut_rates'])],
        index= data['deltaG'].values
    )

    # Only keep the paired bases
    paired = [c in ['(',')'] for c in data['structure'].iloc[0]]
    df = df.loc[:,paired]
    
    # only keep the A and C bases
    idx_AC = [col[0] in ['A','C'] for col in df.columns]
    df = df.loc[:,idx_AC]

    # remove the bases that do not have a value for deltaG == 0.0
    try:
        df = df.loc[:,df.loc[0.0].notna().sum() > 0]
    except KeyError:
        print('No data for deltaG == 0.0 for sample {} and family {}'.format(sample, family))
        return pd.DataFrame({'Kfold':[]})
    
    
    # Change the index to be linear if needed
    if stride == 'child#':    
        df = df.reset_index().rename(columns={'index':'deltaG'})
        for idx, (dG, row) in enumerate(df.groupby('deltaG')):
            df.loc[row.index, 'child#'] = idx
        
        df.set_index('child#', inplace=True)

    # Function to fit
    def sigmoid(x, a, b, c):
        RT = 1.987204258*310/1000
        return a / (1 + b*np.exp(-x/RT)) + c
    
    # Reverse sigmoid function
    def rev_sigmoid(y, a, b, c):
        RT = 1.987204258*310/1000
        return -RT*np.log((a-y+c)/((y-c)*b))

    # Output values 
    base_Kfold = {}

    for base, mut_rate in df.iteritems():
        
        if base == 'deltaG':
            continue
        
        x_data = df.index[~np.isnan(mut_rate)].values
        mut_rate = mut_rate[~np.isnan(mut_rate)].values

        if len(mut_rate) >= 3: # at least 3 points to fit the sigmoid
            
            # Fit the sigmoid
            popt, pcov = curve_fit(sigmoid, x_data, mut_rate, p0=[0.04, 0.02, 0.00], bounds=([0, 0, 0], [0.1, np.inf, 0.05]), max_nfev=1000)
        
            # Compute the sigmoid midpoint 
            LARGE_VALUE = 100   
            midpoint_y = np.mean([sigmoid(LARGE_VALUE, *popt), sigmoid(-LARGE_VALUE, *popt)])  
            midpoint_y = min(max(midpoint_y, min(mut_rate)), max(mut_rate))
            
            midpoint_x = max(min(max(df.index),rev_sigmoid(midpoint_y, *popt)), min(df.index))
            
            # Store the results       
            base_Kfold[base] = {'avg_mr': midpoint_y, 'Kfold': midpoint_x}
            

    df = pd.DataFrame(base_Kfold).T

    # make a gaussian fit to the data
    df['norm'] = norm.pdf(df['Kfold'], np.mean(df['Kfold']), np.std(df['Kfold']))
    
    return df
    