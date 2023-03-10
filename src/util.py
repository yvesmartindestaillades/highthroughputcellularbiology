
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import datetime
from config import version, generate_plots
from scipy import stats
from sklearn.linear_model import LinearRegression

import mixem

def fasta_to_df(fasta_file):
    """Converts a fasta file to a pandas dataframe"""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines]
    names = [l[1:] for l in lines if l.startswith('>')]
    seqs = [l for l in lines if not l.startswith('>')]
    return pd.DataFrame({'reference': names, 'sequence': seqs})

def save_plotly_fig(cwd, filename, fig, format='html'):
    """Save a plotly figure and create the directory if it doesn't exists.

    Args:
        cwd: directory
        filename: name of the html file
        fig: plotly figure
        format: format of the file, either 'html' or 'png'. Default is 'html'
    """

    cwd = os.path.normpath(cwd)
    path = os.path.join('..','..', 'figs', cwd.split('/')[-2], cwd.split('/')[-1].split('-')[0], *(version+'-'+filename).split('/')[:-1])
    path = make_path(path)

    if format == 'png':
        fig.write_image(path+(version+'-'+filename).split('/')[-1]+'.png')
    elif format == 'html':
        fig.write_html(path+(version+'-'+filename).split('/')[-1]+'.html', full_html=True)
    else:
        raise ValueError('format must be either html or png')

def savefig(file:str, close=True)->None:
    """Save a matplotlib figure and create the directory if it doesn't exists.

    Args:
        file: path+title.
        facecolor: color of the background 
    """

    path = make_path('/'.join(file.split('/')[:-1]))
    plt.savefig(path+file.split('/')[-1], bbox_inches='tight')
    if close:
        # Clear the current axes.
        plt.cla() 
        # Clear the current figure.
        plt.clf() 
        # Closes all the figure windows.
        plt.close('all')   
        
def savefig2(cwd, filename):
    cwd = os.path.normpath(cwd)
    path = os.path.join('..','..', 'figs', cwd.split('/')[-2], cwd.split('/')[-1].split('-')[0], *(version+'-'+filename).split('/')[:-1])
    path = make_path(path)
    savefig(path+(version+'-'+filename).split('/')[-1]+'.png', close=generate_plots)

def make_path(path:str)->str:
    """Create directories until path exists on your computer. Turns the keyword 'date' into today's date.

    Args:
        path: series of directories that you want to create.
    
    Returns:
        Updated path with today's date instead of the keyword 'date'  
    """

    path = os.path.normpath(path)
    path=path.split(os.sep)
    try:
        path[path.index('date')] = str(datetime.datetime.now())[:10]
    except:
        'No date in path'
    full_path = ''
    for repo in path:
        full_path = full_path + f"{repo}/"
        if not os.path.exists(full_path):
            os.mkdir(full_path)
    return full_path

def compute_mutation_rates(df):
    """Compute the mutation rate for each column of a dataframe.
    mutation_rate = (# of mutations) / (# of reads). Ignore missing values.
    """
    mut_rates = []
    for col in df.columns:
        col = df[col]
        sub_hist = col[col.notna()].astype(int).sum()
        num_of_reads = col[col.notna()].astype(int).count()
        mut_rate = round(sub_hist / num_of_reads, 6)
        mut_rates.append(mut_rate)
    return mut_rates

def bitvector_to_df(path, *args, **kwargs):
    df = pd.read_csv(path, *args, **kwargs)['Bit_vector']
    df = df.str.split('', expand=True).drop(columns=[0, 171])
    values_map = {
        'A': 1,
        'C': 1,
        'G': 1,
        'T': 1,
        'N': 1,
        '1': 1,
        '0': 0,
        '?': np.nan,
        '.': np.nan
    }
    df = df.applymap(lambda x: values_map[x])
    return df

def findall(p, s):
    '''Yields all the positions of
    the pattern p in the string s.'''
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)

def custom_pearsonr(x, y, percentile=90, diff_filter=0.2):
    x, y = np.array(x), np.array(y)
    min_unique = min(len(np.unique(x)), len(np.unique(y)))
    if min_unique == 0:
        return np.nan
    if len(x) != len(y):
        raise ValueError('x and y must have the same length')
    x , y = x[~np.isnan(x) & ~np.isnan(y)], y[~np.isnan(x) & ~np.isnan(y)]
    if min_unique == 1:
        return np.nan
    if min_unique == 2:
        return stats.pearsonr(x, y)[0]
    else:
        for vx, vy in [[x, y], [y, x]]:
            model = LinearRegression()
            model.fit(vx.reshape(-1,1), vy)
            std_err_per_data_point = np.abs(vy - model.predict(vx.reshape(-1,1)))
            score_all = stats.pearsonr(vx, vy)[0]
            idx_keep = np.where(std_err_per_data_point < np.percentile(std_err_per_data_point, percentile))[0]
            min_unique = min(len(np.unique(vx.take(idx_keep))), len(np.unique(vy.take(idx_keep))))
            if min_unique >=2:
                score_filter_worst_data_point = stats.pearsonr(vx.take(idx_keep), vy.take(idx_keep))[0]
                if score_filter_worst_data_point - score_all > diff_filter:
                    return score_filter_worst_data_point    
        return score_all


class logGMM:

    def __init__(self):
        
        self.weights = [] 
        self.distributions = []
        self.data = []

    def fit_logGMM(self, data):
        data = np.array(data)
        self.data = data[data>0.0]

        # Fit and plot log gaussian mixture
        self.weights, self.distributions, log_likelihood = mixem.em(self.data, [
            mixem.distribution.LogNormalDistribution(mu=-6, sigma=2.4),
            mixem.distribution.LogNormalDistribution(mu=-2.7, sigma=0.1),
            # mixem.distribution.NormalDistribution(mu=300, sigma=50),
        ], initial_weights=[0.8, 0.2], progress_callback=None)

    def get_pdf(self, x_axis):
        return mixem.probability(x_axis, self.weights, self.distributions)

    def find_midpoint(self, interval=[1e-3, 1e-1]):
        x_search = np.linspace(interval[0], interval[1], 10000)

        pdf1 = np.exp(self.distributions[0].log_density(x_search))
        pdf2 = np.exp(self.distributions[1].log_density(x_search))

        x_max_1 = np.argmax(pdf1)
        x_max_2 = np.argmax(pdf2)

        if x_max_1 < x_max_2:
            interval = [x_max_1, x_max_2]
        else:
            interval = [x_max_2, x_max_1]

        pdf_diff = np.abs(pdf1-pdf2)[interval[0]:interval[1]]
        if pdf_diff.size >0:
            return x_search[np.argmin(pdf_diff)]
        else:
            return 0.0

    def get_mode(self, dist_idx):
        return np.exp(self.distributions[dist_idx].mu - self.distributions[dist_idx].sigma**2)


def compute_wilson_interval(p, n, z = 1.96):
    denominator = 1 + z**2/n
    centre_adjusted_probability = p + z*z / (2*n)
    adjusted_standard_deviation = np.sqrt((p*(1 - p) + z*z / (4*n)) / n)
    
    lower_bound = (centre_adjusted_probability - z*adjusted_standard_deviation) / denominator
    upper_bound = (centre_adjusted_probability + z*adjusted_standard_deviation) / denominator
    return (lower_bound, upper_bound)


def family_from_reference(reference):
    return reference.split('=')[1].split('-')[0]

def compute_affine_transformation(s1, s2):
    """Compute the affine transformation that aligns s1 to s2.
    Returns the y = mx + c transformation.
    """
    s1 = np.array(s1)
    s2 = np.array(s2)
    model = LinearRegression()
    model.fit(s1, s2)
    # ignore nan values
    return lambda x: x*model.coef_[0][0] + model.intercept_[0]
