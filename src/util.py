
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import datetime
import config

import mixem

def fasta_to_df(fasta_file):
    """Converts a fasta file to a pandas dataframe"""
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines]
    names = [l[1:] for l in lines if l.startswith('>')]
    seqs = [l for l in lines if not l.startswith('>')]
    return pd.DataFrame({'construct': names, 'sequence': seqs})

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
    path = os.path.join('..','..', 'figs', cwd.split('/')[-2], cwd.split('/')[-1].split('-')[0], *(config.version+'-'+filename).split('/')[:-1])
    path = make_path(path)
    savefig(path+(config.version+'-'+filename).split('/')[-1]+'.png', close=config.generate_plots)

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
