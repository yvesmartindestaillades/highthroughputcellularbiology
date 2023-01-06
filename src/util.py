
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import datetime
import config

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
