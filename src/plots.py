#%load_ext autoreload
#%autoreload 2
import sys
sys.path.append('../../src')
sys.path.append('/Users/ymdt/src/dreem')
from util import *
from config import *
import matplotlib.pyplot as plt
import numpy as np


## Barcodes 
# ---------------------------------------------------------------------------
def mutations_in_barcodes(study, sample):
    bc_reads = study.get_df(sample=sample, section='barcode')['num_of_mutations'].reset_index(drop=True)
    all_mutations = []
    for read in bc_reads:
        all_mutations += read
    
    plt.hist(all_mutations, rwidth=0.3, bins=np.arange(-0.5, 0.5+max(all_mutations), 1))
    plt.yscale('log')
    plt.xticks(np.arange(0, max(all_mutations)))
    plt.grid()
    plt.xlabel('Number of mutations on the barcode')
    plt.ylabel('Count')
    plt.title('Mutations in barcodes - {}'.format(sample))

def num_aligned_reads_per_construct_frequency_distribution(study, sample):
    num_aligned_reads = study.get_df(sample=sample, section='full')['num_aligned'].to_list()
    plt.hist(num_aligned_reads, bins=np.arange(500,max(num_aligned_reads), 1000), rwidth=0.3)
    plt.grid()
    plt.xlabel('Number of reads (binned)')
    plt.ylabel('Number of constructs')
    plt.title('Number of aligned reads per construct freq. distrib. - {}'.format(sample))

def num_aligned_reads_per_construct(study, sample):
    num_aligned_reads = study.get_df(sample=sample, section='full')['num_aligned'].to_list()
    num_aligned_reads.sort()
    num_aligned_reads.reverse()
    plt.bar(np.arange(1,1+len(num_aligned_reads)), num_aligned_reads)
    plt.xlabel('Number of constructs')
    plt.ylabel('Number of reads')
    plt.title('Number of aligned reads per construct - {}'.format(sample))
    plt.grid()
# ---------------------------------------------------------------------------

## DMS-MaPseq 
# ---------------------------------------------------------------------------
def mutations_per_read(study, sample):
    bc_reads = study.get_df(sample=sample, section='full')['num_of_mutations'].reset_index(drop=True)
    all_mutations = []
    for read in bc_reads:
        all_mutations += list(read)
    x = np.arange(0, max(all_mutations), 10)
    plt.hist(all_mutations, rwidth=0.9, bins=x)
    plt.yscale('log')
    plt.xticks(x)
    plt.grid()
    plt.xlabel('Number of mutations per read')
    plt.ylabel('Count')
    plt.title('Mutations per read - {}'.format(sample))

def mutation_identity_at_each_position(study, sample, construct):
    data = study.get_df(sample=sample, construct=construct, section='full').iloc[0]
    df = pd.DataFrame(index = list(data['sequence']))
    for base in ['A','C','G','T']:
        df[base] = np.array(data['mod_bases_'+base])/np.array(data['info_bases'])
    df.plot.bar(stacked=True, figsize= (20,4), color={'A':'r','C':'b','G':'y','T':'g'})
    plt.xticks(rotation=0)
    plt.xlabel('')
    
def mutation_fraction_at_each_position(study, sample, construct):
    data = study.get_df(sample=sample, construct=construct, section='full').iloc[0]
    df = pd.DataFrame(index = list(data['sequence']))
    for base in ['A','C','G','T']:
        df[base] = [mr if b==base  else np.nan for mr, b in zip(data['mut_rates'], data['sequence'])]
    df.plot.bar(stacked=True, figsize= (20,4), color={'A':'r','C':'b','G':'y','T':'g'})
    plt.xticks(rotation=0)

def read_coverage_per_position(study, sample, construct):
    data = study.get_df(sample=sample, construct=construct)
    sections, section_start, section_end = data['section'].unique(), data['section_start'].unique(), data['section_end'].unique()
    idx = np.argsort(section_start)
    sections, section_start, section_end = sections[idx], section_start[idx], section_end[idx]
    fig = plt.figure()
    ax = plt.gca()
    for s, ss, se in zip(sections, section_start, section_end):
        ax.plot(np.arange(ss-1, se), data[data['section']==s]['mut_rates'].values[0], label=s)
    plt.legend()