import sys
sys.path.append('../../src')
sys.path.append('/Users/ymdt/src/dreem') # remove that
from util import *
from config import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import sklearn.metrics
import generate_dataset
import seaborn as sns
# pearsonr = scipy.stats.pearsonr
import scipy.stats

## b / i - Demultiplexing 
# ---------------------------------------------------------------------------
def mutations_in_barcodes(study, sample):
    bc_reads = study.get_df(sample=sample, section='barcode')['num_of_mutations'].reset_index(drop=True).values
    all_mutations = []
    for read in bc_reads:
        all_mutations += list(read)
    
    plt.hist(all_mutations, rwidth=0.9, bins=np.arange(-0.5, 0.5+max(all_mutations), 1))
    plt.yscale('log')
    plt.xticks(np.arange(0, max(all_mutations)))
    plt.grid()
    plt.xlabel('Number of mutations on the barcode')
    plt.ylabel('Count')
    plt.title('Mutations in barcodes - {}'.format(sample))

def num_aligned_reads_per_construct_frequency_distribution(study, sample):
    num_aligned_reads = study.get_df(sample=sample, section='full')['num_aligned'].to_list()
    plt.hist(num_aligned_reads, bins=np.arange(500,max(num_aligned_reads), 1000), rwidth=0.9)
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

## b / ii - DMS-MaPseq
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
    plt.ylabel('Mutation fraction')
    plt.xlabel('Position')
    plt.title('Mutation identity at each position - {} - {}'.format(sample, construct))
    
def mutation_fraction_at_each_position(study, sample, construct):
    data = study.get_df(sample=sample, construct=construct, section='full').iloc[0]
    df = pd.DataFrame(index = list(data['sequence']))
    for base in ['A','C','G','T']:
        df[base] = [mr if b==base  else np.nan for mr, b in zip(data['mut_rates'], data['sequence'])]
    df.plot.bar(stacked=True, figsize= (20,4), color={'A':'r','C':'b','G':'y','T':'g'})
    plt.xticks(rotation=0)
    plt.ylabel('Mutation fraction')
    plt.xlabel('Position')
    plt.title('Mutation fraction at each position - {} - {}'.format(sample, construct))

def read_coverage_per_position(study, sample, construct):
    data = study.get_df(sample=sample, construct=construct)
    sections, section_start, section_end = data['section'].unique(), data['section_start'].unique(), data['section_end'].unique()
    idx = np.argsort(section_start)
    sections, section_start, section_end = sections[idx], section_start[idx], section_end[idx]
    fig = plt.figure()
    ax = plt.gca()
    for s, ss, se in zip(sections, section_start, section_end):
        ax.plot(np.arange(ss-1, se), data[data['section']==s]['cov_bases'].values[0], label=s)
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Read coverage')
    plt.title('Read coverage per position - {} - {}'.format(sample, construct))
    
# ---------------------------------------------------------------------------

## b / iv - Error estimation
# ---------------------------------------------------------------------------
def distribution_of_the_mean(means, residues_set, bounds):
    fig = plt.figure(figsize=(100, 70))
    plt.subplots(sharex=True, sharey=True)
    plt.suptitle('Distribution of the mean per residue')
    ax = plt.gca()
    for i, n in enumerate(residues_set):
        idx_means = n-min(residues_set)
        plt.subplot(3,3,i+1)
        plt.xlim(-0.001,0.015)
        plt.ylim(0, len(means[idx_means]))
        plt.hist(means[idx_means], bins=np.arange(-0.0001,0.015,0.0005))
        plt.plot([bounds[0][idx_means]]*2, (0, len(means[idx_means])), 'r', linewidth=1, label='2.5%:{}'.format(bounds[0][idx_means]))
        plt.plot([bounds[1][idx_means]]*2, (0, len(means[idx_means])), 'r', linewidth=1, label='97.5%:{}'.format(bounds[1][idx_means]))
        plt.title('Residue {}'.format(n+1))
        plt.legend()
        plt.tight_layout()
    plt.tight_layout()
        
# ---------------------------------------------------------------------------

## c / iii - Barcode
# ---------------------------------------------------------------------------
SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES = ['MS2','TC1','ROI','TC2','LAH','buffer']

def barcode_comparison_scatter_plot(study, sample, construct, replicate):
    plt.figure(figsize=(5,6))
    
    # plot scatter plot
    X, Y = [], []
    for section in SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES:
        x = study.get_df(sample = sample, construct = construct, section = section, base_type = ['A','C'])['mut_rates'].iloc[0]
        y = study.get_df(sample = sample, construct = replicate, section = section, base_type = ['A','C'])['mut_rates'].iloc[0]
        plt.plot(x, y, 'x', label = section)
        X, Y = X + list(x), Y + list(y)
        
    # reshape plot
    plt.axis('square')
    (xmin, xmax), (ymin, ymax) = plt.xlim(), plt.ylim()
    plt.xlim(min(xmin, ymin), max(xmax, ymax))
    xmin, xmax = plt.xlim()
    
    # plot y=x line
    plt.plot([xmin, xmax], [xmin, xmax], 'k--', label='y=x')

    # plot regression line
    x, y = np.array(X), np.array(Y)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
    plt.plot([xmin, xmax], [xmin*slope + intercept, xmax*slope + intercept], 'r-', label='Linear regression')
    
    # plot labels
    plt.title('Barcodes comparison - {}'.format(sample))
    plt.text(0.5, 0.9, 'Pearson correlation: {:.2f}'.format(custom_pearsonr(x,y)), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
    plt.text(0.5, 0.85, 'R2: {:.2f}'.format(r_value), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
    plt.xlabel('Barcode: %s' % construct)
    plt.ylabel('Barcode: %s' % replicate)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    
def combined_barcode_replicates_in_sample(study, sample):
    replicates_lists = generate_dataset.generate_barcode_replicates_pairs(study, sample)       
    pearson_scores = generate_dataset.compute_pearson_scores(study, sample, replicates_lists, SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES)
    plt.bar(np.arange(len(pearson_scores)), pearson_scores)
    plt.ylabel('Pearson correlation score average across barcode replicates')
    plt.xlabel('Constructs (barcode vs all replicates)')
    plt.title('Combined barcodes replicates - {}'.format(sample))

# ---------------------------------------------------------------------------

## c / iv - Reproducibility
# ---------------------------------------------------------------------------
def barcode_replicates(study, sample):
    replicates_lists = generate_dataset.generate_barcode_replicates_pairs(study, sample)   
    pearson_scores = generate_dataset.compute_pearson_scores(study, sample, replicates_lists, SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES)
    plt.hist(pearson_scores, rwidth=0.9)
    plt.xlabel('Pearson correlation score average across barcode replicates')
    plt.ylabel('Number of constructs')
    plt.title('Barcodes replicates - {}'.format(sample))
    
# ---------------------------------------------------------------------------

## c / v - Proportional
# ---------------------------------------------------------------------------
def change_in_temp_mut_frac_vs_temperature(study, samples, construct):
    # get data
    data = study.get_df(
        sample = samples,
        construct=construct,
        section='ROI',
        base_type = ['A','C']) 

    data = data[['sequence','mut_rates','index_selected','structure','temperature_k']]
    
    if len(data) == 0:
        print('No data: {}, {}'.format(samples, construct))
        return None
    
    # turn it into a dataframe
    df = pd.DataFrame(
        columns= [base + str(idx+1) for base, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])],
        data = np.array(data['mut_rates'].tolist()),
        index= data['temperature_k'].values
    )
    
    # plot
    paired = [False if residue == '.' else True for residue in data['structure'].iloc[0]]
    df.plot(style= ['-x' if paired[idx] else '-o' for idx in range(len(paired))])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(construct)
    plt.xlabel('Time (secs)')

def mut_rate_across_family_vs_deltaG(study, sample, family):
    
    # get a neat dataframe with the mutation rates for each base at each deltaG
    study.df['frame_shift_ROI'] = generate_dataset.find_frame_shift_ROI(study)
    data = study.get_df(sample=sample, family=family, section='ROI')
    
    data['deltaG'] = data['deltaG'].apply(lambda x: 0 if x == 'void' else x)
    
    # turn it into a dataframe
    df = pd.DataFrame(
        columns= [base + str(idx+1) for base, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])],
        data = [int(offset)*[np.nan] + list(mr) for offset, mr in zip(data['frame_shift_ROI'], data['mut_rates'])],
        index= data['deltaG'].values
    )
    
    # only keep the A and C bases
    idx_AC = [col[0] in ['A','C'] for col in df.columns]
    df = df.loc[:,idx_AC]
    
    # differentiate between paired and unpaired bases
    paired = [residue for idx, residue in enumerate([False if residue == '.' else True for residue in data['structure'].iloc[0]]) if idx_AC[idx]]
    
    # plot
    df.plot(style= ['-x' if paired[idx] else '-o' for idx in range(len(paired))])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title('Changement in temperature (w/ prediction) for ' + family + ' in ' + sample)
    plt.xlabel('ΔG (kcal/mol)')
    plt.ylabel('Mutation rate')
    
# ---------------------------------------------------------------------------

## d / kfold
# ---------------------------------------------------------------------------

def heatmap_across_family_members(study, sample, family):
    study.df['frame_shift_ROI'] = generate_dataset.find_frame_shift_ROI(study)
    data = study.get_df(
        sample=sample,
        family=family,
        section = 'ROI',
        )
    
    sequences_ROI = data[data['section'] == 'ROI'].sort_values('sequence', key=lambda x: x.str.len(), ascending=False)['sequence'].values
    reference = sequences_ROI[0]

    df = pd.DataFrame(
        columns = [base + str(idx + 1) for base, idx in zip(reference, range(len(reference)))],
        data = [int(offset)*[np.nan] + list(mr) for offset, mr in zip(data['frame_shift_ROI'], data['mut_rates'])],
        index = data['construct'].values 
    )
    
    plt.figure(figsize=(data.shape[1]//2, data.shape[0]//2))
    plt.gca().xaxis.tick_top()
    sns.heatmap(df, annot=True, fmt='.2f')
    plt.title('Heatmap across family members with all bases - sample {} - family {}'.format(sample, family))
    plt.ylabel('Construct')

