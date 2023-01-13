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
from sklearn.metrics import r2_score
import plotly.graph_objects as go
import plotly.express as px

## b / i - Demultiplexing 
# ---------------------------------------------------------------------------
def mutations_in_barcodes(study, sample):
    bc_reads = study.get_df(sample=sample, section='barcode')['num_of_mutations'].reset_index(drop=True).values
    all_mutations = []
    for read in bc_reads:
        all_mutations += list(read)

    hist, bin_edges = np.histogram(all_mutations, bins=np.arange(0, max(all_mutations)) )
    return go.Bar( x=bin_edges, y=hist, showlegend=False, marker_color='indianred')

def num_aligned_reads_per_construct_frequency_distribution(study, sample):
    num_aligned_reads = study.get_df(sample=sample, section='full')['num_aligned'].to_list()

    return go.Histogram(x=num_aligned_reads, showlegend=False, marker_color='indianred')

def num_aligned_reads_per_construct(study, sample):
    num_aligned_reads = study.get_df(sample=sample, section='full')['num_aligned'].to_list()
    num_aligned_reads.sort()
    num_aligned_reads.reverse()

    return go.Bar(x=np.arange(1,1+len(num_aligned_reads)), y=num_aligned_reads, showlegend=False, marker_color='indianred')
    
# ---------------------------------------------------------------------------

## b / ii - DMS-MaPseq
# ---------------------------------------------------------------------------
def mutations_per_read(study, sample):
    bc_reads = study.get_df(sample=sample, section='full')['num_of_mutations'].reset_index(drop=True)
    all_mutations = []
    for read in bc_reads:
        all_mutations += list(read)

    hist, bin_edges = np.histogram(all_mutations, bins=np.arange(0, max(all_mutations)) )
    return go.Bar( x=bin_edges, y=hist, showlegend=False, marker_color='indianred')

def mutation_identity_at_each_position(study, sample, construct):
    data = study.get_df(sample=sample, construct=construct, section='full').iloc[0]
    df = pd.DataFrame(index = list(data['sequence']))
    stacked_bar = []
    color_map={'A':'red','C':'blue','G':'yellow','T':'green'}
    for base in ['A','C','G','T']:
        df[base] = np.array(data['mod_bases_'+base])/np.array(data['info_bases'])
        stacked_bar.append( go.Bar(x=np.arange(len(data['sequence'])), y=list(df[base]), marker_color=color_map[base], showlegend=False) )

    return {'fig': stacked_bar, 'data': df}

 
def mutation_fraction_at_each_position(study, sample, construct):
    data = study.get_df(sample=sample, construct=construct, section='full').iloc[0]
    df = pd.DataFrame(index = list(data['sequence']))
    stacked_bar = []
    color_map={'A':'red','C':'blue','G':'yellow','T':'green'}
    for base in ['A','C','G','T']:
        df[base] = [mr if b==base  else np.nan for mr, b in zip(data['mut_rates'], data['sequence'])]
        stacked_bar.append( go.Bar(x=np.arange(len(data['sequence'])), y=list(df[base]), marker_color=color_map[base], showlegend=False) )
    
    return {'fig': stacked_bar, 'data': df}

def read_coverage_per_position(study, sample, construct):
    data = study.get_df(sample=sample, construct=construct)
    sections, section_start, section_end = data['section'].unique(), data['section_start'].unique(), data['section_end'].unique()
    idx = np.argsort(section_start)
    sections, section_start, section_end = sections[idx], section_start[idx], section_end[idx]

    scatters = []
    for s, ss, se in zip(sections, section_start, section_end):

        scatters.append(
            go.Scatter(
                x=np.arange(ss-1, se),
                y=data[data['section']==s]['cov_bases'].values[0],
                name=s, 
                mode='lines',
                showlegend=False) )
    data['section'] = sections
    return {'fig': scatters, 'data': data}
    
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
    """ generate plots for comparing barcode replicates
    
    A scatter plot is generated for each section in SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES
    A line is drawn for the mean of the replicates
    A button is added to the plot to select other constructs and other replicates in study.df
    
    """
    x = study.get_df(
        sample = sample, 
        construct = construct, 
        section = SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES, 
        base_type = ['A','C'])
    
    x = x[['sample','construct','section','mut_rates']]
    
    y = study.get_df(
        sample = sample, 
        construct = replicate,
        section = SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES,
        base_type = ['A','C'])
    
    y = y[['sample','construct','section','mut_rates']]

    data = {
        construct: x,
        replicate: y
    }

    fig = go.Figure()
    for section in SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES:
        fig.add_trace(go.Scatter(
            x = x[x['section']==section]['mut_rates'].values[0].tolist(),
            y = y[y['section']==section]['mut_rates'].values[0].tolist(),
            mode = 'markers',
            name = section
        ))
    
    for trace in __corr_scatter_plot(data):
        fig.add_trace(trace)
    
    fig.layout.update(
        title = 'Barcode comparison - {} '.format(sample),
        xaxis_title = 'Barcode: {}'.format(construct),
        yaxis_title = 'Barcode: {}'.format(replicate),
        showlegend = True
    )
    
    return {'fig': fig, 'data': data}
    
def combined_barcode_replicates_in_sample(study):
    
    data = {}
    samples = study.df['sample'].unique()
    samples.sort()
    for sample in samples:
        replicates_lists = generate_dataset.generate_barcode_replicates_pairs(study, sample)       
        data[sample] = generate_dataset.compute_pearson_scores(study, sample, replicates_lists, SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES)
    
    fig = go.Figure()

    # Add trace for each sample and show only one sample at a time
    for sample in data.keys():
        fig.add_trace(go.Bar(
            y = data[sample],
            name = sample,
            visible = False
        ))
        
    # set the first sample to be visible
    fig.data[0].visible = True
    
    fig.layout.update(
        title = 'Combined barcodes replicates - {}'.format(sample),
        xaxis_title = 'Constructs',
        yaxis_title = 'Pearson correlation score average across barcode replicates'
    )
    
    # add a button to show/hide each sample
    # the button should be on the right of the plot, outside of the plot
    fig.layout.updatemenus = [
        go.layout.Updatemenu(
            active = 0,
            buttons = [
                dict(
                    args = [{'visible': [True if i==j else False for i in range(len(data.keys()))]}],
                    label = sample,
                    method = 'update'
                ) for j, sample in enumerate(data.keys())
            ],
            direction = 'down',
            pad = {'r': 10, 't': 10},
            showactive = True,
            x = 0.7,
            xanchor = 'left',
            y = 1.1,
            yanchor = 'top'
        )
    ]
    
    
    return {'fig': fig, 'data': data}

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
    
def barcode_replicates_per_construct(study, samples, construct):
    data = study.get_df(
        sample= samples,
        construct = construct,
        base_type = ['A','C'],
        section='full'
        )
    
    if not len(data) == 2:
        return data
    
    plt.figure()
    
    plt.plot(data['mut_rates'].iloc[0], data['mut_rates'].iloc[1], 'x')
    
    __correlation_scatter_plot(
        x = data['mut_rates'].iloc[0],
        y = data['mut_rates'].iloc[1]
        )
    
    plt.title('Biological replicates - {} '.format(construct))
    plt.xlabel(data['sample'].iloc[0])
    plt.ylabel(data['sample'].iloc[1])
    
    
def sample_replicates_heatmap_per_family(study, samples, family, section):
    if section in ['LAH','MS2']:
        
        data = study.get_df(
            sample=samples, 
            section=section, 
            family=family, 
            base_type=['A','C'])
        
    if section == 'ROI':

        data = study.get_df(
            sample=samples, 
            section=section, 
            family=family)
        
    data = data[['sample','construct','mut_rates','frame_shift_ROI','sequence']]
    
    if section == 'ROI':
        reference = data[data['sample'] == samples[0]].iloc[0]['sequence']
        data['mut_rates'] = data.apply(lambda row: int(row['frame_shift_ROI'])*[np.nan] + list(row['mut_rates']) + (len(reference) - int(row['frame_shift_ROI']) - len(row['mut_rates']))*[np.nan], axis=1)

    data_sample_0 = data[data['sample'] == samples[0]].reset_index(drop=True)
    data_sample_1 = data[data['sample'] == samples[1]].reset_index(drop=True)

    df = pd.DataFrame(
        index = data_sample_0['construct'],
        columns = data_sample_1['construct'], 
    )
        
    for _, row in data_sample_0.iterrows():
        for _, row2 in data_sample_1.iterrows():
            df.loc[row['construct'], row2['construct']] = custom_pearsonr(row['mut_rates'], row2['mut_rates'])
                                
    fig = px.imshow(
                df, 
                x=df.columns, 
                y=df.index, 
                color_continuous_scale='RdBu_r', 
                zmin=-1, 
                zmax=1, 
                title='Pearson between samples {} and {} (family {})'.format(samples[0], samples[1], family),
                labels=dict(x='Sample {}'.format(samples[1]), y='Sample {}'.format(samples[0])),
                width=800,
                height=800,
                text_auto = '.2f'
                )
    fig.show() 
    
    return {'data': data, 'fig': fig}
        
# ---------------------------------------------------------------------------

## c / v - Proportional
# ---------------------------------------------------------------------------


def change_in_temp_mut_frac_vs_temperature(study, samples, construct):
    __mut_frac_vs_sample_attr(study, samples, construct, 'temperature_k', 'Temperature (K)')
    
def change_in_reaction_time(study, samples, construct):
    __mut_frac_vs_sample_attr(study, samples, construct, 'inc_time_tot_secs', 'DMS incubation time (s)')

def change_in_dms_conc(study, samples, construct):
    __mut_frac_vs_sample_attr(study, samples, construct, 'DMS_conc_mM', 'DMS concentration (mM)')

def mut_rate_across_family_vs_deltaG(study, sample, family):
    
    # get a neat dataframe with the mutation rates for each base at each deltaG
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
    
    plt.subplots(nrows=df.shape[1], figsize=(4,20), sharex=True, sharey=True)
    plt.suptitle('Sample: {}, Family: {}'.format(sample, family))
    for i, (col, pair) in enumerate(zip(df.columns, paired)):
        plt.subplot(df.shape[1], 1, i+1)
        plt.scatter(df[col].index, df[col].values, facecolors='r' if pair else 'none', edgecolors='r')
        plt.title(col)
        plt.xlabel('Predicted free energy (kcal/mol)')
        plt.ylabel('Mutation fraction')
        plt.tight_layout()
    plt.tight_layout()
    return df
    
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



# ---------------------------------------------------------------------------

# utils
# ---------------------------------------------------------------------------

def __correlation_scatter_plot(x,y):
        # reshape plot
        plt.axis('square')
        (xmin, xmax), (ymin, ymax) = plt.xlim(), plt.ylim()
        plt.xlim(min(xmin, ymin), max(xmax, ymax))
        xmin, xmax = plt.xlim()
        
        # plot y=x line
        plt.plot([xmin, xmax], [xmin, xmax], 'k--', label='y=x')

        # plot regression line
        x, y = np.array(x), np.array(y)
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
        
        # compute R2
        r2 =  r2_score(y, x)
        
        plt.plot([xmin, xmax], [xmin*slope + intercept, xmax*slope + intercept], 'r-', label='LR: {:.2f}x + {:.2f}'.format(slope, intercept))
        
        # plot labels
        plt.text(0.5, 0.9, 'Pearson correlation: {:.2f}'.format(custom_pearsonr(x,y)), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
        plt.text(0.5, 0.85, 'R2: {:.2f}'.format(r2), horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

def __mut_frac_vs_sample_attr(study, samples, construct, attr, x_label=None):
    # get data
    data = study.get_df(
        sample = samples,
        construct=construct,
        section='ROI',
        base_type = ['A','C']) 

    data = data[['sample','sequence','mut_rates','index_selected','structure',attr]]
    
    if len(data) == 0:
        print('No data: {}, {}'.format(samples, construct))
        return None
    
    # turn it into a dataframe
    df = pd.DataFrame(
        columns= [base + str(idx+1) for base, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])],
        data = np.array(data['mut_rates'].tolist()),
        index= data[attr].values
    )
    
    # plot
    paired = [False if residue == '.' else True for residue in data['structure'].iloc[0]]
    
    plt.subplots(nrows=df.shape[1], figsize=(4,20), sharex=True, sharey=True)
    plt.suptitle('Construct: {}'.format(construct))
    for i, (col, pair) in enumerate(zip(df.columns, paired)):
        plt.subplot(df.shape[1], 1, i+1)
        plt.scatter(df[col].index, df[col].values, facecolors='r' if pair else 'none', edgecolors='r')
        plt.title(col)
        plt.xlabel(x_label) if x_label else plt.xlabel(attr)
        plt.ylabel('Mutation fraction')
        plt.tight_layout()
    plt.tight_layout()
    return df

def __corr_scatter_plot(data):
    x, y = data.values()
    x = np.concatenate(x['mut_rates'].values)
    y = np.concatenate(y['mut_rates'].values)
    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
    r2 =  r2_score(y, x)
    
    traces = []
    # add the regression line
    traces.append(
        go.Scatter(
            x = x,
            y = intercept + slope * x,
            mode = 'lines',
            name = 'LR: {:.2f}x + {:.2f}'.format(slope, intercept),
            line = dict(
                color = ('rgb(205, 12, 24)'),
                width = 2,
            )
        )
    )
    # Plot x = y
    traces.append(
        go.Scatter(
            x = x,
            y = x,
            mode = 'lines',
            name = 'x=y',
            line = dict(
                color = ('rgb(0, 0, 0)'),
                width = 2,
            )
        )
    )
    
    # add the pearson correlation score as a text
    traces.append(
        go.Scatter(
            x = [np.mean(x)],
            y = [np.max(y)],
            mode = 'text',
            text = ['Pearson correlation score: {:.2f}'.format(r_value)],
            textposition = 'top center'
        ))
    
    # add the R2 score as a text
    traces.append(
        go.Scatter(
            x = [np.mean(x)],
            y = [np.max(y)],
            mode = 'text',
            text = ['R2 score: {:.2f}'.format(r2)],
            textposition = 'bottom center'
        ))
    
    return traces
