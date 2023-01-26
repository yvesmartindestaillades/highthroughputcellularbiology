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
import copy 

from scipy.optimize import curve_fit

## b / i - Demultiplexing 
# ---------------------------------------------------------------------------
def mutations_in_barcodes(study):
    
    fig = go.Figure()

    data = study.df[study.df['section'] == 'barcode']

    for sample in data['sample'].unique():
        hist = np.sum(np.stack(data[data['sample']==sample]['num_of_mutations'].values), axis=0)
        bin_edges = np.arange(0, max(np.argwhere(hist != 0)), 1)
        
        fig.add_trace(
            go.Bar(
                x=bin_edges[:-1],
                y=hist,
                name=sample,
                visible=False,
                hovertemplate='Number of mutations: %{x}<br>Number of reads: %{y}<extra></extra>',
                ))
        
    fig.data[0].visible = True
    
    fig.update_layout(barmode='stack', title='Number of mutations in barcodes - {}'.format(data['sample'].unique()[0]))
    fig.update_layout(
        updatemenus=[
            dict(
                active=0,
                buttons=list([
                    dict(label = sample,
                            method = "update",
                            args = [{"visible": [sample == s for s in data['sample'].unique()]},
                                    {"title": 'Number of mutations in barcodes - {}'.format(sample)}])
                    for sample in data['sample'].unique()
                ]), 
                x = 0.7,
                xanchor = 'left',
                y = 1.1,
                yanchor = 'top'
                )
            ])
    
    return {'fig': fig, 'data': data[['sample','construct','num_of_mutations']]}

def num_aligned_reads_per_construct_frequency_distribution(study, sample):
    num_aligned_reads = study.df[(study.df['sample']==sample) & (study.df['section']=='full')]['num_aligned'].to_list()

    return go.Histogram(x=num_aligned_reads, showlegend=False, marker_color='indianred')

def num_aligned_reads_per_construct(study, sample):
    num_aligned_reads = study.df[(study.df['sample']==sample) & (study.df['section']=='full')]['num_aligned'].to_list()
    num_aligned_reads.sort()
    num_aligned_reads.reverse()

    return go.Bar(x=np.arange(1,1+len(num_aligned_reads)), y=num_aligned_reads, showlegend=False, marker_color='indianred')
    
# ---------------------------------------------------------------------------

## b / ii - DMS-MaPseq
# ---------------------------------------------------------------------------
def mutations_per_read(study, sample):
    data = study.df[(study.df['sample']==sample) & (study.df['section']=='full')]['num_of_mutations'].reset_index(drop=True)
    hist = np.sum(np.stack(data.values), axis=0)
    bin_edges = np.arange(0, max(np.argwhere(hist != 0)), 1)
    return go.Bar( x=bin_edges, y=hist, showlegend=False, marker_color='indianred')

def mutation_identity_at_each_position(study, sample, construct):
    data = study.df[(study.df['sample']==sample) & (study.df['construct']==construct) & (study.df['section']=='full')].iloc[0]
    df = pd.DataFrame(index = list(data['sequence']))
    stacked_bar = []
    color_map={'A':'red','C':'blue','G':'yellow','T':'green'}
    for base in ['A','C','G','T']:
        df[base] = np.array(data['mod_bases_'+base])/np.array(data['info_bases'])
        stacked_bar.append( go.Bar(x=np.arange(len(data['sequence'])), y=list(df[base]), marker_color=color_map[base], showlegend=False) )

    return {'fig': stacked_bar, 'data': df}

 
def mutation_fraction_at_each_position(study, sample, construct):
    data = study.df[(study.df['sample']==sample) & (study.df['construct']==construct) & (study.df['section']=='full')].iloc[0]
    df = pd.DataFrame(index = list(data['sequence']))
    stacked_bar = []
    color_map={'A':'red','C':'blue','G':'yellow','T':'green'}
    for base in ['A','C','G','T']:
        df[base] = [mr if b==base  else np.nan for mr, b in zip(data['mut_rates'], data['sequence'])]
        stacked_bar.append( go.Bar(x=np.arange(len(data['sequence'])), y=list(df[base]), marker_color=color_map[base], showlegend=False) )
    
    return {'fig': stacked_bar, 'data': df}

def read_coverage_per_position(study, sample, construct):
    data = study.df[(study.df['sample']==sample) & (study.df['construct']==construct)]
    sections, section_start, section_end = data['section'].unique(), data['section_start'].unique(), data['section_end'].unique()
    idx = np.argsort(section_start)
    sections, section_start, section_end = sections[idx], section_start[idx], section_end[idx]

    scatters = []
    for i, (s, ss, se) in enumerate(zip(sections, section_start, section_end)):
        y_data = copy.copy(data[data['section']==s]['cov_bases'].values[0])
        if s=='full':
            y_data[min(section_start[1:]-1):max(section_end[1:])] = np.zeros(max(section_end[1:])-min(section_start[1:])+1)
        scatters.append(
            go.Bar(
                x=np.arange(ss-1, se),
                y=y_data,
                name=s, 
                marker_color=px.colors.qualitative.Plotly[i],
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

def barcode_comparison_scatter_plot(study, sample):
    """ generate plots for comparing barcode replicates
    
    A scatter plot is generated for each section in SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES
    A line is drawn for the mean of the replicates
    A button is added to the plot to select other constructs and other replicates in study.df
    
    """
    
    fig = go.Figure()
    
    replicates_lists = generate_dataset.generate_barcode_replicates_pairs(study, sample)       
    assert len(replicates_lists) > 0, 'No barcode replicates found for sample {}'.format(sample)
    
    data = study.get_df(
            sample = sample, 
            section = SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES, 
            base_type = ['A','C'])[['sample','construct','section','mut_rates']]
    
    showed_pairs = []
    uniquepairs = []
    for construct, replicates in replicates_lists.items():
        for replicate in replicates:
            if not (replicate, construct) in showed_pairs:
                
                uniquepairs.append((construct, replicate))
                x = data[data['construct']==construct]
                y = data[data['construct']==replicate]

                df = {
                    construct: x,
                    replicate: y
                }

                for section in SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES:
                    fig.add_trace(go.Scatter(
                        x = x[x['section']==section]['mut_rates'].values[0].tolist(),
                        y = y[y['section']==section]['mut_rates'].values[0].tolist(),
                        mode = 'markers',
                        name = section,#' '.join([construct, replicate]),
                        visible=False,
                            
                    ))
                    showed_pairs.append((construct, replicate))
                
                for trace in __corr_scatter_plot(df):
                    fig.add_trace(trace)
                    showed_pairs.append((construct, replicate))
                
                fig.layout.update(
                    title = 'Barcode comparison - {} '.format(sample),
                    xaxis_title = 'Barcode: {}'.format(construct),
                    yaxis_title = 'Barcode: {}'.format(replicate),
                    showlegend = True,
                )

                
    # Add dropdown menu to select construct and replicate
    fig.update_layout(
        updatemenus=[
            dict(
                active=len(uniquepairs)-1,
                buttons=list([
                    dict(label = ' vs '.join([constructs[0], constructs[1]]),
                            method = 'update',
                            args = [{
                                'visible': [True if constructs == c else False for c in showed_pairs]},
                                {'title': 'Barcode comparison - {} - {} vs {}'.format(sample, constructs[0], constructs[1]),
                                'xaxis':
                                    {'title': 'Barcode: {}'.format(constructs[0])},
                                'yaxis':
                                    {'title': 'Barcode: {}'.format(constructs[1])},
                                }], 
                            )
                    for constructs in uniquepairs
                ]),
            
            direction = 'down',
            pad = {'r': 10, 't': 10},
            showactive = True,
            x = 1,
            xanchor = 'right',
            y = 1,
            yanchor = 'bottom'
            )],
    )
    
    # the x and y axis have the same scale
    fig.update_xaxes(scaleanchor = "y")
    fig.update_yaxes(scaleanchor = "x")
    

    # activate the first trace
    for i in range(10):
        fig.data[i].visible = True
        
    data = study.get_df(
            sample = sample, 
            section = 'full', 
            base_type = ['A','C'])[['sample','construct','mut_rates']]
    
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
            y = data[sample]['scores'],
            name = sample,
            visible = False,
            # add hover info to display the construct names from the list data[sample]['constructs']
            hovertemplate='<b>Construct</b>: %{text}<br><b>Score</b>: %{y:.2f}<extra></extra>',
            text = data[sample]['constructs'],
            # remove the text in the bar
            textposition = 'none'
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
    
    data = pd.concat({sample:pd.DataFrame(data[sample]) for sample in data.keys()}).reset_index().rename(columns={'level_0':'sample'}).drop(columns='level_1')

    return {'fig': fig, 'data': data}

# ---------------------------------------------------------------------------

## c / iv - Reproducibility
# ---------------------------------------------------------------------------
def barcode_replicates(study):
    
    data = {}
    samples = study.df['sample'].unique()
    samples.sort()
    
    for sample in samples:    
        replicates_lists = generate_dataset.generate_barcode_replicates_pairs(study, sample)   
        data[sample] = generate_dataset.compute_pearson_scores(study, sample, replicates_lists, SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES)
    
    fig = go.Figure()
    
    # plot a trace per sample and show only one sample at a time
    for sample in data.keys():
        fig.add_trace(go.Histogram(
            x = data[sample]['scores'],
            name = sample,
            visible = False
        ))
        
    # set the first sample to be visible
    fig.data[0].visible = True
    
    fig.layout.update(
        title = 'Barcodes replicates - {}'.format(sample),
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
    
    data = pd.concat({sample:pd.DataFrame(data[sample]) for sample in data.keys()}).reset_index().rename(columns={'level_0':'sample'}).drop(columns='level_1').rename(columns={'constructs':'construct', 'scores':'score'})

    return {'fig': fig, 'data': data}

    
def bio_replicates_per_construct(study, samples, family):
        
    fig = go.Figure()

    unique_constructs = study.get_df(sample=samples, family=family)['construct'].unique()

    big_data = study.get_df(
        sample= samples,
        construct = unique_constructs,
        base_type = ['A','C'],
        section='full'
        )[['sample','construct', 'mut_rates']]
    
    for construct in unique_constructs:
        # add pearson correlation and r2 in big_data columns
        data = big_data.loc[big_data['construct'] == construct]
        if len(data) == 2:
            big_data.loc[big_data['construct'] == construct, 'pearson'] = custom_pearsonr(data['mut_rates'].iloc[0], data['mut_rates'].iloc[1])
            big_data.loc[big_data['construct'] == construct, 'r2'] = r2_score(data['mut_rates'].iloc[0], data['mut_rates'].iloc[1])
    
    assert len(unique_constructs) > 0, 'No constructs found for the given samples and family'
    
    # For each construct plot three traces: the correlation line, the y=x line and the scatter plot
    for i_c, construct in enumerate(unique_constructs):

        data = big_data.loc[big_data['construct'] == construct]
        
        if not len(data) == 2:
            continue
        
        # Plot the scatter plot of the two replicates
        fig.add_trace(
            go.Scatter(x=data['mut_rates'].iloc[0], y=data['mut_rates'].iloc[1], 
            mode='markers', marker=dict(color='blue'),
            showlegend=False, visible=False,
            ))
        
            # Add Pearson correlation and r2 as annotation

        
        # Plot the correlation line and y=x line
        __correlation_scatter_plot(
            x = data['mut_rates'].iloc[0],
            y = data['mut_rates'].iloc[1],
            fig=fig
            )

    assert len(fig.data) > 0, 'No pairs found for the given samples and family'
    # Show the three traces for the first construct
    for i in range(3):
        fig.data[i].visible = True

    # add a button to show/hide the traces for each construct
    # the button should be on the right of the plot, outside of the plot
    fig.layout.updatemenus = [
        go.layout.Updatemenu(
            active = 0,
            buttons = [
                dict(
                    # show/hide the traces for the construct with the annotation
                    args = [
                        {'visible': [True if i==j or i==j+1 or i==j+2 else False for i in range(len(fig.data))]},
                        # update the title of the plot with r2 and pearson correlation
                        {'title': '{} - {} - r2: {:.2f} - pearson: {:.2f}'.format(family, construct, big_data.loc[big_data['construct'] == construct]['r2'].iloc[0], big_data.loc[big_data['construct'] == construct]['pearson'].iloc[0])}
                        ],

                    label = construct,
                    method = 'update'
                ) for j, construct in enumerate(unique_constructs)
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

    # Set the titles and aspect ratio to 1
    fig.layout.update(
        title = 'Biological replicates - {}'.format(family),
        xaxis_title = samples[0],
        yaxis_title = samples[1],
        xaxis = dict(scaleanchor = "y", scaleratio = 1),
        yaxis = dict(scaleanchor = "x", scaleratio = 1),
        height = 700
    )


    

    #TODO data
    return {'fig': fig, 'data': big_data}

    
    
def sample_replicates_heatmap_per_family(study, samples, family, section):
    if section in ['LAH','MS2']:
        
        data = study.get_df(
            sample=samples, 
            section=section, 
            family=family, 
            base_type=['A','C'])
        
        data = data[['sample','construct','section','mut_rates','sequence']]

        
    if section == 'ROI':

        data = study.get_df(
            sample=samples, 
            section=section, 
            family=family)
        
        data = data[['sample','construct','section','mut_rates','frame_shift_ROI','sequence']]
    
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
    return __mut_frac_vs_sample_attr(study, samples, construct, 'temperature_k', 'Temperature (K)')
    
def change_in_reaction_time(study, samples, construct):
    return __mut_frac_vs_sample_attr(study, samples, construct, 'inc_time_tot_secs', 'DMS incubation time (s)')

def change_in_dms_conc(study, samples, construct):
    return __mut_frac_vs_sample_attr(study, samples, construct, 'DMS_conc_mM', 'DMS concentration (mM)')

def mut_rate_across_family_vs_deltaG(study, sample, family):
    
    # get a neat dataframe with the mutation rates for each base at each deltaG
    data = study.get_df(sample=sample, family=family, section='ROI')
    
    data['deltaG'] = data['deltaG'].apply(lambda x: 0 if x == 'void' else float(x))
    
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
    
    fig = go.Figure()
    
    for col, pair in zip(df.columns, paired):
        fig.add_trace(go.Scatter(
            x = df.index,
            y = df[col],
            mode = 'markers',
            name = col,
            visible = 'legendonly',
            legendgroup=col,
            marker=dict(
                size=20,
                line=dict(
                    color='black',
                    width=0 if pair else 3
                )
            ),
    ))

    # Function to fit
    def sigmoid(x, a, b, c):
        RT = 1.987204258*310/1000
        return a / (1 + b*np.exp(-x/RT)) + c

    # For each column of df, fit a sigmoid, and plot it with the confidance interval
    colors = px.colors.qualitative.Plotly + px.colors.qualitative.D3 + px.colors.qualitative.G10 + px.colors.qualitative.T10
    for i_b, (base, mut_rate) in enumerate(df.iteritems()):
        dG = df.index[~np.isnan(mut_rate)].values
        mut_rate = mut_rate[~np.isnan(mut_rate)].values

        if len(mut_rate) >= 3:
            popt, pcov = curve_fit(sigmoid, dG, mut_rate, p0=[0.04, 0.02, 0.00], bounds=([0, 0, 0], [0.1, np.inf, 0.05]), max_nfev=1000)

            xdata_MC = np.linspace(min(dG), max(dG), 100)
            y_fit = sigmoid(xdata_MC, *popt)

            # Do a Monte Carlo simulation to estimate the uncertainty in the fit parameters using a multinormal distribution
            # with the covariance matrix as the covariance matrix
            N = 1000
            param_samples = np.clip(np.random.multivariate_normal(popt, pcov, N).T.reshape(3, N, 1), 0, np.inf)

            # Compute the sigmoid for each set of parameters for each x value
            y_MC = sigmoid(xdata_MC.reshape(1, -1) , param_samples[0], param_samples[1], param_samples[2])

            fig.add_trace(go.Scatter(x=xdata_MC, y=y_fit, mode='lines', name='fit',
                                    visible = 'legendonly', legendgroup=base, showlegend=False, 
                                    marker_color=colors[i_b]))

            fig.add_trace(go.Scatter(x=np.concatenate((xdata_MC, xdata_MC[::-1])), # x, then x reversed
                                    y=np.concatenate((np.percentile(y_MC, 97.5, axis=0), np.percentile(y_MC, 2.5, axis=0)[::-1])), # upper, then lower reversed
                                    fill='toself',
                                    fillcolor=colors[i_b],
                                    line=dict(color=colors[i_b]),
                                    opacity=0.1,
                                    hoverinfo="skip",
                                    visible = 'legendonly', legendgroup=base, showlegend=False
                                ))

        
    fig.update_layout(
        title = 'Mutation rate across family members for sample {} (family {})'.format(sample, family),
        xaxis_title = 'Predicted free energy (kcal/mol)',
        yaxis_title = 'Mutation rate',
    )
    
    fig.update_layout(dict(updatemenus=[
                        dict(
                            type = "buttons",
                            direction = "left",
                            buttons=list([
                                dict(
                                    args=["visible", "legendonly"],
                                    label="Deselect All",
                                    method="restyle"
                                ),
                                dict(
                                    args=["visible", True],
                                    label="Select All",
                                    method="restyle"
                                )
                            ]),
                            pad={"r": 10, "t": 10},
                            showactive=False,
                            x=1,
                            xanchor="right",
                            y=1.1,
                            yanchor="top"
                        ),
                    ]
              ))

    return {'data': data, 'fig': fig}
    
# ---------------------------------------------------------------------------

## d / kfold
# ---------------------------------------------------------------------------

def heatmap_across_family_members(study, sample):
    
    fig = go.Figure()
    
    big_data = study.df[(study.df['sample'] == sample) & (study.df['section'] == 'ROI')]
    
    # add each family as a plot but show only one at the time
    for family in big_data['family'].unique():
        
        data = big_data[(big_data['family'] == family)]
        
        data_all_samples = study.df[(study.df['family'] == family) & (study.df['section'] == 'ROI')]
        sequences_ROI = data_all_samples.sort_values('sequence', key=lambda x: x.str.len(), ascending=False)['sequence'].values
        reference = sequences_ROI[0]
        
        df = pd.DataFrame(
            columns = [base + str(idx + 1) for base, idx in zip(reference, range(len(reference)))],
            data = [int(offset)*[np.nan] + list(mr) + [np.nan]*(len(reference)-int(offset)-len(mr)) for offset, mr in zip(data['frame_shift_ROI'], data['mut_rates'])],
            index = data['construct'].values 
        ).sort_index(ascending=False)
        
        
        fig.add_trace(go.Heatmap(
            z = df.values,
            x = df.columns,
            y = df.index,
            name = 'Family {}'.format(family),
            visible = False,
            colorscale = 'Viridis',
            colorbar = dict(
                title = 'Mutation fraction',
                titleside = 'right',
                tickmode = 'array',
                tickvals = [0, 0.02, 0.04, 0.06, 0.08],
                ticktext = ['0', '0.02', '0.04', '0.06', '0.08'],
            ),
            showscale = True,
            zmin = 0,
            zmax = 0.08,
            hovertemplate = 'Construct: %{y}<br>Base: %{x}<br>Mutation fraction: %{z:.2f}<extra></extra>',
        ))
    
        
        # place the xticks on top of the plot
        fig.update_xaxes(side="top")
        
            
        
     #   fig.update_xaxes(ticklabelposition="inside top")
        
        # make the background white
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
        )
   
    # add a button to select a single family to display

    # show the first plot
        # set the first sample to be visible
    fig.data[0].visible = True
    
    # add a button to show/hide each sample
    # the button should be on the right of the plot, outside of the plot
    
    fig.layout.updatemenus = [
        go.layout.Updatemenu(
            active = 0,
            buttons = [
                dict(
                    args = [{'visible': [True if i==j else False for i in range(len(big_data['family'].unique()))]}],
                    label = sample,
                    method = 'update'
                ) for j, sample in enumerate(big_data['family'].unique())
            ],
            direction = 'down',
            pad = {'r': 10, 't': 10},
            showactive = True,
            x = 1,
            xanchor = 'left',
            y = 1.2,
            yanchor = 'top'
        )
    ]
        
    fig.layout.title = 'Mutation rate across family members with all bases - sample {}'.format(sample)
    
    return {'data': big_data, 'fig': fig}
    """
    plt.figure(figsize=(data.shape[1]//2, data.shape[0]//2))
    plt.gca().xaxis.tick_top()
    sns.heatmap(df, annot=True, fmt='.2f')
    plt.title('Heatmap across family members with all bases - sample {} - family {}'.format(sample, family))
    plt.ylabel('Construct')
    """


# ---------------------------------------------------------------------------

# utils
# ---------------------------------------------------------------------------

def __correlation_scatter_plot(x, y, fig):
        # reshape plot
        (xmin, xmax) = min(x), max(x)
        (ymin, ymax) = min(y), max(y)


        # Plot a straight line y=x
        fig.add_trace(
            go.Scatter(x=[xmin, xmax], y=[xmin, xmax], mode='lines', name='y=x', line=dict(color='black', dash='dash'),
            showlegend=True, visible=False))

        # plot regression line
        x, y = np.array(x), np.array(y)
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
        
        # compute R2
        r2 =  r2_score(y, x)
        
        # Plot regression line
        fig.add_trace(
            go.Scatter(x=[xmin, xmax], y=[xmin*slope + intercept, xmax*slope + intercept],
            mode='lines', name='LR: {:.2f}x + {:.2f}'.format(slope, intercept), line=dict(color='red'), 
            showlegend=True, visible=False))

        
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
        return {'data': [], 'fig': go.Figure()}
    
    # turn it into a dataframe
    df = pd.DataFrame(
        columns= [base + str(idx+1) for base, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])],
        data = np.array(data['mut_rates'].tolist()),
        index= data[attr].values
    )
    
    # plot
    paired = [False if residue == '.' else True for residue in data['structure'].iloc[0]]
    
    fig = go.Figure()
    
    for col, pair in zip(df.columns, paired):
        fig.add_trace(go.Scatter(
            x = df.index,
            y = df[col],
            mode = 'markers',
            name = col,
            visible = 'legendonly',
            marker=dict(
                size=20,
                line=dict(
                    color='black',
                    width=0 if pair else 3
                )
            ),
    ))
        
    fig.update_layout(
        title='Construct: {}'.format(construct),
        xaxis_title=x_label if x_label else attr,
        yaxis_title='Mutation fraction')

    
    return {'data': df, 'fig': fig}

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
