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
import plotly.graph_objects as go
from scipy.stats import norm
from plotly.subplots import make_subplots
from util import compute_wilson_interval

MARKER_RATIO = 2


## b / i - Demultiplexing 
# ---------------------------------------------------------------------------
def mutations_in_barcodes(study):
    
    fig = go.Figure()

    data = study.df[study.df['section'] == 'barcode']

    for sample in data['sample'].unique():
        hist = np.sum(np.stack(data[data['sample']==sample]['sub_hist'].values), axis=0)
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
    
    return {'fig': fig, 'data': data[['sample','reference','sub_hist']]}

def num_aligned_reads_per_reference_frequency_distribution(study, sample):
    num_aligned_reads = study.df[(study.df['sample']==sample) & (study.df['section']=='full')]['num_aligned'].to_list()

    return go.Histogram(x=num_aligned_reads, showlegend=False, marker_color='indianred')

def num_aligned_reads_per_reference(study, sample):
    num_aligned_reads = study.df[(study.df['sample']==sample) & (study.df['section']=='full')]['num_aligned'].to_list()
    num_aligned_reads.sort()
    num_aligned_reads.reverse()

    return go.Bar(x=np.arange(1,1+len(num_aligned_reads)), y=num_aligned_reads, showlegend=False, marker_color='indianred')
    
# ---------------------------------------------------------------------------

## b / ii - DMS-MaPseq
# ---------------------------------------------------------------------------
def _mutations_per_read_subplot(study, sample):
    data = study.df[(study.df['sample']==sample) & (study.df['section']=='full')]['sub_hist'].reset_index(drop=True)
    hist = np.sum(np.stack(data.values), axis=0)
    bin_edges = np.arange(0, max(np.argwhere(hist != 0)), 1)
    return go.Bar( x=bin_edges, y=hist, showlegend=False, marker_color='indianred')

def mutations_per_read_per_sample(study, unique_samples):
    if not isinstance(unique_samples, list) and not isinstance(unique_samples, tuple) and not isinstance(unique_samples, np.ndarray) and not isinstance(unique_samples, set):
        unique_samples = [unique_samples]
        
    fig = make_subplots(rows=len(unique_samples), cols=1, vertical_spacing=0.4/len(unique_samples),
                        subplot_titles=['Number of mutations per read - {}'.format(sample) for sample in unique_samples])
    for i_s, sample in enumerate(unique_samples):
        fig.add_trace(_mutations_per_read_subplot(study, sample), row=i_s+1, col=1 )
        fig.update_yaxes(title='Count')
        fig.update_xaxes(dtick=10)

    fig.update_layout(autosize=True, height=len(unique_samples)*500, title='Number of mutation per read across samples')
    return {
        'fig':fig,
        'data':study.df[study.df['section']=='full'][['sample','reference','sub_hist']]
        }
    
def mutation_per_read_per_reference(study, sample, reference):
    data = study.df[(study.df['sample']==sample) & (study.df['reference']==reference) & (study.df['section']=='full')]['sub_hist'].reset_index(drop=True)
    if len(data) == 0:
        return None
    data = data.iloc[0]
    MAX_MUTATIONS = 20
    # normalize by the number of reads
    fig = px.bar(x=np.arange(0,MAX_MUTATIONS), y=data[:MAX_MUTATIONS])

    fig.update_layout(barmode='stack')
    fig.update_layout(title='Number of mutations per read - {} - {}'.format(sample, reference))
    fig.update_yaxes(title='Count')
    fig.update_xaxes(title='Number of mutations per read')
    return {
        'fig':fig,
        'data':data
        }

def _mutation_identity_at_each_position_subplot(study, sample, reference, section='full'):
    data = study.df[(study.df['sample']==sample) & (study.df['reference']==reference) & (study.df['section']==section)].iloc[0]
    df, df_err_min, df_err_max = pd.DataFrame(index = list(data['sequence'])), pd.DataFrame(index = list(data['sequence'])), pd.DataFrame(index = list(data['sequence']))
    stacked_bar = []
    color_map={'A':'red','C':'blue','G':'yellow','T':'green'}

    data['err_min'] = [compute_wilson_interval(p, data['num_aligned'])[0] for p in data['sub_rate']]
    data['err_max'] = [compute_wilson_interval(p, data['num_aligned'])[1] for p in data['sub_rate']]

    for base in ['A','C','G','T']:
        df[base] = np.array(data['sub_'+base])/np.array(data['info'])
        stacked_bar.append( go.Bar(x=np.arange(len(data['sequence'])), y=list(df[base]), marker_color=color_map[base], showlegend=False) )
    
    # add error bars to stacked_bar[-1]
    stacked_bar[-1]['error_y'] = dict(
        type='data',
        array= [data['err_max'][i]-data['sub_rate'][i] for i in range(len(data['sequence']))],
        arrayminus = [data['sub_rate'][i]-data['err_min'][i] for i in range(len(data['sequence']))],
        visible=True,
        symmetric=False,
        thickness=1.5,
        width=2,
        color='black'
    )
        
    return {'fig':stacked_bar, 'data':df}

def mutation_identity_at_each_position(study, sample, unique_references, section='full'):

    if not isinstance(unique_references, list) and not isinstance(unique_references, tuple)  and not isinstance(unique_references, np.ndarray)  and not isinstance(unique_references, set):
        unique_references = [unique_references]
    
    fig = make_subplots(rows=len(unique_references), cols=1, vertical_spacing=0.2/len(unique_references),
                    subplot_titles=['Mutation identity at each position - {} - {} reads'.format(cst, reads) for cst, reads in zip(unique_references, study.df[(study.df['sample']==sample)&(study.df['reference'].isin(unique_references))&(study.df['section']=='full')]['num_aligned'])])
    
    for i_c, reference in enumerate(unique_references):
        muts_identity = _mutation_identity_at_each_position_subplot(study, sample, reference, section)

        for bar in muts_identity['fig']:
            fig.add_trace( bar, row=i_c+1, col=1 )
        
        fig.update_xaxes(tickangle=0, 
                tickvals=np.arange(len(muts_identity['data'].index)), ticktext=list(muts_identity['data'].index), tickfont={'size':8},
                row=i_c+1, col=1)
            
    for trace, name in zip(fig["data"][:4], ['A','C','G','T']):
        trace.update(showlegend=True)
        trace["name"] = name

    fig.update_yaxes(title='Mutation fraction')
    fig.update_layout(barmode='stack', height=500*len(unique_references), width=1500)
    plot = {
        'fig':fig,
        'data':study.df[
            (study.df['sample']==sample)&
            (study.df['reference'].isin(unique_references))&
            (study.df['section']=='full')]\
        [['sample','reference','sub_A','sub_C','sub_G','sub_T','num_aligned']]
        }
    return plot

 
def _mutation_fraction_at_each_position_subplot(study, sample, reference, section='full'):
    data = study.df[(study.df['sample']==sample) & (study.df['reference']==reference) & (study.df['section']==section)]
    
    if len(data) >= 1:
        data = data.iloc[0]
    else:
        print('No data for sample {} reference {} section {}'.format(sample, reference, section))
        return {'fig':[], 'data':[]}
    
    df, df_err_min, df_err_max = pd.DataFrame(index = list(data['sequence'])), pd.DataFrame(index = list(data['sequence'])), pd.DataFrame(index = list(data['sequence']))
    stacked_bar = []
    
    data['err_min'] = [compute_wilson_interval(p, data['num_aligned'])[0] for p in  data['sub_rate']]
    data['err_max'] = [compute_wilson_interval(p, data['num_aligned'])[1] for p in data['sub_rate']]

    color_map={'A':'red','C':'blue','G':'yellow','T':'green'}
    for base in ['A','C','G','T']:
        for d, col in zip([df, df_err_min, df_err_max], ['sub_rate', 'err_min', 'err_max']):
            d[base] = [mr if b==base  else np.nan for mr, b in zip(data[col], data['sequence'])]
        stacked_bar.append(
            go.Bar
            (x=np.arange(len(data['sequence'])), 
            y=list(df[base]), 
            marker_color=color_map[base], 
            error_y=dict(type='data', symmetric=False, array=list(df_err_max[base]-df[base]), arrayminus=list(df[base]-df_err_min[base])), 
            showlegend=False)) 
    
    return {'fig': stacked_bar, 'data': df}

def mutation_fraction_at_each_position(study, sample, unique_references, section='full'):

    if not isinstance(unique_references, list) and not isinstance(unique_references, tuple)  and not isinstance(unique_references, np.ndarray)  and not isinstance(unique_references, set):
        unique_references = [unique_references]

    fig = make_subplots(rows=len(unique_references), cols=1, vertical_spacing=0.2/len(unique_references),
                        subplot_titles=['Mutation fraction at each position - {}'.format(cst) for cst in unique_references])
    for i_c, reference in enumerate(unique_references):
        muts_identity = _mutation_fraction_at_each_position_subplot(study, sample, reference, section)

        if len(muts_identity['fig']) == 0:
            continue
        
        for bar in muts_identity['fig']:
            fig.add_trace( bar, row=i_c+1, col=1 )
        
        fig.update_xaxes(tickangle=0, 
                tickvals=np.arange(len(muts_identity['data'].index)), ticktext=list(muts_identity['data'].index), tickfont={'size':8},
                row=i_c+1, col=1)
            
    for trace, name in zip(fig["data"][:4], ['A','C','G','T']):
        trace.update(showlegend=True)
        trace["name"] = name

    fig.update_yaxes(title='Mutation fraction')
    fig.update_layout(barmode='stack', height=500*len(unique_references), width=1500)
    plot = {
        'fig':fig,
        'data':study.df[
            (study.df['sample']==sample)&
            (study.df['reference'].isin(unique_references))&
            (study.df['section']=='full')]\
        [['sample','reference','sub_rate','num_aligned']]
        }
    
    return plot



def _read_coverage_per_position_subplot(study, sample, reference):
    data = study.df[(study.df['sample']==sample) & (study.df['reference']==reference)]
    sections, section_start, section_end = data['section'].unique(), data['section_start'].unique(), data['section_end'].unique()
    idx = np.argsort(section_start)
    sections, section_start, section_end = sections[idx], section_start[idx], section_end[idx]

    scatters = []
    for i, (s, ss, se) in enumerate(zip(sections, section_start, section_end)):
        y_data = copy.copy(data[data['section']==s]['cov'].values[0])
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

def read_coverage_per_position(study, sample, unique_references):
    
    if not isinstance(unique_references, list) and not isinstance(unique_references, tuple)  and not isinstance(unique_references, np.ndarray)  and not isinstance(unique_references, set):
        unique_references = [unique_references]

    fig = make_subplots(rows=len(unique_references), cols=1, vertical_spacing=0.2/len(unique_references),
                        subplot_titles=['Read coverage per position - {}'.format(cst) for cst in unique_references])
    
    for i_c, reference in enumerate(unique_references):
        read_coverage = _read_coverage_per_position_subplot(study, sample, reference)

        for bar in read_coverage['fig']:
            fig.add_trace( bar, row=i_c+1, col=1 )

    # print a legend for each section
    for trace, name in zip(fig["data"][:len(read_coverage['data']['section'])], read_coverage['data']['section']):
        trace.update(showlegend=True)
        trace["name"] = name

    fig.update_yaxes(title='Read coverage')
    fig.update_layout(barmode='stack', height=500*len(unique_references), width=1300)
    plot = {
        'fig':fig,
        'data':study.df[
            (study.df['sample']==sample)&
            (study.df['reference'].isin(unique_references))&
            (study.df['section']=='full')]\
        [['sample','reference','cov']]
        }

    return plot
    
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
    
    
def replicates_fisher_pearson_decorator(function):
    def replicates_fisher_pearson(study, sample):
        replicates_lists, data, title = function(study, sample)
        data['sub_rate'] = data.apply(lambda x: [x['sub_rate'][i] for i in range(len(x['sequence'])) if x['sequence'][i] in ['A','C']], axis=1)
        data.reset_index(inplace=True, drop=True)

        p_values_true_data = []
        p_values_fake_data = []
        pearsonr_true_data = []
        pearsonr_fake_data = []
        
        # Calculate the Fisher's exact test
        def fisher_exact_test(n1, p1, n2, p2):
            # Calculate the odds ratio
            oddsratio, pvalue = stats.fisher_exact([[n1*p1, n1*(1-p1)], [n2*p2, n2*(1-p2)]])
            return oddsratio, pvalue
        
        # Combine the p-values using the Fisher's method
        def combine_pvalues(pvalues):
            return -2 * np.sum(np.log(pvalues))


        def fisher_compute_and_combine_pvalues(n1, mr1, n2, mr2):
            # Calculate the p-values
            pval = []
            for p1, p2 in zip(mr1, mr2):
                oddsratio, pvalue = fisher_exact_test(n1, p1, n2, p2)
                pval.append(pvalue)

            # Combine the p-values
            combined_pvalue = combine_pvalues(pval)

            # Compute the p-value from the chi-squared distribution
            pvalue = 1 - stats.chi2.cdf(combined_pvalue, 2 * len(pval))
            
            return pvalue

        for reference, replicates in replicates_lists.items():
            for replicate in replicates:
                if reference not in data['reference'].values or replicate not in data['reference'].values:
                    print('reference or replicate not found in the data: {} {}'.format(reference, replicate))
                    continue
                
                # using the true data
                n1 = data[data['reference'] == reference]['num_aligned'].values[0]
                n2 = data[data['reference'] == replicate]['num_aligned'].values[0]
                mr1 = data[data['reference'] == reference]['sub_rate'].values[0]
                mr2 = data[data['reference'] == replicate]['sub_rate'].values[0]

                p_values_true_data.append(fisher_compute_and_combine_pvalues(n1, mr1, n2, mr2))
                pearsonr_true_data.append(stats.pearsonr(mr1, mr2)[0])

                # Using the fake data
                mr1 = np.random.binomial(n1, mr1)/n1
                mr2 = np.random.binomial(n2, mr1)/n2
                
                p_values_fake_data.append(fisher_compute_and_combine_pvalues(n1, mr1, n2, mr2))
                pearsonr_fake_data.append(stats.pearsonr(mr1, mr2)[0])
                
                # TODO: average the replicates 
                
        fig = make_subplots(rows=2, cols=1, subplot_titles=('Pearson correlation', 'Combined p-values from Fisher exact test'), vertical_spacing=0.2, shared_xaxes=False)
        # increase height
        fig.update_layout(height=800)


        # group the legend with the names
        fig.add_trace(go.Histogram(x=pearsonr_true_data, histnorm='percent',  marker_color='blue', name='True data', xbins=dict(start=0, end=1, size=0.05), legendgroup='True data'), row=1, col=1)
        fig.add_trace(go.Histogram(x=pearsonr_fake_data, histnorm='percent',  marker_color='red', name='Fake data: split a bitvector', xbins=dict(start=0, end=1, size=0.05), legendgroup='Fake data'), row=1, col=1)

        fig.add_trace(go.Histogram(x=p_values_true_data, histnorm='percent', marker_color='blue',  showlegend=False, xbins=dict(start=0, end=1, size=0.05), hovertemplate='p-value: %{x:.2f}<extra></extra>', legendgroup='True data'), row=2, col=1)
        fig.add_trace(go.Histogram(x=p_values_fake_data, histnorm='percent', marker_color='red', showlegend=False, xbins=dict(start=0, end=1, size=0.05), hovertemplate='p-value: %{x:.2f}<extra></extra>', legendgroup='Fake data'), row=2, col=1)

        fig['layout']['title'] = title
        fig['layout']['xaxis']['title']= 'Pearson correlation'
        fig['layout']['xaxis2']['title']='p-value of the Fisher exact test per residue, combined per reference with the Fisher method'
        fig['layout']['yaxis']['title']='Percent'
        fig['layout']['yaxis2']['title']='Percent'
        
        return {'data': data, 'fig': fig}
    return replicates_fisher_pearson


@replicates_fisher_pearson_decorator
def barcodes_replicates_fisher_pearson(study, sample):
    replicates_lists = generate_dataset.generate_barcode_replicates_pairs(study, sample)       
    data = study.df[(study.df['sample'] == sample) & (study.df['section'] == 'full')][['sample','reference','section','sub_rate','num_aligned','sequence']]
    barcode_bounds = [139,151]

    for col in ['sub_rate','sequence']:
        data[col] = data[col].apply(lambda x: np.array(x[:barcode_bounds[0]]).tolist() + np.array(x[barcode_bounds[1]:]).tolist())
        
    title = 'Similarity test between barcode replicates using Pearson correlation and Fisher exact test for sample {}'.format(sample)

    return replicates_lists, data, title

@replicates_fisher_pearson_decorator
def biological_replicates_fisher_pearson(study, sample):
    data = study.df[(study.df['sample'].isin(sample)) & (study.df['section'] == 'full')][['sample','reference','section','sub_rate','num_aligned','sequence']]
    if sample[0] != sample[1]:
        data = data.groupby('reference').filter(lambda x: len(x) == 2)
    replicates_lists = {reference+'_'+sample[0]: [reference+'_'+sample[1]] for reference in data['reference'].unique()}
    data['reference'] = data.apply(lambda x: x['reference'] + '_' + x['sample'], axis=1)
    
    title = 'Similarity test between biological replicates using Pearson correlation and Fisher exact test for samples {} and {}'.format(sample[0], sample[1])
    
    return replicates_lists, data, title
        
# ---------------------------------------------------------------------------

## c / iii - Barcode
# ---------------------------------------------------------------------------
SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES = ['MS2','TC1','ROI','TC2','LAH','buffer']

def barcode_comparison_scatter_plot(study, sample): #TODO
    """ generate plots for comparing barcode replicates
    
    A scatter plot is generated for each section in SECTION_TO_COMPARE_FOR_BARCODE_REPLICATES
    A line is drawn for the mean of the replicates
    A button is added to the plot to select other references and other replicates in study.df
    
    """
    
    fig = go.Figure()
    
    replicates_lists = generate_dataset.generate_barcode_replicates_pairs(study, sample)       
    assert len(replicates_lists) > 0, 'No barcode replicates found for sample {}'.format(sample)
    
    barcode_bounds = [139,151]
    
    data = study.get_df(
            sample = sample, 
            section = 'full')[['sample','reference','section','sub_rate','sequence']]
    

    for col in ['sub_rate','sequence']:
        data[col] = data[col].apply(lambda x: np.array(x[:barcode_bounds[0]]).tolist() + np.array(x[barcode_bounds[1]:]).tolist())
    
    data['sub_rate'] = data.apply(lambda x: np.array([x['sub_rate'][i] for i in range(len(x['sequence'])) if x['sequence'][i] in ['A','C']]), axis=1)    
    
    showed_pairs = []
    uniquepairs = []
    for reference, replicates in replicates_lists.items():
        for replicate in replicates:
            if not (replicate, reference) in showed_pairs:
                
                uniquepairs.append((reference, replicate))
                x = data[data['reference']==reference]['sub_rate'].values[0]
                y = data[data['reference']==replicate]['sub_rate'].values[0]
                
                assert len(x) == len(y), 'The length of the two replicates are not the same: {} vs {}'.format(len(x), len(y))

                fig.add_trace(go.Scatter(
                    x = x,
                    y = y, 
                    mode = 'markers',
                    name = 'mutation rates',#' '.join([reference, replicate]),
                    visible=False,
                        
                ))
                showed_pairs.append((reference, replicate))
                    
                for trace in __corr_scatter_plot(x, y, visible=False):
                    fig.add_trace(trace)
                    showed_pairs.append((reference, replicate))
                
                fig.layout.update(
                    title = 'Barcode comparison - {} '.format(sample),
                    xaxis_title = 'Barcode: {}'.format(reference),
                    yaxis_title = 'Barcode: {}'.format(replicate),
                    showlegend = True,
                )

                
    # Add dropdown menu to select reference and replicate
    fig.update_layout(
        updatemenus=[
            dict(
                active=len(uniquepairs)-1,
                buttons=list([
                    dict(label = ' vs '.join([references[0], references[1]]),
                            method = 'update',
                            args = [{
                                'visible': [True if references == c else False for c in showed_pairs]},
                                {'title': 'Barcode comparison - {} - {} vs {}'.format(sample, references[0], references[1]),
                                'xaxis':
                                    {'title': 'Barcode: {}'.format(references[0])},
                                'yaxis':
                                    {'title': 'Barcode: {}'.format(references[1])},
                                }], 
                            )
                    for references in uniquepairs
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
    for i in range(5):
        fig.data[i].visible = True
        
    data = study.get_df(
            sample = sample, 
            section = 'full', 
            base_type = ['A','C'])[['sample','reference','sub_rate']]
    
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
            # add hover info to display the reference names from the list data[sample]['references']
            hovertemplate='<b>reference</b>: %{text}<br><b>Score</b>: %{y:.2f}<extra></extra>',
            text = data[sample]['references'],
            # remove the text in the bar
            textposition = 'none'
        ))
        
    # set the first sample to be visible
    fig.data[0].visible = True
    
    fig.layout.update(
        title = 'Combined barcodes replicates - {}'.format(sample),
        xaxis_title = 'references',
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
        xaxis_title = 'references',
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
    
    data = pd.concat({sample:pd.DataFrame(data[sample]) for sample in data.keys()}).reset_index().rename(columns={'level_0':'sample'}).drop(columns='level_1').rename(columns={'references':'reference', 'scores':'score'})

    return {'fig': fig, 'data': data}

    
def bio_replicates_per_reference(study, samples, family, correct_bias=False):
        
    fig = go.Figure()

    unique_references = study.get_df(sample=samples, family=family)['reference'].unique()

    big_data = study.get_df(
        sample= samples,
        reference = unique_references,
        base_type = ['A','C'],
        section='full'
        )[['sample','reference', 'sub_rate']]
    
    if correct_bias:
        sub_rate = {}
        for sample in samples:
            sub_rate[sample] = big_data[big_data['sample']==sample]['sub_rate'].values[0]
    
    for reference in unique_references:
        # add pearson correlation and r2 in big_data columns
        data = big_data.loc[big_data['reference'] == reference]
        if len(data) == 2:
            big_data.loc[big_data['reference'] == reference, 'pearson'] = custom_pearsonr(data['sub_rate'].iloc[0], data['sub_rate'].iloc[1])
            big_data.loc[big_data['reference'] == reference, 'r2'] = r2_score(data['sub_rate'].iloc[0], data['sub_rate'].iloc[1])
    
    assert len(unique_references) > 0, 'No references found for the given samples and family'
    
    trace_id = []
    
    # For each reference plot three traces: the correlation line, the y=x line and the scatter plot
    for i_c, reference in enumerate(unique_references):

        data = big_data.loc[big_data['reference'] == reference]
        
        if not len(data) == 2:
            continue
        
        x=data['sub_rate'].iloc[0]
        y=data['sub_rate'].iloc[1]
        
        # Plot the scatter plot of the two replicates
        fig.add_trace(
            go.Scatter(x=x, y=y, 
            mode='markers', marker=dict(color='blue'),
            showlegend=False, visible=False,
            ))
        
        
            # Add Pearson correlation and r2 as annotation
        # add the pearson correlation score as a text
        
        fig.add_trace(
            go.Scatter(
                x = [np.mean(x)],
                y = [np.max(y)],
                mode = 'text',
                text = ['Pearson correlation score: {:.2f}'.format(custom_pearsonr(x,y))],
                textposition = 'top center',
                visible=False,
                showlegend=False
            ))
        
        # add the R2 score as a text
        fig.add_trace(
            go.Scatter(
                x = [np.mean(x)],
                y = [np.max(y)],
                mode = 'text',
                text = ['R2 score: {:.2f}'.format(r2_score(x,y))],
                textposition = 'bottom center',
                visible=False,
                showlegend=False
            ))
            
        # Plot the correlation line and y=x line
        __correlation_scatter_plot(
            x = data['sub_rate'].iloc[0],
            y = data['sub_rate'].iloc[1],
            fig=fig
            )

        trace_id += [i_c]*5
        
    assert len(fig.data) > 0, 'No pairs found for the given samples and family'
    # Show the three traces for the first reference
    for i in range(5):
        fig.data[i].visible = True

    # add a button to show/hide the traces for each reference
    # the button should be on the right of the plot, outside of the plot
    fig.layout.updatemenus = [
        go.layout.Updatemenu(
            active = 0,
            buttons = [
                dict(
                    # show/hide the traces for the reference with the annotation
                    args = [
                        {'visible': [True if i==j else False for i in trace_id]},
                        # update the title of the plot with r2 and pearson correlation
                        {'title': '{} - {} - r2: {:.2f} - pearson: {:.2f}'.format(family, reference, big_data.loc[big_data['reference'] == reference]['r2'].iloc[0], big_data.loc[big_data['reference'] == reference]['pearson'].iloc[0])}
                        ],

                    label = reference,
                    method = 'update'
                ) for j, reference in enumerate(unique_references)
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
        
        data = data[['sample','reference','section','sub_rate','sequence']]

        
    if section == 'ROI':

        data = study.get_df(
            sample=samples, 
            section=section, 
            family=family)
        
        data = data[['sample','reference','section','sub_rate','frame_shift_ROI','sequence']]
    
    if section == 'ROI':
        reference = data[data['sample'] == samples[0]].iloc[0]['sequence']
        data['sub_rate'] = data.apply(lambda row: int(row['frame_shift_ROI'])*[np.nan] + list(row['sub_rate']) + (len(reference) - int(row['frame_shift_ROI']) - len(row['sub_rate']))*[np.nan], axis=1)

    data_sample_0 = data[data['sample'] == samples[0]].reset_index(drop=True)
    data_sample_1 = data[data['sample'] == samples[1]].reset_index(drop=True)

    df = pd.DataFrame(
        index = data_sample_0['reference'],
        columns = data_sample_1['reference'], 
    )
        
    for _, row in data_sample_0.iterrows():
        for _, row2 in data_sample_1.iterrows():
            df.loc[row['reference'], row2['reference']] = custom_pearsonr(row['sub_rate'], row2['sub_rate'])
                                
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
                text_auto = '.2f',
                )
    
    return {'data': df, 'fig': fig}


def biological_replicates_histogram(study, bio_replicates_samples):
    data_out = pd.DataFrame()
    fig = go.Figure()

    for samples in bio_replicates_samples:
        data = {}
        for s in samples:
            data[s] = study.df[(study.df['sample'] == s) & (study.df['section'] == 'full')][['sample','reference','sub_rate']].set_index('reference')

        # merge per reference, use samples as suffix for columns
        data = pd.concat([data[s].rename(columns={'sub_rate': 'sub_rate_{}'.format(s)}) for s in samples], axis=1)
        data.drop(columns=['sample'], inplace=True)

        def pearson_score_average(list_of_replicates):
            """
            Calculate the average of the person correlation scores of all pairs of replicates.
            """
            list_of_replicates = list_of_replicates.dropna()
            scores = []
            for i in range(len(list_of_replicates)):
                for j in range(i+1, len(list_of_replicates)):
                    scores.append(custom_pearsonr(list_of_replicates[i], list_of_replicates[j]))
            if len(scores) == 0:
                return np.nan
            return np.mean(scores)

        data['pearson'] = data.apply(lambda x: pearson_score_average(x), axis=1)
        data_out = pd.concat([data_out, data], axis=1)
        
        # plot pearson as a histogram
        fig.add_trace(go.Histogram(x=data['pearson'], name='{}-{}'.format(samples[0], samples[1]), visible=False))
        
    fig.update_layout(
        title='Pearson correlation between replicates',
        xaxis_title='Pearson correlation',
        yaxis_title='Count',
    )

    # add a dropdown menu
    fig.update_layout(
        updatemenus=[
            dict(
                active=0,
                buttons=[
                    dict(
                        label= ' - '.join(samples)[:50],
                        method="update",
                        args = [
                            {'visible' : [True if i == j else False for i in range(len(bio_replicates_samples))]},
                            {'title' : 'Pearson correlation between replicates '+ ' - '.join(samples)}
                        ],
                    ) for j, samples in enumerate(bio_replicates_samples)
                    ],
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=1,
                xanchor="left",
                y=1.2,
                yanchor="top"
            )
        ]
    )
    return {'data': data_out, 'fig': fig}
        
# ---------------------------------------------------------------------------

## c / v - Proportional
# ---------------------------------------------------------------------------

def change_in_temp_mut_frac_vs_temperature(study, samples, reference):
    return __mut_frac_vs_sample_attr(study, samples, reference, 'temperature_k', 'Temperature (K)')
    
def change_in_reaction_time(study, samples, reference):
    return __mut_frac_vs_sample_attr(study, samples, reference, 'inc_time_tot_secs', 'DMS incubation time (s)')

def change_in_dms_conc(study, samples, reference):
    return __mut_frac_vs_sample_attr(study, samples, reference, 'DMS_conc_mM', 'DMS concentration (mM)')

def mut_rate_across_family_vs_deltaG(study, sample, family):
    
    # get a neat dataframe with the mutation rates for each base at each deltaG
    data = study.get_df(sample=sample, family=family, section='ROI', index_selected=True)
    
    data['deltaG'] = data['deltaG'].apply(lambda x: 0 if x == 'void' else float(x))
    
    # turn it into a dataframe
    df = pd.DataFrame(
        columns= [base + str(idx+1) for base, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])],
        data = [int(offset)*[np.nan] + list(mr) for offset, mr in zip(data['frame_shift_ROI'], data['sub_rate'])],
        index= data['deltaG'].values
    )
    
    # only keep the A and C bases
    idx_AC = [col[0] in ['A','C'] for col in df.columns]
    df = df.loc[:,idx_AC]
    
    # differentiate between paired and unpaired bases
    paired = [residue for idx, residue in enumerate([False if residue == '.' else True for residue in data['structure'].iloc[0]]) if idx_AC[idx]]
    aligned = data['num_aligned']
    
    fig = go.Figure()
    
    for col, pair in zip(df.columns, paired):
        fig.add_trace(go.Scatter(
            x = df.index,
            y = df[col],
            mode = 'markers+text',
            name = col,
            visible = 'legendonly',
            legendgroup=col,
            marker=dict(
                size = np.sqrt(data['num_aligned'].values)/np.mean(np.sqrt(data['num_aligned'].values))*20,
                line=dict(
                    color='black',
                    width=0 if pair else 3
                )
            ),
            text= [str(int(x)) for x in data['num_aligned'].values],
            textposition='bottom center',
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
            try:
                param_samples = np.clip(np.random.multivariate_normal(popt, pcov, N).T.reshape(3, N, 1), 0, np.inf)
            except:
                continue
            
            # Compute the sigmoid for each set of parameters for each x value
            y_MC = sigmoid(xdata_MC.reshape(1, -1) , param_samples[0], param_samples[1], param_samples[2])

            fig.add_trace(go.Scatter(x=xdata_MC, y=y_fit, mode='lines', name='fit',
                                    visible = 'legendonly', legendgroup=base, showlegend=False, 
                                    marker_color=colors[i_b]))

            fig.add_trace(go.Scatter(x=np.concatenate((xdata_MC, xdata_MC[::-1])), # x, then x reversed
                                    y= np.clip(np.concatenate((np.percentile(y_MC, 97.5, axis=0), np.percentile(y_MC, 2.5, axis=0)[::-1])), 0, 1), # upper, then lower reversed
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
    
    big_data = study.df[(study.df['sample'] == sample) & (study.df['section'] == 'ROI')][['sample','family', 'reference', 'frame_shift_ROI', 'sub_rate','sequence']].reset_index(drop=True)
    # reference is the longest sequence of the family
    for family in big_data['family'].unique():
        reference = big_data[(big_data['family'] == family)].sort_values('sequence', key=lambda x: x.str.len(), ascending=False)['sequence'].values[0]
        big_data.loc[big_data['family']==family, 'sub_rate'] = big_data[big_data['family']==family].apply(lambda x: np.array(int(x['frame_shift_ROI'])*[np.nan] + list(x['sub_rate']) + [np.nan]*(len(reference)-int(x['frame_shift_ROI'])-len(x['sub_rate']))), axis=1)
    
    # add each family as a plot but show only one at the time
    for family in big_data['family'].unique():
        
        data = big_data[(big_data['family'] == family)]
        
        data_all_samples = study.df[(study.df['family'] == family) & (study.df['section'] == 'ROI')]
        sequences_ROI = data_all_samples.sort_values('sequence', key=lambda x: x.str.len(), ascending=False)['sequence'].values
        reference = sequences_ROI[0]
                
        df = pd.DataFrame(
            columns = [base + str(idx + 1) for base, idx in zip(reference, range(len(reference)))],
            data = np.vstack(data['sub_rate'].values),
            index = data['reference'].values 
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
            hovertemplate = 'reference: %{y}<br>Base: %{x}<br>Mutation fraction: %{z:.2f}<extra></extra>',
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
    plt.ylabel('reference')
    """

def sub_rate_for_kfold(study, samples, family, stride = 'turner'):
    
    assert stride in ['turner', 'child#'], 'stride must be either "turner" or "child#"'
    
    if isinstance(samples, str):
        samples = [samples]
    
    df = {}
    for sample in samples:
        df[sample] = generate_dataset.select_data_for_kfold(study, sample, family, stride)
    
    df = pd.concat(df, axis=0).reset_index().rename(columns={'level_0': 'sample', 'level_1': 'deltaG'})
    
    # Make affine transformation to have the same scale for all samples
    affine_tranformations = generate_dataset.compute_affine_transformation_for_replicates(study, samples)
    
    for sample in samples[1:]:
        for col in df.columns[2:]:
            df.loc[df['sample']==sample, col] = affine_tranformations[sample](df[df['sample']==sample][col].values)

    fig = go.Figure()
    
    # Function to fit
    def sigmoid(x, a, b, c):
        RT = 1.987204258*310/1000
        return a / (1 + b*np.exp(-x/RT)) + c
    
    # Make a color generator for the different samples
    colors = {}
    for i, sample in enumerate(samples):
        colors[sample] = px.colors.qualitative.Plotly[i]

    visible_list = []

    for sample in samples:
        for base in df.columns[2:]:
            x, y = df[df['sample']==sample]['deltaG'].values, df[df['sample']==sample][base].values
            # remove nans
            x, y = x[~np.isnan(y)], y[~np.isnan(y)]
            
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode='markers',
                name=sample,
                marker=dict(size=10),
                line=dict(color=colors[sample], width=2),
                legendgroup=samples.index(sample),
                visible=False,
            ))
            visible_list.append(base)
            # fit the data
            try:
                popt, pcov = curve_fit(sigmoid, x, y, p0=[0.04, 0.02, 0.00], bounds=([0, 0, 0], [0.1, np.inf, 0.05]), max_nfev=1000)
            except:
                print('Could not fit data for sample {} and base {}'.format(sample, base))
                print(x,y)
                continue
            x_fit = np.linspace(x.min(), x.max(), 100)
            y_fit = sigmoid(x_fit, *popt)
            
            # plot the fit
            fig.add_trace(go.Scatter(
                x=x_fit,
                y=y_fit,
                mode='lines',
                name=sample,
                showlegend=False,
                line=dict(color=colors[sample], width=2),
                legendgroup=samples.index(sample),
                visible=False,
            ))
            visible_list.append(base)
            
            # draw a vertical line 

    # add title
    fig.update_layout(
        title_text="Mutation rates for bases used in kfold across samples",
        xaxis_title="deltaG",
        yaxis_title="mutation rate (corrected)",
        legend_title_text='Samples')

    # increase title and labels font size 
    fig.update_layout(title_font_size=20)
    fig.update_layout(xaxis_title_font_size=16)
    fig.update_layout(yaxis_title_font_size=16)

    # window y to 0-0.15
    fig.update_yaxes(range=[0, 0.15])
    
    # add a button to select a single base to display
    fig.update_layout(
        updatemenus=[
            go.layout.Updatemenu(
                active=0,
                buttons=list([
                    dict(label = base,
                         method = 'update',
                         args = [{'visible': [True if base in visible_list[i] else False for i in range(len(visible_list))]},
                                 {'title': 'Mutation rates for base {} of family {} across {}'.format(base, family, samples)}])
                    for j, base in enumerate(df.columns[2:])
                ]),
                direction = 'down',
                pad = {'r': 10, 't': 10},
                showactive = True,
                x = 0.1,
                xanchor = 'left',
                y = 1,
                yanchor = 'top'
                
            )
        ])

    # make the first base visible
    for i, b in enumerate(visible_list):
        fig.data[i].visible = b == df.columns[2]
        

    return {'data': df, 'fig': fig}

def kfold_per_family(study, sample, family, stride = 'turner'):
    
    assert stride in ['turner', 'child#'], 'stride must be either "turner" or "child#"'
    
    # get a neat dataframe with the mutation rates for each base at each deltaG
    df = generate_dataset.compute_k_fold_fit(study, sample, family, stride)
    
    # make a gaussian fit to the data
    x = np.linspace(1.2*min(df['Kfold']), 0.8*max(df['Kfold']), 100) if stride == 'turner' else np.linspace(0.8*min(df['Kfold']), 1.2*max(df['Kfold']), 100)
    y = norm.pdf(x, np.mean(df['Kfold']), np.std(df['Kfold']))
    
    fig = go.Figure()
    # legend the data "experimental data"
    fig.add_trace(
        go.Scatter(
            x=df['Kfold'],
            y = df['norm'],
            mode='markers+text',
            text=df.index,
            name='Experimental data',
            textposition='bottom center',
        )
    )
    fig.add_trace(
        go.Line(
            x=x,
            y=y,
            mode='lines',
            name='Gaussian fit',
            line=dict(color='red', width=2)
        )
    )

    # add peak at the mean and show it in the legend
    fig.add_shape(
        type="line",
        x0=np.mean(df['Kfold']),
        y0=0,
        x1=np.mean(df['Kfold']),
        y1=norm.pdf(np.mean(df['Kfold']), np.mean(df['Kfold']), np.std(df['Kfold'])),
        line=dict(
            color="red",
            width=2,
        ),
        name='Fit mean')
    
    # add the standard deviation
    fig.add_annotation(
        x=x[33],
        y=y[33],
        xref="x",
        yref="y",
        text="Std: {:.4f}".format(np.std(df['Kfold'])),
        showarrow=True,
        arrowhead=3,
        arrowcolor="red",
        ax=50,
        ay=0,
        font=dict(
            size=18
        )
    )

    # show the mean in the x axis
    fig.add_annotation(
        x=np.mean(df['Kfold']),
        y=0,
        xref="x",
        yref="y",
        text="Mean: {:.4f}".format(np.mean(df['Kfold'])),
        showarrow=True,
        arrowhead=3,
        arrowcolor="red",
        ax=0,
        ay=-40,
        font=dict(
            size=18
        )
    )


    fig.update_layout(
        title='Kfold distribution for sample {} and family {}'.format(sample, family),
        xaxis_title='Kfold' if stride == 'turner' else 'child number',
        yaxis_title='Probability density',
        legend_title='Legend',
        font=dict(
            size=18
        )
    )

    if stride == 'child#':
        df.rename(columns={'Kfold': 'child#'}, inplace=True)

    return {'fig': fig, 'data': df}


def kfold_per_samples(study, samples, family, stride='turner'):

    assert stride in ['turner', 'child#'], 'stride must be either "turner" or "child#"'

    data = {}
    for sample in samples:
        data[sample] = {}
        for family in study.df[study.df['sample']==sample].family.unique():
            df = generate_dataset.compute_k_fold_fit(study, sample, family, stride)
            data[sample][family] = df['Kfold'].mean()
    
    df = pd.DataFrame(data)

    # build a barplot of the Kfold values
    fig = go.Figure()
    for sample in samples:
        fig.add_trace(go.Bar(x=df.index, y=df[sample], name=sample))
    fig.update_layout(
        title='Kfold Turner-based estimation for {}'.format(' and '.join(samples)) if stride == 'turner' else 'Kfold child-number-based estimation for {}'.format(' and '.join(samples)),
        xaxis_title='Family',
        yaxis_title='Kfold' if stride == 'turner' else 'child number',
        bargap=0.2,
        bargroupgap=0.1
    )

    # reverse the y axis
    if stride == 'turner':
        fig.update_yaxes(autorange="reversed")
        
    elif stride == 'child#':
        df.rename(columns={'Kfold': 'child#'}, inplace=True)

    return {'fig': fig, 'data': df}


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
        pearson = custom_pearsonr(x,y)
        
        # compute R2
        r2 =  r2_score(y, x)
        
        # Plot regression line
        fig.add_trace(
            go.Scatter(x=[xmin, xmax], y=[xmin*slope + intercept, xmax*slope + intercept],
            mode='lines', name='LR: {:.2f}x + {:.2f}'.format(slope, intercept), line=dict(color='red'), 
            showlegend=True, visible=False))

        
def __mut_frac_vs_sample_attr(study, samples, reference, attr, x_label=None):
    # get data
    data = study.get_df(
        sample = samples,
        reference=reference,
        section='ROI',
        base_type = ['A','C']) 

    data = data[['sample','sequence','sub_rate','index_selected','structure','num_aligned',attr]]
    
    if len(data) == 0:
        print('No data: samples:{}, reference:{}'.format(samples, reference))
        return {'data': [], 'fig': go.Figure()}
    
    # turn it into a dataframe
    df = pd.DataFrame(
        columns= [base + str(idx+1) for base, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])],
        data = np.array(data['sub_rate'].tolist()),
        index= data[attr].values
    ).sort_index().rename_axis(index=attr)
    
    # plot
    paired = [False if residue == '.' else True for residue in data['structure'].iloc[0]]
    aligned = data['num_aligned'].values
    
    
    fig = go.Figure()
    
    for col, pair in zip(df.columns, paired):
        fig.add_trace(go.Scatter(
            x = df.index,
            y = df[col],
            mode = 'markers+text',
            name = col,
            visible = 'legendonly',
            marker=dict(
                size = np.sqrt(data['num_aligned'].values)/np.mean(np.sqrt(data['num_aligned'].values))*20,
                line=dict(
                    color='black',
                    width=0 if pair else 3
                )
            ),
            text = aligned,
            textposition = 'bottom center',
    ))
        
    fig.update_layout(
        title='reference: {}'.format(reference),
        xaxis_title=x_label if x_label else attr,
        yaxis_title='Mutation fraction')

    
    return {'data': df, 'fig': fig}

def __corr_scatter_plot(x, y, visible=True):

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
            ),
            visible=visible
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
            ),
            visible=visible
        )
    )
    
    # add the pearson correlation score as a text
    traces.append(
        go.Scatter(
            x = [np.mean(x)],
            y = [np.max(y)],
            mode = 'text',
            text = ['Pearson correlation score: {:.2f}'.format(custom_pearsonr(x,y))],
            textposition = 'top center',
            showlegend=False,
            visible=visible
        ))
    
    # add the R2 score as a text
    traces.append(
        go.Scatter(
            x = [np.mean(x)],
            y = [np.max(y)],
            mode = 'text',
            text = ['R2 score: {:.2f}'.format(r2)],
            textposition = 'bottom center',
            showlegend=False,
            visible=visible
        ))
    
    return traces
