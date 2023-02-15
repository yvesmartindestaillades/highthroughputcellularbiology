import util, plots
import io
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from study_gen import study

matplotlib.use('agg')

def get_avg_pearson_score_for_construct_in_family(study, sample, construct, section):
    family = util.family_from_construct(construct)
    data = plots.sample_replicates_heatmap_per_family(study, [sample,sample], family, section)['data']
    score = data[construct].loc[~data.index.str.contains(construct)].mean()
    return score

DEF_TEMP_FILE = 'temp.png'
def read_and_return_plot(func):
    def read_image(*args):
        func(*args)
        im = plt.imread(DEF_TEMP_FILE)
        os.remove(DEF_TEMP_FILE)
        return im
    return read_image

@read_and_return_plot
def convert_matplotlib_to_image():
    plt.gcf()
    plt.savefig(DEF_TEMP_FILE, bbox_inches='tight', dpi=300)


@read_and_return_plot
def convert_pyplot_to_image(fig):
    fig.write_image(DEF_TEMP_FILE, scale=2)


def nreads_color_factory(n_reads):
    # 1000 - 2000 reads: dark red
    # 2000 - 3000 reads: light red
    # 3000 - 4000 reads: orange
    # 4000 - 5000 reads: yellow
    # 5000 - 6000 reads: light green
    # 6000+ reads: dark green
    
    if n_reads < 2000:
        return 'darkred'
    elif n_reads < 3000:
        return 'red'
    elif n_reads < 4000:
        return 'orange'
    elif n_reads < 5000:
        return 'yellow'
    elif n_reads < 6000:
        return 'lightgreen'
    else:
        return 'darkgreen'

def pearson_color_factory(pearson):
    # 0.0 - 0.3: dark red
    # 0.3 - 0.6: light red
    # 0.6 - 0.7: orange
    # 0.7 - 0.8: yellow
    # 0.8 - 0.9: light green
    # 0.9 - 1.0: dark green

    if pearson < 0.3:
        return 'darkred'
    elif pearson < 0.6:
        return 'red'
    elif pearson < 0.7:
        return 'orange'
    elif pearson < 0.8:
        return 'yellow'
    elif pearson < 0.9:
        return 'lightgreen'
    else:
        return 'darkgreen'
    

def generate_one_pager_construct(study, sample, construct, output_dir):
    
    family = util.family_from_construct(construct)

    # Get data
    data = study.get_df(sample=sample, section='ROI', construct=construct)
    assert len(data) <= 1, 'More than one sequence found for {} - {}'.format(sample, construct)
    assert len(data) > 0, 'No sequence found for {} - {}'.format(sample, construct)
    data = data.iloc[0]

    # Construct info
    construct = construct
    sequence = data['sequence']
    type = 'Canonical base pair'
    gc_content = '{:.2f}%'.format(100*np.array([base in ['G','C'] for base in sequence]).mean())
    deltaG = data['deltaG']
    print('DeltaG: {}'.format(deltaG))
    print('Construct info \n--------------')
    print('Construct: {}'.format(construct))
    print('Sequence: {}'.format(sequence))
    print('Type: {}'.format(type))
    print('GC content: {}'.format(gc_content))
    print('')

    # Sample info
    library = 'TODO' # should be in samples.csv
    exp_env = data['exp_env']
    exp_env_name = 'cell_line' if exp_env == 'in_vivo' else 'buffer'
    exp_env_var = data['cell_line'] if exp_env == 'in_vivo' else data['buffer']
    DMS = data['DMS_conc_mM']
    reaction_time = data['inc_time_tot_secs']
    temperature = data['temperature_k']

    print('Sample info \n-----------')
    print('Library: {}'.format(library))
    print('Cell line: {}'.format(exp_env_var)) if exp_env == 'in_vivo' else print('Buffer: {}'.format(exp_env_var))
    print('DMS concentration: {} mM'.format(DMS))
    print('Reaction time: {} secs'.format(reaction_time))
    print('Temperature: {} K'.format(temperature))
    print('')

    # Quality control
    num_reads = data['num_aligned']
    pearson_R_5_hp = get_avg_pearson_score_for_construct_in_family(study, sample, construct, 'MS2')
    pearson_R_3_hp = get_avg_pearson_score_for_construct_in_family(study, sample, construct, 'LAH')
    print('Quality control \n---------------')
    print('Number of reads: {}'.format(num_reads))
    print('Pearson R 5 hp: {}'.format(pearson_R_5_hp))
    print('Pearson R 3 hp: {}'.format(pearson_R_3_hp))
    print('')

    ## Plots
    figs = {}
    images = {}
    figs['read_coverage_per_position'] = plots.read_coverage_per_position(study, sample, construct)['fig']
    figs['mutation_fraction_at_each_position'] = plots.mutation_fraction_at_each_position(study, sample, construct)['fig']
    figs['mutation_identity_at_each_position'] = plots.mutation_identity_at_each_position(study, sample, construct)['fig']
    figs['mutation_per_read_per_construct'] = plots.mutation_per_read_per_construct(study, sample, construct)['fig']
    for fig_name, fig in figs.items():
        if fig_name != 'mutation_per_read_per_construct':
            fig['layout']['annotations'][0]['font']['size'] = 35
        else:
            fig['layout']['title']['font']['size'] = 35
            fig['layout']['title']['text'] = '# of mutations per read'
        images[fig_name] = convert_pyplot_to_image(fig)


    plot_to_letter = {
        'read_coverage_per_position': 'A',
        'mutation_fraction_at_each_position': 'D',
        'mutation_identity_at_each_position': 'C',
        'mutation_per_read_per_construct': 'B'}


    tab = plt.table(
        cellText=[['Construct', construct],
                    ['Sequence', sequence],
                    ['Type', type],
                    ['GC content', gc_content],
                    ['ΔG', str(deltaG) + ' kcal/mol'],
                    ],
        cellLoc='left',
        loc='bottom',
        bbox=[0, 0, 1, 1],
        colWidths=[0.25, 0.75],
        )

    tab.auto_set_font_size(False)
    tab.set_fontsize(12)

    # make the left column bold
    for (row, col), cell in tab.get_celld().items():
        if (row == 0) or (col == 0):
            cell.set_text_props(weight='bold')    

    plt.title('Construct info', fontsize=20, fontweight='bold')
    plt.axis('off')
    plot_to_letter['construct_info'] = 'K'
    images['construct_info'] = convert_matplotlib_to_image()

    tab = plt.table(
        cellText=[['Library', library],
                    [exp_env_name, exp_env_var] if exp_env == 'in_vivo' else ['Buffer', exp_env_var],
                    ['%DMS', str(DMS)+' mM'],
                    ['Reaction time', str(reaction_time)+' secs'],
                    ['Temperature', str(temperature)+' K (' + str(temperature-273.15) + ' °C)' ],
                    ],

        cellLoc='left',
        loc='bottom',
        bbox=[0, 0, 1, 1],
        colWidths=[0.3, 0.7],
        )

    tab.auto_set_font_size(False)
    tab.set_fontsize(12)

    # make the left column bold
    for (row, col), cell in tab.get_celld().items():
        if col == 0:
            cell.set_text_props(weight='bold')

    plt.title('Sample info', fontsize=20, fontweight='bold')
    plt.axis('off')
    plot_to_letter['sample_info'] = 'M'
    images['sample_info'] = convert_matplotlib_to_image()
    
    # plot quality control
    # plot a balck circle with the score written in white inside

    # 1st row twice larger than the second row
    fig, axs = plt.subplots(3,3, figsize=(7.5,3), gridspec_kw={'height_ratios': [1, 10, 1], 'width_ratios': [1, 1, 1]})
    # plot a circle



    axs[0,1].text(0.5, 0.5, 'Quality Control', horizontalalignment='center', verticalalignment='center', fontsize=20, fontweight='bold')
    axs[1,0].add_patch(plt.Circle((0.5, 0.5), 0.45, color=nreads_color_factory(num_reads)))
    axs[1,0].text(0.5, 0.5, '{:.1f}k'.format(num_reads/1000), horizontalalignment='center', verticalalignment='center', fontsize=30, fontweight='bold')
    axs[2,0].text(0.5, 0.5, '# reads', horizontalalignment='center', verticalalignment='center', fontsize=20)
    axs[1,1].add_patch(plt.Circle((0.5, 0.5), 0.45, color=pearson_color_factory(pearson_R_5_hp)))
    axs[1,1].text(0.5, 0.5, '{:.2f}'.format(pearson_R_5_hp), horizontalalignment='center', verticalalignment='center', fontsize=30, fontweight='bold')
    axs[2,1].text(0.5, 0.5, 'Pearson R\n5\' hp', horizontalalignment='center', verticalalignment='center', fontsize=20)
    axs[1,2].add_patch(plt.Circle((0.5, 0.5), 0.45, color=pearson_color_factory(pearson_R_3_hp)))
    axs[1,2].text(0.5, 0.5, '{:.2f}'.format(pearson_R_3_hp), horizontalalignment='center', verticalalignment='center', fontsize=30, fontweight='bold')
    axs[2,2].text(0.5, 0.5, 'Pearson R\n3\' hp', horizontalalignment='center', verticalalignment='center', fontsize=20)

    for ax in axs.flat:
        ax.axis('off')

    # save the figure

    plot_to_letter['quality_control'] = 'L'
    images['quality_control'] = convert_matplotlib_to_image()

    axd = plt.figure(constrained_layout=True, figsize=(25,40)).subplot_mosaic(
        """
        KKLL
        MMBB
        AAAA
        CCCC
        DDDD
        """
    )    

    for (figname, fig) in images.items():
        axd[plot_to_letter[figname]].imshow(images[figname])
        axd[plot_to_letter[figname]].axis('off')
        plt.tight_layout()
    plt.show()
    plt.tight_layout()
    out_file = os.path.join(output_dir, '{}_{}.pdf'.format(construct, sample))
    if os.path.isfile(out_file):
        os.system('rm {}'.format(out_file))
    plt.savefig(out_file, bbox_inches='tight', pad_inches=0)

if __name__ == '__main__':
    
        
    ###############
    # Parameters
    ###############
    sample = 'lauren470_S1'
    construct = '3042-O-flank_1=hp1-DB'
    output_dir = ''
    ###############

    generate_one_pager_construct(
        study=study,
        sample=sample,
        construct=construct,
        output_dir=output_dir
    )