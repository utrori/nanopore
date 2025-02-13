import matplotlib.pyplot as plt
import utilities
import config
import os
import numpy as np
import methylation_analysis
import fast5_loader

plt.rcParams['font.family']= 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'


def signle_read_plot_structure(f5reader: object, split_align_res: np.array, savedir, filename="", mouse='no', met='yes', pdf='no', adenine='no'):
    """Plot read structure."""
    if not os.path.exists(savedir):
        os.mkdir(savedir)

    loadeddata = f5reader.read_fast5()
    read_id = loadeddata['read_id']
    seq = loadeddata['sequence']
    quality = loadeddata['quality_scores']
    methylation = loadeddata['methylation']
    fig = plt.figure()
    plt.subplots_adjust(left=0.2)
    ax = fig.add_subplot()
    met_window = 200
    if not filename:
        filename = read_id

    #read structure plot
    lc = utilities.get_read_structure(read_id, config.SPLIT_LENGTH, split_align_res, mouse=mouse)
    ax.add_collection(lc)

    #quality plot
    bin1 = 200
    x = []
    y = []
    for n in range(0, len(seq) - bin1, bin1):
        x.append(n)
        y.append(np.mean(quality[n:n+bin1]) * 200 - 10000)
    ax.plot(x, y)

    #methylation plot
    if met == 'yes':
        x, y = methylation_analysis.cpg_methylation_binning(met_window, seq, methylation, method='threshold')
        ax.bar(x, y, width=met_window, color = 'magenta', zorder=0)
    elif met == 'rerio':
        x, y = methylation_analysis.cpg_methylation_binning(met_window, seq, methylation, method='threshold')
        ax.bar(x, y, width=met_window, color = 'magenta', zorder=0, alpha=0.3)
    if adenine == 'yes':
        x, y = methylation_analysis.cpg_methylation_binning(met_window, seq, methylation, method='threshold')
        ax.bar(x, y, width=met_window, color = 'blue', zorder=0, alpha=0.3)

    if mouse == 'yes':
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000, 50000))
        ax.set_yticklabels(('UM', '0', '10k', '20k', '30k', '40k', '50k'), fontsize=20)
    else:
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(('UM', '0', '10k', '20k', '30k', '40k'), fontsize=20)
    ax.set_xticks((0, 10000, 20000, 30000, 40000, 50000))
    ax.set_xticklabels(('0', '10k', '20k', '30k', '40k', '50k'), fontsize=20)
    ax.set_title('Mean quality: {:.1f}'.format(np.mean(quality)), fontsize=15)
    ax.autoscale()
    ax.set_ylim([-12000, 47000])
    if savedir[-1] == '/':
        savedir = savedir[:-1]
    if pdf == 'yes':
        plt.savefig(f'{savedir}/{read_id}.pdf')
    else:
        plt.savefig(f'{savedir}/{read_id}.png', dpi=300)
    plt.close()