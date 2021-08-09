
from collections import Counter
from order import order_muts, order_to_plot_dbs, order_to_plot_indel
import matplotlib.pyplot as plt

# split into even chunks
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

        
# function to plot the SNV processes
def plot_snvs(sig, title, outpath):

    config_params(12)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 3.3),
                            gridspec_kw={'height_ratios': [1, 9]})
    order_plot = order_muts("snv")

    vals = []
    colors = []
    colors_mut = ['#1ebff0', '#050708', '#e62725',
                  '#cbcacb', '#a1cf64', '#edc8c5']
    bot = -0.5
    for ix, c in enumerate(chunks(sig, 16)):
        colors.extend([colors_mut[ix] for s in c])
        axs[0].barh(1, 16, left=bot, color=colors_mut[ix])
        bot += 16
        vals.extend(c)

    axs[0].set_xlim(-1, 96)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=96, lw = 0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.1, xmin=-1, xmax=96, lw = 0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.15, xmin=-1, xmax=96, lw = 0.6, color='grey', alpha=0.2)

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(['{}{}{}'.format(a[0], a[2], a[-1]) for a in order_plot],
                           verticalalignment="center", ha='center', rotation=90, fontsize=12,
                           color='grey')

    plt.tight_layout()
    plt.xlim(-1, 96)

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Mutation count', fontsize = 12)


    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)
        
    prev_pos = 6
    for count_lab, lab in enumerate(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']):
        color = 'black'
        if count_lab == 1:
            color = 'white'
        axs[0].text(prev_pos, 0.85, lab, color=color, weight='bold')
        prev_pos = prev_pos + 16
        
    axs[1].xaxis.set_tick_params(pad=10)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.title(title)
    plt.savefig('{}/{}_snvs.svg'.format(outpath, title))
    plt.savefig('{}/{}_snvs.png'.format(outpath, title), dpi = 600)

    plt.show()

    
# function to plot the DBS processes
def plot_dbs(sig, title, outpath):
    
    config_params(12)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 3.3),
                            gridspec_kw={'height_ratios': [1, 9]})

    vals = []
    colors_mut = ['#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5']

    dbs_color = {
        'AC': '#a6cee3', 'AT': '#1f78b4', 'CC': '#b2df8a', 'CG': '#33a02c', 'CT': '#fb9a99',
        'GC': '#e3211d', 'TA': '#fdbf6f', 'TC': '#ff7f00', 'TG': '#cab2d6',
        'TT': '#6a3d9a',
    }

    order_color = [
        '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
        '#e3211d', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'
    ]

    order_dbs_list = order_to_plot_dbs()

    vals = sig

    colors = [dbs_color[db.split('_')[0]] for db in order_dbs_list]
    counter_colors = Counter(colors)

    bot = -0.5

    for c in order_color:
        axs[0].barh(1, counter_colors[c], left=bot, color=c, align='center',)
        bot += counter_colors[c]

    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    axs[0].set_xlim(-1, 78)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=96, lw=0.3, color='grey', alpha=0.2)
    axs[1].axhline(y=0.1, xmin=-1, xmax=96, lw=0.3, color='grey', alpha=0.2)
    axs[1].axhline(y=0.15, xmin=-1, xmax=96, lw=0.3, color='grey', alpha=0.2)

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center', alpha=1)
    axs[1].set_xticks(x)

    axs[1].set_xticklabels(['{}{}'.format(a[-2], a[-1],) for a in order_dbs_list],
                           rotation=90, fontsize=12, verticalalignment="center", ha='center',
                           color='grey')

    axs[1].set_xlim(-1, 78)

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Mutation count')

    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=10)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    plt.tight_layout()
    plt.savefig('{}/{}_dbs.svg'.format(outpath, title))
    plt.savefig('{}/{}_dbs.png'.format(outpath, title), dpi = 600)
    plt.close()

def plot_indel(sig, title, outpath):

    config_params(12)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(20, 3.3),
                            gridspec_kw={'height_ratios': [1, 9]})

    vals = sig

    colors = ['#fdbe6f'] * 6 + ['#ff8002'] * 6 + ['#b0dd8b'] * 6 + ['#36a12e'] * 6 + \
             ['#fdcab5'] * 6 + ['#fc8a6a'] * 6 + ['#f14432'] * 6 + ['#bc191a'] * 6 + \
             ['#d0e1f2'] * 6 + ['#94c4df'] * 6 + ['#4a98c9'] * 6 + ['#1764ab'] * 6 + \
             ['#e1e1ef'] * 1 + ['#b6b6d8'] * 2 + ['#8683bd'] * 3 + ['#62409b'] * 5

    order_colors = ['#fdbe6f', '#ff8002', '#b0dd8b', '#36a12e', '#fdcab5', '#fc8a6a', '#f14432',
                    '#bc191a', '#d0e1f2', '#94c4df', '#4a98c9', '#1764ab', '#e1e1ef', '#b6b6d8', '#8683bd',
                    '#62409b']

    counter_colors = Counter(colors)

    bot = -0.5

    for c in order_colors:
        axs[0].barh(1, counter_colors[c], left=bot, color=c)
        bot += counter_colors[c]

    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    axs[0].set_xlim(-1, 83)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=83, lw=0.3, color='grey', alpha=0.3)
    axs[1].axhline(y=0.1, xmin=-1, xmax=83, lw=0.3, color='grey', alpha=0.3)
    axs[1].axhline(y=0.15, xmin=-1, xmax=83, lw=0.3,  color='grey', alpha=0.3)
    axs[1].bar(x, vals, color=colors, width=0.7, linewidth=0, align='center', alpha=1)
    axs[1].set_xticks(x)
    axs[1].set_xticklabels([i.split('_')[-1] for i in order_to_plot_indel()], fontsize=12,
                           verticalalignment="center", ha='center', color='grey')
    axs[1].set_xlim(-1, 83)

    plt.tight_layout()

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Mutation count')

    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=10)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    plt.savefig('{}/{}_indels.svg'.format(outpath, title))
    plt.savefig('{}/{}_indels.png'.format(outpath, title), dpi = 600)
    plt.close()