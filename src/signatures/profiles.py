import matplotlib.pyplot as plt
import sys

sys.path.append("../config/")
from config_plot import config_plot_params


def order_muts():
    order = []
    first = ['A', 'C', 'G', 'T']
    pyr = ['C', 'T']
    for p in pyr:
        for mut in first:
            if mut != p:
                for f in first:
                    for f2 in first:
                        comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                        order.append(comb)
    return order


# Split into even chunks
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


# Function to plot the SNV processes
def plot_snvs(sig, title, outpath):
    config_plot_params(3)
    fig, axs = plt.subplots(
        nrows=2, ncols=1, figsize=(3.2, 1), gridspec_kw={'height_ratios': [1, 9]}
    )
    order_plot = order_muts()

    vals = []
    colors = []
    colors_mut = [
        '#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'
    ]
    bot = -0.5
    for ix, c in enumerate(chunks(sig, 16)):
        colors.extend([colors_mut[ix] for _ in c])
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

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(
        ['{}{}{}'.format(a[0], a[2], a[-1]) for a in order_plot],
        verticalalignment="center", ha='center', rotation=90, fontsize=2,
        color='grey'
    )
    axs[1].set_ylabel('Relative Probability')

    plt.tight_layout()
    plt.xlim(-1, 96)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    prev_pos = 6
    for count_lab, lab in enumerate(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']):
        color = 'black'
        if count_lab == 1:
            color = 'white'
        axs[0].text(prev_pos, 0.85, lab, color=color, weight='bold')
        prev_pos = prev_pos + 16

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    plt.suptitle(title, y=1.05)

    plt.savefig(outpath + "/plot_snvs.svg")

    plt.show()
