from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sys.path.append("../config/")
from config_plot import config_plot_params


def return_clon_dict(df_ml, sig_cols, treatment_wanted, signature_treatment):
    dic_clon = defaultdict(list)

    for sample, data in df_ml.groupby(by='SAMPLE'):

        clonal = data[data['Clonal1'] == 'clonal']
        subclonal = data[data['Clonal1'] == 'subclonal']

        for signature in sig_cols:
            if (clonal[signature].sum() > 0) & (subclonal[signature].sum() > 0):

                clonal_v = clonal[signature].sum() / len(clonal)
                subclonal_v = subclonal[signature].sum() / len(subclonal)
                foldchange = clonal_v / subclonal_v

                if signature == signature_treatment:
                    if data['PATIENT'].tolist()[0] in treatment_wanted:
                        dic_clon[signature].append(foldchange)
                else:
                    dic_clon[signature].append(foldchange)
    return dic_clon


def return_clon_dict_meta(df_ml, sig_cols, treatment_wanted, signature_treatment):
    dic_clon = defaultdict(list)

    for sample, data in df_ml.groupby(by='SAMPLE'):

        # The clonality is defined by MutationalTime script
        clonal = data[data['timing_class'] != 'subclonal']
        subclonal = data[data['timing_class'] == 'subclonal']

        for signature in sig_cols:

            sub_data_clonal = clonal[clonal['ML'] == signature]
            sub_data_subclonal = subclonal[subclonal['ML'] == signature]

            if (len(sub_data_clonal) > 50) & (len(sub_data_subclonal) > 50):

                clonal_v = len(sub_data_clonal) / len(clonal)
                subclonal_v = len(sub_data_subclonal) / len(subclonal)
                foldchange = clonal_v / subclonal_v

                if signature == signature_treatment:
                    if data['SAMPLE'].tolist()[0] in treatment_wanted:
                        dic_clon[signature].append(foldchange)
                else:
                    dic_clon[signature].append(foldchange)
    return dic_clon


def load_clonality_treatment(file_ML, sig_treated, all_treat):
    df = pd.read_csv(file_ML, sep='\t')
    sig_cols = df['ML'].unique().tolist()
    dic_clon = return_clon_dict_meta(df, sig_cols, all_treat, sig_treated)

    return dic_clon[sig_treated]


def boxplot_clonality(healthy, ovary_t, urinary_t, outpath):
    toplot = [healthy, ovary_t, urinary_t]
    labels = ['t-AML', 'Ovary', 'Urinary-tract']

    color_ttype = {
        't-AML': 'brown',
        'Ovary': 'orange',
        'Urinary-tract': 'darkgreen',
    }

    np.random.seed(12345)
    config_plot_params(5)

    fig, ax = plt.subplots(1, 1, figsize=(1, 1.3))

    ax = sns.boxplot(data=toplot, showfliers=False, color='white',
                     medianprops={'color': 'black'}, linewidth=0.4, ax=ax)

    for ix, l in enumerate(toplot):
        sign = labels[ix]
        ax.scatter([ix + np.random.uniform(-0.2, 0.2, 1)[0] for i in range(len(l))],
                   l, s=1, alpha=0.75, color=color_ttype[sign])

    plt.ylabel('foldchange clonal/subclonal')
    ax.set_xticklabels(['{}'.format(l) for ixl, l in enumerate(labels)], rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-0.5, len(toplot) + 0.5)

    ax.hlines(1, 0, len(toplot), alpha=0.5)

    plt.savefig('{}/clonal_subclonal_sigs_all_merged.svg'.format(outpath))
    plt.show()
