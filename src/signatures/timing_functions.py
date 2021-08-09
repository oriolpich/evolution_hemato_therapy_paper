import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import mannwhitneyu, fisher_exact
from statannot import add_stat_annotation
from collections import Counter
from biopsy import return_biopsies
import sys

sys.path.append("../config")
from config import HARTWIG_METADATA, FIG_2
from config_plot import config_plot_params

ttype_conversion = {
    "Colon-Rectum": "CR",
    "Ovary": "O",
    "Urinary-tract": "UT",
    "Breast": "BRCA"
}


def plot_agreement(agree_exposed, agree_not_exposed, ttype, drug):
    config_plot_params(6)

    fig, ax = plt.subplots(1, 1, figsize=(0.25, 1))
    total = len(agree_exposed) + len(agree_not_exposed)
    bot = len(agree_exposed) / total
    plt.bar(0, bot, width=0.8, color='#F3B53C')
    plt.bar(0, 1 - bot, bottom=bot, width=0.8, color='#173E57')
    plt.ylabel('Proportion of\ntreated samples')
    plt.xticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlim(-0.8, 0.8)
    plt.xticks([0], [f'{drug} {ttype_conversion[ttype]}\n({total})'], rotation=90)

    plt.savefig('../figures/{}/{}_{}.agreement.svg'.format(FIG_2, ttype, drug), va='top', ha='right')
    plt.show()


def create_subtypes(treatment_sample_level, ttype, drug, msigdf, nmfdf, qval, sig_found=True, sig_fix=None):
    if drug == 'Capecitabine':
        ttype_drug_treat1 = treatment_sample_level[ttype][drug]
        ttype_drug_treat2 = treatment_sample_level[ttype].get('Fluorouracil', set())
        ttype_drug_treat3 = treatment_sample_level[ttype].get('Tegafur', set())
        treated = set(ttype_drug_treat1) | set(ttype_drug_treat2) | set(ttype_drug_treat3)
    if drug == 'Platinum':
        ttype_drug_treat1 = treatment_sample_level[ttype]['Carboplatin']
        ttype_drug_treat2 = treatment_sample_level[ttype]['Cisplatin']
        treated = set(ttype_drug_treat1) | set(ttype_drug_treat2)
    else:
        treated = treatment_sample_level[ttype][drug]

    mSigAct_sig = set(msigdf[msigdf['apval'] < qval].index.tolist())
    mSigAct_notsig = set(msigdf[msigdf['apval'] >= qval].index.tolist())

    treat_msig_sig = [c for c in mSigAct_sig if c in treated]
    treat_msig_nosig = [c for c in mSigAct_notsig if c in treated]

    if sig_found:
        nmfdfcol = nmfdf.columns
        nmfdf.columns = [x.replace('-', '.') for x in nmfdfcol]
        sig_col = msigdf.columns.tolist()[-1]
        if sig_fix is None:
            if sig_col[0] == 'X':
                sig_col = sig_col[1:]

        high_exposed = set(nmfdf[nmfdf[sig_col] > 0].index.tolist())
        not_exposed = set(nmfdf[nmfdf[sig_col] == 0].index.tolist())
        treat_nmf_sig = [c for c in high_exposed if c in treated]
        treat_nmf_nosig = [c for c in not_exposed if c in treated]
        agree_exposed = list(set(treat_msig_sig) & set(treat_nmf_sig))
        agree_not_exposed = list(set(treat_msig_nosig) & set(treat_nmf_nosig))
    else:
        agree_exposed = treat_msig_sig
        agree_not_exposed = treat_msig_nosig

    plot_agreement(agree_exposed, agree_not_exposed, ttype, drug)

    return agree_exposed, agree_not_exposed


def classify_exposition(drug, treatment_days, treatment, agree_exposed, agree_not_exposed):
    exposed = [np.max(treatment[s][drug]) for s in agree_exposed if 'unknown' not in treatment_days[s][drug]]
    not_exposed = [np.max(treatment[s][drug]) for s in agree_not_exposed if 'unknown' not in treatment_days[s][drug]]
    return exposed, not_exposed


def get_timing_stats(agree_exposed, agree_not_exposed, drug, treatment_sample_days_treatment,
                     treatment_sample_days_biopsy_since_start_treatment,
                     treatment_sample_days_biopsy_since_end_treatment):
    # Days of treatment
    exposed_days = [np.sum(treatment_sample_days_treatment[s][drug]) for s in agree_exposed if
                    'unknown' not in treatment_sample_days_treatment[s][drug]]
    notexposed_days = [np.sum(treatment_sample_days_treatment[s][drug]) for s in agree_not_exposed if
                       'unknown' not in treatment_sample_days_treatment[s][drug]]

    # Days since start treatment
    exposed_days_since_start, notexposed_days_since_start = classify_exposition(
        drug, treatment_sample_days_treatment, treatment_sample_days_biopsy_since_start_treatment, agree_exposed,
        agree_not_exposed)

    # Days since end treatment

    exposed_days_since_end, notexposed_days_since_end = classify_exposition(
        drug, treatment_sample_days_treatment, treatment_sample_days_biopsy_since_end_treatment, agree_exposed,
        agree_not_exposed)

    return exposed_days, notexposed_days, exposed_days_since_start, notexposed_days_since_start, exposed_days_since_end, notexposed_days_since_end


def boxplot_two_list(l1, l2, ax, label, alternative, cut=1000):
    sns.violinplot(data=[l1, l2], palette=['white', 'white'],
                   ax=ax, showfliers=False, alpha=0.7, zorder=0)

    sns.boxplot(data=[l1, l2], palette=['#173E57', '#F3B53C'],
                ax=ax, showfliers=False, width=0.8, zorder=10,
                boxprops={'zorder': 10, 'alpha': 0.3},
                whiskerprops={"zorder": 10})

    ax.scatter([0 + np.random.uniform(-0.2, 0.2, 1)[0] for i in range(len(l1))],
               [v if v < cut else cut for v in l1], color='#173E57', s=1, alpha=0.8)

    ax.scatter([1 + np.random.uniform(-0.2, 0.2, 1)[0] for i in range(len(l2))],
               [v if v < cut else cut for v in l2], color='#F3B53C', s=1, alpha=0.8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    stat, pval = mannwhitneyu(l1, l2, alternative=alternative)
    ax.set_title(pval)
    ax.set_xticklabels(['Treated unexposed\n({})'.format(len(l1)), 'Treated exposed\n({})'.format(len(l2)),
                        ], rotation=90)
    ax.set_ylabel(label)


def plot_timing_distributions(exposed_days, notexposed_days, exposed_days_since_end, notexposed_days_since_end,
                              drug, ttype, cut1=500, cut2=4000, alternative='two-sided'):
    config_plot_params(font_size=7)
    fig, ax = plt.subplots(1, 2, figsize=(1.25, 1))

    boxplot_two_list(notexposed_days, exposed_days, ax[0], 'Days of treatment', alternative, cut=cut1)
    boxplot_two_list(notexposed_days_since_end, exposed_days_since_end, ax[1], 'Days since end treatment',
                     alternative, cut=cut2)
    plt.suptitle(f'{ttype}-{drug}')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    ax[0].set_ylim(0, cut1 + 100)
    ax[1].set_ylim(0, cut2 + 100)

    plt.savefig(f'../figures/{FIG_2}/{drug}_{ttype}_timing.svg')

    plt.show()


def distant_seeding(agree_exposed, agree_not_exposed, drug, ttype, distant, close, dic_prop_ttypes,
                    treatment_sample_days_treatment, treatment_sample_days_biopsy_since_end_treatment, toregs):
    metadata = pd.read_csv(HARTWIG_METADATA, sep='\t', encoding='latin-1')

    bio_dic = return_biopsies()
    metadata['BIOPSY'] = metadata['biopsySite'].apply(lambda x: bio_dic.get(x, x))
    dict_biopsy = dict(zip(metadata['sampleId'], metadata['BIOPSY']))

    v1 = Counter([dict_biopsy[s] for s in agree_exposed])
    v2 = Counter([dict_biopsy[s] for s in agree_not_exposed])

    config_plot_params(7)

    fig, ax = plt.subplots(1, 1, figsize=(0.5, 1))

    agree_out = [evaluate_sample_time(s, drug, close, distant, dict_biopsy, 1, treatment_sample_days_treatment,
                                      treatment_sample_days_biopsy_since_end_treatment) for s in agree_exposed]

    agree_not_out = [evaluate_sample_time(s, drug, close, distant, dict_biopsy, 0, treatment_sample_days_treatment,
                                          treatment_sample_days_biopsy_since_end_treatment) for s in agree_not_exposed]

    toreg = agree_out + agree_not_out

    # Save for regression
    testing = pd.DataFrame(toreg)
    testing.dropna(inplace=True)
    testing_2 = testing.copy()
    testing_2[ttype] = 1
    testing_2[drug] = 1
    toregs.append(testing_2)

    distant_exposed = np.sum([v for k, v in v1.items() if k in distant])
    distant_notexposed = np.sum([v for k, v in v2.items() if k in distant])
    close_exposed = np.sum([v for k, v in v1.items() if k in close])
    close_notexposed = np.sum([v for k, v in v2.items() if k in close])

    dic_prop_ttypes[f'{ttype}_{drug}'] = [close_exposed, close_notexposed,
                                          distant_exposed, distant_notexposed]

    plt.bar(1, distant_exposed / (distant_exposed + close_exposed), width=0.75, color='#5D9631')
    plt.bar(0, distant_notexposed / (distant_notexposed + close_notexposed), width=0.75, color='grey')

    odds_ratio, pvalue = fisher_exact([[close_notexposed, close_exposed],
                                       [distant_notexposed, distant_exposed]])
    print(odds_ratio, pvalue)

    plt.xticks([0, 1], ['Not Exposed', 'Exposed', ], rotation=90)
    plt.xlim(-0.5, 1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel('Proportion of\ndistal metastasis')

    plt.show()

    return dic_prop_ttypes, toregs


def evaluate_sample_time(s, drug, prox, dist, dict_biopsy, exp, treatment_sample_days_treatment,
                         treatment_sample_days_biopsy_since_end_treatment):
    try:
        days_treat = np.sum(treatment_sample_days_treatment[s][drug]) / 200
    except:
        days_treat = np.nan
    try:
        days_since_end = np.max(treatment_sample_days_biopsy_since_end_treatment[s][drug]) / 200
    except:
        days_since_end = np.nan

    prox_v = 0
    if dict_biopsy[s] in prox:
        prox_v = 0
    elif dict_biopsy[s] in dist:
        prox_v = 1

    return days_treat, days_since_end, prox_v, exp
