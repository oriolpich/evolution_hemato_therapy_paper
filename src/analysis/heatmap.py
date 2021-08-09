from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import seaborn as sns
import json
from scipy.cluster import hierarchy
import json

from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import Table2x2


from scipy.spatial.distance import pdist
import matplotlib as mpl

# FIXME: whole file repeated

# config for matplotlib
def config_params(font_size=7):
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'


def get_orders(subset_drivers, sample_taml_beat):

    dic_sample = defaultdict(dict)
    for sample, data in subset_drivers.groupby(by='SAMPLE'):
        d = data['Hugo_Symbol'].drop_duplicates().tolist()
        for t in d:
            dic_sample[sample][t] = 1
            
    matrix_driver_genes = pd.DataFrame(dic_sample).fillna(0)
    dic_counts = matrix_driver_genes.sum(axis = 1).to_dict()
    sorted_genes = sorted(dic_counts, key=dic_counts.get, reverse = True)
    sorted_genes_nonzero = [s for s in sorted_genes if dic_counts[s]>1]

    wanted_genes = sorted_genes_nonzero
    g = classic_mutual_exclusivity_visualization(matrix_driver_genes.T, wanted_genes, [])
    new_order_view = g.dendrogram_col.reordered_ind
    sample_list = matrix_driver_genes.T.reset_index().loc[new_order_view]['index'].tolist()

    primary = [s for s in sample_list if s not in sample_taml_beat]
    second = [s for s in sample_list if s in sample_taml_beat]
    final_ord = primary+second

    return final_ord, primary, second, wanted_genes, dic_sample


def classic_mutual_exclusivity_visualization(regression_table, sorted_treatments, colors):

    df = regression_table.copy()
    df = df[sorted_treatments]
    X = df.values
    
    weights = [1 / (i + 1) for i, c in enumerate(df.columns)]
    Y = pdist(X, metric='hamming', w=np.array(weights))
    linkage = hierarchy.linkage(Y, method='ward')
    g = sns.clustermap(df.T, col_linkage=linkage, col_cluster=True, 
                       row_cluster=False, figsize=(4, 6),  #col_colors = colors,
                       cmap = 'Blues', xticklabels = True, yticklabels = True )
    plt.close()

    return g

def get_type_alts(subset_drivers, final_ord, primary, second, wanted_genes, dic_sample):

    for sample, data in subset_drivers[subset_drivers['Hugo_Symbol'].isin(wanted_genes)].groupby(by='SAMPLE'):
        for type_mut, data_t in data.groupby(by='Consequence'):
            if type_mut == 'frameshift_variant':
                d = data_t['Hugo_Symbol'].drop_duplicates().tolist()
                for t in d:
                    dic_sample[sample][t] = 3

            else:
                for i, row in data_t.iterrows():
                    if row['Consequence'] == 'missense_variant':
                        dic_sample[sample][row['Hugo_Symbol']] = 1
                    else:
                        dic_sample[sample][row['Hugo_Symbol']] = 2
                    
    matrix_driver_genes = pd.DataFrame(dic_sample).fillna(0)
    matrix_driver_genes = matrix_driver_genes.loc[wanted_genes]
    matrix_driver_genes = matrix_driver_genes.loc[(matrix_driver_genes.sum(axis=1) != 0),
                                                (matrix_driver_genes.sum(axis=0) != 0)]

    final_ord = [s for s in final_ord if s in matrix_driver_genes.columns]
    primary = [s for s in primary if s in matrix_driver_genes.columns]
    second = [s for s in second if s in matrix_driver_genes.columns]
    
    number_drivers = subset_drivers['SAMPLE'].value_counts().to_dict()
    n_drivers_sample = [number_drivers[s] for s in final_ord]

    return final_ord, primary, second, number_drivers, n_drivers_sample, matrix_driver_genes


def get_samples_data(mutations_beatAML, clinical_data, drivers_AML):

    df_AML = pd.read_csv(drivers_AML, sep ='\t')
    driver_AML = list(df_AML['SYMBOL'].unique())

    #cbioportal to get indels
    muts_cbioport = pd.read_csv(mutations_beatAML, sep ='\t')
    muts_cbioport['SAMPLE'] = muts_cbioport['Tumor_Sample_Barcode'].apply(lambda x : x.split('_')[-1])

    # metadata from cbioportal
    df_metadata_beat = pd.read_csv(clinical_data, sep ='\t', skiprows=4)

    C_A = ['15-00056', '14-00527', '16-00118']

    df_metadata_beat['SAMPLE'] = df_metadata_beat['SAMPLE_ID'].apply(lambda x : x.split('_')[-1])

    sample_taml_beat = df_metadata_beat[(df_metadata_beat['SPECIFIC_DIAGNOSIS_AT_INCLUSION']=='Therapy-related myeloid neoplasms')]['SAMPLE'].unique()

    muts_cbio = muts_cbioport[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'Consequence','Variant_Type', 'short_aa_change', 
            'Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele2', 'SAMPLE']].drop_duplicates()


    muts_cbio['tAML'] = muts_cbio['SAMPLE'].apply(lambda x : 'tAML' if x in sample_taml_beat else 'primaryAML')

    subset_drivers = muts_cbio[muts_cbio['Hugo_Symbol'].isin(driver_AML)]

    return subset_drivers, sample_taml_beat, muts_cbio, driver_AML


def do_barplot(mutations_beatAML, clinical_data, drivers_AML, outpath):

    subset_drivers, sample_taml_beat, muts_cbio, driver_AML = get_samples_data(mutations_beatAML, clinical_data, drivers_AML)

    bars = []
    gen = []
    config_params(5)
    fig, ax = plt.subplots(1, 1, figsize = (4.5, 0.75))
    keep_c = {}
    for gene, data in muts_cbio[muts_cbio['Hugo_Symbol'].isin(driver_AML)].groupby(by='Hugo_Symbol'):
        keep_c[gene] = len(data)
        bars.append(len(data))
        gen.append(gene)

    sorted_k = sorted(keep_c, key=keep_c.get, reverse=True)

    plt.bar(np.arange(len(sorted_k)), [keep_c[c] for c in sorted_k], color = '#D4AA00')
    plt.xticks(np.arange(len(bars)), ["$\it{0}$".format(lab) for lab in sorted_k], rotation = 90)
    plt.ylabel('Number of samples\nwith non-silent mutations')
    plt.xlim(-1, len(sorted_k))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Genes with identified\nsignals of positive selection')

    plt.savefig('{}/barplot_beat.svg'.format(outpath))
    
    plt.show()

def do_heatmap(mutations_beatAML, clinical_data, drivers_AML, outpath):

    subset_drivers, sample_taml_beat, muts_cbio, driver_AML = get_samples_data(mutations_beatAML, clinical_data, drivers_AML)
    final_ord, primary, second, wanted_genes, dic_sample = get_orders(subset_drivers, sample_taml_beat)
    final_ord, primary, second, number_drivers, n_drivers_sample, matrix_driver_genes = get_type_alts(subset_drivers, final_ord,
                                                                                        primary, second, wanted_genes, dic_sample)

    config_params(6)

    fig, ax = plt.subplots(3, 2, figsize = (6.2, 6.25), sharex='col',
                        gridspec_kw={'height_ratios': [0.1, 0.02, 1], 
                                    'width_ratios':[6,1]})

    # ----------------
    # AXIS DEFINITION
    # mut rate

    mutrate = ax[0, 0]
    mutrate.set_ylabel('Number of mutations\nin driver genes', fontsize = 6)
    # mut spectra
    spectra_axis = ax[2,0]

    bar_axis1 = ax[1, 0]

    # mut bar
    bar_axis = ax[2,1]
    mutrate.bar(np.arange(len(n_drivers_sample)),n_drivers_sample, color = 'black')

    mutrate.set_xlim(-1, len(final_ord)+1)
    mutrate.spines['top'].set_visible(False)
    mutrate.get_xaxis().set_visible(False)

    bar_axis1.barh(0, len(primary),
            height = 0.2,
            color = '#0088aaff')

    bar_axis1.barh(0, len(second), left = len(primary),
            height = 0.2,
            color = '#c83737ff')

    bar_axis1.get_xaxis().set_visible(False)
    bar_axis1.get_yaxis().set_visible(False)
    bar_axis1.set_xlim(0, len(final_ord))

    new_cmap = LinearSegmentedColormap.from_list("", ["white", "green", "orange", "darkred"])
    sns.heatmap(matrix_driver_genes[final_ord].loc[wanted_genes],cbar=False, 
                square=False, cmap = new_cmap, ax =spectra_axis, alpha = 0.6,zorder = 1,
    )

    spectra_axis.set_xlim(0, len(final_ord))
    spectra_axis.get_xaxis().set_visible(False)
    spectra_axis.set_yticklabels(["$\it{0}$".format(lab) for lab in wanted_genes])

    for ix, gene in enumerate(wanted_genes[::-1], start=1):
        ns = matrix_driver_genes[primary].loc[gene]
        len_prim = len(ns[ns>0])
        ns = matrix_driver_genes[second].loc[gene]
        len_sec = len(ns[ns>0])
        tot = len_prim+len_sec
        height = 0.4
        bar_axis.barh(ix-0.5 - 0.2, len_prim/len(primary),  color = '#0088aaff', height=height) #align='edge',
        bar_axis.barh(ix-0.5 + 0.2, len_sec/len(second), height=height, color = '#c83737ff')
        if len_prim>0:
            bar_axis.text((len_prim/len(primary)) +0.01,
                            ix-1,  len_prim, fontsize = 5.5 )
        if len_sec>0:
            bar_axis.text(( len_sec/len(second))+0.01,
                            ix-0.5,  len_sec, fontsize = 5.5)
    
    mutrate.set_ylim(-3, 10)

    bar_axis.set_ylim(0, len(wanted_genes))

    bar_axis.spines['left'].set_visible(False)
    bar_axis.spines['bottom'].set_visible(False)
    bar_axis.spines['right'].set_visible(False)
    bar_axis.get_yaxis().set_visible(False)

    bar_axis.xaxis.tick_top()
    bar_axis.set_xlabel('Proportion of patients', fontsize = 6)
    bar_axis.xaxis.set_label_position('top') 

    for i in range(0,2):
        ax[i, 1].spines['left'].set_visible(False)
        ax[i, 1].spines['bottom'].set_visible(False)
        ax[i, 1].spines['right'].set_visible(False)
        ax[i, 1].spines['top'].set_visible(False)
        ax[i, 1].get_yaxis().set_visible(False)
        ax[i, 1].get_xaxis().set_visible(False)

    plt.subplots_adjust(wspace=0.05, hspace=0.001)
    plt.savefig('{}/beatAML_heatmap.png'.format(outpath), dpi = 1000)
    plt.show()


def get_ratio(res_mut_drivers):

    total_prim = res_mut_drivers['primaryAML'].sum()
    total_tAML = res_mut_drivers['tAML'].sum()
    res_mut_drivers  = res_mut_drivers[(res_mut_drivers['primaryAML']>1)&(res_mut_drivers['tAML']>1)]

    fig,ax = plt.subplots(1, 1, figsize = (5,5))
    xlab = []
    keep_res = {}
    for ix, (i, row) in enumerate(res_mut_drivers.sort_values(by='tAML', ascending = False).head(25).iterrows()):

        t = Table2x2([[row['tAML'], total_tAML],
                    [row['primaryAML'], total_prim]])

        pval = t.oddsratio_pvalue()
        o = t.oddsratio_confint()
        cmin, cmax = o
        odds = t.oddsratio

        if odds >5:
            odds = 5
        keep_res[i] = (np.log2(odds), np.log2(cmin), np.log2(cmax), pval )

    return keep_res

def odds_ratio_tAMLs(mutations_beatAML, clinical_data, drivers_AML, outpath):

    # load data
    subset_drivers, sample_taml_beat, muts_cbio, driver_AML = get_samples_data(mutations_beatAML, clinical_data, drivers_AML)

    dic_res = defaultdict(dict)
    for aml, data in muts_cbio.groupby(by='tAML'):
        for gene in driver_AML:

            # accomodate reviewers comments
            if gene == "DNMT3A":
                set_dnmt = len(data[(data['Hugo_Symbol']==gene)&(data["short_aa_change"].str.contains("R882"))]['SAMPLE'].unique())
                dic_res[aml]["DNMT3A p.R882"] = set_dnmt
                set_dnmt = len(data[(data['Hugo_Symbol']==gene)&(~data["short_aa_change"].str.contains("R882", na=False) )]['SAMPLE'].unique())
                dic_res[aml]["DNMT3A other"] = set_dnmt
            else:
                set_g = len(data[data['Hugo_Symbol']==gene]['SAMPLE'].unique())
                dic_res[aml][gene] = set_g

    res_mut_drivers = pd.DataFrame(dic_res)

    # get the ratios
    keep_res = get_ratio(res_mut_drivers)

    tt_df = pd.DataFrame(keep_res, index = ['Val', 'LC', 'HC', 'pvalue']).T

    config_params(5)
    res_reg = tt_df.sort_values(by='pvalue', ascending = True)
    qvals = multipletests(res_reg.pvalue, alpha=0.05, method='fdr_bh')
    res_reg['QVAL'] = qvals[1]
    fig, ax = plt.subplots(1, 1, figsize = (1.5,1.8
                                        ))
    labels = []
    plt.title('AML vs t-AML')
    for ix, (i, row) in enumerate(res_reg.sort_values(by='pvalue', ascending = False).iterrows()):
        
        edgecolor = 'white'
        if row['QVAL']<0.05:
            edgecolor = 'black'
        
        color = '#0088aaff'
        if row['Val']>0:
            color = '#c83737ff'
        plt.scatter( row['Val'], ix, color = color, edgecolor = edgecolor)
        plt.plot([row['LC'], row['HC']], [ix, ix], color = color)
        labels.append(i)

    plt.yticks(np.arange(len(res_reg)), ["$\it{0}$".format(lab) for lab in labels] )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.vlines(0, 0, len(res_reg), ls = '--', color = 'grey')
    plt.xlabel('Odds-ratio')
    plt.savefig('{}/odds_tAML_AML.svg'.format(outpath))
    plt.show()