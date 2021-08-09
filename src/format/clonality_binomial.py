
from scipy.stats import binom
import numpy as np
import matplotlib.pyplot as plt

def Binomials(var_counts, depth, P, CNt, CNn):
    vaf_to_ccf_factor = P / ((CNn * (1 - P) + P * CNt) * 100)
    expVAF = [i * vaf_to_ccf_factor for i in range(1, 101)]
    return binom.pmf(var_counts, depth, expVAF)


def SubCint(x, prob):

    xnorm = x / sum(x)  
    xsort = sorted(xnorm, reverse=True)
    xcumLik = np.cumsum(xsort)
    n = sum(xcumLik < prob) + 1
    LikThresh = xsort[n - 1]
    OrderedCCF = np.arange(0.01, 1.01, 0.01)
    cint = []
    name = []
    for i in range(len(xnorm)):
        if xnorm[i] >= LikThresh:
            cint.append(x[i])
            name.append(OrderedCCF[i])

    rt = name[len(name) - 1]
    AbsCCF = name[cint.index(max(cint))]
    probSubclonal = sum(xnorm[1:95])
    probClonal = sum(xnorm[95:101])

    if rt >= 1:
        result1 = 'clonal'
        if probClonal >= 0.75:
            result2 = 'clonal'
        else:
            result2 = 'undefined'
    else:
        result1 = 'subclonal'
        if probSubclonal >= 0.75:
            result2 = 'subclonal'
        else:
            result2 = 'undefined'

    return result1, result2, AbsCCF


def getClonalStatus(df):

    var_counts = df['VAR_COUNTS']
    CNt = df['TOTAL_CN']
    CNn = df['NORMAL_CN']
    P = df['PURITY']
    ref_counts =df['REF_COUNTS']

    depth = ref_counts + var_counts

    try :
        # binomial function
        x = Binomials(var_counts, depth, P, CNt, CNn)

        # Call main function of the method
        c1, c2, ccf = SubCint(x, 0.95)

        df['Clonal1'] = c1
        df['Clonal2'] = c2
        df['CCF'] = ccf

    except:

        df['Clonal1'] = 'NA'
        df['Clonal2'] = 'NA'
        df['CCF'] = 'NA'

    return df


def plot_clonality(merge, sample, outpath):
    ## plot
    config_params(12)
    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    plt.title('INDELS n = {}'.format(len(merge[merge['CLASS'] == 'INDEL'])))
    merge[merge['CLASS'] == 'INDEL']['VAF'].hist(bins=100)
    plt.savefig('{}/{}_CCF_INDEL.png'.format(outpath, sample))
    plt.close()
    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    plt.title('SNVs n = {}'.format(len(merge[merge['CLASS'] == 'SNV'])))

    merge[merge['CLASS'] == 'SNV']['VAF'].hist(bins=100)
    plt.savefig('{}/{}_CCF_SNV.png'.format(outpath, sample))
    plt.close()
