import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from statsmodels.robust.robust_linear_model import RLM

sys.path.append("../config/")
from config_plot import config_plot_params


def run_irls(x_coords, y_coords):
    exog_data = provide_exog(x_coords)

    model_irls = RLM(y_coords, exog_data)
    try:
        res_irls = model_irls.fit()
    except ZeroDivisionError:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    weighted_r = w_corr(x_coords, y_coords, res_irls.weights)
    wmean = w_mean(y_coords, res_irls.weights)
    pvalue = res_irls.pvalues[0]
    slope = res_irls.params[0]
    intercept = res_irls.params[1]
    return pvalue, weighted_r, slope, intercept, wmean


def provide_exog(covariate):
    exog_data = np.array(np.concatenate((np.array(covariate), np.ones_like(covariate))).reshape(2, len(covariate)))
    return np.transpose(exog_data)


def w_mean(x, w):
    """weighted mean"""
    return np.dot(x, w) / np.sum(w)


def w_cov(x, y, w):
    """weighted covariance"""
    return np.sum(w * (x - w_mean(x, w)) * (y - w_mean(y, w))) / np.sum(w)


def w_corr(x, y, w):
    """weighted correlation"""
    return w_cov(x, y, w) / np.sqrt(w_cov(x, x, w) * w_cov(y, y, w))


def do_plot(row, method_of_extraction, high_mut_platinum, outpath, dic_samples_ttypes, dic_aging_sample_level, dic_HSC):
    config_plot_params(5)
    all_cols = row.columns.tolist()

    fig, ax = plt.subplots(1, 2, figsize=(2, 1.3), sharey=True, )

    age_patients_taml = [dic_aging_sample_level[c] for c in dic_samples_ttypes['tAML'] if
                         ((c in all_cols) & (dic_aging_sample_level.get(c, 'no_data') != 'no_data'))]
    mut_patients_taml = [float(row[c]) for c in dic_samples_ttypes['tAML'] if
                         ((c in all_cols) & (dic_aging_sample_level.get(c, 'no_data') != 'no_data'))]

    ax[1].set_title('tAML')
    ax[0].set_title('Primary')

    pval_taml, coef_taml, slope_taml, intercept, wmean = run_irls(age_patients_taml, mut_patients_taml)

    print("tAML", "pval", pval_taml, "R", coef_taml, "Slope", slope_taml)

    sns.regplot(age_patients_taml, mut_patients_taml, color='#c83737ff', ax=ax[1], robust=True,
                scatter_kws={'s': 5})

    plat_cols = [c for c in high_mut_platinum if c in row.columns]

    aging_vals = [dic_aging_sample_level[c] for c in plat_cols if dic_aging_sample_level[c] != 'no_data']
    allowed_samples = [c for c in plat_cols if dic_aging_sample_level[c] != 'no_data']

    ax[1].scatter(aging_vals, [float(row[c]) for c in allowed_samples], color='orange', s=4)

    # Primary AML
    age_patients = [dic_aging_sample_level[c] for c in dic_samples_ttypes['primary'] if
                    ((c in all_cols) & (c not in high_mut_platinum) &
                     (dic_aging_sample_level.get(c, 'no_data') != 'no_data'))]

    mut_patients = [float(row[c]) for c in dic_samples_ttypes['primary'] if
                    ((c in all_cols) & (c not in high_mut_platinum) &
                     (dic_aging_sample_level.get(c, 'no_data') != 'no_data'))]

    p = sns.regplot(age_patients, mut_patients, color='#0088aaff', ax=ax[0], robust=True, scatter_kws={'s': 5})

    pval_aml, coef_aml, slope_aml, intercept, wmean = run_irls(age_patients, mut_patients)
    print("Primary", "pval", pval_aml, "R", coef_aml, "Slope", slope_aml)

    for t in range(2):
        ax[t].set_xlabel('Age')
        ax[t].set_ylabel('HSC sig Mutations')

    # HSC
    age_patients = [dic_aging_sample_level[c] for c in dic_HSC if
                    (c in all_cols) & (c in dic_aging_sample_level)]
    mut_patients = [float(row[c]) for c in dic_HSC if
                    (c in all_cols) & (c in dic_aging_sample_level)]

    p = sns.regplot(age_patients, mut_patients, color='grey', ax=ax[0], robust=True,
                    scatter_kws={'s': 5, 'alpha': 0.5})

    y = p.get_lines()[0].get_ydata()

    pvalue2, r_value2, slope2, intercep2t, wmean2 = run_irls(age_patients, mut_patients)
    print("HSC Boxtel", "pval", pvalue2, "R", r_value2, "Slope", slope2, "intercept", intercep2t)

    x_vals = np.array(ax[0].get_xlim())
    y_vals = intercep2t + slope2 * x_vals

    for t in range(2):
        ax[t].plot(x_vals, y_vals, '--', c='#737373', lw=0.2)
        ax[t].spines['top'].set_visible(False)
    
    plt.savefig('{}/aging_{}_{}.svg'.format(outpath, method_of_extraction,
                                            row.index.tolist()[0]))
    plt.show()
