import pandas as pd
import numpy as np
import tqdm
import gzip
import pickle
import glob

import matplotlib.pyplot as plt


class Exposures:

    def __init__(self, fn):
        self.exposures = pd.read_csv(fn, sep='\t', index_col=0)
        self.signatures = self.exposures.index.tolist()
        self.samples = self.exposures.columns.tolist()

    def active(self, sig):
        func = lambda x: sum(x > 0) / len(x)
        return func(self.exposures.loc[sig].values)

    def gamma_fit(self, sig):
        x = np.log(self.exposures.loc[sig].values + 1)
        x = x[x > 0]
        fit_alpha, fit_loc, fit_beta = stats.gamma.fit(x)
        return fit_alpha, fit_loc, fit_beta

    def mean_log(self, sig):
        x = np.log(self.exposures.loc[sig].values + 1)
        x = x[x > 0]
        return np.mean(x)

    def std_log(self, sig):
        x = np.log(self.exposures.loc[sig].values + 1)
        x = x[x > 0]
        return np.std(x)

    def plot(self, **kwargs):
        sig = np.random.choice(self.signatures)
        x = np.log(self.exposures.loc[sig].values + 1)
        x = x[x > 0]
        plt.hist(x, **kwargs)
        plt.xlabel('exposure')
        plt.ylabel('freq')
        plt.title(sig)
        plt.show()


class Processes:

    def __init__(self, fn):
        self.processes = pd.read_csv(fn, sep='\t', index_col=0)
        self.contexts = list(map(lambda x: x[0] + x[2] + x[6] + x[4], self.processes.index))
        self.processes.index = self.contexts
        self.names = self.processes.columns.tolist()

    def frequency(self, sig, context):
        return self.processes.loc[context, sig]


from scipy import stats


def negbinom_generator(mu, dispersion, size=1):
    """
    Returns a list of random neg binom draws
    """

    _lambdas = stats.gamma.rvs(1 / dispersion, scale=mu * dispersion, size=size)
    return stats.poisson.rvs(_lambdas, size=size)


def get_noisy_catalogue(n, exposures=None, processes=None):
    """
    produces table with:
        rows: contexts (in Rozen notation: f1+ref+f2+alt)
        columns: samples
    """

    signature_list = processes.names

    catalogue = {}  # number of mutations per channel sample-wise
    makeup = {}  # exposures per signature sample-wise

    for sample in range(n):

        makeup_sample = [0] * len(signature_list)

        counts_per_channel = np.zeros(len(processes.contexts), dtype=int)

        # Decide which processes are active

        active_set = [sig for sig in exposures.signatures if stats.binom.rvs(1, exposures.active(sig)) == 1]

        for sig in active_set:

            fit_alpha, fit_loc, fit_beta = exposures.gamma_fit(sig)
            log_exposure = stats.gamma.rvs(fit_alpha, loc=fit_loc, scale=fit_beta, size=1)[0]
            exposure = np.exp(log_exposure)
            makeup_sample[signature_list.index(sig)] = exposure

            for i, context in enumerate(processes.contexts):
                mu = exposure * processes.frequency(sig, context)
                c = negbinom_generator(mu, 0.1, size=1)[0]
                counts_per_channel += c * np.eye(len(processes.contexts), dtype=int)[i]

        catalogue[f'sample_{sample}'] = counts_per_channel
        makeup[f'sample_{sample}'] = makeup_sample

    return catalogue, makeup


fn_exposures = '../../data/reversed.snvs.exposures.tsv'
exposures = Exposures(fn_exposures)

fn_processes = '../../data/reversed.snvs.processes.tsv'
processes = Processes(fn_processes)


class SignatureInjector:
    """
    given a collection of foreign signatures, it has a method to inject
    noisy --negative binomial-- exposures to a mutational catalogue
    """

    def __init__(self, sig_dict):

        self.signatures = sig_dict
        self.names = sig_dict.keys()

    def inject(self, sig_name, exposure, catalogue):

        signature = self.signatures[sig_name]

        if sig_name not in self.names:
            raise ValueError('Unknown signature name')

        new_catalogue = {}

        for sample in catalogue.columns:
            v = catalogue[sample].values.copy()
            if exposure > 0:
                for i, freq in enumerate(signature):
                    mu = exposure * freq
                    c = negbinom_generator(mu, 0.1, size=1)[0]
                    v += c * np.eye(len(processes.contexts), dtype=int)[i]
            new_catalogue[sample] = v

        return pd.DataFrame(new_catalogue, index=catalogue.index)


def foreign_signatures(processes):
    foreign_signatures = {}

    df = pd.read_csv('../../data/Colon-Rectum.snvs.processes.tsv', sep='\t', index_col=0)

    sig = '2_1.0_SBS17b-0.96'
    foreign_signatures['SBS17b'] = df[sig].values.tolist()

    sig = '10_0.99_NA'
    foreign_signatures['oxaliplatin'] = df[sig].values.tolist()

    df = pd.read_csv('../../data/basic_filters_mutect_CN.SBS31.sigs.tsv', sep='\t', index_col=0)

    sig = 'SBS31'
    foreign_signatures[sig] = list(df[sig].values / df[sig].sum())

    df = pd.read_csv('../../data/basic_filters_mutect_CN.SBS35.sigs.tsv', sep='\t', index_col=0)

    sig = 'SBS35'
    foreign_signatures[sig] = list(df[sig].values / df[sig].sum())

    for sig in foreign_signatures:
        foreign = pd.DataFrame({sig: foreign_signatures[sig]}, index=processes.contexts)
        list_df = [processes.processes, foreign]
        df = pd.concat(list_df, axis=1)
        df.index.name = 'index'
        df.to_csv(f'../../data/synthetic_data_basic_filters/signatures_{sig}.tsv', sep='\t')

    return foreign_signatures
