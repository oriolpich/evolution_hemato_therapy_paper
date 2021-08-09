import os
import gzip
import pandas as pd
import pickle
import collections
import functools


def leveled_default_dict(inner, level=1):
    level = level - 1
    if level == 0:
        return collections.defaultdict(inner)
    else:
        f = functools.partial(leveled_default_dict, inner=inner, level=level)
        return collections.defaultdict(f)


# this will create a dictionary with the ML assignment for each signature
def ML_assignment_dictionary(type_mut, process_file, exposures_file):
    rev = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    name_tumor = os.path.basename(process_file).split('.')[0]
    outfile_name = '{}/{}.ML_assign.pckl.gz'.format(os.path.dirname(process_file), name_tumor)

    if not os.path.isfile(outfile_name):

        class_d = leveled_default_dict(str, level=2)
        class_full = leveled_default_dict(dict, level=2)

        if os.path.isfile(exposures_file):

            W = pd.read_csv(process_file, sep='\t', index_col=0)
            H = pd.read_csv(exposures_file, sep='\t', index_col=0)

            H = H.T.loc[:, ~H.T.columns.duplicated()].T
            W = W.loc[:, ~W.columns.duplicated()]

            for sample_index, sample in enumerate(H.columns):

                # multiply exposures by W
                ex = H[sample] * W
                ex['row_sum'] = ex.sum(axis=1)
                new = ex.div(ex['row_sum'], axis=0)
                new.drop('row_sum', axis=1, inplace=True)
                class_full[sample] = new.to_dict(orient='index')

                # get the one with the maximum probability
                res = new.idxmax(axis=1)

                for i, row in (W * H[sample]).iterrows():

                    # assign the signature to this context
                    class_d[sample][i] = res[i]

                    # now the complementary reverse too, for DBS and SNV
                    if type_mut == 'snvs':
                        rev_order = '{}[{}>{}]{}'.format(rev[i[-1]], rev[i[2]], rev[i[4]], rev[i[0]])
                        class_d[sample][rev_order] = res[i]

                    elif type_mut == 'dbs':
                        rev_order = '{}{}_{}{}'.format(rev[i[1]], rev[i[0]], rev[i[-1]], rev[i[-2]])
                        class_d[sample][rev_order] = res[i]

                for i, v in new.to_dict(orient='index').items():

                    if type_mut == 'SNV':
                        rev_order = '{}[{}>{}]{}'.format(rev[i[-1]], rev[i[2]], rev[i[4]], rev[i[0]])
                        class_full[sample][rev_order] = v

                    elif type_mut == 'dbs':
                        rev_order = '{}{}_{}{}'.format(rev[i[1]], rev[i[0]], rev[i[-1]], rev[i[-2]])
                        class_full[sample][rev_order] = v

            pickle.dump(dict(class_d), gzip.open(outfile_name, 'wb'))
            outfile_name = '{}/{}.total_assign.pckl.gz'.format(os.path.dirname(process_file), name_tumor)
            pickle.dump(dict(class_full), gzip.open(outfile_name, 'wb'))
        else:
            print(exposures_file, "not found")
    else:
        class_d = pickle.load(gzip.open(outfile_name))
        class_full = pickle.load(gzip.open(outfile_name.replace('ML_assign', 'total_assign')))

    return class_d, class_full


# retrieve ML
def assing_ML(df, dML):
    return dML[df['SAMPLE']].get(df['VARIANT_CLASS'], 'not_found')


def assign_probability(df, dML_full, signature):
    return dML_full[df['SAMPLE']][df['VARIANT_CLASS']][signature]


def doMLsignatures(path_original, path_processes_file, exposures_file, type_mut, outpath, annotate=False):
    name_ttype = os.path.basename(path_processes_file).split('.')[0]

    # get the dictionary with the MLikelihood
    dML, dML_full = ML_assignment_dictionary(type_mut, path_processes_file, exposures_file)

    if annotate:

        df = pd.read_csv(path_original, sep='\t', low_memory=False)
        df = df[df['CLASS'] == type_mut]
        all_samples_extraction = set(list(dML.keys()))
        all_samples = set(df['SAMPLE'].tolist())

        wanted_samples = [sample for sample in all_samples if sample in all_samples_extraction]
        df = df[df['SAMPLE'].isin(wanted_samples)]
        df['ML'] = df.apply(assing_ML, axis=1, args=(dML,))

        df_exposures = pd.read_csv(exposures_file, sep='\t', index_col=0)
        for signature in df_exposures.index.tolist():
            df[signature] = df.apply(assign_probability, axis=1, args=(dML_full, signature))

        for sig, data in df.groupby(by='ML'):
            outf = '{}/split/{}.{}.gz'.format(outpath, name_ttype, sig)
            data.to_csv(outf, sep='\t', index=False, header=True, compression='gzip')
        # this is the full file, we can also generate the file per signature to parallelize afterwards
        df.to_csv('{}/{}.ML.{}.gz'.format(outpath, name_ttype, type_mut), sep='\t', index=False, compression='gzip')


    else:
        return dML, dML_full
