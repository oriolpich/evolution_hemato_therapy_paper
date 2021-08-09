import io
import os
import sys
import gzip
from glob import glob
from bgreference import hg19
from collections import defaultdict
import numpy as np
from tqdm import tqdm

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

import pandas as pd

# vcf parser
def vcf_reader(path):
    with open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(str.join(os.linesep, lines)),
        sep='\t', low_memory=False,
    ).rename(columns={'#CHROM': 'CHROM'})


# get reads info from VCF file
def get_reads(df, last_column):
    list_format = df['FORMAT'].split(':')
    d = {list_format[i]: i for i in range(len(list_format))}

    colvaf_list = df[last_column].split(':')

    reference = int(colvaf_list[d['AD']].split(',')[0])
    variant = int(colvaf_list[d['AD']].split(',')[1])

    total_reads = reference + variant
    df['total_reads'] = total_reads
    df['ref_reads'] = reference
    df['var_reads'] = variant

    if total_reads > 0:
        df['VAF'] = df['var_reads'] / df['total_reads']
    else:
        df['VAF'] = np.nan

    return df

# create order format similar to SigProfiler
def create_snv_class(df):
    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    x = df['TRIPLET']
    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])
    return out

# generate order of SBS
def snvs_order():
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

def return_ordered_matrix(df):
    samples_dict = defaultdict(dict)
    order = snvs_order()

    for sample, data in df.groupby(by='SAMPLE'):
        dic_count = data['VARIANT_CLASS'].value_counts().to_dict()
        for i in order:
            samples_dict[sample][i] = dic_count.get(i, 0)

    matrix = pd.DataFrame.from_dict(samples_dict)
    matrix = matrix.loc[order]

    return matrix



def process_boxtel_osorio(path):

    path_osorio = glob(path)

    all_osorio = []
    wanted_Cols = ['CHROM', 'pos-1', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
        'total_reads', 'ref_reads', 'var_reads', 'VAF', 'len_alt', 'len_ref', 'TYPE', 
        'TRIPLET', 'EXTENDED', 'CLASS', 'VARIANT_CLASS', 'SAMPLE']

    for f in tqdm(path_osorio):
        
        df = vcf_reader(f)
        
        sample_file = os.path.basename(f).split('_')[0]
        df.reset_index(inplace = True)
        df.columns = df.columns.tolist()[1:] + ['score']
        
        df_reads = df.apply(get_reads, axis=1, args=(sample_file,))
        # select whether we have SNVs or others
        df_reads['len_alt'] = df_reads['ALT'].str.len()

        # number of characters in ref
        df_reads['len_ref'] = df_reads['REF'].str.len()

        # first classification between SNV and others
        df_reads['TYPE'] = df_reads.apply(
            lambda x: 'SNV' if ((x['len_alt'] == 1) and (x['len_ref'] == 1) and (x['ALT'] != '-') and (x['REF'] != '-')) else 'INDEL', axis=1
        )

        df_reads['pos-1'] = df_reads['POS'] - 1

        # get the triplet
        df_reads['TRIPLET'] = df_reads.apply(lambda x: hg19(x['CHROM'], x['pos-1'], 3), axis=1)
        df_reads['EXTENDED'] = df_reads.apply(lambda x: hg19(x['CHROM'], int(x['POS']) - 2, 5), axis=1)

        snv_df = df_reads[df_reads['TYPE'] != 'INDEL']
        snv_df['CLASS'] = 'SNV'
        snv_df['VARIANT_CLASS'] = snv_df.apply(create_snv_class, axis=1)
        snv_df['SAMPLE'] = sample_file
        
        all_osorio.append(snv_df[wanted_Cols])

    merged = pd.concat(all_osorio)
    matrix = return_ordered_matrix(merged)
    
    # path to save Osorio's data
    matrix.to_csv('/workspace/projects/reverse_calling/data/signatures_2020/boxtel.snvs.tsv', 
                sep ='\t', header = True, index = False)

# path to vcf from Osorio et al
process_boxtel_osorio('/workspace/projects/all_aecc/Osorio_data/vcfs/*.vcf')