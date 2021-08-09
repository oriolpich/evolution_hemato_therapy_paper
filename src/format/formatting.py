import os
import sys

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

import pandas as pd
import io
import os
from pybedtools import BedTool
import numpy as np
from bgreference import hg38
import gzip

from indels_class import IndelsClassifier
from clonality_binomial import getClonalStatus, plot_clonality
from filter_gnomad import filter_gnomad_38


# Get mutations only in mappable regions
def get_mappable_regions(bed, path_mappable):
    mapp_bed = BedTool(path_mappable)
    mappable_mut = bed.intersect(mapp_bed, wa=True, v=True)

    return mappable_mut


# VCF parser
def vcf_reader(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        }, sep='\t', low_memory=False,
    ).rename(columns={'#CHROM': 'CHROM'})


# get reads info from VCF file
def return_reads(row):
    list_format = row['FORMAT'].split(':')
    d = {list_format[i]: i for i in range(len(list_format))}
    colvaf_list = row["TUMOR"].split(':')

    reference_v = row['REF']
    alternate_v = row['ALT']
    reference, variant = 0, 0

    if len(row['REF']) == 2:
        reference_v = row['REF'][0]
        alternate_v = row['ALT'][0]

    if row['FCLASS'] == 'SNV':
        reference = int(colvaf_list[d['{}U'.format(reference_v)]].split(',')[0])
        variant = int(colvaf_list[d['{}U'.format(alternate_v)]].split(',')[0])

    elif row['FCLASS'] == 'INDEL':
        total_depth = int(colvaf_list[d['DP']])
        variant = int(colvaf_list[d['TIR']].split(',')[0])
        reference = total_depth - variant

    total_reads = reference + variant
    row['total_reads'] = total_reads
    row['ref_reads'] = reference
    row['var_reads'] = variant

    if total_reads > 0:
        row['VAF'] = row['var_reads'] / row['total_reads']
    else:
        row['VAF'] = np.nan

    return row


# Create order format similar to SigProfiler
def create_snv_class(df):
    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    x = df['TRIPLET']
    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])
    return out


# Get the major CN in the tumor. If it has no change, then we put the normal cn
def get_major_cn(df):
    if df['overlapp'] == 1:
        major = round(np.float(df['MAJOR_CN_TEMP']))

    else:
        major = round(np.float(df['NORMAL_CN']))
    return major


# Get the CN paths
def get_CNA_name(mutation_file):
    path_ascat = '{}/ASCAT/'.format('/'.join(mutation_file.split('/')[:-2]))
    sample_tumor = os.path.basename(mutation_file).split('_vs_')[0].replace('StrelkaBP_', '')

    path_CN = '{}/{}.cnvs.txt'.format(path_ascat, sample_tumor)
    path_purity = '{}/{}.purityploidy.txt'.format(path_ascat, sample_tumor)

    return path_CN, path_purity


# Get DBS
def get_DBS(df):
    forbidden_ids = []
    keep_dbs_records = []

    for chrom, data in df.groupby(by='CHROM'):

        data.sort_values(by='POS', inplace=True)
        index_done = data.index.tolist()
        data.reset_index(inplace=True)
        new_index = data.index.tolist()

        dic_ix = dict(zip(new_index, index_done))

        data['DIST'] = data['POS'].diff()

        closest_list = data[data['DIST'] == 1].index.tolist()
        forbidden_in_chrom = []

        for i in closest_list:
            row = data.loc[i]
            next_pos = i + 1
            if next_pos in closest_list:
                forbidden_in_chrom.append(i)
                forbidden_in_chrom.append(next_pos)
                forbidden_ids.append(dic_ix[i])
                forbidden_ids.append(dic_ix[i - 1])
                forbidden_ids.append(dic_ix[i + 1])

            if i not in forbidden_in_chrom:
                previous_mut = data.loc[i - 1]
                previous_mut['REF'] = '{}{}'.format(previous_mut['REF'], row['REF'])
                previous_mut['ALT'] = '{}{}'.format(previous_mut['ALT'], row['ALT'])
                keep_dbs_records.append(previous_mut)
                forbidden_ids.append(dic_ix[i])
                forbidden_ids.append(dic_ix[i - 1])

    all_ix = df.index.tolist()
    good_ix = [ix for ix in all_ix if ix not in forbidden_ids]
    SNVs_df = df.loc[good_ix]
    dbs_df = pd.DataFrame(keep_dbs_records)
    dbs_df.drop(['DIST', 'index'], axis=1, inplace=True)
    df = pd.concat([SNVs_df, dbs_df], sort=False)

    return df


# Classify indels
def indels_classification(indel_df):
    # Remove those cases where we have "," in the ALT
    indel_df = indel_df[~indel_df['ALT'].str.contains(',')]

    # set up first classification of variants
    indel_df['CLASS'] = indel_df.apply(
        lambda x: 'INS' if (x['len_alt'] > 1) & (x['len_ref'] == 1) else (
            'DEL' if (x['len_ref'] > 1) & (x['len_alt'] == 1) else (
                'DBS' if (x['len_alt'] == 2) & (x['len_ref'] == 2) else (
                    'MNV' if (x['len_alt'] == x['len_ref']) else (
                        'COMPLEX_INDELS')))), axis=1
    )

    complex_indels_df = indel_df[indel_df['CLASS'] == 'COMPLEX_INDELS']
    complex_indels_df['VARIANT_CLASS'] = 'COMPLEX_INDELS'
    mnv_indels_df = indel_df[indel_df['CLASS'] == 'MNV']
    mnv_indels_df['VARIANT_CLASS'] = 'MNV'

    # we won't be subclassifying these types
    toremove_class = ['MNV', 'COMPLEX_INDELS']

    # remove unwanted others
    indel_df = indel_df[~indel_df['CLASS'].isin(toremove_class)]

    # assing each class to the indels and dbs using PCAWG classification
    indel_df['VARIANT_CLASS'] = indel_df.apply(IndelsClassifier, axis=1, args=['hg38'])

    # DBS class
    indel_df1 = indel_df[indel_df['CLASS'] == 'DBS']
    indel_df1['CLASS'] = 'DBS'

    # INDEL class
    indel_df2 = indel_df[indel_df['CLASS'] != 'DBS']
    indel_df2['CLASS'] = 'INDEL'

    # merge both types
    indels = pd.concat([indel_df1, indel_df2, complex_indels_df, mnv_indels_df])

    indels.drop(['len_ref', 'len_alt'], axis=1, inplace=True)

    return indels


def get_total_CN(df, gender):
    if df['overlapp'] == 0:
        if df['CHROM'] not in ['chrY', 'chrX']:
            total_CN = 2
        else:
            if gender == 'MALE':
                total_CN = 1
            else:
                total_CN = 2
    else:
        total_CN = int(df['MAJOR_CN']) + int(df['MINOR_CN'])

    return total_CN


# Main processing script, this only works with the Sarek default output organisation, where indels and snvs are located
# in the same folder release
def process_strelka_calling(mutation_file, gender, outpath, path_mappable, file_gnomad):
    pd.options.mode.chained_assignment = None

    # Get indel files
    indel_file = mutation_file.replace('somatic_snvs.vcf.gz', 'somatic_indels.vcf.gz')

    # Get CNVS and purity from ASCAT
    cnvs_file, purity_file = get_CNA_name(mutation_file)

    sample_file_name = os.path.basename(mutation_file).replace('StrelkaBP_', '').split('_somatic_snvs')[0]

    sample = '{}_{}'.format(sample_file_name.split('_vs_')[0], sample_file_name.split('_vs_')[1])

    patient = mutation_file.split('/')[-7]
    df = vcf_reader(mutation_file)
    df['FCLASS'] = 'SNV'

    df_indel = vcf_reader(indel_file)
    df_indel['FCLASS'] = 'INDEL'

    # get only canonical chromosomes
    wantedchroms = ['chr{}'.format(i) for i in range(1, 23)]
    wantedchroms.append('chrY')
    wantedchroms.append('chrX')

    # select only variants with the PASS filter and in the chromosomes we are interested in
    df = df[df['FILTER'] == 'PASS']
    df = df[df['CHROM'].isin(wantedchroms)]

    df_indel = df_indel[df_indel['FILTER'] == 'PASS']
    df_indel = df_indel[df_indel['CHROM'].isin(wantedchroms)]

    # read ASCAT copy number estimations
    df_cns = pd.read_csv(cnvs_file, sep='\t')
    df_cns['chr'] = df_cns['chr'].apply(lambda x: 'chr{}'.format(x))
    df_cns.to_csv('{}/{}.{}.CN.gz'.format(outpath, patient, sample), sep='\t', index=False, header=True,
                  compression='gzip')

    cnv_bed = BedTool.from_dataframe(df_cns)
    df_purity = pd.read_csv(purity_file, sep='\t')
    df_purity.to_csv('{}/{}.{}.purity.gz'.format(outpath, patient, sample), sep='\t', index=False, header=True,
                     compression='gzip')
    purity_score = np.float(df_purity['AberrantCellFraction'].tolist()[0])

    # get double-base substitutions
    df = get_DBS(df)

    df_merged = pd.concat([df, df_indel])

    df_reads = df_merged.apply(return_reads, axis=1)

    # select whether we have SNVs or others
    df_reads['len_alt'] = df_reads['ALT'].str.len()

    # number of characters in ref
    df_reads['len_ref'] = df_reads['REF'].str.len()

    # first classification between SNV and others
    df_reads['TYPE'] = df_reads.apply(
        lambda x: 'SNV' if ((x['len_alt'] == 1) and (x['len_ref'] == 1) and (x['ALT'] != '-') and (
                x['REF'] != '-')) else 'INDEL', axis=1
    )

    df_reads['pos-1'] = df_reads['POS'] - 1

    # get the triplet
    df_reads['TRIPLET'] = df_reads.apply(lambda x: hg38(x['CHROM'], x['pos-1'], 3), axis=1)
    df_reads['EXTENDED'] = df_reads.apply(lambda x: hg38(x['CHROM'], int(x['POS']) - 2, 5), axis=1)

    snv_df = df_reads[df_reads['TYPE'] != 'INDEL']
    snv_df['CLASS'] = 'SNV'
    snv_df['VARIANT_CLASS'] = snv_df.apply(create_snv_class, axis=1)

    # classify indel0s
    indel_df = df_reads[df_reads['TYPE'] == 'INDEL']
    indels = indels_classification(indel_df)
    columns = indels.columns

    df_reads_merged = pd.concat([snv_df, indels], sort=True)
    df_reads_merged = df_reads_merged[columns]

    df_reads_merged['sample'] = sample

    # create bed file
    mut_bed = BedTool.from_dataframe(df_reads_merged[[
        'CHROM', 'pos-1', 'POS', 'ref_reads', 'var_reads', 'VAF', 'total_reads', 'REF',
        'ALT', 'sample', 'TYPE', 'CLASS', 'VARIANT_CLASS', 'TRIPLET', 'EXTENDED'
    ]])

    mut_bed = get_mappable_regions(mut_bed, path_mappable)

    # intersect with CN data
    out = mut_bed.intersect(cnv_bed, wao=True)

    # merge to dataframe
    merge = out.to_dataframe(names=[
        'CHROM', 'POS-1', 'POS', 'REF_COUNTS', 'VAR_COUNTS', 'VAF', 'TOTAL_READS', 'REF', 'ALT', 'SAMPLE', 'TYPE',
        'CLASS', 'VARIANT_CLASS', 'TRIPLET', 'EXTENDED', 'multiallelism',
        'c1', 'p1', 'p2', 'MAJOR_CN', 'MINOR_CN', 'overlapp'
    ])

    # get the normal copy number values
    sex_chrom = ('Y', 'X')

    # get normal CN in the chromosome
    merge['NORMAL_CN'] = merge['CHROM'].apply(lambda x: 1 if x in sex_chrom and gender == "MALE" else 2)

    # add the purity score we got from PURPLE
    merge['PURITY'] = purity_score
    merge['GENDER'] = gender

    merge['TOTAL_CN'] = merge.apply(get_total_CN, axis=1, args=[gender])

    # gnomAD annotation
    merge = filter_gnomad_38(merge, file_gnomad)

    # Clonality according to McGranaham et al
    merge = merge.apply(getClonalStatus, axis=1)

    merge['PATIENT'] = str(patient)

    # save the file
    merge.to_csv('{}/{}.{}.muts.gz'.format(outpath, patient, sample), sep='\t', index=False, header=True,
                 compression='gzip')
