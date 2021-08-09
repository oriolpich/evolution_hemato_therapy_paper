import os
import gzip
import tabix
import sys


def return_gnomad_stats(chrom, pos, ref, alt, tb):
    AC_count = 0
    # Some queries crash
    try:
        records = tb.query(chrom, pos, pos)
        all_records = [r for r in records]
        if len(all_records) > 0:
            for record in all_records:
                if len(record) > 0:
                    chrom_t = record[0]
                    pos_t = int(record[1])
                    ref_t = record[3]
                    alt_t = record[4]

                    if (chrom == chrom_t) & (pos == pos_t) & (ref == ref_t) & (alt == alt_t):
                        AC_count = record[7].split('AC=')[1].split(';')[0]

    except ValueError:
        pass

    return AC_count


def apply_gnomad_filter(df, tb):
    chrom = df['CHROM'].replace('chr', '')
    pos = df['POS']
    ref = df['REF']
    alt = df['ALT']
    AC_gnomad = return_gnomad_stats(chrom, pos, ref, alt, tb)

    return AC_gnomad


def filter_gnomad_38(df, file_gnomad):
    tb = tabix.open(file_gnomad)
    AC_gnomad_list = df.apply(apply_gnomad_filter, args=(tb,), axis=1)
    df['gnomAD'] = AC_gnomad_list
    return df
