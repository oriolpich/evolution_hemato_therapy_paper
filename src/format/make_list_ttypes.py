import sys
import pandas as pd
import json

sys.path.append("../config")
from config import LIST_TTYPE, SAMPLE_2_PATIENT, GRU_PHENOTYPES, GRU, DS_HEM

df_meta = pd.read_csv(GRU_PHENOTYPES, sep='\t')

prior_neoplasia_in_AML = [str(x) for x in
                          df_meta[df_meta['Prior diagnosis of prior neoplasm?'] == 'Yes']['SUBJID'].tolist()]

# Read dbGAP files
df_hem = pd.read_csv(DS_HEM, sep=',', comment='#')
df_gru = pd.read_csv(GRU, sep=',', comment='#')

hem_gru = pd.concat([df_hem, df_gru], sort=False)
wanted = ['tAMLTumor', 'secondary AML tumor']
secondary_AML = hem_gru[(hem_gru['histological_type'].isin(wanted)) & (hem_gru['Assay Type'] == 'WGS')][
    'submitted_subject_id'].tolist()

full_tAML_dbgap = secondary_AML + prior_neoplasia_in_AML

secondary_full = hem_gru[hem_gru['submitted_subject_id'].isin(full_tAML_dbgap)]

dic_sex = {'female': 'XX', 'male': 'XY', 'nan': 'XX'}
samples_good = set(secondary_full['submitted_subject_id'].tolist())
dic_equiv_patient_sample = {}
set_tAML_samples = set()

for sample, data in hem_gru[hem_gru['submitted_subject_id'].isin(samples_good)].groupby(by='submitted_subject_id'):

    normal = data[data['is_tumor'] == 'No']['Run'].tolist()[0]
    tumoral = data[data['is_tumor'] == 'Yes']['Run'].tolist()[0]

    sex = dic_sex.get(data[data['is_tumor'] == 'No']['sex'].tolist()[0], 'XX')

    set_tAML_samples.add('{}_{}'.format(tumoral, normal))
    dic_equiv_patient_sample['{}_{}'.format(tumoral, normal)] = sample

wanted = ['Primary Tumor', 'primary tumor', 'primary-tumor', 'Tumor']
primary_AML = hem_gru[(hem_gru['histological_type'].isin(wanted)) & (hem_gru['Assay Type'] == 'WGS')]['submitted_subject_id']
primary_full = hem_gru[(hem_gru['submitted_subject_id'].isin(primary_AML)) & (hem_gru['Assay Type'] == 'WGS')]

tumors_to_analyze = set([m for m, v in primary_full['submitted_subject_id'].value_counts().to_dict().items() if v < 4])
primary_full = hem_gru[hem_gru['submitted_subject_id'].isin(tumors_to_analyze)]

other_primaries = {
    '182896': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '288033': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '308955': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Tumor'},
    '400220': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '461282': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '548327': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '573988': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '758168': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '804168': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '831711': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Primary Tumor'},
    '852559': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Tumor'},
}

relapse = {
    # we remove the first because it is WXS
    # '308955':{'healthy':'Solid Tissue Normal', 'tumoral':'Relapse Tumor'},
    '400220': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Relapse Tumor'},
    '573988': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Relapse Tumor'},
    '758168': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Relapse Tumor'},
    '804168': {'healthy': 'Solid Tissue Normal', 'tumoral': 'Relapse Tumor'},

}

to_run = []

primary_set = set()
relapse_set = set()

for sample, data in primary_full[primary_full['Assay Type'] == 'WGS'].groupby(by='submitted_subject_id'):

    if len(data) == 2:
        if len(data[data['is_tumor'] == 'No']) & len(data[data['is_tumor'] == 'Yes']):
            normal = data[data['is_tumor'] == 'No']['Run'].tolist()[0]
            sex = dic_sex.get(data[data['is_tumor'] == 'No']['sex'].tolist()[0], 'XX')
            tumoral = data[data['is_tumor'] == 'Yes']['Run'].tolist()[0]
            primary_set.add('{}_{}'.format(tumoral, normal))
            dic_equiv_patient_sample['{}_{}'.format(tumoral, normal)] = sample
    else:
        if len(data) > 1:
            if sample in relapse:
                normal = data[data['histological_type'] == relapse[sample]['healthy']]['Run'].tolist()[0]
                sex = dic_sex.get(data[data['is_tumor'] == 'No']['sex'].tolist()[0], 'XX')
                tumoral = data[data['histological_type'] == relapse[sample]['tumoral']]['Run'].tolist()[0]
                relapse_set.add('{}_{}'.format(tumoral, normal))
                dic_equiv_patient_sample['{}_{}'.format(tumoral, normal)] = sample

# Label as primary if not in relapse and not in tAML
final_primary = [s for s in primary_set if s not in set_tAML_samples and s not in relapse_set]
# Label as relapse if not in tAML
final_relapse = [s for s in relapse_set if s not in set_tAML_samples]

# Add the samples from other cohorts, first our samples
final_tAML = list(set_tAML_samples) + ['AP0388_AP0389', 'AP0390_AP0395', 'AP0394_AP0391']

# Then the WTSI samples
wtsi_samples = ['PATIENT_WTS_1_TUM1_PATIENT_WTS_1_BUC',
                'PATIENT_WTS_3_TUM1_PATIENT_WTS_3_BUC']

final_tAML = final_tAML + wtsi_samples

dic = {'tAML': final_tAML, 'primary': final_primary, 'relapse': final_relapse}

# Dump the json
json.dump(dic, open(LIST_TTYPE, 'wt'),  indent=4)

dic_equiv_patient_sample['AP0388_AP0389'] = 'SUBJECT_1'
dic_equiv_patient_sample['AP0390_AP0395'] = 'SUBJECT_TUMOR_2'
dic_equiv_patient_sample['AP0394_AP0391'] = 'SUBJECT_TUMOR_3'
dic_equiv_patient_sample['PATIENT_WTS_1_TUM1_PATIENT_WTS_1_BUC'] = 'PATIENT_WTS_1'
dic_equiv_patient_sample['PATIENT_WTS_3_TUM1_PATIENT_WTS_3_BUC'] = 'PATIENT_WTS_3'

# Dump the ID equivalences
json.dump(dic_equiv_patient_sample, open(SAMPLE_2_PATIENT, 'wt'),  indent=4)
