import pandas as pd 
import json

f = '/workspace/datasets/reverse_calling/AML_Sequencing_Project/phs000159/test/GRU_phenoGenotypeFiles/GRU_phenotypes.txt'
df_meta = pd.read_csv(f, sep ='\t')

df_hem = pd.read_csv('/workspace/datasets/reverse_calling/AML_Sequencing_Project/phs000159/test/run_info_table/SraRunTable_DS-HEM.txt', 
                sep =',',  comment='#')

df_gru = pd.read_csv('/workspace/datasets/reverse_calling/AML_Sequencing_Project/phs000159/test/run_info_table/SraRunTable_GRU.txt', 
                sep =',',  comment='#')

df = pd.concat([df_hem, df_gru], sort = False)

dic_timing = dict(zip(df_meta['SUBJID'], df_meta['Age at Initial Banking']))

dic_timing_f = {}
for sample, data in df.groupby(by='submitted_subject_id'):
    try:
        dic_timing_f[sample] = dic_timing[int(sample)]
    except:
        continue
        
# add wong dataset
df_wong = pd.read_csv('/workspace/datasets/reverse_calling/AML_Sequencing_Project/Wong_2014/dataset_wgs.txt', 
                     sep ='\t')
dic_tAML = dict(zip(df_wong['UPN'].astype(str), df_wong['Age at            t-AML/t-MDS Diagnosis']))
dic_timing_f.update(dic_tAML)

# add sant pau dataset
df_stpau= pd.read_csv('/workspace/datasets/reverse_calling/CNAG/dataset_wgs.txt', 
                     sep ='\t')

dic_tAML_stpau = dict(zip(df_stpau['UPN'].astype(str), df_stpau['Age at            t-AML/t-MDS Diagnosis']))
dic_timing_f.update(dic_tAML_stpau)

# add HSC data
meta_hsc = pd.read_csv('/workspace/projects/all_aecc/Osorio_data/metadata_samples.txt', sep ='\t')
dic_HSC = dict(zip(meta_hsc['Identifier'], meta_hsc['Age (years)']))
dic_timing_f.update(dic_HSC)

# add latency
dic_lat = {}
dic_latency = dict(zip(df_wong['UPN'].astype(str), df_wong['Latency (years)']))
dic_latency_stpau = dict(zip(df_stpau['UPN'].astype(str), df_stpau['Latency (years)']))
dic_lat.update(dic_latency)
dic_lat.update(dic_latency_stpau)


dic_equiv_patient_sample = json.load(open('/workspace/projects/reverse_calling/data/samples/sample_2_patient.json'))

dic_aging_sample_level = {}
for sample, patient in dic_equiv_patient_sample.items():
    dic_aging_sample_level[sample] = dic_timing_f.get(patient, 'no_data')
    
dic_latency_sample_level = {}
for sample, patient in dic_equiv_patient_sample.items():
    dic_latency_sample_level[sample] = dic_lat.get(patient, 'no_data')


dic_aging_sample_level.update(dic_HSC)

json.dump(dic_aging_sample_level, open('/home/opich/bg/tAML_code_review/evolution_hemato_therapy/data/samples/sample_2_age.json', 
                                      'wt'))