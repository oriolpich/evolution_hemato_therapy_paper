from collections import defaultdict
import pandas as pd 
import pickle 
import gzip 


composed_drugs = {
    
    'Fulverstrant':['Fulvestrant'],
    
    'Capecitabine, Ixabepilone':[
        'Capecitabine', 'Ixabepilone'
    ],
    'Fluoropyrimidine':['Fluorouracil'],
    
    'Fluorouracil/Irinotecan/Leucovorin/Oxaliplatin' : [
                'Fluorouracil', 'Leucovorin', 'Irinotecan', 'Oxaliplatin'

    ],
    
    'Fluorouracil, Leucovorin Calcium (FU/LV)' : ['Fluorouracil', 'Leucovorim'],
    'Docetaxel, Doxorubicin, Cyclophosphamide (TAC)': [
        'Docetaxel', 'Doxorubicin', 'Cyclophosphamide'
    ],
    'Dexamethasone, High-dose cytarabine, Cisplatin (DAHP) - Ifosfamide, mitoxantrone and etoposide (VIM) - Dexamethasone, High-dose cytarabine, Cisplatin (DAHP)': [
        'Dexamethasone', 'Cytarabine', 'Cisplatin', 'Ifosfamide', 'Mitoxantrone', 'Etoposide'
    ],
    'Fluorouracil, Epirubicin, Cyclophosphamide, Docetaxel (FEC-D)': [
        'Fluorouracil', 'Epirubicine', 'Cyclophosphamide', 'Docetaxel'
    ],
    'Fluorouracil, Epirubicine and Cyclofosfamide (FEC)': [
        'Fluorouracil', 'Epirubicine', 'Cyclofosfamide'
    ],
    'Fluorouracil, Leucovorin Calcium (FU': [
        'Fluorouracil', 'Leucovorin'
    ],
    'Tegafur, Gimeracil, Oteracil (S1)': [
        'Tegafur', 'Gimeracil', 'Oteracil'
    ],
    'Procarbazine, Lomustine, Vincristine (PCV)': [
        'Procarbazine', 'Lomustine', 'Vincristine'
    ],
    'Doxorubicin, Bleomycin, Vinblastine, Dacarbazin (ABVD)': [
        'Doxorubicin', 'Bleomycin', 'Vinblastine', 'Dacarbazine'
    ],
    'Doxorubicin, Cyclophosphamide (AC)': [
        'Doxorubicin', 'Cyclophosphamide'
    ],
    'Pemetrexed, Carboplatin': [
        'Premetrexed', 'Carboplatin'
    ],
    'leucovorin, fluorouracil, irinotecan (FOLFIRI)': [
        'Leucovorin', 'Fluorouracil', 'Irinotecan'
    ],
    'leucovorin, fluorouracil, oxaliplatin (FOLFOX)': [
        'Leucovorin', 'Fluorouracil', 'Oxaliplatin'
    ],
    'mustine, Vincristine Sulfate, Procarbazine Hydrochloride, Prednisone (MOPP)': [
        'Limustine', 'Vincristine', 'Procarbazine', 'Prednisone'
    ],
    'Vincristine Sulfate, Dactinomycin, Cyclophosphamide (VAC)': [
        'Vincristine', 'Dactinomycin', 'Cyclophosphamide'
    ],
    'Vincristine Sulfate, Etoposide Phosphate, Prednisone, Doxorubicin Hydrochloride (OEPA)': [
        'Vincristine', 'Etoposide', 'Prednisone', 'Doxorubicin'
    ],
    'Cabazitaxel, Zoledronic Acid': [
        'Cabazitaxel', 'Zoledronic Acid'
    ],
    'Capecitabine, Oxaliplatin (CAPOX)': [
        'Capecitabine', 'Oxaliplatin'
    ],
    'Carboplatin, Paclitaxel (PC)': [
        'Carboplatin', 'Paclitaxel'
    ],
    'Carmustine, Etoposide, Cytarabine, Melphalan (BEAM)': [
        'Carmustine', 'Etoposide', 'Cytarabine', 'Melphalan'
    ],
    'Cyclophosphamide, Methotrexate, Fluorouracil (CMF)': [
        'Cyclophosphamide', 'Methotrexate', 'Fluorouracil'
    ],
    'Cyclophosphamide, Vincristine Sulfate, Procarbazine Hydrochloride, Prednisone (COPP)': [
        'Cyclophosphamide', 'Vincristine', 'Procarbazine', 'Prednisone',
    ],
    'Rituximab, Cyclophosphamide, Doxorubicin, Vincristine, Prednisolone (R-CHOP)': [
        'Rituximab', 'Cyclophosphamide', 'Doxorubicin', 'Vincristine', 'Prednisolone'
    ],
    'ABVD': [
        'Doxorubicin', 'Bleomycine', 'Vinblastine', 'Dacarbazine'
    ],
    'AC': [
        'Doxorubicin', 'Cyclophosphamide'
    ],
    'Ad[1': [
        'Doxorubicin'
    ],
    'Adriamycin': [
        'Doxorubicin'
    ],
    'BEAM': [
        'Carmustine', 'Etoposide', 'Cytarabine', 'Melphalan'
    ],
    'CMF': [
        'Cyclophosphamide', 'Methotrexate', 'Fluorouracil'
    ],
    'COPP': [
        'Cyclophosphamide', 'Vincristine', 'Procarbazine', 'Prednisone'
    ],
    'Capox': [
        'Capecitabine', 'Oxaliplatin'
    ],
    'DHAP-VIM-DHAP': [
        'Rituximab', 'Dexamethason', 'Cytarabine', 'Cisplatin', 'Etoposide', 'Ifosfamide', 'Methotrexate'
    ],
    'FEC': [
        'Fluorouracil', 'Epirubicine', 'Cyclophosphamide'
    ],
    'FEC-D': [
        'Fluorouracil', 'Epirubicine', 'Cyclophosphamide', 'Docetaxel'
    ],
    'Folfirinox': [
        'Fluorouracil', 'Irinotecan', 'Oxaliplatin'
    ],
    'MOPP': [
        'Nitrogen_mustard', 'Oncovin', 'Procarbazine', 'Prednisone'
    ],
    'OEPA': [
        'Oncovin', 'Etoposide', 'Prednisone', 'Doxorubicin'
    ],
    'PCV': [
        'Procarbazine', 'Oncovin', 'Lomustine'
    ],
    'R-CHOP': [
        'Rituximab', 'Cyclophosphamide', 'Doxorubicin', 'Oncovin', 'Prednisolone'
    ],
    'S1': [
        'Tegafur', 'Gimeracil', 'Oteracil'
    ],  
    'SYD985': [
        'Trastuzumab', 'Duocarmycin',
    ],
    'TAC': [
        'Docetaxel', 'Doxorubicin', 'Cyclophosphamide'
    ],
    'Taxol': [
        'Paclitaxel'
    ],
    'VAC': [
        'Oncovin', 'Dactinomycin', 'Cyclophosphamide'
    ],
    'Zoladex': [
        'Gosereline'
    ],
    'Casodex': [
        'Bicalutamide'
    ],
    'Eligard': [
        'Leuproreline'
    ],
    'Carbotaxol': [
        'Paclitaxel', 'Carboplatin'
    ],
    'LV': [
        'Fluorouracil', 'Folic acid'
    ],
    'Abirateron, Prednison': [
        'Abirateron', 'Prednisone'
    ],
    'Carboplatin, Gemcitabine': [
        'Carboplatin', 'Gemcitabine'
    ],
    'GDC': [
        'Gemcitabine', 'Docetaxel', 'Carboplatin'
    ],
    'GDC-0032': [
        'Gemcitabine', 'Docetaxel', 'Carboplatin'
    ],
    'Letrozol (+ LEE01': [
        'Letrozole', 'LEE01'],
    'Fluorouracil, Epirubicin, Cyclophosphamide (FEC)': [
        'Fluorouracil', 'Epirubicin', 'Cyclophosphamide'
    ],
    'Mechlorethamine Hydrochloride, Vincristine Sulfate, Procarbazine Hydrochloride, Prednisone (MOPP)': [
        'Mechlorethamine', 'Vincristine', 'Procarbazine', 'Prednisone'
     ],
    'Doxorubicin, Bleomycin, Vinblastine, Dacarbazine (ABVD)': [
        'Doxorubicin', 'Bleomycin', 'Vinblastine', 'Dacarbazine'
    ],
    'Leucovorin, fluorouracil, irinotecan (FOLFIRI)': [
        'Leucovorin', 'Fluorouracil', 'Irinotecan'
    ],
    'Abiraterone, Prednisone': [
        'Abiraterone', 'Prednisone'
    ],
    'Letrozole (+ LEE011 or placebo in Monaleesa trial)': [
      'Letrozole', 'Unknown'
    ],
    'Leucovorin, fluorouracil, oxaliplatin (FOLFOX)': [
        'Leucovorin', 'Fluorouracil', 'Oxaliplatin'
    ],
    'Unknown (neoadjuvant chemotherapy)' :['Unknown'], 
    
    'Unknown chemotherapy for breast tumor':['Unknown'], 
    'Taselisib or placebo':['Unknown'], 
    'Somatostatin':['Ocreotide'], 
    'Olaratumab or placebo':['Unknown'], 
    'UNK':['Unknown'],
    'Unknown adjuvant chemotherapy':['Unknown'],
    
    'Ribociclib or placebo':['Unknown'],
 'Dendritic cell therapy (DCVAC) or placebo':['Unknown'],
 'LEE011 or placebo':['Unknown'],
 'GDC or placebo':['Unknown'],
 'Olaparib or placebo':['Unknown'],
 'Atezolizumab or placebo':['Unknown'],
 'Olaratumab or placebo':['Unknown'],
 'Radium-223 or placebo':['Unknown'],
 'Pembrolizumab or placebo':['Unknown'],
 'Nivolumab or placebo':['Unknown'],
 'Orteronel or placebo':['Unknown'],
 'CDK 4/6 inhibitor or placebo in monaleesa study':['Unknown'],
 'GDC-0032 or placebo':['Unknown'],
 'Taselisib or placebo':['Unknown'],
    'Unknown':['Unknown'],
    'Unknown chemotherapy':['Unknown'],
    'Unknown adjuvant chemotherapy':['Unknown'],
    'Hyperthermic Intraperitoneal Chemotherapy (HIPEC) - COLOPEC':['Hyperthermic Intraperitoneal Chemotherapy (HIPEC)'],

    'Dexamethasone (liposomal)':['Dexamethasone'],
    'Palbociclib, Fulvestrant':['Palbociclib', 'Fulverstrant'],
    'Letrozole, Lapatinib':['Letrozole', 'Lapatinib'],
    'Doxorubicin (liposomal)' : ['Doxorubicin'],
    'Insulin-like growth factor (IGF) inhibitor (study)':['Unknown'],  
            
}

# get the date
def create_date(x):
    if (str(x) !=
        'nan') and (str(x) != 'unknown'):
        date1 = date(int(x.split('-')[0]), int(x.split('-')[1]), int(x.split('-')[2]))
    else:
        date1 = 'unknown'
    return date1

def process_dates_treament(startDate,endDate,biopsy_date):
    
    # get the days treated, days since start treatmend and days since end treatment

    if endDate == 'unknown':
        t_days = 'unknown'
        days_treat_2_biopsy = 'unknown'
        days_since_start_treatment = 'unknown'

    # if treatment ends before biopsy
    elif endDate < biopsy_date:
        days_treat = endDate-startDate
        t_days = days_treat.days
        days_treat_2_biopsy = (biopsy_date-endDate).days
        days_since_start_treatment = (biopsy_date-startDate).days

    else:
        days_treat = biopsy_date-startDate
        t_days = days_treat.days
        days_treat_2_biopsy = 0
        days_since_start_treatment = (biopsy_date-startDate).days

    return t_days, days_treat_2_biopsy, days_since_start_treatment


path_metadata = '/workspace/datasets/hartwig/20200117/metadata_update/metadata/metadata.tsv'
path_treatments = '/workspace/datasets/hartwig/20200117/metadata_update/metadata/pre_biopsy_drugs_by_patient.tsv'
treat = pd.read_csv(path_treatments, sep ='\t',  encoding='latin-1')
pat = pd.read_csv(path_metadata, sep='\t', encoding='latin-1')

pat['primaryTumorLocation'] = pat['primaryTumorLocation'].replace('Bone/soft tissue', 'Bone/Soft tissue')

# fix primary location
pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation'].apply(
    lambda x: str(x).replace(' ', '-').replace('/', '-')
)
pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('Head-and-Neck', 'Head-and-neck')
pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('nan', 'Unknown')
pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('CUP', 'Unknown')

 
no_biopsy_date = pat[pat['biopsyDate'].isnull()]
with_biopsy_date = pat[~pat['biopsyDate'].isnull()]


keep_treated = defaultdict(int)
treatment_sample_level = defaultdict(lambda: defaultdict(set))

treatment_sample_days_treatment = defaultdict(lambda: defaultdict(list))
treatment_sample_days_biopsy_since_end_treatment = defaultdict(lambda: defaultdict(list))
treatment_sample_days_biopsy_since_start_treatment = defaultdict(lambda: defaultdict(list))

for patient, data in with_biopsy_date.groupby(by='#patientId'):
    
    treated = treat[treat['#patientId']==patient]
    
    # if the patient has been treated
    if len(treated):
        
        for ix, treatm in treated.iterrows():
            
            startDate = create_date(treatm['startDate'])
            endDate = create_date(treatm['endDate'])
            drug =  treatm['name']
            
            # loop over possible samples treated
            for i, row in data.iterrows():
                sample = row['sampleId']
                ttype = row['primaryTumorLocation_fixed']
                biopsy_date = create_date(row['biopsyDate'])
                
                # if the start of the treatment is before the biopsy date
                if (startDate < biopsy_date):
                    
                    # when subtracting days, the longest date is the one with more days
                    t_days, days_treat_2_biopsy, days_since_start_treatment = process_dates_treament(startDate,endDate,biopsy_date)
                    
                    if drug in composed_drugs:
                        for small_d in composed_drugs[drug]:
                            treatment_sample_level[ttype][small_d].add(sample)
                            treatment_sample_level['PAN'][small_d].add(sample)
                            
                            treatment_sample_days_treatment[sample][small_d].append(t_days)
                            treatment_sample_days_biopsy_since_end_treatment[sample][small_d].append(days_treat_2_biopsy)
                            treatment_sample_days_biopsy_since_start_treatment[sample][small_d].append(days_since_start_treatment)

                    else:

                        treatment_sample_level[ttype][drug].add(sample)
                        treatment_sample_level['PAN'][drug].add(sample)
                        treatment_sample_days_treatment[sample][drug].append(t_days)
                        treatment_sample_days_biopsy_since_end_treatment[sample][drug].append(days_treat_2_biopsy)
                        treatment_sample_days_biopsy_since_start_treatment[sample][drug].append(days_since_start_treatment)

     
pickle.dump(dict(treatment_sample_level), gzip.open('treatment_hartwig_drug_level.pckl.gz', 'wb'))
pickle.dump(dict(treatment_sample_days_treatment), gzip.open(
    '../../data/signatures/hartwig/treatment_hartwig_drug_level_days_treated.pckl.gz', 'wb'))
pickle.dump(dict(treatment_sample_days_biopsy_since_end_treatment), gzip.open(
    '../../data/signatures/hartwig/treatment_hartwig_drug_level_since_end_treat.pckl.gz', 'wb'))
pickle.dump(dict(treatment_sample_days_biopsy_since_start_treatment), gzip.open(
    '../../data/signatures/hartwig/treatment_hartwig_drug_level_since_start_treat.pckl.gz', 'wb'))