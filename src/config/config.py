# Path to all figures
FIG_1 = './fig1'
FIG_2 = './fig2'
FIG_3 = './fig3'
FIG_4 = './fig4'

# Relative paths for data files
DIR_SAMPLES = '../../data/samples'

# These files are obtained from the dbGAP portal when you are granted acces
GRU_PHENOTYPES = '../../data/AML_Seq_phs000159/GRU_phenoGenotypeFiles/GRU_phenotypes.txt'
GRU = '../../data/AML_Seq_phs000159/run_info_table/SraRunTable_GRU.txt'
DS_HEM = '../../data/AML_Seq_phs000159/run_info_table/SraRunTable_DS-HEM.txt'

# These are all the mutations from AMLs after formatting. The file is created using the formatting script
ALL_MUTS = '../../data/samples/AML.all.muts.gz'
NEW_ALL_MUTS = '../../data/samples/new_AML.all.muts.gz'

# Information about the datasets
WSG_WONG = '../../data/dataset_wgs_Wong_2014.txt'
WSG_CNAG = '../../data/dataset_wgs_CNAG.txt'
WSG_EXTRA = '../../data/extra_dataset_wgs_AML.txt'
WSG_WTSI = '../../data/extra_dataset_wgs_wtsi.txt'

# Dictionaries with clinical data created with the repo scripts
SAMPLE_2_PATIENT = '../../data/samples/sample_2_patient.json'
SAMPLE_2_AGE = '../../data/samples/sample_2_age.json'
LIST_TTYPE = '../../data/samples/list_ttype.json'

# Result of SigProfilerJulia
SIG_ALL_AML_EXPOSURES = '../../data/signatures/AML_all_hemato.snvs.exposures.tsv'
SIG_ALL_AML_PROCESSES = '../../data/signatures/AML_all_hemato.snvs.processes.tsv'
SIG_ALL_AML_INPUT = '../../data/signatures/AML_all.snvs.tsv'

# Signatures from COSMIC and SBS
SIG_PROFILER_SBS = '../../data/software/mSigAct_pipeline/sigProfiler_SBS_signatures_2018_03_28.csv'
SIG_PROFILER_SBS_DECONS = '../../data/software/mSigAct_pipeline/SigProfiler_SBS_signatures_2018_03_28.deconstructsigs.tsv'
PATH_MSIG = '../../data/software/mSigAct_pipeline/'
REGRESSION_INPUT = '../analysis/reg_input.tsv'

ML_SIGNATURES = '../../data/ML'

AML_ALL_HEMATO_ML = '../../data/ML/AML_all_hemato.ML.SNV.gz'

TREATMENT_DIC = '../../data/samples/treatment_dic.json'

UNTREATED_SAMPLES = '../../data/untreated_samples.json.gz'

# Hartwig data
HARTWIG_METADATA = '../../data/hartwig_metadata/metadata.tsv'
SIG_TREATMENT_DRUG_LVL = '../../data/signatures/hartwig/treatment_hartwig_drug_level.pckl.gz'
SIG_TREATMENT_DRUG_LVL_DAYS = '../../data/signatures/hartwig/treatment_hartwig_drug_level_days_treated.pckl.gz'
SIG_TREATMENT_DRUG_LVL_END = '../../data/signatures/hartwig/treatment_hartwig_drug_level_since_end_treat.pckl.gz'
SIG_TREATMENT_DRUG_LVL_START = '../../data/signatures/hartwig/treatment_hartwig_drug_level_since_start_treat.pckl.gz'

# Results of signature exposures of several metastatic cohorts from HMF
SIG_EXTRACT_OVARY = '../../data/mutfootprints_2020/hartwig/signatures/extraction/Ovary.snvs.exposures.tsv'
SIG_EXTRACT_URINARY = '../../data/mutfootprints_2020/hartwig/signatures/extraction/Urinary-tract.snvs.exposures.tsv'
SIG_EXTRACT_COLON = '../../data/mutfootprints_2020/hartwig/signatures/extraction/Colon-Rectum.snvs.exposures.tsv'
SIG_EXTRACT_LUNG = '../../data/mutfootprints_2020/hartwig/signatures/extraction/Lung.snvs.exposures.tsv'
SIG_EXTRACT_BREAST = '../../data/mutfootprints_2020/hartwig/signatures/extraction/Breast.snvs.exposures.tsv'

# MSIG ACTS
M_SIG_ACT_OBSERVED_UNTREATED = '../../data/msigact_observed_untreated/'
M_SIG_ACT_CATALOGUE = 'results.catalogue_*_0.tsv.mSigAct.*.tsv'
RESULT_CATALOGUE = 'results.catalogue_*.tsv'
BASIC_FILTERS = '../../data/output-basic-filters'
M_SIG_ACT_BASIC_FILTERS = '../../data/msigact_basic_filters/'
RESULTS_BASIC_FILTERS = 'results.basic_filters.muts.tsv.mSigAct.*.tsv'

SIG_RESULT_CR_CAPECITABINE = '../../data/mutfootprints_2020/hartwig/signatures/mSigAct_results/Colon-Rectum' \
                             '-Capecitabine/Colon-Rectum_muts.tsv/results.Colon-Rectum_muts.tsv.mSigAct.X2_1.0_SBS17b' \
                             '.0.96.tsv'
SIG_RESULT_CR_OXALIPLATIN = '../../data/mutfootprints_2020/hartwig/signatures/mSigAct_results/Colon-Rectum' \
                            '-Oxaliplatin/Colon-Rectum_muts.tsv/results.Colon-Rectum_muts.tsv.mSigAct.X10_0.99_NA.tsv'
SIG_RESULT_O_CAPECITABINE = '../../data/mutfootprints_2020/hartwig/signatures/mSigAct_results/Ovary-Carboplatin' \
                            '/Ovary_muts.tsv/results.Ovary_muts.tsv.mSigAct.X6_1.0_SBS31.0.92.tsv'

SIG_RESULT_UT_PLATINUM = '../../data/mutfootprints_2020/hartwig/signatures/mSigAct_results/Urinary-tract-Platinum' \
                         '/Urinary-tract_muts.tsv/results.Urinary-tract_muts.tsv.mSigAct.X7_1.0_SBS31.0.97.tsv'

SIG_RESULT_BR_CAPECITABINE = '../../data/mutfootprints_2020/hartwig/signatures/mSigAct_results/Breast-Capecitabine' \
                             '/Breast_muts.tsv/results.Breast_muts.tsv.mSigAct.Capecitabine.tsv'

SIG_RESULT_BR_PLATINUM = '../../data/mutfootprints_2020/hartwig/signatures/mSigAct_results/Breast-Platinum' \
                         '/Breast_muts.tsv/results.Breast_muts.tsv.mSigAct.SBS31.tsv'

COLON_MUTS = '../../data/mutfootprints_2020/hartwig/signatures/colon-rectum/Colon-Rectum.snvs.dlm'
COLON = '../../data/mutfootprints_2020/hartwig/signatures/colon-rectum/Colon-Rectum.snvs.processes.tsv'

M_SIG_ACT = '../../data/mutfootprints_2020/hartwig/signatures/mSigAct'

ML_OVARY = '../../data/ML/Ovary.ML.SNV.gz'
ML_UT = '../../data/ML/Urinary-tract.ML.SNV.gz'

# Osorio data
OSORIO_METADATA = '../../data/osorio/metadata_samples.txt'
BOXTEL = '../../data/mutfootprints_2020/hartwig/signatures/boxtel.snvs.tsv'

# BEATAml data
OHSU_DATA_CLINICAL = '../../data/OHSU_2018/data_clinical_sample.txt'
OHSU_DATA_MUTATIONS = '../../data/OHSU_2018/data_mutations_extended.txt'

# IntoGen data
AML_INTOGEN = '../../data/intogen/AML_intogen.tsv'
