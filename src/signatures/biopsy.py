# dictionary with biopsies


def return_biopsies():
    biopsy_dict = {
        'Abdomen right': 'Abdomen',
        'Bijnier Rechts': 'Adrenal',
        'Chest wall': 'Peritoneum',
        'Cutaneous lesions': 'Skin',
        'Fibrosarcoma': 'Sarcoma',
        'Gastro-intestinal': 'GI',
        'Glioblastoma': 'Brain',
        'Ileocecal': 'Colon',
        'Mesothelioma left diaphragm': 'Mesothelioma',
        'Psoas': 'Muscle',
        'RIP': 'unknown',
        'Small pelvis': 'unknown',
        'Soft mass': 'Soft_tissue',
        'Soft tissue mass pelvix': 'Soft_tissue',
        'Soft tissue mass sternum': 'Soft_tissue',
        'Subcutane': 'Subcutaneous',
        'Thoracal wall': 'Peritoneum',
        'Vagina': 'Vagina',
        'abdomen left side': 'Peritoneum',
        'abdominal mass': 'unknown',
        'belly button': 'unknown',
        'between liver and kiney': 'unknown',
        'leasion': 'unknown',
        'mass inguinal': 'unknown',
        'mesenterium': 'Peritoneum',
        'next to spleen': 'Unknown',
        'omentale soft tissue': 'Omentum',
        'ossale vertebrae': 'Bone',
        'other': 'unknown',
        'pelvic': 'unknown',
        'retroperitoneal': 'Peritoneum',
        'rib': 'Bone',
        'right peritoneal towards small basin': 'Peritoneum',
        'sacrum': 'Bone',
        'shoulder blade': 'Bone',
        'skin abdomen': 'Skin',
        'skin chest wand': 'Skin',
        'subcutaan lower leg right': 'Subcutaneous',
        'subcutane, right side just above the umbilicus': 'Subcutaneous',
        'vertebral column': 'Bone',
        'vulva': 'Vulva',
        'Liver': 'Liver',
        'Lymph node': 'Lymph',
        'Lung': 'Lung',
        'soft tissue': 'Soft_tissue',
        'abdominal wall': 'Peritoneum',
        'abdomen': 'Peritoneum',
        'Adrenal gland': 'Adrenal',
        'subcutis': 'Subcutaneous',
        'omentum': 'Omentum',
        'mamma': 'Breast',
        'Mamma': 'Breast',
        'Soft tissue': 'Soft_tissue',
        'brain' : 'Brain',
        'Thorax wall': 'Peritoneum',
        'Breast': 'Breast',
        'skin/subcutaneous': 'Skin',
        'nan': 'nan',
        'Psoas left': 'Muscle',
        'neck right side': 'Neck',
        'caudaal below liver': 'Liver',
        'mass ventrale bladder': 'Bladder',
        'Vaginal': 'Vagina',
        'Parasternal': 'Peritoneum',
        'pericardiale laesie': 'Pericardium',
        'Clavicula': 'Bone',
        'just above a. subclavia left side': 'nan',
        'esophagus': 'Esophagous',
        'Groin left': 'nan',
        'subcutaneous nodus': 'Subcutaneous',
        'gland in the neck': 'Neck',
        'Lesion glandula parotis': 'Parotida',
        'Skin/ Subcutaneous': 'Skin',
        'prostate': 'Prostate',
        'Subcutis': 'Subcutaneous',
        'intraperitoneal mass': 'Peritoneum',
        'small intestines': 'SI',
        'Laesie level costa 9': 'nan',
        'cutaneous': 'Skin',
        'chest wall': 'Peritoneum',
        'neck': 'Neck',
        'Central in abdomen': 'Peritoneum',
        'Muscle': 'Muscle',
        'epigastrium': 'Peritoneum',
        'epicardial': 'Peritoneum',
        'pancreas': 'Pancreas',
        'Parevertebral': 'Bone',
        'tumor near clavicula medial left': 'nan',
        'Thoracic wall': 'Peritoneum',
        'cutane': 'Skin',
        'omentum cake': 'Omentum',
        'urinary bladder': 'Bladder',
        'gastric mucosa': 'Stomach',
        'Mandible': 'Bone',
        'right paracolisch right paracolisch': 'Colon',
        'Omentum': 'Omentum',
        'Left mamma': 'Breast',
        'ovary': 'Ovary',
        'Mass abdomen': 'Peritoneum',
        'Mamma right': 'Breast',
        'sinus': 'Bone',
        'Oesophagus': 'Esophagous',
        'Mestastasis': 'nan',
        'brain left temporal GBM': 'Brain',
        'axillary': 'nan',
        'thoracic wall': 'Peritoneum',
        'Small pelvic': 'nan',
        'skin, medial left upper leg': 'Skin',
        'peri-anal': 'nan',
        'iliac': 'Muscle',
        'peritoneale': 'Peritoneum',
        'thorax wall': 'Peritoneum',
        'thorax': 'Peritoneum',
        'presternal': 'Bone',
        'subcutane': 'Subcutaneous',
        'Gastric antrum': 'Stomach',
        'intra costal': 'Bone',
        'parasternaal': 'nan',
        'upperarm': 'nan',
        'Pericardium': 'Pericardium',
        'left mandibula': 'Bone',
        'parafalciene lesion': 'Brain',
        'vaginal': 'Vagina',
        'Uterus': 'Uterus',
        'melanonoma': 'Skin',
        'upper abdomen': 'Peritoneum',
        'mediastinale mass': 'Peritoneum',
        'ovarium': 'Ovary',
        'subcutane laesion': 'Skin',
        'paravertebral': 'Bone',
        'nasal cavity': 'nan',
        'left axilla': 'nan',
        'mediastinum': 'Peritoneum',
        'Intra abdominal mass': 'Peritoneum',
        'SKIN': 'Skin',
        'cardia': 'Stomach',
        'Bottom': 'nan',
        'peritoneal deposition': 'Peritoneum',
        'skin leasion': 'Skin',
        'shoulder': 'Skin',
        'Thyroid': 'Thyroid',
        'vagina': 'Vagina',
        'small intestine': 'SI',
        'Para-vertebral mass': 'nan',
        'Pleura': 'Pleura',
        'Skin/Subcutaneous': 'Skin',
        'major duodenal papilla': 'SI',
        'Omental cake': 'Peritoneum',
        'Methastase': 'nan',
        'supraclaviculair': 'nan',
        'subcutane meta': 'Subcutaneous',
        'left retroperitoneal lesion': 'Peritoneum',
        'sternum': 'Bone',
        'axilla': 'nan',
        'Musculus Pectoralis Major': 'Muscle',
        'sub cutaneous': 'Skin',
        'retroperitoneal mass': 'Peritoneum',
        'Subcutane metastasis.': 'Subcutaneous',
        'axillair': 'nan',
        'abdominal': 'Peritoneum',
        'Amputated upper leg (Leg amputation 16/06/2016)': 'nan',
        'Muscle back': 'Muscle',
        'Neck lesion': 'Neck',
        'Muscle dorsum': 'Muscle',
        'colon': 'Colon',
        'ingrinal': 'nan',
        'Skin/subcutaneous': 'Skin',
        'infraclaviculare': 'nan',
        'duodenum': 'SI',
        'skin, subcutaneous': 'Skin',
        'oesophagus': 'Esophagous',
        'Ovarium': 'Ovary',
        'ventral abdomen': 'Peritoneum',
        'Subcutaneous': 'Subcutaneous',
        'sarcoma': 'Subcutaneous',
        'Left para vertebral mass': 'nan',
        'bladder': 'Bladder',
        'musculus rectus abdominis': 'Muscle',
        'small pelvis': 'nan',
        'skin lesion': 'Skin',
        'inguinal': 'nan',
        'hypo-echogenic mass chest wall': 'Peritoneum',
        'Skin/ subcutaneous': 'Skin',
        'Bone': 'Bone',
        'Neck': 'Neck',
        'Abdomen muscle': 'Muscle',
        'gall bladder': 'Biliary',
        'squamous cell carcinoma': 'Skin',
        'Subcutaneus leasion': 'Subcutaneous',
        'right side': 'nan',
        'Intra-abdominal mass': 'Peritoneum',
        'breast/mamma': 'Breast',
        'sigmoid': 'Colon',
        'right shoulder': 'nan',
        'Cutane lesion': 'Skin',
        'left mamma': 'Breast',
        'groin': 'nan',
        'ABDOMINAL WALL': 'Peritoneum',
        'Subcutus': 'Subcutaneous',
        'retroperitoneum': 'Peritoneum',
        'trachea': 'Trachea',
        'stomach': 'Stomach',
        'cervix': 'Cervix',
        'Epigastrium': 'Stomach',
        'bronchus': 'Lung',
        'Left scapula': 'Bone',
        'Abdominal wall': 'Peritoneum',
        'Vulva': 'Vulva',
        'Pleural fluid': 'Pleura',
        'Brain': 'Brain',
        'Kidney': 'Kidney',
        'breast left': 'Breast',
        'soft tissue mass sacrum': 'Soft_tissue',
        'armpit': 'nan',
        'Cutaneous': 'Skin',
        'Musculus Erector Truncus Spinae': 'Muscle',
        'Sarcoma': 'Subcutaneous',
        'Urinary tract': 'Urinary',
        'Penis': 'Penis',
        'Back': 'nan',
        'Left maximus gluteus': 'Subcutaneous',
        'unknown': 'Unknown',
        'Soft tissue mass': 'Soft_tissue',
        'tongue': 'Muscle',
        'Below left nipple': 'Breast',
        'Rectum': 'Colon',
        'subcutane: right side just above the umbilicus': 'Subcutaneous',
        'lesion subcutaneous': 'Subcutaneous',
        'm. rectus abdominus': 'Muscle',
        'Abdominal mass': 'Peritoneum',
        'cutis': 'Skin',
        'parasternal': 'Peritoneum',
        'Nates': 'Brain',
        'peritonial lesion': 'Peritoneum',
        'muscle': 'Muscle',
        'subcatenous': 'Subcutaneous',
        'Woundbed necrotic cavity perineal': 'Peritoneum',
        'neck right': 'Neck',
        'skinlesion': 'Skin',
        'esophagus - small intestine': 'Esophagous',
        'Paravertebral': 'Paravertebral',
        'right breast': 'Breast',
        'Endonasal': 'nan',
        'paravertebrale massa': 'Paravertebral',
        'Skin/ cutaneous': 'Skin',
        'pleural': 'Pleura',
        'Abdominal wall cutaneous': 'Peritoneum',
        'CNS': 'Brain',
        'Peritoneum': 'Peritoneum',
        'Head': 'Neck',
        'muscle tissue': 'Muscle',
        'lacrimal sac': 'nan',
        'pleural lesion': 'Pleura',
        'mamma right parasternal': 'Breast',
        'adnex': 'Ovary',
        'retro-peritoneal': 'Peritoneum',
        'Occiput': 'Bone',
        'psoas': 'Muscle',
        'subcutaneous mass': 'Peritoneum',
        'swelling back': 'nan',
        'mass anterieur': 'Peritoneum',
        'Scapula': 'Bone',
        'glandular package in aorta': 'nan',
        'Ascites': 'nan',
        'Salivary gland': 'Salivary_Gland',
        'Sternum': 'Bone',
        'Abdominal': 'Peritoneum',
        'neck lymph node': 'Lymph',
        'right axilla': 'nan',
        'Adnex': 'Ovary',
        'Primary': 'Primary',
        'Mouth': 'Mouth',
        'Abdominal metastatic': 'Peritoneum',
        'soft tissue mass near sternum': 'Soft_tissue',
        'subcutaneous metastasis': 'Subcutaneous',
        'In lower pelvis': 'nan',
        'pelvis': 'nan',
        'Sacrum': 'Bone',
        'bone': 'Bone',
        'Thorax': 'Peritoneum',
        'Leg': 'nan',
        'Palatum Molle': 'Mouth',
        'Skin': 'Skin',
        'Unknown': 'Unknown',
        'Intra abdominal': 'Peritoneum',
        'Stomach': 'Stomach',
        'subcutane left abdominal': 'Subcutaneous',
        'rectum': 'Colon',
        'pleura': 'Pleura',
        'mass axillary': 'nan',
        'Ulcus neck': 'SI',
        'ilium': 'SI',
        'mass Brain': 'Brain',
        'peritoneum': 'Peritoneum',
        'subcateneous': 'Subcutaneous',
        'lower abdomen': 'Peritoneum',
        'intestins': 'SI',
        'Neck gland': 'Neck',
        'Para iliacal mass': 'Peritoneum',
        'Mass head': 'Head',
        'Cutane metastasis abdomen': 'nan',
        'Metastasis Pancreascarcinoma': 'Pancreas',
        'breast': 'Breast',
        'Axilla': 'nan',
        'Abdomen': 'Peritoneum',
        'subcutaneous': 'Subcutaneous',
        'para-aortal abdominal': 'Peritoneum',
        'esophagus/stomach': 'Esophagous',
        'skin': 'Skin',
        'cutan': 'Skin',
        'subcutane right gluteal': 'Subcutaneous',
        'Skin subcutaneous': 'Subcutaneous',
        'NaN': 'nan',
        '(Sub)cutaneous': 'Subcutaneous',
        'Abdominal wall paramedian': 'nan',
        'Axillair right gland': 'nan',
        'Bladder': 'Bladder',
        'Intra peritoneal': 'Peritoneum',
        'Intramuscular': 'Muscle',
        'Lesion sacral': 'nan',
        'Mass': 'nan',
        'Mass right mamma': 'Breast',
        'Metastasis': 'nan',
        'Ovarium CA': 'Ovary',
        'Ovary': 'Ovary',
        'Pectoralis major': 'Muscle',
        'Pelvic metastasis': 'nan',
        'Perineal': 'nan',
        'Skull': 'Bone',
        'Soft Tissue': 'Soft_tissue',
        'Spleen/Pancreas': 'Pancreas',
        'Subcutaneous lesion on the back': 'Subcutaneous',
        'belly': 'Stomach',
        'cavum nasi': 'nan',
        'costa metastase': 'nan',
        'doorgroei naar subcutaan linker flank': 'Subcutaneous',
        'dorsal from vena cava inferior/medial right kidney': 'nan',
        'glioma': 'Brain',
        'gluteus maximus': 'Muscle',
        'hemipelvis': 'nan',
        'interpulmonal': 'Lung',
        'intestine': 'Colon',
        'intra-abdominale mass': 'Peritoneum',
        'intramuscular': 'Muscle',
        'intraperitoneale deposition': 'Peritoneum',
        'large mass Starting from the uterus.': 'Uterus',
        'leasion abdomen': 'nan',
        'lesion median in the lower abdomen.': 'nan',
        'mamma right': 'Breast',
        'mass mesorectum': 'Colon',
        'mesenteriale nodus': 'nan',
        'mesentery': 'nan',
        'metastasis right leg': 'nan',
        'musculus': 'Muscle',
        'occipital': 'Bone',
        'omental cake': 'Omentum',
        'pleural mass': 'Pleura',
        'scalp': 'nan',
        'scin': 'nan',
        'skull skin': 'nan',
        'soft tissue inside os ilium': 'Soft_tissue',
        'subcutaneously': 'Subcutaneous',
        'subcutanous': 'Subcutaneous',
        'cerebral metastasis':'Brain', 
        'abdominal sampling/debulking':'Peritoneum', 
        'unknown':'nan', 
        'Femur':'Bone', 
         'os ilium':'Bone', 
        'costa':'nan', 
        'Unknown':'nan', 
        'mammae':'Breast', 
        'scapula':'Bone', 
         'Prostaat':'Prostate', 
         'lymfadenopathy left':'Lymph', 
        'Mammae':'Breast', 
    }
    return biopsy_dict
