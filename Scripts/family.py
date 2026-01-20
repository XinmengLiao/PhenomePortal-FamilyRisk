#%%Cell 1 Define the path for parameters and databases
# ==============================================================================
import sys
import pandas as pd
import numpy as np
import os
import gzip
import paramiko
import re
from scp import SCPClient
import time
from itertools import islice
from datetime import datetime
import pickle
from collections import defaultdict
from pathlib import Path

all_start_time = time.time()
# ====== 1) setting up paths and running parameters ======
path = 'sysmed'  # local or sysmed

# User input filtration parameters

# Default setting 
# user_af_clinvar = 1
# user_af_predict = 0.05
# user_ada_score = 0.6
# user_rf_score = 0.6
# user_revel_score=0.75
# user_spliceai_al=0.5
# user_spliceai_dg=0.5
# user_spliceai_dl=0.5
# user_spliceai_ag=0.5
# user_bayesdel_addaf_score=0.0692655
# user_bayesdel_noaf_score=-0.0570105
# user_am_classification=['likely_pathogenic','ambigous']
# user_am_pathogenicity=0.564
# user_clinvar= ["Pathogenic","Likely_pathogenic","Uncertain_significance","Conflicting_classifications_of_pathogenicity"]
# user_acmg_classification = ["Pathogenic","Likely_pathogenic","Uncertain_significance","Benign","Likely_benign"]

# config file of setting paths
config = {
    "local": {
        "fileName": sys.argv[1],
        "output_file": sys.argv[2],
        "ped_file": sys.argv[3],
        'genedb_file': sys.argv[19],
        'customized_genedb_file': sys.argv[20],
        'screening_list': '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/Preset_screening_list.txt'
    },
    "sysmed": {
        "fileName": sys.argv[1],
        "output_file": sys.argv[2],
        "ped_file": sys.argv[3],
        'genedb_file': sys.argv[19],
        'customized_genedb_file': sys.argv[20],
        'screening_list': '/mnt/nas/Genomics/Genome/FamilyRisk/Datasets/Preset_screening_list.txt'
    }
}

user_af_clinvar = sys.argv[4]
user_af_clinvar = float(user_af_clinvar)
user_af_predict = sys.argv[5]
user_af_predict = float(user_af_predict)
user_ada_score = sys.argv[6]
user_ada_score = float(user_ada_score)
user_rf_score = sys.argv[7]
user_rf_score = float(user_rf_score)
user_revel_score = sys.argv[8]
user_revel_score = float(user_revel_score)
user_spliceai_al = sys.argv[9]
user_spliceai_al = float(user_spliceai_al)
user_spliceai_dg = sys.argv[10]
user_spliceai_dg = float(user_spliceai_dg)
user_spliceai_dl = sys.argv[11]
user_spliceai_dl = float(user_spliceai_dl)
user_spliceai_ag = sys.argv[12]
user_spliceai_ag = float(user_spliceai_ag)
user_bayesdel_addaf_score = sys.argv[13]
user_bayesdel_addaf_score = float(user_bayesdel_addaf_score)
user_bayesdel_noaf_score = sys.argv[14]
user_bayesdel_noaf_score = float(user_bayesdel_noaf_score)
user_am_classification = sys.argv[15]
user_am_pathogenicity = sys.argv[16]
user_am_pathogenicity = float(user_am_pathogenicity)
user_clinvar= sys.argv[17].split(",")
user_acmg_classification = sys.argv[18].split(",")

if path not in config:
    raise ValueError(f"unknown path: {path}")

cfg = config[path]

# ====== 3) Reading necessary files ======

def read_db_file(filepath, encoding="ISO-8859-1", sep="\t", fillna_str="No info", drop_allna_cols=False):
    """help to formally read files"""
    df = pd.read_csv(filepath, sep=sep, encoding=encoding)
    if drop_allna_cols:
        df = df.dropna(axis=1, how='all')
    df = df.replace(np.nan, fillna_str)
    return df

# Read the ped file
ped_df = pd.read_csv(cfg['ped_file'], sep="\t", header=None)
ped_df.columns = ['FamilyID', 'IndividualID', 'FatherID', 'MotherID', 'Sex', 'Phenotype']
fam_df = ped_df[(ped_df['FatherID'] != '0') & (ped_df['MotherID'] != '0')]

# GeneDB list
if cfg["genedb_file"] != "" and cfg["customized_genedb_file"] == "":
    print("Useing predefined genedb option.")
    project_list = cfg["genedb_file"].split(',')
    genedb = pd.read_csv(cfg["screening_list"], sep="\t")
    genedb = genedb[genedb['Project'].isin(project_list)]
elif cfg["customized_genedb_file"] != "" and cfg["genedb_file"] == "":
    print("Using customized genedb file.")
    genedb = pd.read_csv(cfg["customized_genedb_file"], sep="\t")
elif cfg["genedb_file"] == "" and cfg["customized_genedb_file"] == "":
    print("Using default Newborn Screening List.")
    genedb = pd.read_csv(cfg["screening_list"], sep="\t")
    genedb = genedb[genedb['Project'].isin(['NBScreening'])]
else:
    raise ValueError("Please provide either a predefined genedb option or a customized genedb file.")

#%%Cell 2 start the analysis program for the vep annotated vcf files, generated report A B C and D
# ==============================================================================
start_time = time.time()
# col_map is used to build a mapping once colNames_CSQ is obtained
col_map = {}  # Will be populated after parsing the "ID=CSQ" line
headings = []

## For REVEL scores refinement from 0.460&0.460&.&. to 0.460
def extract_decimal_from_string(s):
    if not s or not isinstance(s, str):
        return None
    matches = re.findall(r"\d+\.\d+", s)
    if not matches:
        return None
    return matches[0]


# Reading vcf.gz file 
file = gzip.open(cfg["fileName"],'rt')
tLine = file.readline()
i = 0
reportA,reportgwas, reporteqtl = [], [], []

while tLine:
    # remove the newline character
    tLine = tLine.rstrip('\n')
    # split the current line
    iContent = tLine.split('\t')
    i += 1
    ##get the content from VCF annotation header
    if tLine.startswith('#'):
        if 'ID=CSQ' in tLine:
            annoText = iContent[0].split('Format: ')
            colNames_CSQ = annoText[1].replace('">','')
            colNames_CSQ = colNames_CSQ.split('|')
            # construct col_map for all use
            col_map = { name: idx for idx, name in enumerate(colNames_CSQ) }
        elif tLine.startswith('#CHROM'):
            headings = iContent
            print(headings)
        # directly goes into next line
        tLine = file.readline()
        #print(tLine)
        continue
    if not headings:
        tLine = file.readline()
        #print(tLine)
        continue
    
    iText = [s for s in iContent[headings.index('INFO')].split(';') if 'CSQ=' in s]
    iText = iText[0].replace('CSQ=','').split(',')
    
    saveFlag1, saveFlagPGx = False, False
    
    for j in range(0,len(iText)):
        jText = iText[j].split('|')
        # fixing REVEL score 
        if 'REVEL' in col_map and 'REVEL_score' in col_map:
            revel_idx = col_map['REVEL']
            revel_score_idx = col_map['REVEL_score']
            if jText[revel_idx] == '':
                parsed_value = extract_decimal_from_string(jText[revel_score_idx])
                if parsed_value is not None:
                    jText[revel_idx] = parsed_value

        ichr = iContent[headings.index('#CHROM')].split('_')[0]
        ipos = iContent[headings.index('POS')]
        iref = iContent[headings.index('REF')]
        ialt = iContent[headings.index('ALT')]
        ivariation3 = f"{ipos}_{iref}_{ialt}"
        ivariation4 = f"{ichr}_{ipos}_{iref}_{ialt}"
                        
        # 1) ClinVar and prediction filtered based on users' needs 
        if 'MAX_AF' in col_map:
            max_af_val = jText[col_map['MAX_AF']]
            try:
                max_af_numeric = float(max_af_val) if max_af_val != '' else 0.0
            except (ValueError, TypeError):
                max_af_numeric = 0.0
            jText[col_map['MAX_AF']] = max_af_numeric

        if jText[col_map['ClinVar_CLNSIG']] != "" and max_af_numeric < user_af_clinvar:
            saveFlag1 = "Keep"
        
        elif jText[col_map['ClinVar_CLNSIG']] == "" and max_af_numeric < user_af_predict:   
            predicted_impact = (
                    ('IMPACT' in col_map and jText[col_map['IMPACT']] == 'HIGH')
                    or ('ada_score' in col_map and jText[col_map['ada_score']] != '' and float(jText[col_map['ada_score']]) > user_ada_score)
                    or ('rf_score' in col_map and jText[col_map['rf_score']] != '' and float(jText[col_map['rf_score']]) > user_rf_score)
                    or ('REVEL' in col_map and jText[col_map['REVEL']] != '' and float(jText[col_map['REVEL']]) > user_revel_score)
                    or ('SpliceAI_pred_DS_AL' in col_map and jText[col_map['SpliceAI_pred_DS_AL']] != '' and float(jText[col_map['SpliceAI_pred_DS_AL']])>user_spliceai_al)
                    or ('SpliceAI_pred_DS_DG' in col_map and jText[col_map['SpliceAI_pred_DS_DG']] != '' and float(jText[col_map['SpliceAI_pred_DS_DG']])>user_spliceai_dg)
                    or ('SpliceAI_pred_DS_DL' in col_map and jText[col_map['SpliceAI_pred_DS_DL']] != '' and float(jText[col_map['SpliceAI_pred_DS_DL']])>user_spliceai_dl)
                    or ('SpliceAI_pred_DS_AG' in col_map and jText[col_map['SpliceAI_pred_DS_AG']] != '' and float(jText[col_map['SpliceAI_pred_DS_AG']])>user_spliceai_ag)
                    or ('BayesDel_addAF_score' in col_map and jText[col_map['BayesDel_addAF_score']] != '' and float(jText[col_map['BayesDel_addAF_score']])>user_bayesdel_addaf_score)
                    or ('BayesDel_noAF_score' in col_map and jText[col_map['BayesDel_noAF_score']] != '' and float(jText[col_map['BayesDel_noAF_score']])>user_bayesdel_noaf_score)
                    or ('am_class' in col_map and jText[col_map['am_class']] in user_am_classification
                        and 'am_pathogenicity' in col_map
                        and jText[col_map['am_pathogenicity']] != ''
                        and float(jText[col_map['am_pathogenicity']])>user_am_pathogenicity))
            if predicted_impact:
                saveFlag1 = "Keep"
        
    # after for j in range(len(iText)) loop, if saveFlag1/2/3/4/5 has value, append the line to respective report
    if saveFlag1:
        reportA.append(tLine)


    # print progress every 1000000 lines
    if i % 1000000 == 0:
        print(f"{i} lines processed!")

    # read the next line
    tLine = file.readline()


file.close()
total_row = i
print('Cell 2 VEP annotated File processing done! Now start to map GeneDB and DiseaseDB')
end_time = time.time()
print("Total processing time: {:.2f} seconds".format(end_time - start_time))

# Manage text file into a dataframe
base_vcf_columns = headings
#base_vcf_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','Sample_Info']

csq_columns = [''] * len(col_map)  
for col_name, col_index in col_map.items():
    csq_columns[col_index] = col_name

all_output_columns = base_vcf_columns + csq_columns

def process_vcf_line_to_expanded_format(vcf_line, col_map, csq_columns):
    expanded_rows = []

    fields = vcf_line.split('\t')
    info_field = fields[7]  # INFO
    csq_info = [s for s in info_field.split(';') if 'CSQ=' in s]
    
    if not csq_info:
        return expanded_rows
    
    csq_data = csq_info[0].replace('CSQ=', '').split(',')
    
    for transcript_annotation in csq_data:
        transcript_fields = transcript_annotation.split('|')
        
        while len(transcript_fields) < len(csq_columns):
            transcript_fields.append('')
        
        expanded_row = fields + transcript_fields[:len(csq_columns)]
        expanded_rows.append(expanded_row)
    
    return expanded_rows

expanded_reportA = []
for i, vcf_line in enumerate(reportA):
    expanded_rows = process_vcf_line_to_expanded_format(vcf_line, col_map, csq_columns)
    expanded_reportA.extend(expanded_rows)
if expanded_reportA:
    expanded_reportA = pd.DataFrame(expanded_reportA, columns=all_output_columns)

# Keeps the rows with Users defined prediction and ClinVar criteria
def contains_any_or_empty(x, keywords):
    if str(x).strip() == "":
        return True
    return any(kw in str(x) for kw in keywords)

expanded_reportA = expanded_reportA[expanded_reportA["ClinVar_CLNSIG"].apply(lambda x: contains_any_or_empty(x, user_clinvar))]
expanded_reportA = expanded_reportA[~expanded_reportA["#CHROM"].astype(str).str.contains("_")].copy()

# match with genedb 
expanded_reportA = pd.merge(expanded_reportA, genedb, how="left", left_on="SYMBOL", right_on="Genes")
expanded_reportA = expanded_reportA[
    (expanded_reportA['Disease'].notna()) & 
    (expanded_reportA['Disease'] != '') & 
    (expanded_reportA['Disease'] != ' ')
]

#%% Cell 3 GeneBe ACMG classification ======================================
import genebe as gnb

small_df = expanded_reportA.loc[:,["#CHROM","POS","REF","ALT"]]
small_df = small_df.rename(columns={"#CHROM":"chr","POS":"pos","REF":"ref","ALT":"alt"})
unique_small_df = small_df.drop_duplicates()

try:
    annotated_df = gnb.annotate(
        unique_small_df,
        genome='hg38',
        use_ensembl=False,
        use_refseq=True,
        flatten_consequences=True,
        output_format="dataframe"
    )
except Exception as e:
    print(f"GeneBe annotation failed: {str(e)}")
    # create an empty DataFrame to keep the structure
    annotated_df = pd.DataFrame(columns=[
        'chr', 'pos', 'ref', 'alt', 'gene_symbol', 
        'acmg_score', 'acmg_classification', 'acmg_criteria'
    ])
annotated_df = annotated_df.rename(columns={"chr":"#CHROM","pos":"POS","ref":"REF","alt":"ALT"})

required_columns = ["#CHROM","POS","REF","ALT","gene_symbol",
                   "acmg_score","acmg_classification",'acmg_criteria']
for col in required_columns:
    if col not in annotated_df.columns:
        annotated_df[col] = None  # add missing columns

small_annotate_all = annotated_df[required_columns].rename(columns={'gene_symbol': 'SYMBOL'})
# merge automatically handles missing values
expanded_reportA = pd.merge(
    expanded_reportA, 
    small_annotate_all, 
    how="left", 
    on=["#CHROM","POS","REF","ALT","SYMBOL"]
)


expanded_reportA['rsID'] = expanded_reportA['Existing_variation'].str.extract(r'(rs\d+)', expand=False)
expanded_reportA['variant_info'] = expanded_reportA[['#CHROM', 'POS', 'REF', 'ALT']].fillna('').astype(str).agg('_'.join, axis=1)
expanded_reportA = expanded_reportA[expanded_reportA['acmg_classification'].isin(user_acmg_classification)]
expanded_reportA = expanded_reportA.drop_duplicates()


#%% Cell 4 Parse kid genotype ======================================
kidID = []
kidSex = []
for i, row in fam_df.iterrows():
        fatherID = row["FatherID"]
        motherID = row["MotherID"]
        kidID.append(row["IndividualID"])
        kidSex.append(str(row["Sex"]))

def extract_gt(val):
    if isinstance(val, str):
        return val.split(":")[0]
    return None

expanded_reportA["FatherGenotype"] = expanded_reportA[fatherID].apply(extract_gt)
expanded_reportA["MotherGenotype"] = expanded_reportA[motherID].apply(extract_gt)

for i in range(len(kidID)):
    kid = kidID[i]
    sex = str(kidSex[i])
    if kid in expanded_reportA.columns:
        varID = f'Kid{i}Genotype_{kid}'
        varSex = f'Kid{i}Gender_{kid}'
        expanded_reportA[varID] = expanded_reportA[kid].apply(extract_gt)
        expanded_reportA[varSex] = ["Male" if sex == "1" else "Female"] * len(expanded_reportA)
                
kid_genotype_cols = [col for col in expanded_reportA.columns if col.startswith("Kid") and "Genotype" in col]
columns_to_replace = ["FatherGenotype", "MotherGenotype"] + kid_genotype_cols

def replace_missing_gt(gt):
    if isinstance(gt, str):
        gt = gt.replace("./.", "0/0").replace(".|.", "0|0")
        gt = re.sub(r"\.", "0", gt)
    return gt

for col in columns_to_replace:
    expanded_reportA[col] = expanded_reportA[col].apply(replace_missing_gt)

# expanded_reportA['keep'] = False
# for col in kid_genotype_cols:
#     expanded_reportA['keep'] |= (
#         (expanded_reportA[col] != "./.") &
#         (expanded_reportA["FatherGenotype"] != "./.") &
#         (expanded_reportA["MotherGenotype"] != "./.")
#     )

# expanded_reportA = expanded_reportA[expanded_reportA['keep']].copy()
# expanded_reportA.drop(columns=['keep'], inplace=True)

#%% Cell 5 Parse kid pattern ======================================

for kid_gt_col in kid_genotype_cols:
    pattern_col = kid_gt_col.replace("Genotype_","Pattern_")

    def classify_variant(row):
        kid = row[kid_gt_col]
        father = row["FatherGenotype"]
        mother = row["MotherGenotype"]
        inh = row['Inheritance']
        chrom = row['#CHROM']

        if any(gt == "./." for gt in [kid, father, mother]):
            return "NA"

        # De novo: 
        if father in ["0/0", "0|0"] and mother in ["0/0", "0|0"] and kid in ["0/1", "0|1" ,"1/1", "1|1", "1/0", "1|0"]:
            return "de novo"

        # Recessive: 
        elif father in ["0/1", "0|1", "1/0", "1|0"] and mother in ["0/0","0|0"] and kid in ["0/1", "0|1", "1/0", "1|0"] and inh in ["AR","XLR"]:
            return "recessive"
        elif mother in ["0/1", "0|1", "1/0", "1|0"] and father in ["0/0","0|0"] and kid in ["0/1", "0|1", "1/0", "1|0"] and inh in ["AR","XLR"]:
            return "recessive"
        elif father in ["1"] and mother in ["0"] and kid in ["1"] and chrom == "chrX" and inh in ["AR","XLR"]:
            return "recessive"
        elif mother in ["1"] and father in ["0"] and kid in ["1"] and chrom == "chrX" and inh in ["AR","XLR"]:
            return "recessive"
        
        # Dominant
        elif kid in ["1/1", "1|1", "0/1", "0|1", "1/0", "1|0","1"] and inh in ["AD","Multiple and/or complex pattern","AR; AD","AD; AR","X-linked multiple and/or complex pattern","XLD"]:
            return "AD/XLD dominant"
        elif kid in ["1/1", "1|1","1"] and inh in ["AR","XLR"] :
            return "AR/XLR dominant"
        else:
            return "NA"

    expanded_reportA[pattern_col] = expanded_reportA.apply(classify_variant, axis=1)

#%% Cell 6 Output python managed file and extract unique genes ======================================
print("\nPython output Statistic")
print(f"Original rows: {total_row}")
print(f"Filtered and output rows: {len(expanded_reportA)}")

expanded_reportA.to_csv(cfg["output_file"], sep="\t", index=False, quoting=3)