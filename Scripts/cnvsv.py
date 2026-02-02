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

# config file of setting paths
config = {
    "local": {
        "fileName": sys.argv[1],
        "output_file": sys.argv[2],
        'genedb_file': sys.argv[3],
        'customized_genedb_file': sys.argv[4],
        'compound_het': sys.argv[5],
        'data_type': sys.argv[6],
        'screening_list': '/Users/xinmengliao/Documents/Project/20250710_FamilyRisk/Datasets/genelists/Preset_screening_list20251125.txt'
    },
    "sysmed": {
        "fileName": sys.argv[1],
        "output_file": sys.argv[2],
        'genedb_file': sys.argv[3],
        'customized_genedb_file': sys.argv[4],
        'compound_het': sys.argv[5],
        'data_type': sys.argv[6],
        'screening_list': '/mnt/nas/Genomics/Genome/FamilyRisk/Datasets/Preset_screening_list20251125.txt'
    }
}

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
    genedb = genedb[genedb['Project'].isin(['BabySeq_GroupB', 'BabySeq_GroupA'])]
else:
    raise ValueError("Please provide either a predefined genedb option or a customized genedb file.")

#%%Cell 2 Filter the genes that only included in the screening list
start_time = time.time()
resultdf = pd.read_csv(cfg["fileName"], sep="\t")
resultdf_filter = resultdf[resultdf["Gene_name"].isin(genedb["Genes"])]

if cfg["compound_het"].lower() == "yes":
    if cfg["data_type"].lower() == "cnv":
        output_filename=cfg["output_file"].replace(".cnvsv.txt", ".cnv.with_compound_het.tsv")
    if cfg["data_type"].lower() == "sv":
        output_filename=cfg["output_file"].replace(".cnvsv.txt", ".sv.with_compound_het.tsv")
else:
    if cfg["data_type"].lower() == "cnv":
        output_filename=cfg["output_file"].replace(".cnvsv.txt", ".cnv.tsv")
    if cfg["data_type"].lower() == "sv":
        output_filename=cfg["output_file"].replace(".cnvsv.txt", ".sv.tsv")
    output_filename=cfg["output_file"]
    
resultdf_filter.to_csv(output_filename, sep="\t", index=False, quoting=3)
