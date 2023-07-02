#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import math
import os

import warnings
warnings.filterwarnings("ignore")


def cad_genotype_code(x1, x2, x3, y):
    if y == x1:
        return 2  #2
    elif y == x2:
        return 1  #1
    elif y == x3:
        return 0  #0
    else :
        return np.NaN


xls = pd.ExcelFile('CAD PRS_Calci.xlsx')
df = pd.read_excel(xls, 'Sheet3')

PRS_samples_df = pd.read_csv("PRS _ demo data for four samples.csv")

df = df.dropna()
col_list = ['Band', 'rs-number', 'Gene(s)', 'Risk allele', 'Risk allele frequency', 'OR ', '95% CI)', 'V3md calls',
            'Genotype Call Dose 2', 'Genotype call  Dose 1', 'Genotype call Dose 0']


template_df = df[col_list]

for sample_id in list(PRS_samples_df.columns):
    print(sample_id)
    temp_df = pd.concat([template_df, PRS_samples_df[sample_id]], axis=1).dropna()
    temp_df["CAD Genotype code_cal"] = temp_df[['Genotype Call Dose 2', 'Genotype call  Dose 1', 'Genotype call Dose 0', sample_id]].apply(lambda x: cad_genotype_code(x["Genotype Call Dose 2"], x["Genotype call  Dose 1"], x["Genotype call Dose 0"], x[sample_id]), axis=1)
    temp_df["Beta_cal"]              = temp_df["OR "].apply(lambda x: math.log(x))
    temp_df["Population score_cal"]  = temp_df[["Beta_cal", 'Risk allele frequency']].apply(lambda x: (x['Beta_cal'] * x['Risk allele frequency']), axis = 1)
    temp_df["Zero center score_cal"] = temp_df[["Beta_cal", "CAD Genotype code_cal", 'Population score_cal']].apply(lambda x: ((x["Beta_cal"] * x["CAD Genotype code_cal"])-x['Population score_cal']), axis = 1)
    temp_df["z score_cal"] = temp_df["Zero center score_cal"]/(temp_df["Population score_cal"].std())
    temp_df["Population score_std"] = temp_df["Population score_cal"].std()
    temp_df["z score avg"] = temp_df["z score_cal"].mean()
    
    path = "PRS_Outputs"
    if not os.path.exists(path):
        os.makedirs(path)
    sample_file_name = path + "/" + "PRS_" + sample_id + ".csv"
    temp_df.to_csv(sample_file_name, index= False)