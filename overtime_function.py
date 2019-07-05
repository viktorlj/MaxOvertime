import pandas as pd
import numpy as np
import io
import os
import seaborn as sns
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
import base64

#Function to read and process deleterious variants file
def read_var_file(path):
    deleterious_variants = pd.read_csv(path, sep='\t')
    del_var = deleterious_variants.iloc[: ,[0,1,3,4,5,6]]
    
    proper_transcript = []
    variant_type = []
    reference_variant = []
    alternate_variant = []

    for i in del_var['Transcript Variant']:
        for j in i.replace(" ", "").split(";"):        
            #Pick first entry with c and exit
            if j[0] == "c":
                proper_transcript.append(j)
                if len(j.split(">")) > 1:
                    variant_type.append("SNV")
                    reference_variant.append(j.split(">")[0][-1])
                    alternate_variant.append(j.split(">")[1])
                elif len(j.split("del")) > 1:
                    variant_type.append("del")
                    reference_variant.append(j.split("del")[1])
                    alternate_variant.append("-")
                elif len(j.split("dup")) > 1:  
                    variant_type.append("ins")
                    reference_variant.append("-")
                    alternate_variant.append(j.split("dup")[1]) 
                elif len(j.split("ins")) > 1:
                    variant_type.append("ins")
                    reference_variant.append("-")
                    alternate_variant.append(j.split("ins")[1])
                else:
                    variant_type.append("-")
                    reference_variant.append("-")
                    alternate_variant.append("-")
                break
    
    #Add lists as new columns
    del_var = del_var.assign(Proper_Transcript=proper_transcript)
    del_var = del_var.assign(Variant_Type=variant_type)
    del_var = del_var.assign(Reference_Variant=reference_variant)
    del_var = del_var.assign(Alternate_Variant=alternate_variant)

    #Deletions need to have position -1
    del_var.loc[del_var['Variant_Type'] == 'del', 'Position'] -= 1
    return del_var
    
#Function to read vcf file    
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

#Function to read gVCDF and select variants
def select_genomic_variants(path, concat_del_vars, date):
    genomic_vcf = read_vcf(path)
    genomic_vcf['CHROM'] = genomic_vcf['CHROM'].str.replace('chr', '')

    #Select variants from genomic VCF file
    selected_variants = pd.merge(concat_del_vars, genomic_vcf,  how='left', left_on=['Chromosome','Position'], right_on = ['CHROM','POS'])

    #Select wanted columns
    selected_variants = selected_variants.iloc[: ,[0,1,2,3,4,5,6,7,8,-3,-1]]

    #Rename sample unique column to generic name, add Unique_Symbol column, split info field
    selected_variants = selected_variants.rename(columns={ selected_variants.columns[-1]: "Sample_Info" })
    selected_variants['Unique_Symbol'] = selected_variants['Gene Symbol']+"_"+selected_variants['Proper_Transcript']
    selected_variants['VAF'] = selected_variants['Sample_Info'].str.split(":", expand = True)[3].astype('float64').round(4)
    selected_variants['INFO'] = selected_variants['INFO'].str.replace('DP=', '')
    selected_variants['Depth'] = selected_variants['INFO'].str.split(";", expand = True)[0]
    selected_variants['Date'] = date
    
    #Sort by VAF in order to remove duplicates and only keep the variant with largets VAF
    selected_variants = selected_variants.sort_values(by='VAF').drop_duplicates(subset=['Chromosome', 'Position'], keep = 'last')

    return selected_variants