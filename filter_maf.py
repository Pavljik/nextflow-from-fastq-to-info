import sys
import os
import contextlib

import pandas as pd
import numpy as np
from myvariant import MyVariantInfo
from tqdm import tqdm

sample = sys.argv[1]

var_info = MyVariantInfo()

cols = [
    'Hugo_Symbol', 'Chromosome',
    'Start_Position', 'End_Position', 
    'EXON', 'Strand', 'Variant_Classification', 'Feature_type',
    'Variant_Type', 'HGVSc', 'HGVSp', 
    'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
    'dbSNP_RS', 'gnomADe_AF', 'CLIN_SIG', 'PolyPhen', 'SIFT', 'PUBMED'
]

def get_pathogenic_rcvs(rcvs, rs: str) -> pd.Series:
    
    if isinstance(rcvs, dict):
        rcvs = [rcvs]

    df = pd.DataFrame()
    for rcv in rcvs:
        try:
            if 'pathogen' in rcv.get('clinical_significance', '').lower():
                df = pd.concat(
                    [
                        df,
                        pd.Series(
                            {
                                'RCV': rcv.get('accession', 'NA'), 
                                'Clinical_Significance': rcv.get('clinical_significance', 'NA'), 
                                'Disease': rcv.get('conditions', {}).get('name', 'NA'),
                            },
                            name = rs
                        )
                    ], axis=1
                )
        except:
            print(rcv)
            continue
            
    return df.T.reset_index().drop('index', axis=1)



table1 = pd.read_csv(f'./result/{sample}.maf', sep='\t', comment='#', 
                        header=0, index_col=None)
table1 = table1.dropna(subset=['CLIN_SIG']).query("CLIN_SIG != 'benign'")[cols]

table2 = pd.DataFrame()

for _, row in tqdm(table1.sort_values('Hugo_Symbol').iterrows(), desc='Mutations'):
    
    with contextlib.redirect_stderr(open(os.devnull, 'w')):
        rs = var_info.getvariant(row.dbSNP_RS, fields='clinvar')

    rs_list = rs if isinstance(rs, list) else [rs] if isinstance(rs, dict) else []

    for record in rs_list:
        clinvar_data = record.get('clinvar')
        if clinvar_data:
            df = get_pathogenic_rcvs(clinvar_data.get('rcv', []), row.dbSNP_RS)
            df = pd.concat(
                [
                    pd.DataFrame(
                        {
                            'Hugo_Symbol': row.Hugo_Symbol,
                            'Variant_Classification': row.Variant_Classification, 
                            'HGVSp': row.HGVSp,
                            'Allele1': row.Tumor_Seq_Allele1,
                            'Allele2': row.Tumor_Seq_Allele2, 
                            'dbSNP_RS': row.dbSNP_RS, 
                            'PolyPhen': row.PolyPhen, 
                            'SIFT': row.SIFT,
                            'gnomADe_AF': row.gnomADe_AF,
                        }, index=df.index.values),
                    df
                ], axis=1
            )
            
            table2 = pd.concat(
                [
                    table2,
                    df
                ], axis=0
            )

table2 = table2.reset_index().drop('index', axis=1)

table2.to_csv(f"./result/{sample}-filtered.maf", sep='\t')
