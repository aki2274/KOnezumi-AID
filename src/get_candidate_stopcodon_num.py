import re
from typing import Tuple
import numpy as np
import pandas as pd

def set_presetnum(transcripts_name:str,gene_df,gene_seq_data:dict):
    data = gene_df[gene_df['name']==transcripts_name].reset_index()
    chrom = str(data.loc[0,'chrom'])
    txStart = str(data.loc[0,'txStart'])
    txEnd = str(data.loc[0,'txEnd'])
    query = f'{transcripts_name}::{chrom}:{txStart}-{txEnd}'
    orf_seq = gene_seq_data[query]
    cdsStart =int(data['cdsStart'])
    start = data['exonStarts']
    end = data['exonEnds']
    exon_start_list = list(start.split(','))
    exon_end_list = list(end.split(','))
    return orf_seq,data,txStart,cdsStart,exon_start_list,exon_end_list

def get_exon_seq(transcripts_name,gene_df,gene_seq_data):
    orf_seq,data,txStart,cdsStart,exon_start_list,exon_end_list = set_presetnum(transcripts_name,gene_df,gene_seq_data)
    exon_seq_list =[] 
    for s in range(len(exon_start_list)):
        start_num = int(exon_start_list[s])- txStart
        end_num= int(exon_end_list[s]) - txStart
        exon_seq =orf_seq[start_num:end_num]
        exon_seq_list.append(exon_seq)
    return ''.join(exon_seq_list)

