import re
from typing import Tuple
import numpy as np
import pandas as pd

def setup(transcripts_name:str,gene_df,gene_seq_data:dict):
    data = gene_df[gene_df['name']==transcripts_name].reset_index()
    chrom = str(data.loc[0,'chrom'])
    txStart = str(data.loc[0,'txStart'])
    txEnd = str(data.loc[0,'txEnd'])
    query = f'{transcripts_name}::{chrom}:{txStart}-{txEnd}'
    orf_seq = gene_seq_data[query]
    cdsStart =int(data['cdsStart'])
    start = data['exonStarts']
    end = data['exonEnds']
    exon_start_list = [int(x) for x in list(start.split(','))]
    exon_end_list = [int(x) for x in list(end.split(','))]
    return orf_seq,data,int(txStart),cdsStart,exon_start_list,exon_end_list

def get_exon_seq(transcripts_name,gene_df,gene_seq_data):
    orf_seq,data,txStart,cdsStart,exon_start_list,exon_end_list = setup(transcripts_name,gene_df,gene_seq_data)
    exon_seq_list =[] 
    for s in range(len(exon_start_list)):
        start_num = exon_start_list[s]- txStart
        end_num= exon_end_list[s] - txStart
        exon_seq =orf_seq[start_num:end_num]
        exon_seq_list.append(exon_seq)
    return ''.join(exon_seq_list)

def get_startcodon_exonindex(transcripts_name,gene_df,gene_seq_data)->int:
    orf_seq,data,txStart,cdsStart,exon_start_list,exon_end_list = setup(transcripts_name,gene_df,gene_seq_data)
    for s in range(len(exon_start_list)):
        if (exon_start_list[s] <= cdsStart<=exon_end_list[s]):
             exon_num =s          
    return exon_num

def get_stopcodon_exonindex(transcripts_name,gene_df,gene_seq_data)->int:
    orf_seq,data,txStart,cdsStart,exon_start_list,exon_end_list = setup(transcripts_name,gene_df,gene_seq_data)
    cdsEnd = int(data['cdsEnd'])
    for s in range(len(exon_start_list)):
        if (exon_start_list[s] <= cdsEnd<= exon_end_list[s]):
             exon_num =s          
    return exon_num

def get_cdsseq(transcripts_name,gene_df,gene_seq_data,exon_seq,startcodon_exonindex,endcodon_exonindex):
    orf_seq,data,txStart,cdsStart,exon_start_list,exon_end_list = setup(transcripts_name,gene_df,gene_seq_data)
    exonCount = int(data['exonCount'])
    cdsEnd = int(data['cdsEns'])

    if endcodon_exonindex-exonCount+1 == 0:
            end_index = exon_end_list[endcodon_exonindex]-cdsEnd
            cds_end = exon_seq[:-end_index]    
    else:     
        end_index = exon_end_list[endcodon_exonindex]-cdsEnd

        for s in range(endcodon_exonindex-exonCount+1):
            end_index += exon_end_list[endcodon_exonindex-s]-exon_start_list[endcodon_exonindex-s]
        cds_end = exon_seq[:-end_index]

    
    if startcodon_exonindex == 0:
            start_index = cdsStart-txStart#exon_start_list[0]の意味のはず
            cds_seq = cds_end[start_index:]    
    else:     
        start_index = cdsStart- exon_start_list[startcodon_exonindex]

        for s in range(startcodon_exonindex):
            start_index += exon_end_list[s]-exon_start_list[s]
        cds_seq = cds_end[start_index:]
    return cds_seq

def get_stopcodon_num(cds_seq)->list:
    matches = re.finditer(r'(?=(CAA)|(?=(CAG))|(?=(CGA))|(?=(TGG)))', cds_seq)
    result = [match.start() for match in matches if (match.start() % 3) == 0]#1の可能性もある
    return result
