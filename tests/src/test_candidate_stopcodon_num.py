'''
import pytest
import pandas as pd
import src.get_candidate_stopcodon_index as gcsi

#################################
# find candidate stopcodon index
#################################

def make_dataframe(name:list,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds):
    gene_data =[]
    for a in zip(name,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds):     
        gene_data.append(['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds'])      

    col_names =['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds'] 
    gene_df=pd.DataFrame(gene_data,columns=col_names)

sample_name =['A']

sample_seq = ''

expected =''

sample_df =''


@pytest.mark.parametrize("transcripts_name,gene_df,gene_seq_data,expected", zip(sample_name,sample_df,sample_seq, expected))
def test_get_sxon_seq(transcripts_name,gene_df,gene_seq_data,expected):
    assert gcsi.get_exon_seq(transcripts_name,gene_df,gene_seq_data) == excepted
'''