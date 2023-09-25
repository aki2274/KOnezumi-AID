import re
from typing import Tuple

def read_fasta(path : str) -> str:
    sequences = {}
    seq_id = None
    current_seq = ""
    
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):  
                if seq_id is not None:
                    sequences[seq_id] = current_seq
                    current_seq = ""
                seq_id = line[1:]  
            else:
                current_seq += line

        if seq_id is not None and current_seq != "":
            sequences[seq_id] = current_seq

    return sequences

def split_by_space(text: str):
    return re.split('\s+', text)

def sort_dataset(text_list: list) -> str:
    contents =[]
    for i in range(len(text_list)-1):
        contents.append(text_list[i])
    return contents

def read_refFlat(path: str) -> str:
    with open(path, "r", encoding="utf-8") as file:
        lines =file.readlines()
    result = []
    for line in lines:
        line =split_by_space(line)
        contents=sort_dataset(line)
        result.append(contents)
    return result
def make_datasets(path,pathfa)->Tuple[any,dict[str:str]]:
    gene_data =read_refFlat(path)
    col_names =['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
    gene_df=pd.DataFrame(gene_data,columns=col_names)
    
    base_data = read_fasta(path_fa)
    return gene_df,base_data