import re
from typing import Tuple
import numpy as np
import pandas as pd

def get_re_comp(seq:str)->str:
    complementary = {'A':'T','T':'A','G':'C','C':'G'}
    rev_seq = seq[::-1]
    result = []
    result = [complementary[base.upper()] for base in rev_seq]
    return ''.join(result)

def find_ct_target(seq:str)->list:
    match = re.findall(r'(?=(\w{21}GG))', seq)
    target_list=[]
    target_list = [True if target[1:4] in ('CAA', 'CAG', 'CGA') or target[2:5] in ('CAA', 'CAG', 'CGA') or target[3:6] in ('CAA', 'CAG', 'CGA') else False for target in match]      
    result =[match[s] for s in range(len(match)) if target_list[s]]
    result = list(dict.fromkeys(result))
    return result

def get_ct_target_num(target_list,que)->list:#CXXのCの位置を出す
    result = []
    for i in target_list:
        matches = re.finditer(i, que)
        result = [match.start() for match in matches]
    result = list(dict.fromkeys(result))
    return result

def find_ag_target(seq:str)->list:
    match = re.findall(r'(?=(CC\w{18,21}))', seq)
    target_list=[]
    target_list =[True if target[20:23] =='TGG' or target[19:22] =='TGG'or target[18:21] =='TGG' or target[17:20] =='TGG' else False for target in match]
    result =[match[s] for s in range(len(match)) if target_list[s]]
    result = list(dict.fromkeys(result))
    return result


def get_ag_target_num(target_list,que)->list:#TGGのTの位置を出す
    result = []
    for i in target_list:
        matches = re.finditer(i, que)
        for match in matches:
            rev_match = match.group()[::-1]
            rev_seq = re.search('GGT',rev_match)
            add_num =rev_seq.start()
            start_index = match.start()
            result.append(start_index+len(rev_match)-3-add_num)
    result = list(dict.fromkeys(result))
    return result
'''
def get_targetstart_num(name,target,gene_df)->list:
    data = gene_df[gene_df['name']==name]
    num = data['txStart'].astype(int).to_numpy()
    target = np.array(target)
    result = target+num
    return result    
'''