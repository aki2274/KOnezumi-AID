import re
from typing import Tuple
import numpy as np
import pandas as pd

def set_presetnum(name:str):
    data = gene_df[gene_df['name']==name]
    cdsStart =int(data['cdsStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start_list = list(start.split(','))
    end_list = list(end.split(','))
    return data,cdsStart,start_list,end_list

def target_gene_start(name,target):
    data = gene_df[gene_df['name']==name]
    num = data['txStart'].astype(int).to_numpy()
    target = np.array(target)
    result = target+num
    return result    

def get_exon_seq(name,sequence_data:str)->str:
    data,cdsStart,start_list,end_list = set_presetnum(name)
    set_num =int(data['txStart'])  
    exon_seq_list =[] 
    for s in range(len(start_list)-1):#最後が空白なので-1している
        start_num = int(start_list[s])- set_num
        end_num= int(end_list[s]) - set_num
        exon_seq =sequence_data[start_num:end_num]
        exon_seq_list.append(exon_seq)
    return ''.join(exon_seq_list)

def where_codonstart(name)->int:
    data,cdsStart,start_list,end_list = set_presetnum(name)
    set_num =int(data['cdsStart'].to_list()[0])-int(data['txStart'].to_list()[0])
    for s in range(len(start_list)-1):
        if (int(start_list[s]) <= cdsStart<=int(end_list[s])):
             exon_num =s
            
    return exon_num

def make_cds_seq(name,exon_data,exon_num):
    data,cdsStart,start_list,end_list = set_presetnum(name)
    set_num =int(data['cdsStart'].to_list()[0])-int(data['txStart'].to_list()[0])
    CDs =''
    cds_list=[]
    if exon_num == 0:
            cds_list.append(exon_data[set_num:])
            
    else:     
        start_num = int(data['cdsStart'].to_list()[0])- int(start_list[exon_num])

        for s in range(exon_num):
            start_num += int(end_list[s])-int(start_list[s])
        cds_list.append(exon_data[start_num:])
        
    return ''.join(cds_list)

def get_stopcodon_num(cds_seq)->list:
    matches = re.finditer('(?=(CAA)|(?=(CAG))|(?=(CGA))|(?=(TGG)))', cds_seq)
    result=[]
    for match in matches:
        start_pos = match.start() 
        
        if (start_pos % 3) == 0:
            result.append(start_pos)
    return result

def add_num_list(name)->list:
    data,cdsStart,start_list,end_list = set_presetnum(name)
    set_num =int(data['cdsStart'].to_list()[0])
    add_list=[]
    num = 0
    for s in range(len(start_list)-1):
        element_list=[]
        element_list.append(num)
        
        num += int(end_list[s])-int(start_list[s])
        
        element_list.append(num)
        add_list.append(element_list)
        
    return add_list

def correct_added_num(name,number_list, interval_list,exon_num):#エクソンのはじめから何番目であるかに変換する
    data,cdsStart,start_list,end_list = set_presetnum(name)
    set_num =int(data['cdsStart'].to_list()[0])
    add_result=[]
    
    number_list = np.array(number_list)
    start_num=int(interval_list[exon_num][0])
    number_list = number_list + (set_num-int(start_list[exon_num])+start_num)
    number_list = number_list.tolist()#ここでエクソン上で何番目かになってる（cdsからエクソン上に）
    
    
    indices = []
    for number in number_list:
        for index, interval in enumerate(interval_list):
           
            if int(interval[0]) <= int(number) < int(interval[1]):
                indices.append(index)
                break
    n =0        
    for i in indices:
        
        if i == 0:
            add_result.append(set_num)
            n +=1
        else:
            add_result.append(int(number_list[n])-int(interval_list[i][0]))
            n += 1

    return add_result

def add_num(name,number_list, interval_list,exon_num):
    data,cdsStart,start_list,end_list = set_presetnum(name)
    set_num =int(data['cdsStart'].to_list()[0])
    
    number_list = np.array(number_list)
    start_num=int(interval_list[exon_num][0])
    number_list = number_list + (set_num-int(start_list[exon_num])+start_num)
    number_list = number_list.tolist()#ここでエクソン上で何番目かになってる（cdsからエクソン上に）

    
    
    indices = []
    for number in number_list:
        for index, interval in enumerate(interval_list):
          
            if int(interval[0]) <= int(number) < int(interval[1]):
                indices.append(index)
                break
    add_result = []     
    for i in indices:
        if i == 0:
            add_result.append(set_num)
        else:
            add_result.append(int(start_list[i]))

    return add_result

def where_codonend(name):
    data,cdsStart,start_list,end_list = set_presetnum(name)
    set_num =int(data['cdsEnd'].to_list()[0])-int(data['txStart'].to_list()[0])
    cdsEnd = int(data['cdsEnd'].to_list()[0])
    for s in range(len(start_list)-1):
        if (int(start_list[s]) <= cdsEnd<=int(end_list[s])):
             exon_num =s
            
    return exon_num