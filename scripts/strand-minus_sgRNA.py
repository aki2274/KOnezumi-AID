#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
from pandas import DataFrame, Series

path = "refFlat.txt"
path_fa = "result.fa"


# In[ ]:


def read_fasta(file_path):
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


# In[ ]:


import re

def split_by_space(text):
    return re.split('\s+', text)


# In[ ]:


def gene_data(text_list):
    contents =[]
    for i in range(len(text_list)-1):
        contents.append(text_list[i])
    return contents


# In[ ]:


def result(path):
    with open(path, "r", encoding="utf-8") as file:
        lines =file.readlines()
    result = []
    for line in lines:
        line =split_by_space(line)
        contents=gene_data(line)
        result.append(contents)
    return result
    


# In[ ]:


gene_data =result(path)
col_names =['geneName','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds']
gene_df=pd.DataFrame(gene_data,columns=col_names)


# In[ ]:


base_data = read_fasta(path_fa)


# In[ ]:





# In[ ]:


def target_base_CT_minus(seq):
    match = re.findall(r'(?=(CC\w{20}))', seq)#ここ不安
    target_list=[]
    for i in range(len(match)):
        target_1 = match[i][19:22]
        target_2 = match[i][18:21]
        target_3 = match[i][17:20]
        if (target_1 =='TCG'or target_1 =='CTG'or target_1 =='TTG')or(target_2 =='TCG'or target_2 =='CTG'or target_2 =='TTG')or(target_3 =='TCG'or target_3 =='CTG'or target_3 =='TTG'):
            target_list.append(True)
        
        else:
            target_list.append(False)
    result =[match[s] for s in range(len(match)) if target_list[s]]
    result = list(dict.fromkeys(result))
    return result


# In[ ]:


def target_base_AG_minus(seq):
    match = re.findall('(?=(\w{18,21}GG))', seq)
    target_list=[]
    for i in range(len(match)):
        target = match[i][0:3]
        
        if target =='CCA':
            target_list.append(True)
        
        else:
            target_list.append(False)
    result =[match[s] for s in range(len(match)) if target_list[s]]
    result = list(dict.fromkeys(result))
    return result


# In[ ]:


def CT_target_start_minus(target_list,qua):    
    result = []
    for i in target_list:
        matches = re.finditer(i, qua)
        for match in matches:
            seq = match.group()[::-1]
            add_seq = re.search(('GCT|GTC|GTT'),seq)
            add_num =add_seq.start()
            start_index = match.start()
            result.append(start_index+21-add_num)
    return result


# In[ ]:


def AG_target_start_minus(target_list,qua):    
    result = []
    for i in target_list:
        matches = re.finditer(i, qua)
        for match in matches:
            start_index = match.start()
            
            result.append(start_index+2)#もしかしたら+3?
    return result


# In[ ]:


def target_gene_start(name,target):
    data = gene_df[gene_df['name']==name]
    num = data['txStart'].astype(int).to_numpy()
    target = np.array(target)
    result = target+num
    #ここ-1が必要かもしれない？
    return result
    


# In[ ]:





# In[ ]:


def transcription(seq):
    complementary = {'A':'T','T':'A','G':'C','C':'G'}
    rev_seq = seq[::-1]
    result = []
    for base in rev_seq:
        base =base.upper()
        result.append(complementary[base])
    
    return ''.join(result)


# In[ ]:


def where_start_codon(name):
    data = gene_df[gene_df['name']==name]
    cdsStart =int(data['cdsEnd'].to_list()[0])#strand-minusなのでcdsEndが開始コドン
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    for s in range(len(start)-1):
        if (int(start[s]) <= cdsStart<=int(end[s])):
            exon_num =s
    return exon_num


# In[ ]:


def transcription_exon_start(name,sequence_data):#transcription名乗ってるけど翻訳される前のエクソンを出してる。
 
    
    
    data = gene_df[gene_df['name']==name]
    
    set_num =int(data['txStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    exon_seq_list = []
    start_list = [num for num in start.split(",")]
    end_list = [num for num in end.split(",")]
    
    for s in range(len(start_list)-1):
        start_num = int(start_list[s])- set_num
        
        end_num= int(end_list[s]) - set_num
    
        
        exon_seq =sequence_data[start_num:end_num]
        exon_seq_list.append(exon_seq)

    result = ''.join(exon_seq_list)
    
    exon_num=where_start_codon(name)
    
    exon_count = data['exonCount']
    exon_start = int(exon_count)-1-exon_num
    
    
    
    if exon_start == 0:
        result = result[:- (int(end_list[len(start_list)-2])-int(data['cdsEnd'].to_list()[0]))]
        

    else:
        exon_start_num =0
        for s in range(exon_start):
            exon_start_num += int(end_list[-2-s])-int(start_list[-2-s])
            
        exon_start_num += int(end_list[-2-exon_start])-int(data['cdsEnd'].to_list()[0])
        result = result[:-(exon_start_num)]
        
    
    return result
    
    


# In[ ]:


def end_codon(cds_seq):#seqはtranscription
    matches = re.finditer('(?=(CAA)|(?=(CAG))|(?=(CGA))|(?=(TGG)))', cds_seq)
    result=[]
    for match in matches:
        start_pos = match.start() 
        
        if (start_pos % 3) == 0:
            result.append(start_pos)
    return result


# In[ ]:


def add_num_list(name):
    data = gene_df[gene_df['name']==name]
    set_num =int(data['cdsStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    add_list=[]
    num = 0
    for s in range(len(start)-1):
        element_list=[]
        element_list.append(num)
        
        num += int(end[s])-int(start[s])
        
        element_list.append(num)
        add_list.append(element_list)
        
    return add_list


# In[ ]:





# In[ ]:


def correct_codon_num(name,end_codon_data,range_list):#txStartからの距離に変換
    
    
    data = gene_df[gene_df['name']==name]
    exonCount = int(data['exonCount'])
    set_num =int(data['txStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    exon_seq_list = []
    start_list = [num for num in start.split(",")]
    end_list = [num for num in end.split(",")]
    
    exon_num=where_start_codon(name)
    
    exon_count = data['exonCount']
    exon_start = int(exon_count)-1-exon_num
    

    if exon_start == 0:#確認してない
        add_num =  (int(end_list[len(start_list)-2])-int(data['cdsEnd'].to_list()[0]))
        

    else:#あってそう
        add_num =0
        for s in range(exon_start): 
            add_num += int(end_list[-2-s])-int(start_list[-2-s])  

        
        add_num += int(end_list[-2-exon_start])-int(data['cdsEnd'].to_list()[0])
     
   
    end_codon_data = np.array(end_codon_data)   
    end_codon_data = end_codon_data + add_num
    
    end_codon_data = abs(end_codon_data - range_list[exonCount-1][1]+1)#+1するかは微妙なライン、この後の処理に合わせて変えるべき
    
    result = end_codon_data.tolist()
    
    return result


# In[ ]:


def codon_in_exon(name,correct_num,range_list):

    data = gene_df[gene_df['name']==name]
    set_num =int(data['txStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    add_result=[]
    
    
    
    indices = []
    for number in correct_num:
        for index, interval in enumerate(range_list):
            if int(interval[0]) <= int(number) < int(interval[1]):
                indices.append(index)
                break
            
    for i in indices:
        if i == 0:
            add_result.append(set_num)
        else:
            add_result.append(int(start[i]))

    return add_result
    


# In[ ]:


def add_num_correct(name,correct_num, range_list):
    data = gene_df[gene_df['name']==name]
    set_num =int(data['cdsStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    add_result=[]
         
     
    indices = []
    for number in correct_num:
        for index, interval in enumerate(range_list):
            if int(interval[0]) <= int(number) < int(interval[1]):
                indices.append(index)
                break
    n =0        
    for i in indices:
        add_result.append(int(correct_num[n])-int(range_list[i][0]))
        n += 1

    return add_result


# In[ ]:


def where_end_codon(name):
    data = gene_df[gene_df['name']==name]
    cdsStart =int(data['cdsStart'].to_list()[0])
    set_num =int(data['cdsStart'].to_list()[0])-int(data['txStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    for s in range(len(start)-1):
        if (int(start[s]) <= cdsStart<=int(end[s])):
             exon_num =s
            
    return exon_num


# In[ ]:


def farfrom_last_exon(name,target_start):
    target = target_start.tolist()
    
    data = gene_df[gene_df['name']==name]
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]
    start = list(start.split(','))
    end = list(end.split(','))
    exonCount = where_end_codon(name)
    count = int(data['exonCount'])
    start.pop()
    end.pop()
    target_exon_num=[]
    for i in range(len(target)):
        target_num = target[i]
        for s in range(len(start)):
            if (int(start[s]) <= target_num <=int(end[s])):
                target_exon_num.append(s)
    
    bool_list = []
    if count != 1:
        for t in range(len(target_exon_num)):
        
            if target_exon_num[t]==exonCount:
                bool_list.append(False)

            elif target_exon_num[t]==exonCount+1:
                if target[t]-int(start[exonCount-1])<50:
                    bool_list.append(False)
                else:
                    bool_list.append(True)
                
            else:
                bool_list.append(True)
            
    else:
        for t in range(len(target_exon_num)):
            bool_list.append(True)  
    
    
            

    result =[target_start[s] for s in range(len(target_start)) if bool_list[s]]

    return result
    


# In[ ]:





# In[ ]:


def CT_target_base_strand_minus(name):
    name_df=gene_df[gene_df['name']==name].reset_index()
    chrom = str(name_df.loc[0,'chrom'])
    txStart = str(name_df.loc[0,'txStart'])
    txEnd = str(name_df.loc[0,'txEnd'])
    
    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    available_base_num =target_gene_start(name,CT_target_start_minus(target_base_CT_minus(seq),seq))#配列的に狙えるほう
    
    possibility_end_codon= end_codon(transcription(transcription_exon_start(name,seq)))
    a= codon_in_exon(name,correct_codon_num(name,possibility_end_codon,add_num_list(name)),add_num_list(name))
    b= add_num_correct(name,correct_codon_num(name,possibility_end_codon,add_num_list(name)),add_num_list(name))
    num_codon = np.array(a)+np.array(b)                                               
                                                   
    num_codon = [x for x in num_codon if x >= int(name_df.loc[0,'cdsStart'])]
    result = np.intersect1d(available_base_num,num_codon)
    
    result = np.array(farfrom_last_exon(name,result))
    
    return  result


# In[ ]:


def AG_target_base_strand_minus(name):
    name_df=gene_df[gene_df['name']==name].reset_index()
    chrom = str(name_df.loc[0,'chrom'])
    txStart = str(name_df.loc[0,'txStart'])
    txEnd = str(name_df.loc[0,'txEnd'])
    
    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    available_base_num =target_gene_start(name,AG_target_start_minus(target_base_AG_minus(seq),seq))#配列的
    
    possibility_end_codon= end_codon(transcription(transcription_exon_start(name,seq)))
    a= codon_in_exon(name,correct_codon_num(name,possibility_end_codon,add_num_list(name)),add_num_list(name))
    b= add_num_correct(name,correct_codon_num(name,possibility_end_codon,add_num_list(name)),add_num_list(name))
    num_codon = np.array(a)+np.array(b)                                               
                                                   
    num_codon = [x for x in num_codon if x >= int(name_df.loc[0,'cdsStart'])]
    result = np.intersect1d(available_base_num,num_codon)
    
    result = np.array(farfrom_last_exon(name,result))
    
    return  result


# In[ ]:


def CT_target_minus(name,target_start):
    name_df=gene_df[gene_df['name']==name].reset_index()
    chrom = str(name_df.loc[0,'chrom'])
    txStart = str(name_df.loc[0,'txStart'])
    txEnd = str(name_df.loc[0,'txEnd'])
    
    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    
    target = (target_start-int(txStart)).tolist()
    result=[]
    for i in range(len(target)):
        s = target[i]
        pro_sgRNA = seq[s-21:s+2]
        matches = re.finditer('(?=(CC))',pro_sgRNA)
        for add_seq in matches:
            add_num =add_seq.start()
            if add_num<=2:
                sgRNA = seq[s-21+add_num:s+2+add_num]
                result.append(sgRNA)
    
    
    return result


# In[ ]:





# In[ ]:


def AG_target_minus(name,target_start):
    name_df=gene_df[gene_df['name']==name].reset_index()
    chrom = str(name_df.loc[0,'chrom'])
    txStart = str(name_df.loc[0,'txStart'])
    txEnd = str(name_df.loc[0,'txEnd'])
    
    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    
    target = (target_start-int(txStart)).tolist()
    result=[]
    for i in range(len(target)):
        s = target[i]
        pro_sgRNA = seq[s-2:s+21]
        r_pro_sgRNA = pro_sgRNA[::-1]
        matches = re.finditer('(?=(GG))',r_pro_sgRNA)
        for add_seq in matches:
            add_num =add_seq.start()
            if 0<= add_num <= 3:
                sgRNA = seq[s-2-add_num:s+21-add_num]
                result.append(sgRNA)
    
    
    return result


# In[ ]:


def target_from_genename(genename):
    que_df = gene_df[gene_df['geneName'] == genename].reset_index()
    transcription_num = len(que_df)

    if transcription_num == 1:
        name = que_df.loc[0, 'name']

        CT = CT_target_minus(name, CT_target_base_strand_minus(name))
        AG = AG_target_minus(name, AG_target_base_strand_minus(name))

        return f"CT-editingは'{CT}'. AG-editingは'{AG}'."

    else:
        name = que_df.loc[0, 'name']

        CT_num = CT_target_base_strand_minus(name)
        AG_num = AG_target_base_strand_minus(name)

        for i in range(1, len(que_df)):
            name = que_df.loc[i, 'name']

            pCT_num = CT_target_base_strand_minus(name)
            pAG_num = AG_target_base_strand_minus(name)

            CT_num = np.intersect1d(CT_num, pCT_num)
            AG_num = np.intersect1d(AG_num, pAG_num)

        CT = CT_target_minus(name, CT_num)
        AG = AG_target_minus(name, AG_num)

        return f"CT-editingは'{CT}'. AG-editingは'{AG}'."


# In[ ]:


target_from_genename('Gk2')


# In[ ]:


target_from_genename('Xkr4')

