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


def target_base_CT_plus(seq):
    match = re.findall(r'(?=(\w{18,20}GG))', seq)#ここ不安
    target_list=[]
    for i in range(len(match)):
        target = match[i][0:3]
        if target =='CAA'or target =='CAG'or target =='CGA':
            target_list.append(True)
        
        else:
            target_list.append(False)
    result =[match[s] for s in range(len(match)) if target_list[s]]
    result = list(dict.fromkeys(result))
    return result


# In[ ]:


def target_base_AG_plus(seq):
    match = re.findall('(?=(CC\w{21}))', seq)
    target_list=[]
    for i in range(len(match)):
        target_1 = match[i][20:23]
        target_2 = match[i][19:22]
        target_3 = match[i][18:21]
        target_4 = match[i][17:20]
        if target_1 =='TGG'or target_2 =='TGG'or target_3 =='TGG'or target_4 =='TGG':#kここはCCA?
            target_list.append(True)
        
        else:
            target_list.append(False)
    result =[match[s] for s in range(len(match)) if target_list[s]]
    result = list(dict.fromkeys(result))
    return result


# In[ ]:


def CT_target_start(target_list,que):    
    result = []
    for i in target_list:
        matches = re.finditer(i, que)
        for match in matches:
            start_index = match.start()
            result.append(start_index)
    return result


# In[ ]:


def AG_target_start(target_list,que):    
    result = []
    for i in target_list:
        matches = re.finditer(i, que)
        for match in matches:
            seq = match.group()[::-1]
            add_seq = re.search('GGT',seq)
            add_num =add_seq.start()
            
            start_index = match.start()
            result.append(start_index+20-add_num)
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


def exon_data(name,sequence_data):
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
        
    
    
    return ''.join(exon_seq_list)


# In[ ]:


def where_start_codon(name):
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


def start_codon(name,exon_data,exon_num):
    data = gene_df[gene_df['name']==name]
    set_num =int(data['cdsStart'].to_list()[0])-int(data['txStart'].to_list()[0])
    CDs =''
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    cds_list=[]
    if exon_num == 0:
            cds_list.append(exon_data[set_num:])
            
    
    
    else:     
        start_num = int(data['cdsStart'].to_list()[0])- int(start[exon_num])

        for s in range(exon_num):
            start_num += int(end[s])-int(start[s])
        cds_list.append(exon_data[start_num:])
        
    return ''.join(cds_list)


# In[ ]:


def end_codon(cds_seq):
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


def add_num_correct(name,number_list, interval_list,exon_num):
    data = gene_df[gene_df['name']==name]
    set_num =int(data['cdsStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    add_result=[]
    
    number_list = np.array(number_list)
    start_num=int(interval_list[exon_num][0])
    number_list = number_list + (set_num-int(start[exon_num])+start_num)
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


# In[ ]:


def add_num(name,number_list, interval_list,exon_num):

    data = gene_df[gene_df['name']==name]
    set_num =int(data['cdsStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    add_result=[]
    
    number_list = np.array(number_list)
    start_num=int(interval_list[exon_num][0])
    number_list = number_list + (set_num-int(start[exon_num])+start_num)
    number_list = number_list.tolist()#ここでエクソン上で何番目かになってる（cdsからエクソン上に）

    
    
    indices = []
    for number in number_list:
        for index, interval in enumerate(interval_list):
          
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


def where_end_codon(name):
    data = gene_df[gene_df['name']==name]
    cdsStart =int(data['cdsEnd'].to_list()[0])
    set_num =int(data['cdsEnd'].to_list()[0])-int(data['txStart'].to_list()[0])
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]#最後に空白があるからその対処
    start = list(start.split(','))
    end = list(end.split(','))
    for s in range(len(start)-1):
        if (int(start[s]) <= cdsStart<=int(end[s])):
             exon_num =s
            
    return exon_num


# In[ ]:


def farfrom_last_exon(name,target_start):#exonが1つの場合うまくいかないはず
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
            
            elif target_exon_num[t]==exonCount-1:
                if int(end[exonCount-1])-target[t]<50:
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


def CT_target_base_strand_plus(name):
    name_df=gene_df[gene_df['name']==name].reset_index()
    chrom = str(name_df.loc[0,'chrom'])
    txStart = str(name_df.loc[0,'txStart'])
    txEnd = str(name_df.loc[0,'txEnd'])
    
    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    available_base_num =target_gene_start(name,CT_target_start(target_base_CT_plus(seq),seq))#配列的に狙えるほう
    
    
    
    start_exon_num= where_start_codon(name)
    end_codon_num =end_codon(start_codon(name,exon_data(name,seq),start_exon_num))   
    add_list =add_num(name,end_codon_num,add_num_list(name),start_exon_num)
 
    end_codon_num=add_num_correct(name,end_codon_num,add_num_list(name),start_exon_num)
    
    available_base_num =np.array(available_base_num)
    end_codon_num = np.array(end_codon_num)
    add_list = np.array(add_list)
 
    
    num_able = available_base_num
    num_codon = end_codon_num+add_list

    num_codon = [x for x in num_codon if x <= int(name_df.loc[0,'cdsEnd'])]
    result = np.intersect1d(num_able,num_codon)

    
    result = np.array(farfrom_last_exon(name,result))
    
    return result


# In[ ]:


def AG_target_base_strand_plus(name):
    name_df=gene_df[gene_df['name']==name].reset_index()
    chrom = str(name_df.loc[0,'chrom'])
    txStart = str(name_df.loc[0,'txStart'])
    txEnd = str(name_df.loc[0,'txEnd'])
    
    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    available_base_num =target_gene_start(name,AG_target_start(target_base_AG_plus(seq),seq))#配列的に狙えるほう

    
    start_exon_num= where_start_codon(name)
    end_codon_num =end_codon(start_codon(name,exon_data(name,seq),start_exon_num))   
    add_list =add_num(name,end_codon_num,add_num_list(name),start_exon_num)
    
  
    end_codon_num=add_num_correct(name,end_codon_num,add_num_list(name),start_exon_num)
    
    available_base_num =np.array(available_base_num)
    end_codon_num = np.array(end_codon_num)
    add_list = np.array(add_list)

    num_able = available_base_num
    num_codon = end_codon_num+add_list
    
    num_codon = [x for x in num_codon if x <= int(name_df.loc[0,'cdsEnd'])]
    
    
    result = np.intersect1d(num_able,num_codon)
    
    result =  np.array(farfrom_last_exon(name,result))
    
    return result


# In[ ]:





# In[ ]:


def CT_target_plus(name,target_start):
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
        pro_sgRNA = seq[s-1:s+24]
        r_pro_sgRNA = pro_sgRNA[::-1]
        matches = re.finditer('(?=(GG))',r_pro_sgRNA)
       
        for add_seq in matches:
            add_num =add_seq.start()
            if 2 <= int(add_num) <= 4:
                sgRNA = seq[s+1-add_num:s+24-add_num]
                result.append(sgRNA)
    
    
    return result


# In[ ]:


def AG_target_plus(name,target_start):
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
        pro_sgRNA = seq[s-20:s+3]
        matches = re.finditer('(?=(CC))',pro_sgRNA)
        
        for add_seq in matches:
            add_num =add_seq.start()
            if 0 <= int(add_num) <= 3:
                sgRNA = seq[s-20+add_num:s+3+add_num]
                result.append(sgRNA)
        
        
      
    
    
    return result


# In[ ]:


def target_from_genename(genename):#出力は要検討
    que_df=gene_df[gene_df['geneName']==genename].reset_index()
    transcription_num = len(que_df)
    
    if transcription_num==1:
        name =que_df.loc[0,'name']
        
        CT = CT_target_plus(name,CT_target_base_strand_plus(name))
        AG = AG_target_plus(name,AG_target_base_strand_plus(name))

        
        return f"CT-editingは'{CT}'.AG-editingは'{AG}'."
        
    else:
        name =que_df.loc[0,'name']
        
        CT_num = CT_target_base_strand_plus(name)
        AG_num = AG_target_base_strand_plus(name)
        
        for i in range(1,len(que_df)):
            name =que_df.loc[i,'name']
            
            pCT_num = CT_target_base_strand_plus(name)
            pAG_num = AG_target_base_strand_plus(name)
                                     
            CT_num = np.intersect1d(CT_num,pCT_num)
            AG_num = np.intersect1d(AG_num,pAG_num)
        
        CT = CT_target_plus(name,CT_num)
        AG = AG_target_plus(name,AG_num)
        
        return f"CT-editingは'{CT}'.AG-editingは'{AG}'."
        
        
        
        
    
    


# In[ ]:


target_from_genename('Mul1')


# In[ ]:


target_from_genename('Sulf1')


# In[ ]:


target_from_genename('Trp53')

