#datasetのimportの方法は？
import plus_get_target_num
import plus_transrate_codonnum_to_genenum
AID = plus_get_target_num
conv_sc = plus_transrate_codonnum_to_genenum

#name = input()

def make_preset(name):
    name_df=gene_df[gene_df['name']==name].reset_index()
    chrom = str(name_df.loc[0,'chrom'])
    txStart = str(name_df.loc[0,'txStart'])
    txEnd = str(name_df.loc[0,'txEnd'])
    return name_df, chrom, txStart, txEnd

def where_codonend(name):
    data,cdsStart,start_list,end_list = set_presetnum(name)
    set_num =int(data['cdsEnd'].to_list()[0])-int(data['txStart'].to_list()[0])
    cdsEnd = int(data['cdsEnd'].to_list()[0])
    for s in range(len(start_list)-1):
        if (int(start_list[s]) <= cdsEnd<=int(end_list[s])):
             exon_num =s
            
    return exon_num

def farfrom_last_exon(name,target_start):#exonが1つの場合50%ruleを適応させるべき
    target = target_start.tolist()
    
    data = gene_df[gene_df['name']==name]
    start = data['exonStarts'].to_list()[0]
    end = data['exonEnds'].to_list()[0]
    start = list(start.split(','))
    end = list(end.split(','))
    exonCount = where_codonend(name)
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

def get_CtoT_target_seq(name):
    name_df, chrom, txStart, txEnd = make_preset(name)

    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    available_base_num =AID.get_targetstart_num(name,AID.CtoT_target_start(AID.get_CtoT_target(seq),seq))#配列的に狙えるほう
    
    
    
    start_exon_num= conv_sc.where_codonstart(name)
    end_codon_num =conv_sc.get_stopcodon_num(conv_sc.make_cds_seq(name,conv_sc.get_exon_seq(name,seq),start_exon_num))   
    add_list =conv_sc.add_num(name,end_codon_num,conv_sc.add_num_list(name),start_exon_num)
 
    end_codon_num = conv_sc.correct_added_num(name,end_codon_num,conv_sc.add_num_list(name),start_exon_num)
    
    available_base_num =np.array(available_base_num)
    end_codon_num = np.array(end_codon_num)
    add_list = np.array(add_list)
 
    
    num_able = available_base_num
    num_codon = end_codon_num+add_list

    num_codon = [x for x in num_codon if x <= int(name_df.loc[0,'cdsEnd'])]
    result = np.intersect1d(num_able,num_codon)

    
    result = np.array(farfrom_last_exon(name,result))
    
    return result

def get_AtoG_target_seq(name):
    name_df, chrom, txStart, txEnd = make_preset(name)

    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    available_base_num =AID.get_targetstart_num(name,AID.AtoG_target_start(AID.get_AtoG_target(seq),seq))#配列的に狙えるほう
    
    
    
    start_exon_num= conv_sc.where_codonstart(name)
    end_codon_num =conv_sc.get_stopcodon_num(conv_sc.make_cds_seq(name,conv_sc.get_exon_seq(name,seq),start_exon_num))   
    add_list =conv_sc.add_num(name,end_codon_num,conv_sc.add_num_list(name),start_exon_num)
 
    end_codon_num = conv_sc.correct_added_num(name,end_codon_num,conv_sc.add_num_list(name),start_exon_num)
    
    available_base_num =np.array(available_base_num)
    end_codon_num = np.array(end_codon_num)
    add_list = np.array(add_list)
 
    
    num_able = available_base_num
    num_codon = end_codon_num+add_list

    num_codon = [x for x in num_codon if x <= int(name_df.loc[0,'cdsEnd'])]
    result = np.intersect1d(num_able,num_codon)

    
    result = np.array(farfrom_last_exon(name,result))
    
    return result

