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
    return name_df, crhom, txStart, txEnd

def get_CtoT_target_seq(name):
    name_df, crhom, txStart, txEnd = make_preset(name)

    query = f'{name}::{chrom}:{txStart}-{txEnd}'
    seq = base_data[query]
    available_base_num =AID.get_targetstart_num(name,AID.CtoT_target_start(AID.get_CtoT_target(seq),seq))#配列的に狙えるほう
    
    
    
    start_exon_num= conv_sc.where_codonstart(name)
    end_codon_num =conv_sc.where_codonend(conv_sc.make_cds_seq(name,conv_sc.get_exon_seq(name,seq),start_exon_num))   
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

