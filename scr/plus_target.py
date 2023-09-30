import re

def get_CtoT_target(seq)->list:
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

def CtoT_target_start(target_list,que)->list:#CXXのCの位置を出す
    result = []
    for i in target_list:
        matches = re.finditer(i, que)
        for match in matches:
            start_index = match.start()
            result.append(start_index)
    result = list(dict.fromkeys(result))
    return result

def get_AtoG_target(seq)->list:
    match = re.findall('(?=(CC\w{21}))', seq)
    target_list=[]
    for i in range(len(match)):
        target_1 = match[i][20:23]
        target_2 = match[i][19:22]
        target_3 = match[i][18:21]
        target_4 = match[i][17:20]
        
        if target_1 =='TGG'or target_2 =='TGG'or target_3 =='TGG'or target_4 =='TGG':
            target_list.append(True)
        else:
            target_list.append(False)
    
    result =[match[s] for s in range(len(match)) if target_list[s]]
    result = list(dict.fromkeys(result))
    return result

def AtoG_target_start(target_list,que):#TGGのTの位置を出す
    result = []
    for i in target_list:
        matches = re.finditer(i, que)
        for match in matches:
            seq = match.group()[::-1]
            add_seq = re.search('GGT',seq)
            add_num =add_seq.start()
            start_index = match.start()
            result.append(start_index+20-add_num)
    result = list(dict.fromkeys(result))
    return result

def get_targetstart_num(name,target)->list:
    data = gene_df[gene_df['name']==name]
    num = data['txStart'].astype(int).to_numpy()
    target = np.array(target)
    result = target+num
    return result    