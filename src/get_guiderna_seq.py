from __future__ import annotations
import re
from typing import Tuple
import numpy as np
import pandas as pd


complementary = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

def get_revcomp(seq: str) -> str:
    return "".join(complementary[base.upper()] for base in seq[::-1])


def find_ct_target_seq(seq: str) -> list[str]:
    matches = re.findall(r"(?=(\w{21}GG))", seq)
    # Check for the presence of "CAA", "CAG", or "CGA" in the next 3 positions.
    targets = [
        match for match in matches
        if any(match[i:i+3] in {"CAA", "CAG", "CGA"} for i in range(1, 4))
    ]
    # Remove duplicates while preserving order
    return list(dict.fromkeys(targets))


def get_ct_target_num(que, target_list) -> list:  # CXXのCの位置を出す
    result = []
    for i in target_list:
        matches = re.finditer(i, que)
        num = [
            re.search(r"(CAA|CAG|CGA)", str(match.group())).start() + match.start()
            for match in matches
        ]
        result.append(num[0])
    result = list(dict.fromkeys(result))
    return result


def find_ag_target(seq: str) -> list:
    match = re.findall(r"(?=(CC\w{18,21}))", seq)
    target_list = []
    target_list = [
        True
        if target[20:23] == "TGG"
        or target[19:22] == "TGG"
        or target[18:21] == "TGG"
        or target[17:20] == "TGG"
        else False
        for target in match
    ]
    result = [match[s] for s in range(len(match)) if target_list[s]]
    result = list(dict.fromkeys(result))
    return result


def get_ag_target_num(que, target_list) -> list:  # TGGのTの位置を出す
    result = []
    for i in target_list:
        matches = re.finditer(i, que)
        for match in matches:
            rev_match = match.group()[::-1]
            rev_seq = re.search("GGT", rev_match)
            add_num = rev_seq.start()
            start_index = match.start()
            result.append(start_index + len(rev_match) - 3 - add_num)
    result = list(dict.fromkeys(result))
    return result


"""
def get_targetstart_num(name,target,gene_df)->list:
    data = gene_df[gene_df['name']==name]
    num = data['txStart'].astype(int).to_numpy()
    target = np.array(target)
    result = target+num
    return result    
"""
