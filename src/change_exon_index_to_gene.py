from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass
class DataClass:
    orf_seq: str
    data: pd.DataFrame
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


#####
# Change candidate stopcodon index in exon to index in genome
#####


# 　エクソンごとの範囲を取得、candidateが何番目のエクソンか知る、
def get_range_of_exon(ds: DataClass, cdsStart_exon_index: int) -> int:
    # エクソンごとの範囲を取得
    exon_range_list = []
    num = 0
    if cdsStart_exon_index == 0:
        for s in range(len(ds.exon_start_list)):
            element_list = []
            element_list.append(num)
            num += ds.exon_end_list[s] - ds.exon_start_list[s]
            element_list.append(num)
            exon_range_list.append(element_list)
    else:
        for s in range(len(ds.exon_start_list) - cdsStart_exon_index):
            element_list = []
            element_list.append(num)
            num += (
                ds.exon_end_list[s + cdsStart_exon_index]
                - ds.exon_start_list[s + cdsStart_exon_index]
            )
            element_list.append(num)
            exon_range_list.append(element_list)
    return exon_range_list


def get_candidate_codon_exon_index(
    candidate_codon_index_list: list[int], exon_range_list
):
    # 何番目のエクソンであるかのリストを作成
    exon_index_list = []
    for t in range(len(candidate_codon_index_list)):
        for s in range(len(exon_range_list)):
            if (
                exon_range_list[s][0]
                <= candidate_codon_index_list[t]
                <= exon_range_list[s][1]
            ):
                exon_index_list.append(s)
    return exon_index_list


# 足すべき数値を取得、足す
def add_num_to_orf_index(
    ds,
    candidate_codon_index_list,
    exon_range_list,
    exon_index_list,
    startcodon_exon_index,
) -> int:
    add_num_list = [
        (
            ds.exon_start_list[i + startcodon_exon_index]
            - exon_range_list[i][0]
            - ds.txStart
        )
        for i in exon_index_list
    ]
    candidate_stopcodon_index = np.array(candidate_codon_index_list) + np.array(
        add_num_list
    )
    return candidate_stopcodon_index.tolist()
