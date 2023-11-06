from __future__ import annotations
import re
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
# Seaech candedate stopcodon index
#####


def get_exon_seq(ds: DataClass) -> str:
    # make spliced seq (NOT cds)
    exon_seq_list = []
    for s in range(len(ds.exon_start_list)):
        start_num = ds.exon_start_list[s] - ds.txStart
        end_num = ds.exon_end_list[s] - ds.txStart
        exon_seq = ds.orf_seq[start_num:end_num]
        exon_seq_list.append(exon_seq)
    return "".join(exon_seq_list)


def get_startcodon_exonindex(ds: DataClass) -> int:
    # search the position of startcodon in exon number
    for s in range(len(ds.exon_start_list)):
        if ds.exon_start_list[s] <= ds.cdsStart <= ds.exon_end_list[s]:
            exon_index = s
    return exon_index


def get_stopcodon_exonindex(ds: DataClass) -> int:
    # search the position of stopcodon in exon number
    for s in range(len(ds.exon_start_list)):
        if ds.exon_start_list[s] <= ds.cdsEnd <= ds.exon_end_list[s]:
            exon_index = s
    return exon_index


def get_cdsseq(
    ds: DataClass,
    exon_seq: str,
    startcodon_exonindex: int,
    endcodon_exonindex: int,
):
    if endcodon_exonindex - ds.exonCount + 1 == 0:
        end_index = ds.exon_end_list[endcodon_exonindex] - ds.cdsEnd
        cds_end = exon_seq[:-end_index]
    else:
        end_index = ds.exon_end_list[endcodon_exonindex] - ds.cdsEnd

        for s in range(endcodon_exonindex - ds.exonCount + 1):
            end_index += (
                ds.exon_end_list[endcodon_exonindex - s]
                - ds.exon_start_list[endcodon_exonindex - s]
            )
        cds_end = exon_seq[:-end_index]

    if startcodon_exonindex == 0:
        start_index = ds.cdsStart - ds.exon_start_list[0]
        cds_seq = cds_end[start_index:]
    else:
        start_index = ds.cdsStart - ds.exon_start_list[startcodon_exonindex]

        for s in range(startcodon_exonindex):
            start_index += ds.exon_end_list[s] - ds.exon_start_list[s]
        cds_seq = cds_end[start_index:]
    return cds_seq


def get_candidate_stopcodon_num(cds_seq: str) -> list:
    matches = re.finditer(r"(?=(CAA)|(?=(CAG))|(?=(CGA))|(?=(TGG)))", cds_seq)
    candidate_codon_index_list = [
        match.start() for match in matches if (match.start() % 3) == 0
    ]
    return candidate_codon_index_list


#####
# Change candidate stopcodon index in exon to index in genome
#####


# 　エクソンごとの範囲を取得、candidateが何番目のエクソンか知る、
def get_number_of_exon(
    ds: DataClass, candidate_codon_index_list: list[int], cdsStart_exon_index: int
) -> int:
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
    return exon_range_list, exon_index_list


# 足すべき数値を取得、足す
def a(
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
