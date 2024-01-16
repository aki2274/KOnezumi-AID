from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class GeneData:
    orf_seq: str
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
# 　エクソンごとの範囲を取得、candidate stopcodon が何番目のエクソンか知る、


def get_exon_range(ds: GeneData):
    # エクソンごとの範囲を取得
    exon_range_list = []
    num = 0
    for s in range(len(ds.exon_start_list)):
        element_list = []
        element_list.append(num)
        num += ds.exon_end_list[s] - ds.exon_start_list[s]
        element_list.append(num)
        exon_range_list.append(element_list)

    return exon_range_list


# cdsでのindexをエクソンのindexに変換
def translate_incds_index_to_exon(
    ds: GeneData,
    candidate_stopcodon_cds_index_list: list[int],
    exon_range_list: list[int],
    cdsStart_exon_index: int,
) -> list[int]:
    candidate_index_array = np.array(candidate_stopcodon_cds_index_list)
    # エクソン上でのエクソンのはじめからstartcodonが含まれるエクソンの先頭までの距離
    dist_to_startcodon_exon = exon_range_list[cdsStart_exon_index][0]
    # エクソン上でのエクソンのはじめからstartcodon(cdsStart)までの距離。
    dist_to_startcodon = (
        ds.cdsStart
        - int(ds.exon_start_list[cdsStart_exon_index])
        + dist_to_startcodon_exon
    )
    # ここでエクソン上で何番目かになってる（cdsからエクソン上に）
    result = (candidate_index_array + dist_to_startcodon).tolist()
    return result


def get_candidate_stopcodon_exon_num(
    candidate_stopcodon_index_inexon: list[int], exon_range_list
):
    # 何番目のエクソンであるかのリストを作成
    exon_index_list = []
    for t in range(len(candidate_stopcodon_index_inexon)):
        for s in range(len(exon_range_list)):
            if (
                exon_range_list[s][0]
                <= candidate_stopcodon_index_inexon[t]
                <= exon_range_list[s][1]
            ):
                exon_index_list.append(s)
    return exon_index_list


# 足すべき数値を取得、足す
def translate_idex_in_exon_to_orf(
    ds: GeneData,
    candidate_codon_index_inexon_list: list[int],
    exon_range_list: list[int],
    exon_index_list: list[int],
) -> int:
    # genome上でのcandidatestopcodonを含むexonの先頭までの距離を取得、exon分を2度足しているので引く
    add_num_to_genome_index = [
        (ds.exon_start_list[i] - ds.txStart - exon_range_list[i][0])
        for i in exon_index_list
    ]
    candidate_stopcodon_index = np.array(candidate_codon_index_inexon_list) + np.array(
        add_num_to_genome_index
    )
    return candidate_stopcodon_index.tolist()
