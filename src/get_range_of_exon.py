from __future__ import annotations
from dataclasses import dataclass


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
