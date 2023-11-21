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


#####
# Seaech candedate stopcodon index
#####


def get_exon_seq(ds: GeneData) -> str:
    # make seq spliced by exons (NOT cds)
    exon_seq_list = []
    for s in range(len(ds.exon_start_list)):
        exon_start = ds.exon_start_list[s] - ds.txStart
        exon_end = ds.exon_end_list[s] - ds.txStart
        exon_seq = ds.orf_seq[exon_start:exon_end]
        exon_seq_list.append(exon_seq)
    return "".join(exon_seq_list)


def get_startcodon_exon_num(ds: GeneData) -> int:
    # search the position of startcodon by exon number
    for s in range(len(ds.exon_start_list)):
        if ds.exon_start_list[s] <= ds.cdsStart <= ds.exon_end_list[s]:
            exon_index = s
    return exon_index


def get_stopcodon_exon_num(ds: GeneData) -> int:
    # search the position of stopcodon by exon number
    for s in range(len(ds.exon_start_list)):
        if ds.exon_start_list[s] <= ds.cdsEnd <= ds.exon_end_list[s]:
            exon_index = s
    return exon_index


def get_cdsseq(
    ds: GeneData,
    exon_seq: str,
    startcodon_exon_num: int,
    endcodon_exon_num: int,
) -> str:
    # return cds seq
    stopcodon_index = ds.exon_end_list[endcodon_exon_num] - ds.cdsEnd
    if endcodon_exon_num - ds.exonCount + 1 == 0:  # if stopcodon is in the last exon
        exon_seq_to_stopcodon = exon_seq[:-stopcodon_index]
    else:  # if stopcodon is not in the last exon. (add the len of exons after stopcodon)
        for s in range(endcodon_exon_num - ds.exonCount + 1):
            stopcodon_index += (
                ds.exon_end_list[endcodon_exon_num - s]
                - ds.exon_start_list[endcodon_exon_num - s]
            )
        exon_seq_to_stopcodon = exon_seq[
            :-stopcodon_index
        ]  # get seq from exon_seq to stopcodon

    startcodon_index = ds.cdsStart - ds.exon_start_list[startcodon_exon_num]
    if startcodon_exon_num == 0:  # if startcodon is in the first exon
        cds_seq = exon_seq_to_stopcodon[startcodon_index:]
    else:  # if startcodon is not in the first exon. (add the len of exons before startcodon)
        for s in range(startcodon_exon_num):
            startcodon_index += ds.exon_end_list[s] - ds.exon_start_list[s]
        cds_seq = exon_seq_to_stopcodon[startcodon_index:]
    return cds_seq
