from __future__ import annotations
from src.create_gene_dataclass import GeneData


def generate_exon_seq(ds: GeneData) -> str:
    # make seq spliced by exons (NOT cds)
    exon_seq_list = [
        ds.orf_seq[start - ds.txStart : end - ds.txStart]
        for start, end in zip(ds.exon_start_list, ds.exon_end_list)
    ]
    return "".join(exon_seq_list)


def get_startcodon_exon_num(ds: GeneData) -> int:
    # search the exon number of startcodon
    for s in range(len(ds.exon_start_list)):
        if ds.exon_start_list[s] <= ds.cdsStart <= ds.exon_end_list[s]:
            exon_index = s
    return exon_index


def get_stopcodon_exon_num(ds: GeneData) -> int:
    # search the exon number of stopcodon
    for s in range(len(ds.exon_start_list)):
        if ds.exon_start_list[s] <= ds.cdsEnd <= ds.exon_end_list[s]:
            exon_index = s
    return exon_index


def generate_cdsseq(ds: GeneData) -> str:
    exon_seq = generate_exon_seq(ds)
    startcodon_exon_num = get_startcodon_exon_num(ds)
    endcodon_exon_num = get_stopcodon_exon_num(ds)

    stopcodon_index = ds.exon_end_list[endcodon_exon_num] - ds.cdsEnd
    if endcodon_exon_num - ds.exonCount + 1 == 0:  # if stopcodon is in the last exon
        exon_seq_to_stopcodon = exon_seq[:-stopcodon_index]
    else:  # if stopcodon is not in the last exon. (add the len of exons after stopcodon)
        for s in range(1, ds.exonCount - endcodon_exon_num):
            stopcodon_index += (
                ds.exon_end_list[endcodon_exon_num + s]
                - ds.exon_start_list[endcodon_exon_num + s]
            )
        exon_seq_to_stopcodon = exon_seq[
            :-stopcodon_index
        ]  # get seq from exon_seq start to stopcodon

    startcodon_index = ds.cdsStart - ds.exon_start_list[startcodon_exon_num]
    if startcodon_exon_num == 0:  # if startcodon is in the first exon
        cds_seq = exon_seq_to_stopcodon[startcodon_index:]
    else:  # if startcodon is not in the first exon. (add the len of exons before startcodon)
        for s in range(startcodon_exon_num):
            startcodon_index += ds.exon_end_list[s] - ds.exon_start_list[s]
        cds_seq = exon_seq_to_stopcodon[startcodon_index:]
    return cds_seq
