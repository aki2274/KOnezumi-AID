from __future__ import annotations
from src.create_gene_dataclass import GeneData


def get_exon_range(ds: GeneData):
    # get exon range
    exon_range_list = []
    num = 0
    for s in range(len(ds.exon_start_list)):
        num += ds.exon_end_list[s] - ds.exon_start_list[s]
        exon_range_list.append(
            [num - (ds.exon_end_list[s] - ds.exon_start_list[s]), num - 1]
        )

    return exon_range_list
