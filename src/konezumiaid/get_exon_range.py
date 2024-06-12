from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData


def get_exon_range(transcript_record: GeneData):
    # get exon range
    exon_range = []
    sum_exon_length = 0
    for start, end in zip(
        transcript_record.exon_start_list, transcript_record.exon_end_list
    ):
        exon_length = end - start
        exon_range.append([sum_exon_length, sum_exon_length + exon_length - 1])
        sum_exon_length += exon_length

    return exon_range
