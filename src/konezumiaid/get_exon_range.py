from __future__ import annotations
from konezumiaid.create_gene_dataclass import TranscriptRecord


def get_exon_range(transcript_record: TranscriptRecord):
    # get exon range
    exon_range = []
    sum_exon_length = 0
    for start, end in zip(transcript_record.exon_start_positions, transcript_record.exon_end_positions):
        exon_length = end - start
        exon_range.append([sum_exon_length, sum_exon_length + exon_length - 1])
        sum_exon_length += exon_length

    return exon_range
