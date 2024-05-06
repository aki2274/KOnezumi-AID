from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData


def generate_exon_seq(transcript_recod: GeneData) -> str:
    """Generate exon sequence (concatenate all exons' sequence in the transcript)"""
    exon_seq_list = [
        transcript_recod.orf_seq[exon_start:exon_end]
        for exon_start, exon_end in zip(
            transcript_recod.exon_start_list, transcript_recod.exon_end_list
        )
    ]
    return "".join(exon_seq_list)


def get_startcodon_exon_num(transcript_recod: GeneData) -> int:
    # search the exon number of startcodon
    for s in range(len(transcript_recod.exon_start_list)):
        if (
            transcript_recod.exon_start_list[s]
            <= transcript_recod.cdsStart
            <= transcript_recod.exon_end_list[s]
        ):
            exon_index = s
    return exon_index


def get_stopcodon_exon_num(transcript_recod: GeneData) -> int:
    # search the exon number of stopcodon
    for s in range(len(transcript_recod.exon_start_list)):
        if (
            transcript_recod.exon_start_list[s]
            <= transcript_recod.cdsEnd
            <= transcript_recod.exon_end_list[s]
        ):
            exon_index = s
    return exon_index


def generate_cdsseq(transcript_recod: GeneData) -> str:
    exon_seq = generate_exon_seq(transcript_recod)
    startcodon_exon_num = get_startcodon_exon_num(transcript_recod)
    endcodon_exon_num = get_stopcodon_exon_num(transcript_recod)

    stopcodon_index = (
        transcript_recod.exon_end_list[endcodon_exon_num] - transcript_recod.cdsEnd
    )
    if stopcodon_index == 0:  # if stopcodon is in the last exon
        exon_seq_to_stopcodon = exon_seq
    elif (
        endcodon_exon_num - transcript_recod.exonCount + 1 == 0
    ):  # if stopcodon is in the last exon
        exon_seq_to_stopcodon = exon_seq[:-stopcodon_index]
    else:  # if stopcodon is not in the last exon. (add the len of exons after stopcodon)
        for s in range(1, transcript_recod.exonCount - endcodon_exon_num):
            stopcodon_index += (
                transcript_recod.exon_end_list[endcodon_exon_num + s]
                - transcript_recod.exon_start_list[endcodon_exon_num + s]
            )
        exon_seq_to_stopcodon = exon_seq[
            :-stopcodon_index
        ]  # get seq from exon_seq start to stopcodon
    startcodon_index = (
        transcript_recod.cdsStart
        - transcript_recod.exon_start_list[startcodon_exon_num]
    )
    if startcodon_exon_num == 0:  # if startcodon is in the first exon
        cds_seq = exon_seq_to_stopcodon[startcodon_index:]
    else:  # if startcodon is not in the first exon. (add the len of exons before startcodon)
        for s in range(startcodon_exon_num):
            startcodon_index += (
                transcript_recod.exon_end_list[s] - transcript_recod.exon_start_list[s]
            )
        cds_seq = exon_seq_to_stopcodon[startcodon_index:]
    return cds_seq
