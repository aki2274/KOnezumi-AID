from __future__ import annotations
from konezumiaid.create_gene_dataclass import TranscriptRecord


def generate_exon_seq(transcript_record: TranscriptRecord) -> str:
    """Generate exon sequence (concatenate all exons' sequence in the transcript)"""
    exon_seq_list = [
        transcript_record.transcript_seq[exon_start:exon_end]
        for exon_start, exon_end in zip(transcript_record.exon_start_positions, transcript_record.exon_end_positions)
    ]
    return "".join(exon_seq_list)


def get_startcodon_exon_index(transcript_record: TranscriptRecord) -> int:
    """Search the exon index of the exon has startcodon"""
    startcodon_position = transcript_record.cds_start
    exon_index = None
    for n, (start, end) in enumerate(zip(transcript_record.exon_start_positions, transcript_record.exon_end_positions)):
        if start <= startcodon_position <= end:
            exon_index = n
            break
    return exon_index


def get_stopcodon_exon_index(transcript_record: TranscriptRecord) -> int:
    """Search the exon index of the exon has stopcodon"""
    stopcodon_position = transcript_record.cds_end
    exon_index = None
    for n, (start, end) in enumerate(zip(transcript_record.exon_start_positions, transcript_record.exon_end_positions)):
        if start <= stopcodon_position <= end:
            exon_index = n
            break
    return exon_index


def generate_cdsseq(transcript_record: TranscriptRecord) -> str:
    """Generate CDS sequence(CDS sequence is the sequence from startcodon to stopcodon)"""
    exon_seq = generate_exon_seq(transcript_record)
    startcodon_index = get_startcodon_exon_index(transcript_record)
    stopcodon_index = get_stopcodon_exon_index(transcript_record)

    stopcodon_distance_from_exonend = transcript_record.exon_end_positions[stopcodon_index] - transcript_record.cds_end
    if (
        stopcodon_distance_from_exonend == 0 and stopcodon_index - transcript_record.exon_count + 1 == 0
    ):  # the case of not having 3'UTR
        seq_to_stopcodon = exon_seq
    else:  # the case of having 3'UTR
        if stopcodon_index - transcript_record.exon_count + 1 == 0:  # the case of stopcodon is in the last exon
            seq_to_stopcodon = exon_seq[:-stopcodon_distance_from_exonend]
        else:  # if stopcodon is not in the last exon. (add the len of exons after the exon that has stopcodon)
            for s in range(1, transcript_record.exon_count - stopcodon_index):
                # add the len of exons after the exon that has stopcodon
                stopcodon_distance_from_last_exonend = stopcodon_distance_from_exonend
                stopcodon_distance_from_last_exonend += (
                    transcript_record.exon_end_positions[stopcodon_index + s]
                    - transcript_record.exon_start_positions[stopcodon_index + s]
                )
            seq_to_stopcodon = exon_seq[:-stopcodon_distance_from_last_exonend]

    startcodon_distance_from_exonstart = (
        transcript_record.cds_start - transcript_record.exon_start_positions[startcodon_index]
    )
    if startcodon_index == 0:  # the case of startcodon is in the first exon
        cds_seq = seq_to_stopcodon[startcodon_distance_from_exonstart:]
    else:  # if startcodon is not in the first exon. (add the len of exons before startcodon)
        startcodon_distance_from_first_exonstart = startcodon_distance_from_exonstart
        for s in range(startcodon_index):
            startcodon_distance_from_first_exonstart += (
                transcript_record.exon_end_positions[s] - transcript_record.exon_start_positions[s]
            )
        cds_seq = seq_to_stopcodon[startcodon_distance_from_first_exonstart:]
    return cds_seq
