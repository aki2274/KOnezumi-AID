from __future__ import annotations
import numpy as np
from konezumiaid.get_exon_range import get_exon_range
from konezumiaid.create_gene_dataclass import GeneData


def translate_cds_position_to_exon(
    transcript_record: GeneData,
    candidate_codon_position_in_cds: list[int],
    startcodon_exon_index: int,
    exon_range_list: list[tuple[int, int]],  # list of tuple of exon start and end
) -> list[int]:
    # the distance from the exon start which has startcodon(on splicedexon)
    dist_to_exon_has_startcodon_splicedexon = exon_range_list[startcodon_exon_index][0]
    # the distance from the exon start which has startcodon(on orf)
    dist_to_exon_has_startcodon_from_txstart = int(
        transcript_record.exon_start_list[startcodon_exon_index]
    )
    dist_5UTR = (
        transcript_record.cdsStart
        - dist_to_exon_has_startcodon_from_txstart
        + dist_to_exon_has_startcodon_splicedexon
    )
    # add 5'UTR length to the candidate codon position in splicedexon
    candidate_codon_position_in_splicedexon = [
        posistion + dist_5UTR for posistion in candidate_codon_position_in_cds
    ]
    return candidate_codon_position_in_splicedexon


def get_exonindex_in_cand_codon(
    candidate_codon_position_in_splicedexon: list[int],
    exon_range_list: list[tuple[int, int]],
) -> list[int]:
    # create a list of exon number of candidate codon is in
    exon_indexes = [
        i
        for position in candidate_codon_position_in_splicedexon
        for i, (start, end) in enumerate(exon_range_list)
        if start <= position <= end
    ]
    return exon_indexes


def translate_position_in_splicedexon_to_orf(
    transcript_record: GeneData,
    candidate_codon_positions_in_cds: list[int],
    startcodon_exon_index: int,
) -> list[int]:
    # get the sum of intron length from the orf start to the exon which has candidate codon
    exon_range_list = get_exon_range(transcript_record)
    candidate_codon_positions_in_splicedexon = translate_cds_position_to_exon(
        transcript_record,
        candidate_codon_positions_in_cds,
        startcodon_exon_index,
        exon_range_list,
    )
    exon_indexes_in_ptc = get_exonindex_in_cand_codon(
        candidate_codon_positions_in_splicedexon, exon_range_list
    )

    intron_len_to_exon_has_ptc = [
        (transcript_record.exon_start_list[i] - exon_range_list[i][0])
        for i in exon_indexes_in_ptc
    ]
    # add intron length (exon length + the length from exon which has candidate codon to it)
    candidate_ptc_codon_positions = np.array(
        candidate_codon_positions_in_splicedexon
    ) + np.array(intron_len_to_exon_has_ptc)

    return candidate_ptc_codon_positions.tolist()
