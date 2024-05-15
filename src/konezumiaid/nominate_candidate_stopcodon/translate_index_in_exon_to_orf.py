from __future__ import annotations
import numpy as np
from konezumiaid.get_exon_range import get_exon_range
from konezumiaid.create_gene_dataclass import GeneData


#####
# Change candidate stopcodon index in exon to index in genome
#####
# 　エクソンごとの範囲を取得、candidate stopcodon が何番目のエクソンか知る、


# the candidate indexes in cds are convereted to indexes in exon
def translate_incds_index_to_exon(
    ds: GeneData,
    candidate_stopcodon_cds_index_list: list[int],
    cdsStart_exon_index: int,
) -> list[int]:
    # get the exon range
    exon_range_list = get_exon_range(ds)
    # the distance from the start of orf to the exon which has startcodon
    dist_to_startcodon_exon = exon_range_list[cdsStart_exon_index][0]
    # the distance from the start of the exon(has startcodn) to the startcodon(cdsStart)
    dist_to_startcodon = (
        ds.cdsStart
        - int(ds.exon_start_list[cdsStart_exon_index])
        + dist_to_startcodon_exon
    )
    # calculate candidate codon index in exon
    candidate_codon_index_inexon = [
        index + dist_to_startcodon for index in candidate_stopcodon_cds_index_list
    ]
    return candidate_codon_index_inexon


def get_exonnum_of_candidate(
    candidate_codon_index_inexon: list[int], exon_range_list: list[tuple[int, int]]
) -> list[int]:
    # create a list of exon number of candidate codon is in
    exon_index_list = [
        s
        for t in candidate_codon_index_inexon
        for s, (start, end) in enumerate(exon_range_list)
        if start <= t <= end
    ]
    return exon_index_list


def translate_index_in_exon_to_orf(
    ds: GeneData,
    candidate_stopcodon_cds_index_list: list[int],
    cdsStart_exon_index: int,
) -> int:
    # get the sum of intron length from the orf start to the exon which has candidate codon
    exon_range_list = get_exon_range(ds)
    candidate_codon_index_inexon_list = translate_incds_index_to_exon(
        ds, candidate_stopcodon_cds_index_list, cdsStart_exon_index
    )
    exon_index_list = get_exonnum_of_candidate(
        candidate_codon_index_inexon_list, exon_range_list
    )

    add_num_to_genome_index = [
        (ds.exon_start_list[i] - ds.txStart - exon_range_list[i][0])
        for i in exon_index_list
    ]
    # add intron length abd (exon length + the length from exon which has candidate codon to candidate) to get the index of candidate codon in genome
    candidate_stopcodon_index = np.array(candidate_codon_index_inexon_list) + np.array(
        add_num_to_genome_index
    )
    return candidate_stopcodon_index.tolist()
