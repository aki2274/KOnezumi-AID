from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.nominate_ptc_guide.find_candidate_seq import (
    search_c_to_t_guide_seq,
    search_g_to_a_guide_seq,
)
from konezumiaid.nominate_ptc_guide.transform_candidateseq_to_position import (
    transform_c_to_t_guideseq_to_position,
    transform_g_to_a_guideseq_to_position,
)
from konezumiaid.nominate_ptc_guide.generate_cds_seq import (
    get_startcodon_exon_index,
    generate_cdsseq,
)

from konezumiaid.nominate_ptc_guide.search_candidate_codon import (
    search_candidate_codon_position,
)

from konezumiaid.nominate_ptc_guide.translate_codon_position_to_inorf import (
    translate_position_in_splicedexon_to_orf,
)

from konezumiaid.nominate_ptc_guide.make_grna_from_position import (
    extract_c_to_t_grna_from_position,
    extract_g_to_a_grna_from_position,
)


def nominate_candidate_stopcodon(
    transcript_recod: GeneData,
) -> list[dict]:
    """
    Export candidate gRNA, that can induce PTC by using Target-AID.

    Args:
        transcript_recod (GeneData): A GeneData object representing the transcript record.
    Returns:
        list[dict]: the candidate gRNA.
    """
    # 1. generate cds seq
    cds_seq = generate_cdsseq(transcript_recod)
    # 2. search candidate index
    candidate_ptc_index = search_candidate_codon_position(cds_seq)
    # 3. get startcodon exon number
    exon_index = get_startcodon_exon_index(transcript_recod)
    # 4. translate candidate index in exon to index in genome
    candidate_ptc_codon_positions = translate_position_in_splicedexon_to_orf(
        transcript_recod, candidate_ptc_index, exon_index
    )
    # 5,find ct target seq
    all_ct_guide_seq = search_c_to_t_guide_seq(transcript_recod.orf_seq)
    # 6. find ga target seq
    all_ga_guide_seq = search_g_to_a_guide_seq(transcript_recod.orf_seq)
    # 7. transform ct guideseq to index
    editable_c_positions = transform_c_to_t_guideseq_to_position(
        transcript_recod.orf_seq, all_ct_guide_seq
    )
    # 8. transform ga guideseq to index
    editable_codon_T_positions = transform_g_to_a_guideseq_to_position(
        transcript_recod.orf_seq, all_ga_guide_seq
    )
    # 9. compare candidate index and guideseq index,then export candidate if the same index
    ct_candidate_positions = [
        index
        for index in candidate_ptc_codon_positions
        if index in editable_c_positions
    ]
    ga_candidate_positions = [
        index
        for index in candidate_ptc_codon_positions
        if index in editable_codon_T_positions
    ]
    # 10. transform candidate index to candidate gRNA
    ct_candidate_grna = extract_c_to_t_grna_from_position(
        transcript_recod, ct_candidate_positions
    )
    ga_candidate_grna = extract_g_to_a_grna_from_position(
        transcript_recod, ga_candidate_positions
    )

    return ct_candidate_grna, ga_candidate_grna
