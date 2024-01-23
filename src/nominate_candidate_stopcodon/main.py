from __future__ import annotations
from pathlib import Path
import pickle
from src.create_gene_dataclass import GeneData, create_dataclass
from src.nominate_candidate_stopcodon.find_candidate_seq import (
    find_ct_target_seq,
    find_ga_target_seq,
)
from src.nominate_candidate_stopcodon.transform_candidateseq_to_index import (
    transform_ct_guideseq_to_index,
    transform_ga_guideseq_to_index,
)
from src.nominate_candidate_stopcodon.generate_cds_seq import (
    get_startcodon_exon_num,
    generate_cdsseq,
)
from src.nominate_candidate_stopcodon.search_candidate_codon import (
    search_candidate_index,
)
from src.nominate_candidate_stopcodon.translate_index_in_exon_to_orf import (
    translate_index_in_exon_to_orf,
)
from src.nominate_candidate_stopcodon.make_grna_from_index import (
    convert_ct_grna,
    convert_ga_grna,
)


def nominate_candidate_stopcodon(
    transcript_name: str, refflat_path: Path, seq_path: Path
) -> list[dict]:
    """
    Export candidate gRNA.

    Args:
        transcript_name (str): transcript name.
    Returns:
        list[dict]: the candidate gRNA.

    Example:
        >>> transcript_name = "NM_001011874"
        >>> main(transcript_name)
        ['CCAGTATCAGCCCACTTGGTAGG', 'GACAACTATGTAAAGAGACTTGG'], ['CCACTTATTCTCAGATTTTGGGG', 'CCACTTCATTGATGCTACTGGTT']
    """
    # 0. get gene data
    refflat = pickle.load(open(refflat_path, "rb"))
    seq_dict = pickle.load(open(seq_path, "rb"))
    ds = create_dataclass(transcript_name, refflat, seq_dict)
    # 1. generate cds seq
    cds_seq = generate_cdsseq(ds)
    # 2. search candidate index
    candidate_ptc_index = search_candidate_index(cds_seq)
    # 3. get startcodon exon number
    cdsStart_exon_index = get_startcodon_exon_num(ds)
    # 4. translate candidate index in exon to index in genome
    candidate_ptc_index_in_orf = translate_index_in_exon_to_orf(
        ds, candidate_ptc_index, cdsStart_exon_index
    )
    # 5,find ct target seq
    ct_target_seq = find_ct_target_seq(ds.orf_seq)
    # 6. find ga target seq
    ga_target_seq = find_ga_target_seq(ds.orf_seq)
    # 7. transform ct guideseq to index
    ct_guideseq_index = transform_ct_guideseq_to_index(ds.orf_seq, ct_target_seq)
    # 8. transform ga guideseq to index
    ga_guideseq_index = transform_ga_guideseq_to_index(ds.orf_seq, ga_target_seq)
    # 9. compare candidate index and guideseq index,then export candidate if the same index
    ct_candidate_index = [
        index for index in candidate_ptc_index_in_orf if index in ct_guideseq_index
    ]
    ga_candidate_index = [
        index for index in candidate_ptc_index_in_orf if index in ga_guideseq_index
    ]
    # 10. transform candidate index to candidate gRNA
    ct_candidate_grna = convert_ct_grna(ds, ct_candidate_index)
    ga_candidate_grna = convert_ga_grna(ds, ga_candidate_index)

    return ct_candidate_grna, ga_candidate_grna
