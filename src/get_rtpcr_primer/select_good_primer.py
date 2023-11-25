from __future__ import annotations
from src.get_rtpcr_primer.check_primer_quality import (
    check_exon_junction,
    rewrite_primer_pair_intron_len,
)
from src.get_rtpcr_primer.make_rtpcr_primer import (
    export_candidate_primer,
    get_candidate_primer_pairs,
)
from src.get_rtpcr_primer.label_primers_quality import (
    label_primers_quality,
)

#####
# Need primer uniqueness check
#####


def get_good_primer(
    exon_seq: str,
    exon_range: list[list[int, int]],
) -> list[dict]:
    # get good primer pairs
    candidate = export_candidate_primer(exon_seq)
    candidate_primer_info = get_candidate_primer_pairs(exon_seq, candidate, exon_range)
    # check the primer quality
    exon_junction_bool = check_exon_junction(candidate_primer_info, exon_range)
    candidate_primer_info = rewrite_primer_pair_intron_len(
        candidate_primer_info, exon_range
    )
    # label the primer quality and return
    return label_primers_quality(candidate_primer_info, exon_junction_bool)
