from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.get_reverse_complement import get_revcomp
from konezumiaid.apply_nmd_rules.evaluate_grna import (
    create_candidates_list_dict,
    label_in_start_150bp,
    eliminate_in_back_half,
    eliminate_in_last_exon,
    label_in_50bp_from_LEJ,
)


def apply_nmd_rules(
    transcript_record: GeneData, ct_cand: list[str], ga_cand: list[str]
) -> list[dict]:
    """
    Apply and select gRNA candidates according to NMD rules.

    This function takes a GeneData object representing a transcript record, a list of candidate C to T conversion gRNA (ct_cand), and a list of candidate G to A conversion gRNA (ga_cand). It applies a series of rules to filter the gRNA candidates based on NMD (Nonsense-Mediated Decay) rules.

    Args:
        transcript_record (GeneData): A GeneData object representing the transcript record.
        ct_cand (list[str]): A list of candidate C to T conversion gRNA.
        ga_cand (list[str]): A list of candidate G to A conversion gRNA.

    Returns:
        list[dict]: A list of candidate gRNA that have been filtered based on NMD rules.
    """
    # 1. create gRNA list
    gRNA_list = create_candidates_list_dict(ct_cand, ga_cand)

    # 2. label in start 150bp
    gRNA_list = label_in_start_150bp(gRNA_list, transcript_record)

    if transcript_record.exonCount == 1:  # The case of single exon
        # 3. label in front half
        gRNA_list = eliminate_in_back_half(gRNA_list, transcript_record)

    else:  # The case of multi exon
        # 3. eliminate in last exon
        gRNA_list = eliminate_in_last_exon(gRNA_list, transcript_record)
        # 4. label in 50bp from LEJ
        gRNA_list = label_in_50bp_from_LEJ(gRNA_list, transcript_record)

    # 4 or 5. combine ct_cand and ga_cand to create cand_seq
    candidates = []
    for d in gRNA_list:
        tmp_dict = {}
        for k, v in d.items():
            if k == "ct_seq":
                tmp_dict["seq"] = v
            elif k == "ga_seq":
                tmp_dict["seq"] = get_revcomp(v)
            else:
                tmp_dict[k] = v
        candidates.append(tmp_dict)
    return candidates
