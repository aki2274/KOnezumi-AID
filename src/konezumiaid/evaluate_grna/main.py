from __future__ import annotations
from konezumiaid.create_gene_dataclass import TranscriptRecord
from konezumiaid.get_reverse_complement import get_revcomp
from konezumiaid.evaluate_grna.apply_nmd_rules import (
    create_candidates_list_dict,
    label_in_start_150bp,
    label_in_more_than_400bp_exon,
    eliminate_in_back_half,
    eliminate_in_last_exon,
    label_in_50bp_from_LEJ,
)
from konezumiaid.evaluate_grna.add_grna_info import link_to_crisperdirect


def apply_nmd_rules(transcript_record: TranscriptRecord, ct_cand: list[str], ga_cand: list[str]) -> list[dict]:
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

    

    if transcript_record.exon_count == 1:  # The case of single exon
        # 2. label in front half
        gRNA_list = eliminate_in_back_half(gRNA_list, transcript_record)

    else:  # The case of multi exon
        # 2. label in start 150bp
        gRNA_list = label_in_start_150bp(gRNA_list, transcript_record)
        # 3. eliminate in last exon
        gRNA_list = eliminate_in_last_exon(gRNA_list, transcript_record)
        # 4. label in 50bp from LEJ
        gRNA_list = label_in_50bp_from_LEJ(gRNA_list, transcript_record)
        # 5. label in more than 400bp exon
        gRNA_list = label_in_more_than_400bp_exon(gRNA_list, transcript_record)

    # combine ct_cand and ga_cand to create cand_seq
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

    # link to CRISPRdirect
    candidates = link_to_crisperdirect(candidates)
    return candidates
