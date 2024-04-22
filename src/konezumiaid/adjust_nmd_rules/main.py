from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.adjust_nmd_rules.rete_position import (
    create_candidates_list_dict,
    label_in_start_150bp,
    eliminate_in_front_half,
    eliminate_in_last_exon,
    label_in_50bp_from_LEJ,
)


def adjust_nmd_rules(
    ds: GeneData, ct_cand: list[str], ga_cand: list[str]
) -> list[dict]:
    """
    Filter gRNA candidates based on NMD rules.

    Args:
        ds (GeneData): the DataClass named GeneData .
        ct_cand (list[str]): candidate C to T conversion gRNA.
        ga_cand (list[str]): candidate G to A conversion gRNA.
    Returns:
        list[dict]: list of candidate gRMA that filtered NMD rules.
    """
    # 1. create gRNA list
    gRNA_list = create_candidates_list_dict(ct_cand, ga_cand)

    # 2. label in start 150bp
    gRNA_list = label_in_start_150bp(gRNA_list, ds)
    # check exon count
    if ds.exonCount == 1:# The case of single exon 
        # 3. label in front half
        gRNA_list = eliminate_in_front_half(gRNA_list, ds)
        # 4. combine ct_cand and ga_cand to create cand_seq
        candidates = []
        for d in gRNA_list:
            tmp_dict = {}
            for k, v in d.items():
                if k == "ct_seq" or k == "ga_seq":
                    tmp_dict["candidate"] = v
                else:
                    tmp_dict[k] = v
            candidates.append(tmp_dict)
        return candidates
    else: # The case of multi exon
        # 3. eliminate in last exon
        gRNA_list = eliminate_in_last_exon(gRNA_list, ds)
        # 4. label in 50bp from LEJ
        gRNA_list = label_in_50bp_from_LEJ(gRNA_list, ds)
        # 5. combine ct_cand and ga_cand to create cand_seq
        candidates = []
        for d in gRNA_list:
            tmp_dict = {}
            for k, v in d.items():
                if k == "ct_seq" or k == "ga_seq":
                    tmp_dict["candidate"] = v
                else:
                    tmp_dict[k] = v
            candidates.append(tmp_dict)
        return candidates
