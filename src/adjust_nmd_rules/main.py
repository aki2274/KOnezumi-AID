from __future__ import annotations
from src.create_gene_dataclass import GeneData
from src.adjust_nmd_rules.create_grna_lsit import create_list
from src.adjust_nmd_rules.rete_position import (
    label_in_start_150bp,
    label_in_front_half,
    eliminate_in_last_exon,
    label_in_50bp_from_LEJ,
)


def adjust_nmd_rules(
    ds: GeneData, ct_cand: list[str], ga_cand: list[str]
) -> list[dict]:
    """
    Adjust NMD rules.

    Args:
        ds (GeneData): GeneData class.
        ct_cand (list[str]): candidate ct gRNA.
        ga_cand (list[str]): candidate ga gRNA.
    Returns:
        list[dict]: adjusted gRNA.
    """
    # 1. create gRNA list
    gRNA_list = create_list(ct_cand, ga_cand)

    # 2. label in start 150bp
    gRNA_list = label_in_start_150bp(gRNA_list, ds)
    # check exon count
    if ds.exonCount == 1:
        # 3. label in front half
        gRNA_list = label_in_front_half(gRNA_list, ds)
        return gRNA_list
    else:
        # 4. eliminate in last exon
        gRNA_list = eliminate_in_last_exon(gRNA_list, ds)
        # 5. label in 50bp from LEJ
        gRNA_list = label_in_50bp_from_LEJ(gRNA_list, ds)
        return gRNA_list
