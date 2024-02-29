from __future__ import annotations
from src.create_gene_dataclass import Gene
from src.adjust_nmd_rules.create_grna_lsit import create_list
from src.adjust_nmd_rules.rete_position import eliminate_in_last_exon


def adjust_nmd_rules(gene: Gene) -> Gene:
    gene = create_list(gene)
    gene = eliminate_in_last_exon(gene)
    return gene
