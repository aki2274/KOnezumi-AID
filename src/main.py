from __future__ import annotations
import pandas as pd
from src.create_gene_dataclass import GeneData
from src.nominate_candidate_stopcodon.main import nominate_candidate_stopcodon
from src.nominate_spliceside_guide.search_candidate import search_splice_candidate
from src.adjust_nmd_rules.main import adjust_nmd_rules
from src.get_rtpcr_primer.main import export_primers

ds = None


def show_table(
    adj_gRNA: list[dict], acand: list[dict], dcand: list[dict], cprimer: list[dict]
) -> None:
    adjusted_gRNA_df = pd.DataFrame(adj_gRNA)
    acceptor_cand_df = pd.DataFrame(acand)
    donor_cand_df = pd.DataFrame(dcand)
    candidate_primers_df = pd.DataFrame(cprimer)
    print(adjusted_gRNA_df.to_string())
    print(acceptor_cand_df.to_string())
    print(donor_cand_df.to_string())
    print(candidate_primers_df.to_string())


def konezumi_main() -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    ct_acand, ga_acand = nominate_candidate_stopcodon(ds)
    adjusted_gRNA = adjust_nmd_rules(ds, ct_acand, ga_acand)
    acceptor_cand, donor_cand = search_splice_candidate(ds)
    candidate_primers = export_primers(ds)


def excecute_main(name: str) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    return None
