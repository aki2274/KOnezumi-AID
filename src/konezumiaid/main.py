from __future__ import annotations
import pandas as pd
import pickle
import sys
from pathlib import Path
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.create_gene_dataclass import create_dataclass
from konezumiaid.nominate_ptc_guide.main import nominate_candidate_stopcodon
from konezumiaid.nominate_splicesite_guide.search_candidate import search_site_candidate
from konezumiaid.apply_nmd_rules.main import apply_nmd_rules

# from konezumiaid.get_rtpcr_primer.main import export_primers


def show_table(
    df_ptc_gRNA: pd.DataFrame,
    df_acceptor_cand: pd.DataFrame,
    df_donor_cand: pd.DataFrame,
) -> None:
    print("PTC gRNA")
    print(df_ptc_gRNA.to_string())
    print("Acceptor gRNA")
    print(df_acceptor_cand.to_string())
    print("Donor gRNA")
    print(df_donor_cand.to_string())


def export_csv(
    name: str,
    df_ptc_gRNA: pd.DataFrame,
    df_acceptor_cand: pd.DataFrame,
    df_donor_cand: pd.DataFrame,
) -> None:
    df_ptc_gRNA.to_csv(f"{name}_ptc_gRNA.csv", index=False)
    df_acceptor_cand.to_csv(f"{name}_acceptor_cand.csv", index=False)
    df_donor_cand.to_csv(f"{name}_donor_cand.csv", index=False)


def konezumiaid_main(
    transcript_record: GeneData,
) -> tuple[list[dict], list[dict], list[dict]]:
    ct_cand, ga_cand = nominate_candidate_stopcodon(transcript_record)
    applied_nmd_rules_gRNA = apply_nmd_rules(ct_cand, ga_cand)
    acceptor_cand, donor_cand = search_site_candidate(transcript_record)
    return applied_nmd_rules_gRNA, acceptor_cand, donor_cand


def excecute(name: str) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
    seq_path = Path("data", "sorted_seq_dict.pkl")
    df_refflat = pickle.load(open(refflat_path, "rb"))
    seq_dict = pickle.load(open(seq_path, "rb"))
    if name.startswith("NM_"):
        transcript_record = create_dataclass(name, df_refflat, seq_dict)
        ptc_cand, acceptor_cand, donor_cand = konezumiaid_main(transcript_record)
        df_ptcp_cand = pd.DataFrame(ptc_cand)
        df_acceptor_cand = pd.DataFrame(acceptor_cand)
        df_donor_cand = pd.DataFrame(donor_cand)
        show_table(df_ptcp_cand, df_acceptor_cand, df_donor_cand)
        export_csv(name, df_ptcp_cand, df_acceptor_cand, df_donor_cand)


def main():
    if len(sys.argv) != 2:
        raise ValueError("Please provide a gene name as an argument.")
    gene_name = sys.argv[1]
    excecute(gene_name)
