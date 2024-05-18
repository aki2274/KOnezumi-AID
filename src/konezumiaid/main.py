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
    adjusted_gRNA_df: pd.DataFrame,
    acceptor_cand_df: pd.DataFrame,
    donor_cand_df: pd.DataFrame,
) -> None:
    print("PTC gRNA")
    print(adjusted_gRNA_df.to_string())
    print("Acceptor gRNA")
    print(acceptor_cand_df.to_string())
    print("Donor gRNA")
    print(donor_cand_df.to_string())


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
    refflat = pickle.load(open(refflat_path, "rb"))
    seq_dict = pickle.load(open(seq_path, "rb"))
    if name.startswith("NM_"):
        transcript_record = create_dataclass(name, refflat, seq_dict)
        ptc_cand, acceptor_cand, donor_cand = konezumiaid_main(transcript_record)
        ptcp_cand_df = pd.DataFrame(ptc_cand)
        acceptor_cand_df = pd.DataFrame(acceptor_cand)
        donor_cand_df = pd.DataFrame(donor_cand)
        show_table(ptcp_cand_df, acceptor_cand_df, donor_cand_df)

    else:
        transcript_names = [d["name"] for d in refflat if d["geneName"] == name]
        symbol_cand = {}
        for transcript in transcript_names:
            transcript_record = create_dataclass(transcript, refflat, seq_dict)
            ptc_cand, acceptor_cand, donor_cand = konezumiaid_main(transcript_record)
            ptcp_cand_df = pd.DataFrame(ptc_cand)
            acceptor_cand_df = pd.DataFrame(acceptor_cand)
            donor_cand_df = pd.DataFrame(donor_cand)
            symbol_cand[transcript] = {
                "ptc_cand": ptcp_cand_df,
                "acceptor_cand": acceptor_cand_df,
                "donor_cand": donor_cand_df,
            }
        for k, v in symbol_cand.items():
            print(k)
            show_table(
                v["ptc_cand"],
                v["acceptor_cand"],
                v["donor_cand"],
            )


def main():
    if len(sys.argv) != 2:
        raise ValueError("Please provide a gene name as an argument.")
    gene_name = sys.argv[1]
    excecute(gene_name)
