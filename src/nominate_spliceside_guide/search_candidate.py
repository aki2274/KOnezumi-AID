from __future__ import annotations
from src.create_gene_dataclass import GeneData
from typing import List, Dict


def search_candidate(ds: GeneData) -> List[Dict[str, int]]:
    acceptor_candidates = [
        {"seq": ds.orf_seq[s - 22: s + 3][ds.orf_seq[s - 22: s + 3].index("GG"): ds.orf_seq[s - 22: s + 3].index("GG") + 23],
         "exon_num": i + 2}
        for i, s in enumerate(ds.exon_start_list[1:])
        if "GG" in ds.orf_seq[s - 22: s + 3][:4]
    ]

    donor_candidates = [
        {"seq": ds.orf_seq[e - 21: e + 4][ds.orf_seq[e - 21: e + 4].index("GG"): ds.orf_seq[e - 21: e + 4].index("GG") + 23],
         "exon_num": i + 1}
        for i, e in enumerate(ds.exon_end_list[:-1])
        if "GG" in ds.orf_seq[e - 21: e + 4][:4]
    ]

    return acceptor_candidates, donor_candidates
