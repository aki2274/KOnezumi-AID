from __future__ import annotations
from src.create_gene_dataclass import GeneData


def search_candidate(ds: GeneData):
    acceptor_cand_rna = []
    donor_cand_rna = []
    for s in ds.exon_start_list:
        search_seq = ds.orf_seq[s - 22 : s + 3]
        # if "GG" in 0~5 position, it is acceptor candidate
        if "GG" in search_seq[:4]:
            acceptor_cand_rna.append(search_seq)

    for e in ds.exon_end_list[:-1]:
        search_seq = ds.orf_seq[e - 22 : e + 3]
        if "GG" in search_seq[:4]:
            donor_cand_rna.append(search_seq)
    return acceptor_cand_rna, donor_cand_rna
