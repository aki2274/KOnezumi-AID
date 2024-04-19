from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData


def search_candidate(ds: GeneData) -> tuple[list[dict[int, str], list[dict[int, str]]]]:
    acceptor_cands = [
    {
        "seq": ds.orf_seq[s - 22 : s + 3][cc_idx : cc_idx + 23],
        "exon_num": i + 2
    }
    for i,s in enumerate(ds.exon_start_list[1:])
    if "CC" in ds.orf_seq[s - 22 : s + 3][:4] and "AG" in ds.orf_seq[s - 2 : s]
    for cc_idx in [idx for idx in range(3) if ds.orf_seq[s - 22 : s + 3][idx : idx + 2] == "CC"]
]


    donor_cands = [
        {
            "seq": ds.orf_seq[e - 21 : e + 4][cc_idx : cc_idx + 23],
            "exon_num": i + 1,
        }
        for i, e in enumerate(ds.exon_end_list[:-1])
        if "CC" in ds.orf_seq[e - 21 : e + 4][:4] and "GT" in ds.orf_seq[e : e + 2]
        for cc_idx in [idx for idx in range(3) if ds.orf_seq[e - 21 : e + 4][idx : idx + 2] == "CC"]
    ]

    acceptor_candidates = [
        candidate for candidate in acceptor_cands if "TTTT" not in candidate["seq"]
    ]
    donor_candidates = [
        candidate for candidate in donor_cands if "TTTT" not in candidate["seq"]
    ]

    return acceptor_candidates, donor_candidates
