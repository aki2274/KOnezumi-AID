from __future__ import annotations
import re
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.nominate_candidate_stopcodon.generate_cds_seq import generate_cdsseq

#####
# check grna position
#####


def create_candidates_list_dict(ct_cand: list[str], ga_cand: list[str]) -> list[dict]:
    candidates_list = []
    for grna in ct_cand:
        candidates_list.append({"ct_seq": grna})
    for grna in ga_cand:
        candidates_list.append({"ga_seq": grna})
    return candidates_list


def label_in_start_150bp(cand_grna: list[dict], transcript_record: GeneData) -> list[dict]:
    candidate_ = cand_grna.copy()
    for grna in candidate_:
        if "ct_seq" in grna:
            #Add +1 because the position of PTC is +1~3 from the guide start.
            if re.search(grna["ct_seq"], transcript_record.orf_seq).start() + 1 <= transcript_record.cdsStart + 150:
                grna["in_start_150bp"] = True
            else:
                grna["in_start_150bp"] = False
        elif "ga_seq" in grna:
            #Add +16 for the same reason as above.
            if re.search(grna["ga_seq"], transcript_record.orf_seq).start() + 16 <= transcript_record.cdsStart + 150:
                grna["in_start_150bp"] = True
            else:
                grna["in_start_150bp"] = False
    return candidate_


# the case of single exon
def eliminate_in_front_half(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    # check if the PTC is not in the front half of the CDS
    removed = []
    cds = generate_cdsseq(ds)
    for grna in cand_grna:
        if "ct_seq" in grna:
            ptc_index = re.search(grna["ct_seq"], cds)
            if ptc_index is not None and ptc_index.start() + 3 < len(cds) / 2:
                removed.append(grna)
        elif "ga_seq" in grna:
            ptc_index = re.search(grna["ga_seq"], cds)
            if ptc_index is not None and ptc_index.start() + 19 < len(cds) / 2:
                removed.append(grna)
    return removed


# the case of multi exon
def eliminate_in_last_exon(cand_grna: list[str], ds: GeneData) -> list[str]:
    removed = []
    for grna in cand_grna:
        if "ct_seq" in grna:
            if (
                re.search(grna["ct_seq"], ds.orf_seq).start() + 1
                <= ds.exon_start_list[-1]
            ):
                removed.append(grna)
        elif "ga_seq" in grna:
            if (
                re.search(grna["ga_seq"], ds.orf_seq).start() + 16
                <= ds.exon_start_list[-1]
            ):
                removed.append(grna)
    return removed


def label_in_50bp_from_LEJ(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    candidate_ = cand_grna.copy()
    for grna in candidate_:
        if "ct_seq" in grna:
            if (
                re.search(grna["ct_seq"], ds.orf_seq).start() + 1
                >= ds.exon_end_list[-2] - 50
            ):
                grna["in_50bp_from_LEJ"] = True
            else:
                grna["in_50bp_from_LEJ"] = False
        elif "ga_seq" in grna:
            if (
                re.search(grna["ga_seq"], ds.orf_seq).start() + 16
                >= ds.exon_end_list[-2] - 50
            ):
                grna["in_50bp_from_LEJ"] = True
            else:
                grna["in_50bp_from_LEJ"] = False
    return candidate_
