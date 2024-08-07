from __future__ import annotations
import re
from konezumiaid.create_gene_dataclass import TranscriptRecord
from konezumiaid.nominate_ptc_guide.generate_cds_seq import generate_cdsseq


def create_candidates_list_dict(ct_cand: list[dict], ga_cand: list[dict]) -> list[dict]:
    candidates_list = []
    for cand in ct_cand:
        candidates_list.append({"aminoacid": cand["aminoacid"], "ct_seq": cand["seq"]})
    for cand in ga_cand:
        candidates_list.append({"aminoacid": cand["aminoacid"], "ga_seq": cand["seq"]})
    return candidates_list


def label_in_start_150bp(candidate: list[dict], transcript_record: TranscriptRecord) -> list[dict]:
    candidate_ = candidate.copy()
    for grna in candidate_:
        if "ct_seq" in grna:
            # Add +1 because the position of PTC is +1~3 from the guide start.
            if (
                re.search(grna["ct_seq"], transcript_record.transcript_seq).start() + 1
                <= transcript_record.cds_start + 150
            ):
                grna["in_start_150bp"] = True
            else:
                grna["in_start_150bp"] = False
        elif "ga_seq" in grna:
            # Add +16 for the same reason as above.
            if (
                re.search(grna["ga_seq"], transcript_record.transcript_seq).start() + 16
                <= transcript_record.cds_start + 150
            ):
                grna["in_start_150bp"] = True
            else:
                grna["in_start_150bp"] = False
    return candidate_


def label_in_more_than_400bp_exon(candidate: list[dict], transcript_record: TranscriptRecord) -> list[dict]:
    candidate_ = candidate.copy()
    exon_lengths = [
        end - start for start, end in zip(transcript_record.exon_start_positions, transcript_record.exon_end_positions)
    ]

    for grna in candidate_:
        if "ct_seq" in grna:
            edit_position = re.search(grna["ct_seq"], transcript_record.transcript_seq).start() + 1
        elif "ga_seq" in grna:
            edit_position = re.search(grna["ga_seq"], transcript_record.transcript_seq).start() + 16
        else:
            continue

        exon_num = next(
            (i for i, exon_end in enumerate(transcript_record.exon_end_positions) if edit_position <= exon_end), None
        )

        if exon_num is not None and exon_lengths[exon_num] > 400:
            grna["in_more_than_400bp_exon"] = True
        else:
            grna["in_more_than_400bp_exon"] = False

    return candidate_


# the case of single exon
def eliminate_in_back_half(candidate: list[dict], transcript_record: TranscriptRecord) -> list[dict]:
    # If the transcript is single exon, the PTC should be in the front half of the CDS.
    # so eliminate the candidates in the back half of the CDS.
    removed = []
    cds = generate_cdsseq(transcript_record)
    for grna in candidate:
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
def eliminate_in_last_exon(candidate: list[str], transcript_record: TranscriptRecord) -> list[str]:
    removed = []
    last_exon_start = transcript_record.exon_start_positions[-1]
    for grna in candidate:
        if "ct_seq" in grna:
            if re.search(grna["ct_seq"], transcript_record.transcript_seq).start() + 1 <= last_exon_start:
                removed.append(grna)
        elif "ga_seq" in grna:
            if re.search(grna["ga_seq"], transcript_record.transcript_seq).start() + 16 <= last_exon_start:
                removed.append(grna)
    return removed


# LEJ: Last Exon Junction. The junction between the last exon and the penultimate exon.
def label_in_50bp_from_LEJ(candidate: list[dict], transcript_record: TranscriptRecord) -> list[dict]:
    candidate_ = candidate.copy()
    last_2nd_exon_end = transcript_record.exon_end_positions[-2]
    for grna in candidate_:
        if "ct_seq" in grna:
            if re.search(grna["ct_seq"], transcript_record.transcript_seq).start() + 1 >= last_2nd_exon_end - 50:
                grna["in_50bp_from_LEJ"] = True
            else:
                grna["in_50bp_from_LEJ"] = False
        elif "ga_seq" in grna:
            if re.search(grna["ga_seq"], transcript_record.transcript_seq).start() + 16 >= last_2nd_exon_end - 50:
                grna["in_50bp_from_LEJ"] = True
            else:
                grna["in_50bp_from_LEJ"] = False
    return candidate_
