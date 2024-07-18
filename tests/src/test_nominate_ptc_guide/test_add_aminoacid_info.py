from src.konezumiaid.nominate_ptc_guide.add_aminoacid_info import (
    transrate_codon_to_aminoacid,
    link_position_and_aminoacid,
    add_aminoacid_info,
)


def test_transrate_codon_to_aminoacid():
    seq = "TTTTTCTTATTGTCTTCCTCATCGTATTACTAATAGTGTTGCTGATGGCTTCTCCTACTGCCTCCCCCACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGAAGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG"
    excepted = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    assert transrate_codon_to_aminoacid(seq) == excepted


def test_link_position_and_aminoacid():
    positions = [1, 2, 3]
    cds_position = [0, 3, 6]
    seq = ["ATG", "ATC", "GGA"]
    expected = [
        {"position": 1, "aminoacid": "1M"},
        {"position": 2, "aminoacid": "2I"},
        {"position": 3, "aminoacid": "3G"},
    ]
    assert link_position_and_aminoacid(positions, cds_position, seq) == expected


def test_add_aminoacid_info():
    candidate_grna = [
        {"position": 1, "seq": "ATG"},
        {"position": 2, "seq": "ATC"},
        {"position": 3, "seq": "GGA"},
    ]
    aminoacid = [
        {"position": 1, "aminoacid": "M"},
        {"position": 2, "aminoacid": "I"},
        {"position": 3, "aminoacid": "G"},
    ]
    expected = [
        {"position": 1, "seq": "ATG", "aminoacid": "M"},
        {"position": 2, "seq": "ATC", "aminoacid": "I"},
        {"position": 3, "seq": "GGA", "aminoacid": "G"},
    ]
    assert add_aminoacid_info(candidate_grna, aminoacid) == expected
