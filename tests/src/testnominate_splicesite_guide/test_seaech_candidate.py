from src.konezumiaid.nominate_splicesite_guide.search_candidate import (
    search_site_candidate,
)
from src.konezumiaid.create_gene_dataclass import TranscriptRecord


def test_no_candidate():
    data = TranscriptRecord(
        "NNNNGGCCNNNNNNNAAAANNNNNNGTAGNCCNNNNNNNNNNNNNNNNNNGANNNNNNNCCNNNNNNNNNNNNNNNNAGNNNNNNNNNCCNNNNNNNNNNNNNNNNNNAGNNNNNNNNCCNNNNNNNNNNNNNNNNGTNNNNNNNNNNNN",
        # PAM isn't CC,there is TTTT(AAAA, because this gRNA is revers compliment), Donor isn't GT , Acceptor isn't AG,
        # the exon length is a multiple of 3 and skipping exon has only 3'UTR
        0,
        150,
        0,
        130,
        4,
        [0, 30, 80, 110, 146],
        [25, 50, 100, 137, 149],
    )

    result = search_site_candidate(data)
    excepted = ([], [])
    assert result == excepted


def test_has_candidate():
    data = TranscriptRecord(
        "AGNNCCNNNNNNNNNNNNNNNNNNNGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCNNNNNNNNNNNNNNNNNAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        # PAM isn't CC,there is TTTT(AAAA, because this gRNA is revers compliment), Donor isn't GT , Acceptor isn't AG,
        # the exon length is a multiple of 3 and skipping exon has only 3'UTR
        0,
        150,
        0,
        149,
        5,
        [2, 30, 80, 110, 146],
        [25, 50, 100, 137, 149],
    )
    result = search_site_candidate(data)
    excepted = (
        [
            {
                "seq": "NNCTNNNNNNNNNNNNNNNNNGG",
                "exon_index": 3,
                "link_to_crisprdirect": "https://crispr.dbcls.jp/?userseq=NNCTNNNNNNNNNNNNNNNNNGG&pam=NGG&db=mm39",
            }
        ],
        [
            {
                "seq": "ACNNNNNNNNNNNNNNNNNNNGG",
                "exon_index": 1,
                "link_to_crisprdirect": "https://crispr.dbcls.jp/?userseq=ACNNNNNNNNNNNNNNNNNNNGG&pam=NGG&db=mm39",
            }
        ],
    )
    assert result == excepted
