from scripts.example import read_fasta


def test_read_fasta():
    result = read_fasta("tests/data/example.fasta")
    expected = {"test1": "ACGTACGTACGTACGTACGT", "test2": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"}
    assert result == expected
