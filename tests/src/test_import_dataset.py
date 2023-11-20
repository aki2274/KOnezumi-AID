import pytest
from src import import_dataset
import pickle
from tempfile import NamedTemporaryFile


def test_read_fasta():
    # Create a temporary FASTA file for testing
    fasta_content = ">gene1\nATGC\n>gene2\nTTGGCC\n>gene3\nAAATTT\nCCCGGG"
    with NamedTemporaryFile(mode="w", delete=False) as temp_file:
        temp_file.write(fasta_content)
        temp_file_name = temp_file.name

    # Expected result
    expected = {"gene1": "ATGC", "gene2": "TTGGCC", "gene3": "AAATTTCCCGGG"}

    # Test the function
    assert import_dataset.read_fasta(temp_file_name) == expected


def test_read_pkl():
    # Create a temporary pkl file for testing
    pkl_content = [
        {"id": "1", "name": "A", "value": "30"},
        {"id": "2", "name": "B", "value": "25"},
    ]
    with NamedTemporaryFile(mode="wb", delete=False, suffix=".pkl") as temp_file:
        pickle.dump(pkl_content, temp_file)
        temp_file_name = temp_file.name
    # Expected result
    expected = [
        {"id": "1", "name": "A", "value": "30"},
        {"id": "2", "name": "B", "value": "25"},
    ]

    # Test the function
    assert import_dataset.read_pkl(temp_file_name) == expected
