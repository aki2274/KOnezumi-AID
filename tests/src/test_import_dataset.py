import pytest
from src import import_dataset
from tempfile import NamedTemporaryFile

def test_read_fasta():
    # Create a temporary FASTA file for testing
    fasta_content = ">gene1\nATGC\n>gene2\nTTGGCC\n>gene3\nAAATTT\nCCCGGG"
    with NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file.write(fasta_content)
        temp_file_name = temp_file.name

    # Expected result
    expected = {
        'gene1': 'ATGC',
        'gene2': 'TTGGCC',
        'gene3': 'AAATTTCCCGGG'
    }

    # Test the function
    assert import_dataset.read_fasta(temp_file_name) == expected


def test_read_csv():
    # Create a temporary CSV file for testing
    csv_content = "id,name,value\n1,A,30\n2,B,25"
    with NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as temp_file:
        temp_file.write(csv_content)
        temp_file_name = temp_file.name

    # Expected result
    expected = [
        {'id': '1', 'name': 'A', 'value': '30'},
        {'id': '2', 'name': 'B', 'value': '25'}
    ]

    # Test the function
    assert import_dataset.read_csv(temp_file_name) == expected
