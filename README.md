[![License](https://img.shields.io/badge/License-MIT-9cf.svg)](https://choosealicense.com/licenses/mit/)
[![DOI](https://zenodo.org/badge/673151657.svg)](https://zenodo.org/badge/latestdoi/673151657)

# KOnezumi-AID
`KOnezumi-AID` is the command-line tool to automate the gRNA design for multiplex KO mouse using Target-AID

## Installation
### Prerequisits
- Python 3.9 or later
- Unix-like environment (Linux, macOS, WSL2, etc.)

### Installation🔨
#### From [Bioconda](https://anaconda.org/bioconda/konezumiaid) (Recommended)

`conda install -c conda-forge -c bioconda konezumiaid`

#### From [PyPI](https://libraries.io/pypi/KOnezumiAID):

`pip install KOnezumiAID`

### Required Packages (Not needed if installed via conda)
- bedtools


Follow the [official instruction](https://bedtools.readthedocs.io/en/latest/content/installation.html)

- bowtie

Follow the [official instruction](https://bowtie-bio.sourceforge.net/manual.shtml#:~:text=is%20future%20work.-,Obtaining%20Bowtie,-You%20may%20download)

### Input data set (e.g. Mus musculus mm39)
#### Locus information
`refFlat.txt.gz` from [UCSC](
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/
)

#### genomic sequence
`mm39.fa.gz` from [UCSC](
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/
)

#### Download scripts (bash)

```bash
mkdir -p data
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refFlat.txt.gz |
    gzip -dc > data/refFlat.txt
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz |
    gzip -dc > data/mm39.fa
```



## Usage
### Create data set for KOnezumi-AID

`konezumiaid preprocess <your refFlat.txt Path> <your mm39.fa Path>`

### Examples

`konezumiaid preprocess data/refFlat.txt data/mm39.fa`

### Search candidate by gene symbol or transcript name (Refseq id)

KOnezumi-AID accepts a gene symbol or a transcript name.

- Search by gene symbol

You can obtain the gRNAs that are present in all transcript variants.

- Search by transcript name

You can obtain the transcript's gRNAs and access more information about the gRNAs.


`konezumiaid <-n | --name> <gene symbol | transcript name>`

### Examples
- Search candidate by the gene symbol
```
$konezumiaid Rp1
PTC gRNA
                       seq
0  ACAGGTTATGCAGTGTCCTGTGG
1  ACGACACAGCATCACCAGGCTGG
2  GCCAGGGCCGAGGGCGCCTGCGG
3  ACAGTTTGGCGGCGTTCGGGTGG
4  CCAGGGCCGAGGGCGCCTGCGGG
Acceptor gRNA
No Acceptor gRNA found.
Donor gRNA
No Donor gRNA found.
```

- Search candidate by the transcript name

```
$ konezumiaid NM_001370921
PTC gRNA
                        seq  in_start_150bp  in_50bp_from_LEJ
0   ACAGTTTGGCGGCGTTCGGGTGG            True             False
1   ACGACACAGCATCACCAGGCTGG           False             False
2   ACAGGTTATGCAGTGTCCTGTGG           False             False
3   ACAACCTGTCCTTCCAGGTAAGG           False             False
4   ACCAATCAGAACAATCCCACTGG           False             False
5   ACGAATGTATCTGAGGATTAAGG           False             False
6   TCAGGCCAATGTCACATTGTGGG           False             False
7   CTCAGGCCAATGTCACATTGTGG           False             False
8   CCAGGGCCGAGGGCGCCTGCGGG           False             False
9   GCCAGGGCCGAGGGCGCCTGCGG           False             False
10  TCCAGTGGGATTGTTCTGATTGG           False             False
11  CCAGTACTGGGATTTGTCACTGG           False             False
Acceptor gRNA
                       seq  exon_index
0  ACCTGGGATTGAAAGGAACAAGG          20
1  TCTGTTGGAGAAAAGCCCCATGG          22
2  ACCTGAAGAAAATGGAAAACAGG          23
Donor gRNA
                       seq  exon_index
0  TACCTTGCCCAAGTCCATCATGG           8
1  TTACCTCTCACAGGTGAAGATGG          22
```
