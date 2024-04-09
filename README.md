# KOnezumi-AID
`KOnezumi-AID` is the command-line tool to automate the design of multiplex KO mouse using Target-AID

## Installation
### Prerequisits
- Python 3.9 or later
- Unix-like environment (Linux, macOS, WSL, etc.)

### Installation
Clone the repository

`$ git clone https://github.com/aki2274/KOnezumi-AID.git`

Build 

## Input data set
### genomic sequence
`mm39.fa.gz` from [UCSC](
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/
)


### Locus information
`refFlat.text` from [UCSC](
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/
)

### Download scripts (bash)

```bash
mkdir -p data
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz |
    gzip -dc > data/mm39.fa
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refFlat.txt.gz |
    gzip -dc > data/refFlat.txt
```

## Usage
### KOnezumi-AID accepts a gene symbol or a transcript name.
`KOnezumiAID ${gene symbol or transcript name}`
### Examples
- Search candidate by the gene symbol

`KOnezumiAID Xkr4`

- Search candidate by the transcript name

`KOnezumiAID NM_001011874`