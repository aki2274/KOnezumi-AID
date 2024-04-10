# KOnezumi-AID
`KOnezumi-AID` is the command-line tool to automate the design of multiplex KO mouse using Target-AID

## Installation
### Prerequisits
- Python 3.9 or later
- Unix-like environment (Linux, macOS, WSL, etc.)

### Installation
- Clone the repository

`git clone https://github.com/aki2274/KOnezumi-AID.git`

- Move to the directory and Build package with poetry

`poetry build`

- Install KOnezumiAID!

`pip install ${package_name}`



## Input data set
### Locus information
`refFlat.text` from [UCSC](
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/
)

### genomic sequence
`mm39.fa.gz` from [UCSC](
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/
)


### Download scripts (bash)

```bash
mkdir -p data
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refFlat.txt.gz |
    gzip -dc > data/refFlat.txt
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz |
    gzip -dc > data/mm39.fa
```

## Usage
### Create data set for KOnezumi-AID

`createdata ${refFlat.txt Path} ${mm39.fa Path}`

### KOnezumi-AID accepts a gene symbol or a transcript name.
`konezumiaid ${gene symbol or transcript name}`
### Examples
- Search candidate by the gene symbol

`konezumiaid Xkr4`

- Search candidate by the transcript name

`konezumiaid NM_001011874`