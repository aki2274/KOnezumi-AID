# KOnezumi-AID
`KOnezumi-AID` is the command-line tool to automate the gRNA design for multiplex KO mouse using Target-AID

## Installation
### Prerequisits
- Python 3.9 or later
- Unix-like environment (Linux, macOS, WSL, etc.)

### Input data set
#### Locus information
`refFlat.text` from [UCSC](
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/
)

### genomic sequence
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

### Installation🔨
#### Follow these steps to get started with KOnezumiAID using git clone
- Clone the repository

`git clone https://github.com/aki2274/KOnezumi-AID.git`

- Move into the cloned directory

`cd KOnezumi-AID`

- Build package with poetry

`poetry build`

- Install KOnezumiAID from the generated package

`pip install <package_name>.whl`

Now you're ready to use KOnezumiAID in your Python projects.

## Usage
### Create data set for KOnezumi-AID

`konezumiaid createdata <your refFlat.txt Path> <your mm39.fa Path>`
### Examples
`konezumiaid createdata data/refFlat.txt data/mm39.fa`

### Search candidate by gene symbol or transcript name
KOnezumi-AID accepts a gene symbol or a transcript name.

`konezumiaid <gene symbol or transcript name>`

### Examples
- Search candidate by the gene symbol

`konezumiaid Xkr4`

- Search candidate by the transcript name

`konezumiaid NM_001011874`