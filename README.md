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

#### Required Packages (Not needed if installed via Bioconda)
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
curl https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refFlat.txt.gz |
    gzip -dc > data/refFlat.txt
curl https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz |
    gzip -dc > data/mm39.fa
```



## Usage
### KOnezumi-AID's output   
- gRNAs to generate PTC (premature termination codon)
- gRNAs to disrupt splice acceptor site
- gRNAs to disrupt splice donor site

KOnezumi-AID provides these gRNAs in standerd output and CSV format.   
The CSV file is located in `data/output` directory and named as `<gene symbol/transcript name>_ptc_gRNA.csv` or `<gene symbol/transcript name>_splice_gRNA.csv`.

### Create data set for KOnezumi-AID

`konezumiaid preprocess <your refFlat.txt Path> <your mm39.fa Path>`

### Example

`konezumiaid preprocess data/refFlat.txt data/mm39.fa`

### Search candidate by gene symbol or transcript name (Refseq id)

KOnezumi-AID accepts a gene symbol or a transcript name.  
`konezumiaid <-n | --name> <gene symbol | transcript name>`

#### Search by gene symbol
You can obtain the gRNAs that are present in all transcript variants.  
> [!NOTE]  
> If the gene has one transcript, the result is the same as searching by the transcript name

#### Search by transcript name
You can obtain the transcript's gRNAs and access more information about the gRNAs.  
- in_start_150bp: The gRNA is located in the first 150bp of the transcript or not.
- in_50bp_from_LEJ: The gRNA is located in the 50bp from the last exon junction or not.
- exon_index: the index of the exon where the gRNA is located.

### Examples
#### Search candidate by the gene symbol (gene symbol with multiple transcripts)
```
$konezumiaid -n Rp1     
Processing NM_001370921...
Processing NM_001195662...
Processing NM_011283...
List of gRNAs to generate PTC (premature termination codon)
  Target sequence (20mer + PAM) Target amino acid                                                      link to CRISPRdirect
0       ACAGTTTGGCGGCGTTCGGGTGG                 Q  https://crispr.dbcls.jp/?userseq=ACAGTTTGGCGGCGTTCGGGTGG&pam=NGG&db=mm39
1       ACGACACAGCATCACCAGGCTGG                 R  https://crispr.dbcls.jp/?userseq=ACGACACAGCATCACCAGGCTGG&pam=NGG&db=mm39
2       ACAGGTTATGCAGTGTCCTGTGG                 Q  https://crispr.dbcls.jp/?userseq=ACAGGTTATGCAGTGTCCTGTGG&pam=NGG&db=mm39
3       CCAGGGCCGAGGGCGCCTGCGGG                 W  https://crispr.dbcls.jp/?userseq=CCAGGGCCGAGGGCGCCTGCGGG&pam=NGG&db=mm39
4       GCCAGGGCCGAGGGCGCCTGCGG                 W  https://crispr.dbcls.jp/?userseq=GCCAGGGCCGAGGGCGCCTGCGG&pam=NGG&db=mm39
List of gRNAs to disrupt splice acceptor site
No gRNA found.
List of gRNAs to disrupt splice donor site
No gRNA found.
```
#### Search candidate by the gene symbol (gene symbol with single transcript)  
```
$konezumiaid -n Mafa    
List of gRNAs to generate PTC (premature termination codon)
  Target sequence (20mer + PAM)  Recommended Target amino acid                                                      link to CRISPRdirect
0       CTCAGGCCGGGGGCGCCCCGGGG         True               87Q  https://crispr.dbcls.jp/?userseq=CTCAGGCCGGGGGCGCCCCGGGG&pam=NGG&db=mm39
1       GCTCAGGCCGGGGGCGCCCCGGG         True               87Q  https://crispr.dbcls.jp/?userseq=GCTCAGGCCGGGGGCGCCCCGGG&pam=NGG&db=mm39
2       CCAGCACCACCTGAACCCCGAGG         True              121Q  https://crispr.dbcls.jp/?userseq=CCAGCACCACCTGAACCCCGAGG&pam=NGG&db=mm39
3       GGTCAGAGCTTCGCGGGCGGCGG         True              167Q  https://crispr.dbcls.jp/?userseq=GGTCAGAGCTTCGCGGGCGGCGG&pam=NGG&db=mm39
List of gRNAs to disrupt splice acceptor site
No gRNA found.
List of gRNAs to disrupt splice donor site
No gRNA found.
```

#### Search candidate by the transcript name  
```
$ konezumiaid -n NM_001370921
List of gRNAs to generate PTC (premature termination codon)
   Target sequence (20mer + PAM)  Recommended Target amino acid                                                      link to CRISPRdirect
0        ACAGTTTGGCGGCGTTCGGGTGG        False               46Q  https://crispr.dbcls.jp/?userseq=ACAGTTTGGCGGCGTTCGGGTGG&pam=NGG&db=mm39
1        ACGACACAGCATCACCAGGCTGG         True               88R  https://crispr.dbcls.jp/?userseq=ACGACACAGCATCACCAGGCTGG&pam=NGG&db=mm39
2        ACAGGTTATGCAGTGTCCTGTGG         True              192Q  https://crispr.dbcls.jp/?userseq=ACAGGTTATGCAGTGTCCTGTGG&pam=NGG&db=mm39
3        ACAACCTGTCCTTCCAGGTAAGG         True              389Q  https://crispr.dbcls.jp/?userseq=ACAACCTGTCCTTCCAGGTAAGG&pam=NGG&db=mm39
4        ACCAATCAGAACAATCCCACTGG         True              698Q  https://crispr.dbcls.jp/?userseq=ACCAATCAGAACAATCCCACTGG&pam=NGG&db=mm39
5        ACGAATGTATCTGAGGATTAAGG         True              723R  https://crispr.dbcls.jp/?userseq=ACGAATGTATCTGAGGATTAAGG&pam=NGG&db=mm39
6        TCAGGCCAATGTCACATTGTGGG         True              861Q  https://crispr.dbcls.jp/?userseq=TCAGGCCAATGTCACATTGTGGG&pam=NGG&db=mm39
7        CTCAGGCCAATGTCACATTGTGG         True              861Q  https://crispr.dbcls.jp/?userseq=CTCAGGCCAATGTCACATTGTGG&pam=NGG&db=mm39
8        CCAGGGCCGAGGGCGCCTGCGGG         True              126W  https://crispr.dbcls.jp/?userseq=CCAGGGCCGAGGGCGCCTGCGGG&pam=NGG&db=mm39
9        GCCAGGGCCGAGGGCGCCTGCGG         True              126W  https://crispr.dbcls.jp/?userseq=GCCAGGGCCGAGGGCGCCTGCGG&pam=NGG&db=mm39
10       TCCAGTGGGATTGTTCTGATTGG         True              704W  https://crispr.dbcls.jp/?userseq=TCCAGTGGGATTGTTCTGATTGG&pam=NGG&db=mm39
11       CCAGTACTGGGATTTGTCACTGG         True             1052W  https://crispr.dbcls.jp/?userseq=CCAGTACTGGGATTTGTCACTGG&pam=NGG&db=mm39
List of gRNAs to disrupt splice acceptor site
  Target sequence (20mer + PAM)  Exon index                                                      link to CRISPRdirect
0       ACCTGGGATTGAAAGGAACAAGG          20  https://crispr.dbcls.jp/?userseq=ACCTGGGATTGAAAGGAACAAGG&pam=NGG&db=mm39
1       TCTGTTGGAGAAAAGCCCCATGG          22  https://crispr.dbcls.jp/?userseq=TCTGTTGGAGAAAAGCCCCATGG&pam=NGG&db=mm39
2       ACCTGAAGAAAATGGAAAACAGG          23  https://crispr.dbcls.jp/?userseq=ACCTGAAGAAAATGGAAAACAGG&pam=NGG&db=mm39
List of gRNAs to disrupt splice donor site
  Target sequence (20mer + PAM)  Exon index                                                      link to CRISPRdirect
0       TACCTTGCCCAAGTCCATCATGG           8  https://crispr.dbcls.jp/?userseq=TACCTTGCCCAAGTCCATCATGG&pam=NGG&db=mm39
1       TTACCTCTCACAGGTGAAGATGG          22  https://crispr.dbcls.jp/?userseq=TTACCTCTCACAGGTGAAGATGG&pam=NGG&db=mm39
```