# Hello KOnezumi-AID

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


#### Use as pandas.DataFrame in this code
DataFrame refFlat has 11 columns.
- string  geneName;           "Name of gene as it appears in Genome Browser." 
- string  name;               "Name of gene" 
- string  chrom;              "Chromosome name" 
- char[1] strand;             "+ or - for strand" 
- uint    txStart;            "Transcription start position" 
- uint    txEnd;              "Transcription end position" 
- uint    cdsStart;           "Coding region start" 
- uint    cdsEnd;             "Coding region end" 
- uint    exonCount;          "Number of exons" 
- uint[exonCount] exonStarts; "Exon start positions" 
- uint[exonCount] exonEnds;   "Exon end positions" 
