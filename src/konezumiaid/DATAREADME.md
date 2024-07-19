# Description
This directory contains the output data of the KOnezumi-AID. The data is stored in CSV format.

## Output Description

### `<gene name>_ptc_gRNA.csv`
This file contains the gRNAs designed to induce premature termination codon (PTC) in the specified gene. The file includes the following columns:
- **Target Sequence (20mer + PAM)**: The candidate 20mer gRNA sequence including the PAM sequence (NGG).
- **Recommended**: Indicates whether the target sequence meets the following criteria :
  - The gRNA sequence is not within 150bp downstream of the start codon.
  - For genes with multiple exons:
    - The gRNA sequence is not within 50bp upstream of the last exon junction.
    - The gRNA sequence is not within a long exon (exon length > 400bp).

  **NOTE** This column is present only if the gene has a single transcript.
- **Target Amino Acid**: The amino acid that the gRNA targets. If the gene has one transcript, this column shows the position of the amino acid in the protein.
- **Link to CRISPRdirect**: A link to the CRISPRdirect (mm39).

### `<gene name>_splice_gRNA.csv`
This file contains the gRNAs designed to disrupt the splice site in the specified gene. The file includes the following columns:
- **Target Sequence (20mer + PAM)**: The candidate 20mer gRNA sequence including the PAM sequence (NGG).
- **Exon Index**: The index of the exon where the gRNA is located  
**NOTE** This column is present only if the gene has a single transcript.
- **Link to CRISPRdirect**: A link to the CRISPRdirect (mm39).