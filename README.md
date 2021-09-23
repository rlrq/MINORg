# MINORg
**Mi**nimum **No**n-**R**eference **g**RNA finder
- Finds the minimum gRNA set required to target multiple alignable genes in multiple non-reference genomes

## IMPT
The code in this repository is not complete. Only sections of it are viable. This README file is also not up-to-date. Please refer to slides/PDF in the 'Links' section for execution details for the version on the workstation (accessible only to lab members and guests with accounts).

## Requirements
- BLAST+ suite
- Python 3
   - Biopython
- bedtools
- mafft

## Links
- Detailed overview of steps in the programme: https://tinyurl.com/sr84ae9e (Google slides)
- Flowchart to select & use appropriate input parameters: https://tinyurl.com/jyke76b8 (PDF)

## Overview of steps
1. Identify candidate targets in non-reference genome
   1. Extract user-specified reference gene(s) from a reference genome (.fasta) using GFF3 annotation (.gff)
      - Sequence(s) include introns
      - Optional: User may specify a protein domain (using CDD PSSM-ID) to restrict search
      1. CDS-only regions of user-specified reference gene(s) from a reference genome (.fasta) will be extracted and translated using GFF3 annotation
      2. RPS-BLAST protein sequence(s) to domain database and identify domain range(s)
      3. Extract user-specified reference gene(s) from a reference genome (.fasta) using GFF3 annotation (.gff) and restricted to the corresponding genomic coordinates of the domains
   2. BLASTn reference gene(s) against non-reference genome(s) (.fasta)
   3. Filter hits by minimum % identity (optional)
   4. Merge overlapping hits within specified distance of each other (to accommodate introns/insertions)
   5. Filter merged hits for minimum length and % identity into target sequences
   6. Filter target sequences for those with best alignment to target genes(s) (optional)
      - Ensures that genes that are similar but not part of the set of user-specified target gene(s) will not be targeted
2. Identify candidate gRNA in non-reference targets
   1. Restricted by user-specified PAM and gRNA length
3. Screen candidate gRNA
   1. Eliminate candidate gRNA with off-target hits
      1. Mask targets in non-reference genome(s) (.fasta)
         - Only regions the length of targets with 100% identity to targets will be masked
         - All non-reference genomes provided will be screened simultaneously so all candidate gRNA that pass this screening test should not have off-targets in any of the non-reference genomes provided
         - User may also provide sequences to check against
      2. BLASTn candidate gRNA against masked non-reference and reference genome(s)
         - Optional: Screen reference genome also
      4. Eliminate candidate gRNA that align with masked non-reference genome(s) and fail maximum match/gaps criteria
   2. Eliminate candidate gRNA that do not align within the CDS of reference genes
      1. Extract CDS-only regions of user-specified reference gene(s) from a reference genome (.fasta) using GFF3 annotation (.gff)
         - If the user specified a domain, the range will be restricted accordingly
      2. Align non-reference target sequences (output of step 1.5) with reference sequences from steps 1.1 (or 1.1.3 if domain is specified) and 2.1
      3. For all candidate gRNA, check their position in the alignment (based on where in each non-reference target they originate) and eliminate any gRNA without **AT LEAST ONE** alignment within the reference CDS regions
4. Find minimum gRNA set that covers all target sequences

## Inputs
- Step 1
   - Data:
      - Reference genome (--ref xxx.fasta)
      - Reference GFF3 annotation (--gff xxx.gff)
      - Non-reference sequences/genome (--nonref xxx.fasta)
   - Parameters:
      - Gene IDs (--gene)
         - Used with:
            - Accession/individual (-a) OR
            - Query fasta file (-q xxx.fasta)
      - Target sequences (--target xxx.fasta)
   - Optional parameters:
      - Minimum hit % identity (--minid 85 (%))
      - Minimum candidate target length (--minlen 0 (bp))
      - Maximum merge buffer (--buffer 100 (bp))
   - Optional for domain restriction:
      - PSSM-ID and rpsblast+ database
- Step 2
   - Parameters:
      - PAM (--pam)
      - gRNA length (--length)
- Step 3
   - Optional parameters:
      - Minimum off-target gaps (--gaps 0)
      - Minimum off-target mismatch (--mismatch 1 (bp))
   - Optional data:
      - Background sequences (--background xxx.fasta)
- Step 4
   - Optional paramters:
      - Number of sets to output (--set 1)
      - Minimum set algorithm (--algo LAR)
