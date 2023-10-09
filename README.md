# MINORg
**Mi**nimum **No**n-**R**eference **g**RNA finder
- Finds the minimum gRNA set required to target multiple alignable genes in multiple non-reference genomes
- Available as both command line application and Python package

Preprint: https://www.biorxiv.org/content/10.1101/2022.03.10.481891

Publication (Nucleic Acids Research): https://doi.org/10.1093/nar/gkad142

Zenodo: https://doi.org/10.5281/zenodo.7644871

## Availability
Some dependencies are not available for Windows. Windows users should use a Linux emulator to run MINORg.

## Installation
A version of MINORg is available at pypi.org (minorg) as well as Docker (rlrq/minorg OR rlrq/minorg-lite). You may follow this [guide](https://rlrq.github.io/MINORg/build/html/installation.html) to install MINORg and its dependencies.

## Requirements
- [Python 3](https://www.python.org/)
- [BLAST+ suite](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [BEDTools](https://bedtools.readthedocs.io/en/latest/index.html)

## Links
- Tutorial, example, and documentation: https://rlrq.github.io/MINORg
- Detailed overview of steps in the programme: https://tinyurl.com/sr84ae9e (Google slides) (not up to date)
- Flowchart to select & use appropriate input parameters: https://tinyurl.com/jyke76b8 (PDF) (not up to date)

## IMPT
Please refer to slides/PDF in the 'Links' section for execution details for the version on the workstation (accessible only to lab members and guests with accounts).

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
      3. Eliminate candidate gRNA with hits outside masked regions in non-reference genome(s) and fail maximum match/gaps criteria
   2. Eliminate candidate gRNA that do not align within the CDS of reference genes
      1. Extract CDS-only regions of user-specified reference gene(s) from a reference genome (.fasta) using GFF3 annotation (.gff)
         - If the user specified a domain, the range will be restricted accordingly
      2. Align non-reference target sequences (output of step 1.5) with reference sequences from steps 1.1 (or 1.1.3 if domain is specified) and 2.1
      3. For all candidate gRNA, check their position in the alignment (based on where in each non-reference target they originate) and eliminate any gRNA that do not align within **AT LEAST ONE** reference gene's desired feature
         - Users may change the desired feature (default is CDS)
4. Find minimum gRNA set that covers all target sequences

## Inputs
- Step 1
   - Data:
      - Reference genome (--assembly xxx.fasta; used with --annotation)
      - Reference GFF3 annotation (--annotation xxx.gff3; used with --assembly)
      - Non-reference sequences/genome (--query xxx.fasta)
      - Target sequences (--target xxx.fasta)
   - Parameters:
      - Gene IDs (--gene)
         - Used with:
            - Reference genome (--assembly xxx.fasta --annotation xxx.gff3)
            - Query fasta file (--query xxx.fasta)
   - Optional parameters:
      - Minimum hit % identity (--minid 85 (%))
      - Minimum candidate target length (--minlen 0 (bp))
      - Maximum merge buffer (--buffer 100 (bp))
   - Optional for domain restriction:
      - PSSM-ID (--domain) and rpsblast+ database (--db)
- Step 2
   - Parameters:
      - PAM (--pam SpCas9)
      - gRNA length (--length 20 (bp))
- Step 3
   - Optional parameters:
      - Minimum off-target gaps (--ot-gap 0)
      - Minimum off-target mismatch (--ot-mismatch 1 (bp))
      - Position-dependent off-target thresholds (--ot-pattern <pattern>; overrides --ot-gap and --ot-mismatch)
   - Optional data:
      - Background sequences (--background xxx.fasta)
- Step 4
   - Optional paramters:
      - Number of sets to output (--set 1)
      - Manually approve each gRNA set (--manual)
