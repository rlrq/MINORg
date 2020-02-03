# miNoR
**mi**nimum **no**n-**r**eference gRNA finder
Finds the minimum gRNA set required to target multiple alignable genes in multiple non-reference genomes

## Main steps
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
   5. Filter merged hits for minimum length and % identity
2. Identify candidate gRNA in non-reference targets
   1. Restricted by user-specified PAM and gRNA length
3. Screen candidate gRNA
   1. Eliminate candidate gRNA with off-target hits
      1. Mask non-reference targets in non-reference genome(s) (.fasta)
         - Only regions the length of targets with 100% identity to targets will be masked
         - All non-reference genomes provided will be screened simultaneously so all candidate gRNA that pass this screening test should not have off-targets in any of the non-reference genomes provided
         - WIP: screening against reference genome
      2. BLASTn candidate gRNA against masked non-reference genome(s)
      3. Eliminate candidate gRNA that align with masked non-reference genome(s) and fail maximum match/gaps criteria
   2. Eliminate candidate gRNA that do not align within the CDS of reference genes
      1. Extract CDS-only regions of user-specified reference gene(s) from a reference genome (.fasta) using GFF3 annotation (.gff)
         - If the user specified a domain, the range will be restricted accordingly
      2. Align non-reference target sequences (output of step 1.5) with reference sequences from steps 1.1 (or 1.1.3 if domain is specified) and 2.1
      3. For all candidate gRNA, check their position in the alignment (based on where in each non-reference target they originate) and eliminate any gRNA without **AT LEAST ONE** alignment within the reference CDS regions
4. Find minimum gRNA set that covers all target sequences
