Algorithms
==========

Non-reference homologue inference
---------------------------------

Inferrence of homologues is a BLAST-based step that is controlled by 7 parameters: ``--minid``, ``--minlen``, ``--mincdslen``, ``--merge-within``, ``--check-id-before-merge``, ``--check-recip``, ``--relax-recip``

It consists broadly of 3 steps:

#. BLASTN of target reference genes to non-reference genomes

   * ``--check-id-before-merge`` (flag) and ``--minid``: remove hits below percentage ID specified using ``--minid``
   * Output: BLAST hits
    
#. Merge BLAST hits into candidate targets

   * ``--minlen``: remove candidate targets shorter than specified length (bp)
   * ``--mincdslen``: remove candidate targets with fewer than the specified bases covered by hits to CDS of reference genes
   * ``--minid``: remove candidate targets where none of the hits that were merged into it have a percentage ID equal to or greater than the specified value
   * Output: Candidate targets
    
#. Optional: reciprocal BLAST of candidate targets back to reference genome(s)

   * Filter candidate targets by whether they are more similar to target reference genes than any other non-target reference genes
   * Particularly useful if your target genes have very similar homologues
   * Use one of the following flags to turn ON this option:
     
     * ``--check-recip`` (flag): Remove candidate targets where the hit with the highest bitscore is to a non-target reference genes EVEN if  the hit overlaps with a target reference gene
     * ``--relax-recip`` (flag): Remove candidate targets where the hit with the highest bitscore is to a non-target reference genes ONLY if the hit DOES NOT overlap with a target reference gene
       
   * Output: Filtered candidate targets


Off-target assessment
---------------------

Off-target assessment is a BLAST-based step controlled by 2 parameters: ``--background``, ``--screen-ref``, ``--mask-gene``, ``--unmask-gene``, ``--ot-gap``, ``--ot-mismatch``

It consist broadly of 3 steps

#. Mask targets (and/or target reference genes) in background sequences

   * Masking allows gRNA with hits to the masked region to NOT be considered off-target
   * Only regions in background sequences with 100% identity across the FULL LENGTH of a target/gene to mask will be masked
     
     * If a sequence comprises only part of a target/gene to mask, it will NOT be masked
       
   * ``--mask-gene``: reference gene(s) to mask. By default, all genes passed to ``--gene`` are masked.
   * ``--unmask-gene``: reference gene(s) to unmask. By default, no genes passed to ``--gene`` are unmasked.
     
     * ``--unmask-gene`` has priority. If a gene appears in both ``--unmask-gene`` and ``--mask-gene``, it will NOT be masked.
   
   * Output: Ranges of masked regions in background sequences
   
#. BLASTN of gRNA to background sequences

   * Any files passed to ``--query`` will be included in this screen
   * ``--screen-ref`` (flag): include reference genome(s)
   * ``--background``: file(s) containing additional background sequences
   * Output: BLAST hits
   
#. Identify potentially problematic hits

   * Hits are subjected to the following tests in order:
     
      1. Hits fully in masked regions will be considered non-problematic
      2. Hits that have a minimum number of gaps OR mismatches specified by the user will be considered non-problematic. Users may increase the thresholds for more stringent filtering.
         
         * ``--ot-gap``: minimum number of gaps (default=1; minimum is 1)
         * ``--ot-mismatch``: minimum number of mismatches (default=1; minimum is 1)
           
      3. Hits where a gRNA is unaligned for a length that is greater than '(``--ot-gap`` - 1) + (``--ot-mismatch`` - 1) - <gaps in hit> - <mismatches in hit>' will be considered non-problematic
         
         * If ``--ot-gap 1 --ot-mismatch 1``, then a hit without gaps or mismatches must be perfectly aligned across the full length of a gRNA to be considered problematic
         * If ``--ot-gap 1 --ot-mismatch 2``, then a hit without gaps or mismatches must be perfectly aligned across at least <gRNA length>-1 bp to be considered problematic
      4. Optional: Hits that do not have a PAM pattern nearby will be considered non-problematic
         
         * ``--ot-pamless`` (flag): use this flag to turn this option OFF
   
   * Hits that are not considered non-problematic by AT LEAST ONE of the above test will be considered problematic
   * Output: Problematic BLAST hits
      
#. gRNA with hits that are problematic are considered to have failed off-target assessment
   

Within-feature inference
------------------------

MINORg aligns unannotated targets to annotated reference genes (supplied using ``--gene`` (CLI) or ``genes`` (Python)) in order to infer gene feature positions.

.. image:: images/minorg_within_feature.png

An alignment of an unannotated target sequence with 2 homologous reference genes is shown in the figure above. In this example, the desired feature in which to generate gRNA is the coding region (CDS).

* An effective feature (CDS in this case) range is generated for each target-reference pair separately
* Where there is an insertion relative to a reference gene, an effective feature (CDS in this case) range is only continuous across it if the insertion is smaller than a user-specified max_insertion length (``--max-insertion``/\ ``max_insertion``, default: 15 bp)
  
  * Using a max_insertion of 15 bp, the insertions marked by dashed (smaller than 15 bp) and dotted (larger than 15 bp) lines are included and excluded from the effective CDS range respectively
    
* The minimum requirement for a gRNA to pass this check is to fall entirely within at least one gene's effective feature (CDS in this case) range
* Users may adjust the threshold for minimum number/fraction of effective feature (CDS in this case) ranges a gRNA is required to be contained within to pass
  
  * If using ``--min-within-n``/\ ``min_within_n`` of 1
    
    * gRNA are required to fall entirely within only 1 gene's effective feature (CDS in this case) range
    * Both gRNA 2 and gRNA 4 pass
      
  * If using ``--min-within-fraction``/\ ``min_within_fraction`` of 1
    
    * gRNA are required to fall entirely within ALL genes' effective feature (CDS in this case) ranges
    * Only gRNA 2 passes
    * If parts of your genes are freqently pseudogenised, you may wish to set this value fairly high in order to ensure that most, if not all, gRNA are in conserved coding regions
