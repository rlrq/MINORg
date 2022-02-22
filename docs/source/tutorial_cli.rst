Tutorial (Command line)
=======================

In all the following tutorial, the current directory/working directory is presumed to contain all files in https://github.com/rlrq/MINORg/tree/master/examples. If you have not downloaded the files, please do so and navigate to the directory that contains them.

Note that all command line code snippets in the following tutorial are for **bash terminal**. You may have to adapt them according to your operating system.

Defining target sequences
~~~~~~~~~~~~~~~~~~~~~~~~~

User-provided targets
+++++++++++++++++++++

Let us begin with the simplest MINORg execution:

.. code-block:: bash
   
   $ minorg --directory ./example_00_target --target ./sample_CDS.fasta
   Final gRNA sequence(s) have been written to minorg_gRNA_final.fasta
   Final gRNA sequence ID(s), gRNA sequence(s), and target(s) have been written to minorg_gRNA_final.map
   
   1 mutually exclusive gRNA set(s) requested. 1 set(s) found.
   Output files have been generated in /path/to/current/directory/example00

The above combination of arguments tells MINORg to generate gRNA from targets in a user-provided FASTA file (``--target ./sample_CDS.fasta``) and to output files into directory ``--directory ./example00``. By default, MINORg generates 20 bp gRNA using NGG PAM.


Reference gene(s) as targets
++++++++++++++++++++++++++++

Both examples below specify a reference assembly (``--assembly ./subset_ref_TAIR10.fasta``) and annotation (``--annotation ./subset_ref_TAIR10.gff``) file, allowing MINORg to retrieve the gene sequence as target.

Single gene
^^^^^^^^^^^

.. code-block:: bash
   
   $ minorg --directory ./example_01_singlegene \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff
   Extracting reference sequences
   Finding homologues
   Max acceptable insertion length: 15
   Final gRNA sequence(s) have been written to minorg_gRNA_final.fasta
   Final gRNA sequence ID(s), gRNA sequence(s), and target(s) have been written to minorg_gRNA_final.map

   1 mutually exclusive gRNA set(s) requested. 1 set(s) found.
   Output files have been generated in /path/to/current/directory/example01

In the above example, ``--indv ref`` tells MINORg to generate gRNA for reference gene(s), and ``--gene AT5G45050`` tells MINORg that the target gene is AT5G45050.

Multiple genes
^^^^^^^^^^^^^^

Using ``--gene``
****************

See also: :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

To specify multiple genes, simply use ``--gene`` with comma-separated gene IDs, or ``--gene`` multiple times

.. code-block:: bash
                
   $ minorg --directory ./example_02_multigene \
            --indv ref --gene AT5G45050,AT5G45060,AT5G45200,AT5G45210,AT5G45220,AT5G45230,AT5G45240,AT5G45250 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff

OR

.. code-block:: bash
                
   $ minorg --directory ./example_02_multigene \
            --indv ref --gene AT5G45050 --gene AT5G45060 --gene AT5G45200 \
            --gene AT5G45210 --gene AT5G45220 --gene AT5G45230 --gene AT5G45240 --gene AT5G45250 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff


Using ``--cluster``
*******************
See also: :ref:`Configuration:2-level lookup`, :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

MINORg can also accept preset combinations of genes using ``--cluster`` and ``--cluster-set``. ``--cluster-set`` accepts a tab-separated lookup file that maps alias(es) to a combinations of genes (see :ref:`Configuration:cluster` for format). ``--cluster`` is used to specify the alias of a combination of genes in that lookup file.

.. code-block:: bash
                
   $ minorg --directory ./example_03_cluster \
            --indv ref --cluster RPS6 --cluster-set ./subset_cluster_mapping.txt
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff

The above code snippet is effectively identical to the examples in :ref:`Tutorial_cli:Multiple genes`.

Like ``--gene``, multiple combinations of genes can be specified to ``--cluster``. However, unlike ``--gene``, each combination will be processed separately (i.e. minimum sets will be separately generated for each combination).

.. code-block:: bash
                
   $ minorg --directory ./example_03_cluster \
            --indv ref --cluster RPS6,TTR1 --cluster-set ./subset_cluster_mapping.txt
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff


Multiple and non-standard reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See: :ref:`Tutorial_cli:Defining reference genomes`

Multiple reference genomes may be useful when generating gRNA across species boundaries. See :ref:`Tutorial_cli:Multiple reference genomes` for how to specify and use multiple reference genomes.

Some reference genomes may require non-standard genetic code (applicable only with the use of ``--domain``) or have unusual attribute field names in their GFF3 annotation files. See :ref:`Tutorial_cli:Non-standard genetic code` for how to specify non-standard genetic codes and :ref:`Tutorial_cli:Non-standard GFF3 attribute field names` for how to specify mapping of unusual GFF3 attribute field names to standard field names.


Non-reference gene(s) as targets
++++++++++++++++++++++++++++++++

Annotated genes
^^^^^^^^^^^^^^^

If your target genes have been annotated in their non-reference genomes (i.e. you have a GFF3 file containing annotations of your targets), you can use :ref:`Tutorial_cli:Reference gene(s) as targets` if you have a single non-reference genome, or :ref:`Tutorial_cli:Multiple reference genomes` if you have multiple non-reference genomes. In either case, you may treat your non-reference genome the same way you would a reference genome.


Unannotated genes
^^^^^^^^^^^^^^^^^

Using ``--extend-gene`` and ``--extend-cds``
********************************************
See also: :ref:`Parameters:Extended genome`

If you have both genomic and CDS-only sequences of your target genes but not a GFF3 annotation file, MINORg can infer coding regions (CDS) for your target genes using ``--extend-gene`` and ``--extend-cds``.

.. code-block:: bash

   $ minorg --directory ./example_04_ext \
            --indv ref --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

Note that ``--extend-gene`` and ``--extend-cds`` effectively add new genes to the reference genome, so they can be used just like any reference gene. Therefore, they can also be used in combination with ``--query`` or ``--indv``.

Using ``--query``
*****************

See also: :ref:`Algorithms:Non-reference homologue inference`, :ref:`Parameters:Multi-argument (CLI)`

If you would like MINORg to infer homologues genes in non-reference genomes, you can use ``--query`` to specify the FASTA files of those non-reference genomes. You may provide multiple non-reference genomes by using ``--query`` multiple times.

.. code-block:: bash

   $ minorg --directory ./example_05_query \
            --query ./subset_9654.fasta --query ./subset_9655.fasta \
            --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

``--query`` can be used in combination with ``--indv``. For inference parameters, see :ref:`Algorithms:Non-reference homologue inference`.


Using ``--indv``
****************

See also: :ref:`Algorithms:Non-reference homologue inference`, :ref:`Configuration:2-level lookup`, :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

You can also use ``--indv`` to ask MINORg to infer homologues genes in non-reference genomes. Similar to ``--clusters``, MINORg accepts a lookup file for non-reference genomes using ``--genome-set`` (see :ref:`Configuration:genome` for format) and one or more non-reference genome alias using ``--indv``.

.. code-block:: bash

   $ minorg --directory ./example_06_indv \
            --indv 9654,9655 --genome-set ./subset_genome_mapping.txt \
            --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

The above code snippet is effectively identical to the example in :ref:`Tutorial_cli:Using \`\`--query\`\``.

``--indv`` can be used in combination with ``--query``. For inference parameters, see :ref:`Algorithms:Non-reference homologue inference`.


Domain as targets
+++++++++++++++++

MINORg allows users to specify the identifier of an RPS-BLAST position-specific scoring matrix (PSSM-Id) to further restrict the target sequence to a given domain associated with the PSSM-Id. This could be particularly useful when designing gRNA for genes that do not share conserved domain structures but do share a domain that you wish to knock out.

Local database
^^^^^^^^^^^^^^

.. code-block:: bash

   $ minorg --directory ./example_07_domain \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --rpsblast /path/to/rpsblast/executable --db /path/to/rpsblast/db \
            --domain 214815

In the above example, gRNA will be generated for the WRKY domain (PSSM-Id 214815 as of CDD database v3.18) of the gene AT5G45050. Users are responsible for providing the PSSM-Id of a domain that exists in the gene. Unlike other examples, the database (``--db``) is not provided as part of the example files. You will have to download it yourself. See :ref:`Parameters:RPS-BLAST local database` for more information.

Remote database
^^^^^^^^^^^^^^^

While it is in theory possible to use the remote CDD database & servers instead of local ones, the ``--remote`` option for the 'rpsblast'/'rpsblast+' command from the BLAST+ package has never worked for me. In any case, if your version of local rpsblast is able to access the remote database, you can use ``--remote-rps`` instead of ``--db /path/to/rpsblast/db``.

.. code-block:: bash

   $ minorg --directory ./example_07_domain \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --rpsblast /path/to/rpsblast/executable --remote-rps \
            --domain 214815


Feature as targets
++++++++++++++++++

You may specify gene features to restrict gRNA to. By default, MINORg generates gRNA in coding regions (CDS). However, so long as a feature type is valid in the GFF3 annotation file provided to MINORg, gRNA can theoretically be designed for any feature type.

.. code-block:: bash

   $ minorg --directory ./example_08_feature \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --feature three_prime_UTR

The above example will generate gRNA in the 100 bp 3' UTR of AT5G45050.


Defining gRNA
~~~~~~~~~~~~~

See also: :ref:`Parameters:PAM`

By default, MINORg generates 20 bp gRNA using SpCas9's NGG PAM. You may specify other gRNA length using ``--length`` and other PAM using ``--pam``.

.. code-block:: bash

   $ minorg --directory ./example_09_grna \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --length 19 --pam Cas12a

In the example above, MINORg will generate 19 bp gRNA (``--length 19``) using Cas12a's unusual 5' PAM pattern (TTTV<gRNA>) (``--pam Cas12a``). MINORg has several built-in PAMs (see :ref:`Parameters:Preset PAM patterns` for options), and also supports customisable PAM patterns using ambiguous bases and regular expressions (see :ref:`Parameters:PAM` for format).


Filtering gRNA
~~~~~~~~~~~~~~

MINORg supports 3 different gRNA filtering options, all of which can be used together.

Filter by GC content
++++++++++++++++++++

.. code-block:: bash

   $ minorg --directory ./example_10_gc \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --gc-min 0.2 --gc-max 0.8

In the above example, MINORg will exclude gRNA with less than 20% (``--gc-min 0.2``) or greater than 80% (``--gc-max 0.8``) GC content. By default, minimum GC content is 30% and maximum is 70%.


Filter by off-target
++++++++++++++++++++

See: :ref:`Algorithms:Off-target assessment`

.. code-block:: bash

   $ minorg --directory ./example_11_ot_ref \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --screen-reference \
            --background ./subset_ref_Araly2.fasta --background ./subset_ref_Araha1.fasta \
            --ot-indv 9654,9655 --genome-set ./subset_genome_mapping.txt \
            --ot-gap 2 --ot-mismatch 2

In the above example, MINORg will screen gRNA for off-targets in:

* The reference genome (``--screen-reference``)
* Two different FASTA files (``--background ./subset_Araly2.fasta --background ./subset_Araha1.fasta``)
* Two non-reference genomes (``--ot-indv 9654,9655 --genome-set ./subset_genome_mapping.txt``)
  * ``--ot-indv`` functions similarly to ``--indv`` in that it requires ``--genome-set``, except that ``--ot-indv`` specifies non-refernece genomes for off-target assessment

``--ot-gap`` and ``--ot-mismatch`` control the minimum number of gaps or mismatches off-target gRNA hits must have to be considered non-problematic; any gRNA with at least one problematic gRNA hit will be excluded. See :ref:`Algorithms:Off-target assessment` for more on the off-target assessment algorithm.

In this case, ``--screen-reference`` is actually redundant as the genome from which targets are obtained (which, because of ``--indv ref``, is the reference genome) are automatically included for background check. However, in the example below, when the targets are from non-reference genomes, the reference genome is not automatically included for off-target assessment and thus ``--screen-reference`` is NOT redundant. Additionally, do note that the genes passed to ``--gene`` are masked in the reference genome, such that any gRNA hits to them are NOT considered off-target and will NOT be excluded.

.. code-block:: bash

   $ minorg --directory ./example_12_ot_nonref \
            --indv 9654 --genome-set ./subset_genome_mapping.txt \
            --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --screen-ref --background ./subset_ref_Araly2.fasta --background ./subset_ref_Araha1.fasta \
            --ot-indv 9655 \
            --ot-gap 2 --ot-mismatch 2

To skip off-target check entirely, use ``--skip-bg-check``.

.. code-block:: bash

   $ minorg --directory ./example_13_skipbgcheck \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --skip-bg-check


Filter by feature
+++++++++++++++++

See: :ref:`Algorithms:Within-feature inference`

Defining reference genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~

Single reference genome
+++++++++++++++++++++++

See examples in :ref:`Tutorial_cli:Reference gene(s) as targets`.

Multiple reference genomes
++++++++++++++++++++++++++
See also: :ref:`Parameters:Reference`, :ref:`Configuration:2-level lookup`, :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

Similar to ``--clusters`` and ``--indv``, MINORg accepts a lookup file for reference genomes using ``--reference-set`` and one or more reference genome alias using ``--reference``. See :ref:`Parameters:Reference` for a more comprehensive overview and :ref:`Configuration:reference` for lookup file format.

.. code-block:: bash
                
   $ minorg --directory ./example_xx_multiref \
            --indv ref --gene AT1G33560,AL1G47950.v2.1,Araha.3012s0003.v1.1 \
            --reference-set ./arabidopsis_genomes.txt --reference tair10,araly2,araha1

In the example above, MINORg will design gRNA for 3 highly conserved paralogues in 3 different species. Note that you should be careful that any gene IDs you use should either be unique across all reference genomes OR be shared only among your target genes. Otherwise, MINORg will treat any undesired genes with the same gene IDs as targets as well.

Non-standard genetic codes and mapping of non-standard attribute field names for multiple genomes should be specified in the lookup file passed to ``--reference-set``. See :ref:`Configuration:reference` for file format.

Non-standard reference
++++++++++++++++++++++

Non-standard genetic code
^^^^^^^^^^^^^^^^^^^^^^^^^

When using ``--domain``, users should ensure that the correct genetic code is specified, as MINORg has to first translate CDS into peptides for domain search using RPS-BLAST. The default genetic code is the Standard Code. Please refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for genetic code numbers and names.

.. code-block:: bash

   $ minorgpy --directory ./example_xx_geneticcode \
              --indv ref --gene gene-Q0275 \
              --assembly ./subset_ref_yeast_mt.fasta --annotation ./subset_ref_yeast_mt.gff \
              --domain 366140 --genetic-code 3

In the above example, the gene 'gene-Q0275' is a yeast mitochondrial gene, and ``--domain 366140`` specifies the PSSM-Id for the COX3 domain in the Cdd v3.18 RPS-BLAST database. The genetic code number for yeast mitochondrial code is '3'.

As a failsafe, MINORg does not terminate translated peptide sequences at the first stop codon. This ensures that any codons after an incorrectly translated premature stop codon will still be translated. Typically, a handful of mistranslated codons can still result in the correct RPS-BLAST domain hits, although hit scores may be slightly lower. Nevertheless, to ensure maximum accuracy, the correct genetic code is preferred.


Non-standard GFF3 attribute field names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See also: :ref:`Parameters:Attribute modification`

MINORg requires standard attribute field names in GFF3 files in order to properly map subfeatures to their parent features (e.g. map CDS to mRNA, and mRNA to gene). Non-standard field names should be mapped to standard ones using ``--attr-mod`` (for 'attribute modification').

.. code-block:: bash

   $ minorgpy --directory ./example_xx_attrmod \
              --indv ref --gene Os01t0100100 \
              --assembly ./subset_ref_irgsp.fasta --annotation ./subset_ref_irgsp.gff \
              --attr-mod 'mRNA:Parent=Locus_id'

The IRGSP 1.0 reference genome for rice (*Oryza sativa* subsp. Nipponbare) uses a non-standard attribute field name for mRNA entries in their GFF3 file. Instead of 'Parent', which is the standard name of the field used to map a feature to its parent feature, mRNA entries in the IRGSP 1.0 annotation use 'Locus_id'. See :ref:`Parameters:Attribute modification` for more details on how to format the input to ``--attr-mod``.
