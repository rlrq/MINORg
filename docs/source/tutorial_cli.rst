Tutorial (Command line)
=======================

In all the following tutorial, the current directory/working directory is presumed to contain all files in https://github.com/rlrq/MINORg/tree/master/examples. If you have not downloaded the files, please do so and navigate to the directory that contains them.

Note that all command line code snippets in the following tutorial are for **bash terminal**. You may have to adapt them according to your terminal.


Setting up the tutorial
~~~~~~~~~~~~~~~~~~~~~~~

To ensure that the examples in this tutorial work, please replace '/path/to' in the files 'arabidopsis_genomes.txt', 'athaliana_genomes.txt', and 'subset_genome_mapping.txt' with the full path to the directory containing the example files.


IMPT: Note on executables
~~~~~~~~~~~~~~~~~~~~~~~~~

See: :ref:`Parameters:blastn, rpsblast/rpsblast+, MAFFT`, :ref:`Parameters:BEDTools`

If blastn, rpsblast/rpsblast+, mafft, and/or bedtools is/are not in your command-search path, you will have to append one or more of the appropriate parameters below to your command to tell MINORg where they are.

``--blastn <path to blastn executable>``

* Applicable to: ``minorg`` (full programme), ``minorg seq`` (subcommand :ref:`Tutorial_cli:Subcommand \`\`seq\`\``), ``minorg filter`` (:ref:`Tutorial_cli:Subcommand \`\`filter\`\``)

``--rpsblast <path to rpsblast or rpsblast+ executable>``

* Applicable to: ``minorg`` (full programme), ``minorg seq`` (subcommand :ref:`Tutorial_cli:Subcommand \`\`seq\`\``), ``minorg grna`` (subcommand :ref:`Tutorial_cli:Subcommand \`\`grna\`\``), ``minorg filter`` (:ref:`Tutorial_cli:Subcommand \`\`filter\`\``)

``--mafft <path to MAFFT executable>``

* Applicable to: ``minorg`` (full programme), ``minorg seq`` (subcommand :ref:`Tutorial_cli:Subcommand \`\`seq\`\``), ``minorg grna`` (subcommand :ref:`Tutorial_cli:Subcommand \`\`grna\`\``), ``minorg filter`` (:ref:`Tutorial_cli:Subcommand \`\`filter\`\``)

``--bedtools <path to directory containing BEDTools executables>``

* Applicable to: ``minorg`` (full programme), ``minorg seq`` (subcommand :ref:`Tutorial_cli:Subcommand \`\`seq\`\``)

Defining target sequences
~~~~~~~~~~~~~~~~~~~~~~~~~

User-provided targets
+++++++++++++++++++++

Let us begin with the simplest MINORg execution:

.. code-block:: bash
   
   $ minorg --directory ./example_100_target --target ./sample_CDS.fasta
   Final gRNA sequence(s) have been written to minorg_gRNA_final.fasta
   Final gRNA sequence ID(s), gRNA sequence(s), and target(s) have been written to minorg_gRNA_final.map
   
   1 mutually exclusive gRNA set(s) requested. 1 set(s) found.
   Output files have been generated in /path/to/current/directory/example_100_target

The above combination of arguments tells MINORg to generate gRNA from targets in a user-provided FASTA file (``--target ./sample_CDS.fasta``) and to output files into a user-specified directory (``--directory ./example_100_target``). By default, MINORg generates 20 bp gRNA using NGG PAM.


Reference gene(s) as targets
++++++++++++++++++++++++++++

Both examples below specify a reference assembly (``--assembly ./subset_ref_TAIR10.fasta``) and annotation (``--annotation ./subset_ref_TAIR10.gff``) file, allowing MINORg to retrieve the gene sequence as target.

Single gene
^^^^^^^^^^^

.. code-block:: bash
   
   $ minorg --directory ./example_101_singlegene \
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
                
   $ minorg --directory ./example_102_multigene \
            --indv ref --gene AT5G45050,AT5G45060,AT5G45200,AT5G45210,AT5G45220,AT5G45230,AT5G45240,AT5G45250 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff

OR

.. code-block:: bash
                
   $ minorg --directory ./example_102_multigene \
            --indv ref --gene AT5G45050 --gene AT5G45060 --gene AT5G45200 \
            --gene AT5G45210 --gene AT5G45220 --gene AT5G45230 --gene AT5G45240 --gene AT5G45250 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff


Using ``--cluster``
*******************

See also: :ref:`Configuration:2-level lookup`, :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

MINORg can also accept preset combinations of genes using ``--cluster`` and ``--cluster-set``. ``--cluster-set`` accepts a tab-separated lookup file that maps alias(es) to a combinations of genes (see :ref:`Configuration:cluster` for format). ``--cluster`` is used to specify the alias of a combination of genes in that lookup file.

.. code-block:: bash
                
   $ minorg --directory ./example_103_cluster \
            --indv ref --cluster RPS6 --cluster-set ./subset_cluster_mapping.txt \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff

The above code snippet is effectively identical to the examples in :ref:`Tutorial_cli:Multiple genes`.

Like ``--gene``, multiple combinations of genes can be specified to ``--cluster``. However, unlike ``--gene``, each combination will be processed separately (i.e. minimum sets will be separately generated for each combination).

.. code-block:: bash
                
   $ minorg --directory ./example_103_cluster \
            --indv ref --cluster RPS6,TTR1 --cluster-set ./subset_cluster_mapping.txt \
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

If you have both genomic and CDS-only sequences of your target genes but not a GFF3 annotation file, MINORg can infer coding regions (CDS) for your target genes using ``--extend-gene`` and ``--extend-cds``. See :ref:`Parameters:Extended genome` for how to name your sequences to ensure proper mapping of CDS to genes.

.. code-block:: bash

   $ minorg --directory ./example_104_ext \
            --indv ref --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

Note that ``--extend-gene`` and ``--extend-cds`` effectively add new genes to the reference genome, so they can be used just like any reference gene. Therefore, they can also be used in combination with ``--query`` or ``--indv``.

Using ``--query``
*****************

See also: :ref:`Algorithms:Non-reference homologue inference`, :ref:`Parameters:Multi-argument (CLI)`

If you would like MINORg to infer homologues genes in non-reference genomes, you can use ``--query`` to specify the FASTA files of those non-reference genomes. You may provide multiple non-reference genomes by using ``--query`` multiple times.

.. code-block:: bash

   $ minorg --directory ./example_105_query \
            --query ./subset_9654.fasta --query ./subset_9655.fasta \
            --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

``--query`` can be used in combination with ``--indv``. For inference parameters, see :ref:`Algorithms:Non-reference homologue inference`.


Using ``--indv``
****************

See also: :ref:`Algorithms:Non-reference homologue inference`, :ref:`Configuration:2-level lookup`, :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

You can also use ``--indv`` to ask MINORg to infer homologues genes in non-reference genomes. Similar to ``--clusters``, MINORg accepts a lookup file for non-reference genomes using ``--genome-set`` (see :ref:`Configuration:genome` for format) and one or more non-reference genome alias using ``--indv``.

.. code-block:: bash

   $ minorg --directory ./example_106_indv \
            --indv 9654,9655 --genome-set ./subset_genome_mapping.txt \
            --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

The above code snippet is effectively identical to the example in :ref:`Tutorial_cli:Using \`\`--query\`\``.

``--indv`` can be used in combination with ``--query``. For inference parameters, see :ref:`Algorithms:Non-reference homologue inference`.


Domain as targets
+++++++++++++++++

MINORg allows users to specify the identifier of an RPS-BLAST position-specific scoring matrix (PSSM-Id) to further restrict the target sequence to a given domain associated with the PSSM-Id. This could be particularly useful when designing gRNA for genes that do not share conserved domain structures but do share a domain that you wish to knock out. ``--domain`` can also be used with ``--query`` or ``--indv``.

Local database
^^^^^^^^^^^^^^

.. code-block:: bash

   $ minorg --directory ./example_107_domain \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --rpsblast /path/to/rpsblast/executable --db /path/to/rpsblast/db \
            --domain 214815

In the above example, gRNA will be generated for the WRKY domain (PSSM-Id 214815 as of CDD database v3.18) of the gene AT5G45050. Users are responsible for providing the PSSM-Id of a domain that exists in the gene. Unlike other examples, the database (``--db``) is not provided as part of the example files. You will have to download it yourself. See :ref:`Parameters:RPS-BLAST local database` for more information.

Remote database
^^^^^^^^^^^^^^^

While it is in theory possible to use the remote CDD database & servers instead of local ones, the ``--remote`` option for the 'rpsblast'/'rpsblast+' command from the BLAST+ package has never worked for me. In any case, if your version of local rpsblast is able to access the remote database, you can use ``--remote-rps`` instead of ``--db /path/to/rpsblast/db``.

.. code-block:: bash

   $ minorg --directory ./example_107_domain \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --rpsblast /path/to/rpsblast/executable --remote-rps \
            --domain 214815
..
   Feature as targets
   ++++++++++++++++++

   You may specify gene features to restrict gRNA to. By default, MINORg generates gRNA in coding regions (CDS). However, so long as a feature type is valid in the GFF3 annotation file provided to MINORg, gRNA can theoretically be designed for any feature type.

   .. code-block:: bash

      $ minorg --directory ./example_108_feature \
               --indv ref --gene AT5G45050 \
               --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
               --feature three_prime_UTR

   The above example will generate gRNA in the 100 bp 3' UTR of AT5G45050.


Defining gRNA
~~~~~~~~~~~~~

See also: :ref:`Parameters:PAM`

By default, MINORg generates 20 bp gRNA using SpCas9's NGG PAM. You may specify other gRNA length using ``--length`` and other PAM using ``--pam``.

.. code-block:: bash

   $ minorg --directory ./example_108_grna \
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

   $ minorg --directory ./example_109_gc \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --gc-min 0.2 --gc-max 0.8

In the above example, MINORg will exclude gRNA with less than 20% (``--gc-min 0.2``) or greater than 80% (``--gc-max 0.8``) GC content. By default, minimum GC content is 30% and maximum is 70%.


Filter by off-target
++++++++++++++++++++
See: :ref:`Algorithms:Off-target assessment`

Using total mismatch/gap/unaligned
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See: :ref:`Algorithms:Total mismatch/gap/unaligned`

Thresholds for total number of mismatches or gaps (and unaligned positions) required for an off-target gRNA hit to be considered non-problematic are controlled by ``--ot-mismatch`` and ``--ot-gap`` respectively. See :ref:`Algorithms:Total mismatch/gap/unaligned` for more.

.. code-block:: bash

   $ minorg --directory ./example_110_ot_ref \
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
  * Note that any AT5G45050 homologues in these two genomes will NOT be masked. This means that only gRNA that do not target any AT5G45050 homologues in these two non-reference genomes will pass this off-target check.
    * To mask homologues in these genomes, you will need to provide a FASTA file containing the sequences of their homologues using ``--mask <FASTA>``. You may use subcommand ``seq`` (see :ref:`Tutorial_cli:Subcommand \`\`seq\`\``) to identify these homologues.

``--ot-gap`` and ``--ot-mismatch`` control the minimum number of gaps or mismatches off-target gRNA hits must have to be considered non-problematic; any gRNA with at least one problematic gRNA hit will be excluded. See :ref:`Algorithms:Off-target assessment` for more on the off-target assessment algorithm.

In the case above, ``--screen-reference`` is actually redundant as the genome from which targets are obtained (which, because of ``--indv ref``, is the reference genome) are automatically included for background check. However, in the example below, when the targets are from **non-reference genomes**, the reference genome is not automatically included for off-target assessment and thus ``--screen-reference`` is NOT redundant. Additionally, do note that the genes passed to ``--gene`` are masked in the reference genome, such that any gRNA hits to them are NOT considered off-target and will NOT be excluded.

.. code-block:: bash

   $ minorg --directory ./example_111_ot_nonref \
            --indv 9654 --genome-set ./subset_genome_mapping.txt \
            --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --screen-ref --background ./subset_ref_Araly2.fasta --background ./subset_ref_Araha1.fasta \
            --ot-indv 9655 \
            --ot-gap 2 --ot-mismatch 2
            
Using position-specific mismatch/gap/unaligned
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See: :ref:`Algorithms:Position-specific mismatch/gap/unaligned`

Finer control of off-target definition can be achieved using :attr:`~minorg.MINORg.MINORg.ot_pattern`, which allows users to provide a pattern that specifies different thresholds for different positions along a gRNA. Unlike ``--ot-mismatch`` and ``--ot-gap``, which specify the **LOWER-bound of NON-problematic** hits, ``--ot-pattern`` specifies **UPPER-bound of PROBLEMATIC** hits. By default, unaligned positions will be treated as mismatches, but this behaviour can be altered by raising ``--ot-unaligned-as-mismatch-unset``. See :ref:`Parameters:Off-target pattern` for how to build an off-target pattern, and :ref:`Algorithms:Position-specific mismatch/gap/unaligned` for more on how unaligned positions can be counted.

When ``--ot-pattern`` is specified, ``--ot-mismatch`` and ``--ot-gap`` will be ignored.

The following example is identical to the first in :ref:`Tutorial_py:Using total mismatch/gap/unaligned`, except ``--ot-mismatch`` and ``--ot-gap`` are replaced with ``--ot-pattern``.

.. code-block:: bash

   $ minorg --directory ./example_112_ot_ref_pattern \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --screen-reference \
            --background ./subset_ref_Araly2.fasta --background ./subset_ref_Araha1.fasta \
            --ot-indv 9654,9655 --genome-set ./subset_genome_mapping.txt \
            --ot-pattern '2mg-5-,0mg4'

In the above example, ``--ot-pattern '0mg-4,2mg-5-'`` means that MINORg will discard any gRNA with at least one off-target hit where:

* There are no mismatches or gaps between positions -4 and -1, and there are no more than 2 mismatches or gaps from position -5 to the 5' end.

PAM-less off-target check
^^^^^^^^^^^^^^^^^^^^^^^^^

By default, MINORg does NOT check for the presence of PAM sites next to potential off-target hits. You may override this behaviour using ``--ot-pam``. This tells MINORg to mark off-target hits that fail the ``--ot-gap`` or ``--ot-mismatch`` thresholds (or match ``--ot-pattern``) as problematic ONLY IF there is a PAM site nearby.

.. code-block:: bash

   $ minorg --directory ./example_113_ot_pamless \
            --indv 9654 --genome-set ./subset_genome_mapping.txt \
            --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --screen-ref --background ./subset_ref_Araly2.fasta --background ./subset_ref_Araha1.fasta \
            --ot-indv 9655 \
            --ot-gap 2 --ot-mismatch 2 \
            --ot-pam

Skip off-target check
^^^^^^^^^^^^^^^^^^^^^

To skip off-target check entirely, use ``--skip-bg-check``.

.. code-block:: bash

   $ minorg --directory ./example_114_skipbgcheck \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --skip-bg-check


Filter by feature
+++++++++++++++++

See: :ref:`Algorithms:Within-feature inference`

By default, when ``--gene`` is used, MINORg restricts gRNA to coding regions (CDS). For more on how MINORg does this for inferred, unannotated homologues, see :ref:`Algorithms:Within-feature inference`. You may change the feature type in which to design gRNA using ``--feature``. See column 3 of your GFF3 file for valid feature types (see https://en.wikipedia.org/wiki/General_feature_format for more on GFF file format).

.. code-block:: bash
                
   $ minorg --directory ./example_115_withinfeature \
            --indv ref --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --feature three_prime_UTR

Generating minimum gRNA set(s)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Number of sets
++++++++++++++

By default, MINORg outputs a single gRNA set covering all targets. You may request more (mutually exclusive) sets using ``--set``.

.. code-block:: bash
                
   $ minorg --directory ./example_116_set \
            --indv ref --cluster RPS6 --cluster-set ./subset_cluster_mapping.txt \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --set 5


Prioritise non-redundancy
+++++++++++++++++++++++++

By default, MINORg selects gRNA for sets using these criteria in decreasing order of priority:

#. Coverage (of as yet uncovered targets)
#. Proximity to 5' end
#. Non-redundancy

Proximity is only assessed when there is a tie for coverage, and non-redundancy when there is a tie for both coverage and proximity. You may flip the priority of proximity and non-redundancy using ``--prioritise-nr`` if you prefer to minimise multiple edits in a single target when using a single set of gRNA. (The priority of coverage is NOT modifiable.)

.. code-block:: bash

   $ minorg --directory ./example_117_nr \
            --indv ref --cluster RPS6 --cluster-set ./subset_cluster_mapping.txt \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --prioritise-nr

Excluding gRNA
++++++++++++++

You may specify gRNA sequences to exclude from any final gRNA set using ``--exclude``.

.. code-block:: bash

   $ minorg --directory ./example_118_exclude \
            --indv ref --cluster RPS6 --cluster-set ./subset_cluster_mapping.txt \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --exclude ./sample_exclude_RPS6.fasta

The gRNA names in the file passed to ``--exclude`` do not matter. Only the sequences are used when determining whether to exclude a gRNA.

Accepting unknown checks
++++++++++++++++++++++++

Sometimes, not all filtering checks (GC, background, and feature) are set for all sequences. This is not an issue if you use the full programme (i.e. ``minorg <arguments>``), but may be relevant if you are re-generating sets using the 'minimumset' subcommand (i.e. ``minorg minimumset <arguments>``) with a modified mapping file OR a mapping file from the 'filter' subcommand where not all filters have been applied.

Let us take a look at 'sample_custom_check.map', where we've added a custom check called 'my_custom_check' in the last column::

  gRNA id	gRNA sequence	target id	target sense	gRNA strand	start	end	group	background	GC	feature	my_custom_check
  gRNA_001	CTTCATCTTCTTCTCGAAAT	targetA	NA	+	8	27	1	pass	pass	NA	pass
  gRNA_001	CTTCATCTTCTTCTCGAAAT	targetB	NA	+	80	99	1	pass	pass	NA	pass
  gRNA_002	GATGTTTTCTTGAGCTTCAG	targetA	NA	+	37	56	1	pass	pass	NA	NA
  gRNA_002	GATGTTTTCTTGAGCTTCAG	targetB	NA	+	286	305	1	pass	pass	NA	pass
  gRNA_002	GATGTTTTCTTGAGCTTCAG	targetC	NA	+	109	128	1	pass	pass	NA	fail
  gRNA_002	GATGTTTTCTTGAGCTTCAG	targetD	NA	+	110	129	1	pass	pass	NA	fail
  gRNA_003	ATGTTTTCTTGAGCTTCAGA	targetB	NA	+	38	57	1	pass	pass	NA	NA
  gRNA_003	ATGTTTTCTTGAGCTTCAGA	targetC	NA	+	287	306	1	pass	pass	NA	pass
  gRNA_003	ATGTTTTCTTGAGCTTCAGA	targetD	NA	+	110	129	1	pass	pass	NA	pass

There are three possible values for check status: 'pass', 'fail', and 'NA'.

An invalid/unset check is an 'NA'. If a check is unset for all entries (as is the case with the check 'feature' here), it will be ignored (i.e. the check is treated as 'pass' for all entries). However, when a check has been set for some entries but not others (as is the case with the 'my_custom_check' check here), MINORg will treat invalid/unset checks as 'fail' by default. This is because there isn't enough information on whether this constitutes a pass or fail for the check, and MINORg prefers to be conservative when outputting gRNA. You may override this behaviour using the ``--accept-invalid``. By doing so, MINORg will treat 'NA' as 'pass' for all checks.

.. code-block:: bash

   $ minorg minimumset --directory ./example_119_acceptinvalid \
                       --map ./sample_custom_check.map \
                       --accept-invalid

                       
Manually approve gRNA sets
++++++++++++++++++++++++++

You may opt to manually inspect each gRNA set before MINORg write them to file using the ``--manual`` flag.

.. code-block:: bash

   $ minorg --directory ./example_120_manual --target ./sample_CDS.fasta
            --manual
   	ID	sequence (Set 1)
   	gRNA_001	GGAATACAAGAGATTATCGA
   Hit 'x' to continue if you are satisfied with these sequences. Otherwise, enter the sequence ID or sequence of an undesirable gRNA (case-sensitive) and hit the return key to update this list: x
   Final gRNA sequence(s) have been written to minorg_gRNA_final.fasta
   Final gRNA sequence ID(s), gRNA sequence(s), and target(s) have been written to minorg_gRNA_final.map
   
   1 mutually exclusive gRNA set(s) requested. 1 set(s) found.
   Output files have been generated in /path/to/current/directory/example_119_manual


Subcommands
~~~~~~~~~~~

MINORg comprises of four main steps:

#. Target sequence identification
#. Candidate gRNA generation
#. gRNA filtering
#. Minimum gRNA set generation

As users may only wish to execute a subset of these steps instead of the full programme, MINORg also provides four subcommands corresponding to these four steps:

#. ``seq``
#. ``grna``
#. ``filter``
#. ``minimumset``

The subcommands may be useful if you already have a preferred off-target/on-target assessment software. In this case, you may execute subcommands ``seq`` and ``grna``, submit the gRNA output by MINORg for off-target/on-target assessment, update the .map file output by MINORg with the status of each gRNA for that off-target/on-target assessment, and execute ``minimumset`` to obtain a desired number of minimum gRNA sets.
   
Subcommand ``seq``
++++++++++++++++++

The ``seq`` subcommand identifies target sequences, whether by extracting them from a reference genome or inferring homologues in unannotated genomes. All parameters described in :ref:`Tutorial_cli:Defining target sequences` (except ``--target``) and :ref:`Tutorial_cli:Defining reference genomes` apply.

This step will output target sequences into a file ending with '_targets.fasta'.

To use this subcommand, simply replace the command ``minorg`` with ``minorg seq``.

.. code-block:: bash

   $ minorg seq --directory ./example_121_subcmdseq \
                --query ./subset_9654.fasta --query ./subset_9655.fasta \
                --gene AT1G10920 \
                --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

Subcommand ``grna``
+++++++++++++++++++

The ``grna`` subcommand generates gRNA within target sequences. It incorporates parts of the ``seq`` and ``filter`` subcommands in order to provide rudimentary filtering for gRNA within specific GFF3 features (e.g. CDS) for reference genes as well as by GC content. All parameters described in :ref:`Tutorial_cli:Defining target sequences` (except those related to homology discovery in unannotated genomes such as ``query``, ``indv``, and ``genome-set``), :ref:`Tutorial_cli:Defining reference genomes`, :ref:`Tutorial_cli:Defining gRNA`, and :ref:`Tutorial_cli:Filter by GC content`, and :ref:`Tutorial_cli:Filter by feature` apply.

Unlike the full programme or the ``seq`` subcommand, however, ``--indv ref`` is not necessary to specify reference genes as target. As this subcommand does not support homologue discovery, if ``--gene`` or ``--cluster`` is used, ``--indv ref`` will automatically be filled since unannotated genomes are not allowed.

This step will output target sequences into a file ending with '_targets.fasta' if ``--target`` was not used. gRNA sequences will be written into files ending with '_gRNA_all.fasta' (for all candidate gRNA) and '_gRNA_pass.fasta' (for candidate gRNA that pass GC and feature checks). A file ending with '_gRNA_all.map' that maps gRNA to their targets will also be generated. You may optionally specify the location of the FASTA and .map output files using:

* ``--out-fasta``: path to output file that originally ends with '_gRNA_all.fasta'
* ``--out-pass``: path to output file that originally ends with '_gRNA_pass.fasta'
* ``--out-map``: path to output file that originally ends with '_gRNA_all.map'

To use this subcommand, simply replace the command ``minorg`` with ``minorg grna``.

.. code-block:: bash

   $ minorg grna --directory ./example_122_subcmdgrna \
                 --cluster RPS6 --cluster-set ./subset_cluster_mapping.txt \
                 --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
                 --length 19 --pam Cas12a \
                 --feature three_prime_UTR \
                 --gc-min 0.2 --gc-max 0.8 \
                 --out-map ./example_120.map

Subcommand ``filter``
+++++++++++++++++++++

The ``filter`` subcommand takes in a compulsory MINORg .map file (``--map``) and rewrites some/all checks. You should specify the checks you wish to re-assess using some combination of ``--gc-check``, ``--background-check``, and/or ``--feature-check`` flags OR ``--check-all`` to raise all three flags. For in-place modification of the .map file, use ``--in-place``. Otherwise, MINORg will write a new file using the default naming format of '<prefix>_gRNA_all.map' (this may still overwrite the original file if the directory and prefix are identical to what was used to generate the original file).

gRNA sequences will be written into files ending with '_gRNA_all.fasta' (for all candidate gRNA) and '_gRNA_pass.fasta' (for candidate gRNA that pass the updated checks). A file ending with '_gRNA_all.map' that maps gRNA to their targets will also be generated with the updated check statuses. As with subcommand ``grna``, you may optionally specify the location of the FASTA and .map output files using:

* ``--out-fasta``: path to output file that originally ends with '_gRNA_all.fasta'
* ``--out-pass``: path to output file that originally ends with '_gRNA_pass.fasta'
* ``--out-map``: path to output file that originally ends with '_gRNA_all.map'

To use this subcommand, simply replace the command ``minorg`` with ``minorg filter``.

In all cases, you may rename the gRNA using ``--rename <FASTA>``, where the FASTA file contains the gRNA sequences you wish to rename with sequence IDs of the names you wish to rename them to.

GC check
^^^^^^^^

All parameters described in :ref:`Tutorial_cli:Filter by GC content` apply.

.. code-block:: bash

   $ minorg filter --directory ./example_123_subcmdfilter_gc \
                   --map ./sample_custom_check.map \
                   --gc-check --gc-min 0.2 --gc-max 0.8

Background check
^^^^^^^^^^^^^^^^

All parameters described in :ref:`Tutorial_cli:Filter by off-target` apply. Additionally, you should supply target sequences using ``--target`` so that MINORg can mask them (this tells MINORg that any gRNA hits to them is in fact on-target and NOT off-target). Any additional sequences to be masked may be provided using ``--mask <FASTA>``. If you are using ``--screen-ref`` to include reference genome(s) (see :ref:`Tutorial_cli:Multiple reference genomes` for how to specify multiple reference genomes) in the off-target screen, you may specify reference genes to be masked using ``--mask-gene`` or ``--mask-cluster`` (unlike ``--cluster``, all clusters passed to ``--mask-cluster`` will be processed simultaneously; i.e. there will not be separate executions for each cluster).

Let us first generate a .map file for filtering.

.. code-block:: bash

   $ minorg --directory ./example_124_subcmdfilter_bg_pt1 \
            --indv 9654,9655 --genome-set ./subset_genome_mapping.txt \
            --cluster RPS6 --cluster-set ./subset_cluster_mapping.txt \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
            --skip-bg-check

In the code above, we skipped off-target check by raising the ``--skip-bg-check`` flag. But we've changed out mind and would like to screen the reference genome and the non-reference genomes that these targets are from AND we don't want our gRNA to be able to target any genes in 'subset_9944.fasta' and 'subset_9947'. We can do that using the ``filter`` subcommand.

.. code-block:: bash

   $ minorg filter --directory ./example_124_subcmdfilter_bg_pt2 \
                   --map ./example_123_subcmdfilter_bg_pt1/minorg_RPS6/minorg_RPS6_gRNA_all.map \
                   --background-check \
                   --target ./example_123_subcmdfilter_bg_pt1/minorg_RPS6/minorg_RPS6_gene_targets.fasta \
                   --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
                   --screen-ref \
                   --mask-cluster RPS6 --cluster-set ./subset_cluster_mapping.txt \
                   --ot-indv 9654,9655,9944,9947 --genome-set ./subset_genome_mapping.txt

The above code may be a little unwieldy. However, if the target identification step of MINORg takes a while to run (for example when the genome files are large and take forever to process), you may prefer not to re-run the full MINORg programme with updated parameters and instead use the ``filter`` subcommand on files that have already been generated. You should then use the ``minimumset`` subcommand (see :ref:`Tutorial_cli:Subcommand \`\`minimumset\`\``) to regenerate minimum sets using your update .map file.

Feature check
^^^^^^^^^^^^^

All parameters described in :ref:`Tutorial_cli:Filter by feature` apply. Additionally, you will need to provide a FASTA file of target sequences (using ``--target <FASTA>``), reference genome(s) (see :ref:`Tutorial_cli:Defining reference genomes`), and genes (using ``--gene <gene IDs>`` or ``--cluster <cluster alias>``). The specified reference gene(s) will be extracted from the reference genome(s) and aligned with target sequence(s) in order for MINORg to infer feature boundaries in target sequence(s). See :ref:`Algorithms:Within-feature inference` for the algorithm of how feature boundaries are inferred.

Do note that unlike the full programme or the ``seq`` subcommand, all clusters passed to ``--cluster`` will be processed simultaneously (i.e. there will not be separate executions for each cluster).

Let us first generate a .map file for filtering.

.. code-block:: bash

   $ minorg --directory ./example_125_subcmdfilter_feature_pt1 \
            --indv 9654,9655 --genome-set ./subset_genome_mapping.txt \
            --gene AT5G45050 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff

By default, MINORg sets the desired feature to 'CDS'. You can re-assess and overwrite the 'feature' check in the .map file to only allow gRNA in the 3' UTR using ``minorg filter`` with the ``--feature-check`` flag raised.

.. code-block:: bash

   $ minorg filter --directory ./example_125_subcmdfilter_feature_pt2 \
                   --map ./example_124_subcmdfilter_feature_pt1/minorg/minorg_gRNA_all.map \
                   --feature-check \
                   --target ./example_124_subcmdfilter_feature_pt1/minorg/minorg_gene_targets.fasta \
                   --gene AT5G45050 \
                   --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
                   --feature three_prime_UTR

Combination of checks
^^^^^^^^^^^^^^^^^^^^^

You can execute all checks (or some combination of them) in a single ``minorg filter`` command as well, if you wish. Just make sure that you raise the appropriate flag(s).

To use some combination of checks, simply raise the relevant flags (``--gc-check``, ``--background-check``, and/or ``--feature-check``). In the example below, we filter the gRNA generated by full MINORg execution in :ref:`Tutorial_cli:Feature check` by both GC content (``--gc-check``) as well as gene feature (``--feature-check``).

.. code-block:: bash

   $ minorg filter --directory ./example_126_subcmdfilter_gcfeature \
                   --map ./example_124_subcmdfilter_feature_pt1/minorg/minorg_gRNA_all.map \
                   --gc-check \
                   --gc-min 0.2 --gc-max 0.8 \
                   --feature-check \
                   --target ./example_124_subcmdfilter_feature_pt1/minorg/minorg_gene_targets.fasta \
                   --gene AT5G45050 \
                   --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
                   --feature three_prime_UTR

To execute all checks, use ``--check-all``. In the example below, we filter the gRNA generated by full MINORg execution in :ref:`Tutorial_cli:Feature check` by all checks.

.. code-block:: bash

   $ minorg filter --directory ./example_127_subcmdfilter_all \
                   --map ./example_124_subcmdfilter_feature_pt1/minorg/minorg_gRNA_all.map \
                   --check-all \
                   --gc-min 0.2 --gc-max 0.8 \ ## GC
                   --target ./example_124_subcmdfilter_feature_pt1/minorg/minorg_gene_targets.fasta \ ## feature
                   --gene AT5G45050 \
                   --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \
                   --feature three_prime_UTR \
                   --screen-ref --mask-gene AT5G45050 ## off-target

Subcommand ``minimumset``
+++++++++++++++++++++++++

The ``minimumset`` subcommand generates mutually exclusive minimum set(s) of gRNA, where each set is capable of covering all targets. All parameters described in :ref:`Tutorial_cli:Generating minimum gRNA set(s)` apply.

This step will write final gRNA sequences into a file ending with '_gRNA_final.fasta'. A file ending with '_gRNA_final.map' that maps gRNA to their targets will also be generated. You may optionally specify the location of the FASTA and .map output files using:

* ``--out-fasta``: path to output file that originally ends with '_gRNA_final.fasta'
* ``--out-map``: path to output file that originally ends with '_gRNA_final.map'

**NOTE:** Unlike subcommands ``grna`` and ``filter``, ``--out-fasta`` and ``--out-map`` are used to specify output files for **FINAL** gRNA sets, not all candidate gRNA.

To use this subcommand, simply replace the command ``minorg`` with ``minorg grna``.

.. code-block:: bash

   $ minorg minimumset --directory ./example_128_subcmdminimumset
                       --map ./example_105_query/minorg/minorg_gRNA_all.map \
                       --target ./example_105_query/minorg/minorg_gene_targets.fasta \
                       --set 5 --manual --prioritise-nr

In order for MINORg to better assess a gRNA's proximity to the 5' end (of hopefully sense strand) of a target in the event a tie-breaker is necessary, it is strongly suggested that target sequences be provided using ``--target <FASTA>`` so MINORg knows how long a target sequence is. This is especially so if the target sequences are antisense ones (you can check this using the .map file) generated by MINORg's inferences of homologues in unannotated genomes.

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
                
   $ minorg --directory ./example_129_multiref \
            --indv ref --gene AT1G33560,AL1G47950.v2.1,Araha.3012s0003.v1.1 \
            --reference tair10,araly2,araha1 --reference-set ./arabidopsis_genomes.txt

In the example above, MINORg will design gRNA for 3 highly conserved paralogues in 3 different species. Note that you should be careful that any gene IDs you use should either be unique across all reference genomes OR be shared only among your target genes. Otherwise, MINORg will treat any undesired genes with the same gene IDs as targets as well.

Non-standard genetic codes and mapping of non-standard attribute field names for multiple genomes should be specified in the lookup file passed to ``--reference-set``. See :ref:`Configuration:reference` for file format.

Non-standard reference
++++++++++++++++++++++

Non-standard genetic code
^^^^^^^^^^^^^^^^^^^^^^^^^

When using ``--domain``, users should ensure that the correct genetic code is specified, as MINORg has to first translate CDS into peptides for domain search using RPS-BLAST. The default genetic code is the Standard Code. Please refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for genetic code numbers and names.

.. code-block:: bash

   $ minorgpy --directory ./example_129_geneticcode \
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

   $ minorgpy --directory ./example_130_attrmod \
              --indv ref --gene Os01t0100100 \
              --assembly ./subset_ref_irgsp.fasta --annotation ./subset_ref_irgsp.gff \
              --attr-mod 'mRNA:Parent=Locus_id'

The IRGSP 1.0 reference genome for rice (*Oryza sativa* subsp. Nipponbare) uses a non-standard attribute field name for mRNA entries in their GFF3 file. Instead of 'Parent', which is the standard name of the field used to map a feature to its parent feature, mRNA entries in the IRGSP 1.0 annotation use 'Locus_id'. See :ref:`Parameters:Attribute modification` for more details on how to format the input to ``--attr-mod``.

Multithreading
~~~~~~~~~~~~~~

MINORg supports multi-threading in order to process files in parallel. Any excess threads may also be used for BLAST. This is most useful when you are querying multiple genomes (whether using ``--query`` or ``--indv``), have multiple reference genomes (``--reference``), or multiple background sequences (``--background``).

To run MINORg with parallel processing, use ``--thread <number of threads>``.

.. code-block:: bash

   $ minorg --directory ./example_131_thread \
            --query ./subset_9654.fasta --query ./subset_9655.fasta \
            --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta \
            --thread 2
