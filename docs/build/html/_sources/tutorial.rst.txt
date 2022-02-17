Tutorial
========

In all the following tutorial, the current directory/working directory is presumed to contain all files in https://github.com/rlrq/MINORg/tree/master/examples. If you have not downloaded the files, please do so and navigate to the directory that contains them.


Command line
------------

Note that all command line code snippets in the following tutorial are for **bash terminal**. You may have to adapt them according to your operating system.

Defining target sequences
~~~~~~~~~~~~~~~~~~~~~~~~~

User-provided targets
+++++++++++++++++++++

Let us begin with the simplest MINORg execution:

.. code-block:: bash
   
   $ minorg --directory ./example00 --target ./sample_CDS.fasta
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
   
   $ minorg --directory ./example_01 \
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
                
   $ minorg --directory ./example_02 \
            --indv ref --gene AT5G45050,AT5G45060,AT5G45200,AT5G45210,AT5G45220,AT5G45230,AT5G45240,AT5G45250 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff

OR

.. code-block:: bash
                
   $ minorg --directory ./example_02 \
            --indv ref --gene AT5G45050 --gene AT5G45060 --gene AT5G45200 \
            --gene AT5G45210 --gene AT5G45220 --gene AT5G45230 --gene AT5G45240 --gene AT5G45250 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff


Using ``--cluster``
*******************
See also: :ref:`Configuration:2-level lookup`, :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

MINORg can also accept preset combinations of genes using ``--cluster`` and ``--cluster-set``. ``--cluster-set`` accepts a tab-separated lookup file that maps alias(es) to a combinations of genes (see :ref:`Configuration:cluster` for format). ``--cluster`` is used to specify the alias of a combination of genes in that lookup file.

.. code-block:: bash
                
   $ minorg --directory ./example_03 \
            --indv ref --cluster RPS6 --cluster-set ./subset_cluster_mapping.txt
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff

The above code snippet is effectively identical to the examples in :ref:`Tutorial:Multiple genes`.

Like ``--gene``, multiple combinations of genes can be specified to ``--cluster``. However, unlike ``--gene``, each combination will be processed separately (i.e. minimum sets will be separately generated for each combination).

.. code-block:: bash
                
   $ minorg --directory ./example_03 \
            --indv ref --cluster RPS6,TTR1 --cluster-set ./subset_cluster_mapping.txt
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff


Non-reference genes as targets
++++++++++++++++++++++++++++++

Annotated genes
^^^^^^^^^^^^^^^

If your target genes have been annotated in their non-reference genomes (i.e. you have a GFF3 file containing annotations of your targets), you can use :ref:`Tutorial:Reference gene(s) as targets` if you have a single non-reference genome, or :ref:`Tutorial:Multiple reference genomes` if you have multiple non-reference genomes. In either case, you may treat your non-reference genome the same way you would a reference genome.


Unannotated genes
^^^^^^^^^^^^^^^^^

Using ``--extend-gene`` and ``--extend-cds``
********************************************
See also: :ref:`Parameters:Extended genome`

If you have both genomic and CDS-only sequences of your target genes but not a GFF3 annotation file, MINORg can infer coding regions (CDS) for your target genes using ``--extend-gene`` and ``--extend-cds``.

.. code-block:: bash

   $ minorg --directory ./example_04 \
            --indv ref --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

Note that ``--extend-gene`` and ``--extend-cds`` effectively add new genes to the reference genome, so they can be used just like any reference gene. Therefore, they can also be used in combination with ``--query`` or ``--indv``.

Using ``--query``
*****************
See also: :ref:`Parameters:Multi-argument (CLI)`

If you would like MINORg to infer homologues genes in non-reference genomes, you can use ``--query`` to specify the FASTA files of those non-reference genomes. You may provide multiple non-reference genomes by using ``--query`` multiple times.

.. code-block:: bash

   $ minorg --directory ./example_05 \
            --query ./subset_9654.fasta --query ./subset_9655.fasta \
            --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

``--query`` can be used in combination with ``--indv``.

Using ``--indv``
****************
See also: :ref:`Configuration:2-level lookup`, :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

You can also use ``--indv`` to ask MINORg to infer homologues genes in non-reference genomes. Similar to ``--clusters``, MINORg accepts a lookup file for non-reference genomes using ``--genome-set`` (see :ref:`Configuration:genome` for format) and one or more non-reference genome alias using ``--indv``.

.. code-block:: bash

   $ minorg --directory ./example_05 \
            --indv 9654,9655 --genome-set ./subset_genome_mapping.txt \
            --gene AT1G10920 \
            --extend-gene ./sample_gene.fasta --extend-cds ./sample_CDS.fasta

The above code snippet is effectively identical to the example in :ref:`Tutorial:Using \`\`--query\`\``.

``--indv`` can be used in combination with ``--query``.

Multiple reference genomes
++++++++++++++++++++++++++
See also: :ref:`Parameters:Reference`, :ref:`Configuration:2-level lookup`, :ref:`Parameters:Comma-separated (CLI)`, :ref:`Parameters:Multi-argument (CLI)`

Similar to ``--clusters`` and ``--indv``, MINORg accepts a lookup file for reference genomes using ``--reference-set`` and one or more reference genome alias using ``--reference``. See :ref:`Parameters:Reference` for a more comprehensive overview.

.. code-block:: bash
                
   $ minorg --directory ./example07 \
            --indv ref --gene AT1G33560,AL1G47950.v2.1,Araha.3012s0003.v1.1 \
            --reference-set ./arabidopsis_genomes.txt --reference tair10,araly2,araha1

In the example above, MINORg will design gRNA for 3 highly conserved paralogues in 3 different species. Note that you should be careful that any gene IDs you use should either be unique across all reference genomes OR be shared only among your target genes. Otherwise, MINORg will treat any undesired genes with the same gene IDs as targets as well.

Non-standard reference
++++++++++++++++++++++

Non-standard genetic code
^^^^^^^^^^^^^^^^^^^^^^^^^

When using ``--domain``, users should ensure that the correct genetic code is specified, as MINORg has to first translate CDS into peptides for domain search using RPS-BLAST. The default genetic code is the Standard Code. Please refer to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for genetic code numbers and names.

.. code-block:: bash

   $ minorgpy --directory ./example08 \
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

   $ minorgpy --directory ./example09 \
              --indv ref --gene Os01t0100100 \
              --assembly ./subset_ref_irgsp.fasta --annotation ./subset_ref_irgsp.gff \
              --attr-mod 'mRNA:Parent=Locus_id'

The IRGSP 1.0 reference genome for rice (*Oryza sativa* subsp. Nipponbare) uses a non-standard attribute field name for mRNA entries in their GFF3 file. Instead of 'Parent', which is the standard name of the field used to map a feature to its parent feature, mRNA entries in the IRGSP 1.0 annotation use 'Locus_id'. See :ref:`Parameters:Attribute modification` for more details on how to format the input to ``--attr-mod``.


Python package
--------------

Getting started
~~~~~~~~~~~~~~~

To begin, import the :class:`~minorg.MINORg.MINORg` class.

>>> from minorg.MINORg import MINORg

To create a MINORg object:

>>> my_minorg = MINORg(directory = '/path/to/output/directory', prefix = 'prefix')

Both ``directory`` and ``prefix`` are optional. If not provided, they will default to the current directory and 'minorg' respectively.

If you wish to use the default values specified in a config file, use this instead:

>>> my_minorg = MINORg(config = '/path/to/config.ini', directory = '/path/to/output/directory', prefix = 'prefix')

You may now set your parameters using the attributes of your :class:`~minorg.MINORg.MINORg` object. For a table listing the equivalent CLI arguments and :class:`~minorg.MINORg.MINORg` attributes, see :ref:`Parameters:CLI vs Python`. For example, you can specify executables as such:

>>> my_minorg.blastn = '/path/to/blastn/executable'
>>> my_minorg.rpsblast = '/path/to/rpsblast/executable'
>>> my_minorg.mafft = '/path/to/mafft/executable'

Note that, unlike the command line, the Python package does not support aliases even if the config file has been set up appropriately for command line executions. Therefore, there are no true equivalents to ``--cluster``, ``--indv``, or ``--reference``.

To specify cluster genes (analogous to ``--cluster`` and ``--gene``):

>>> my_minorg.cluster = 'RPS6' ## incorrect; this attribute does not exist; does not throw error now but will cause problems later
>>> my_minorg.genes = ['AT5G46260','AT5G46270','AT5G46450','AT5G46470','AT5G46490','AT5G46510','AT5G46520'] ## correct

To specify query FASTA files (analogous to ``--indv`` and ``--query``):

>>> my_minorg.indv = '9654,9655' ## incorrect; this attribute does not exist; does not throw error now but will cause problems later
>>> my_minorg.query = {'9654': '/path/to/subset_9654.fasta', '9655': '/path/to/subset_9655.fasta'} ## correct

To specify reference genomes (analogous to ``--reference``, ``--assembly``, ``--annotation``, ``--attr-mod``, and ``--genetic-code``; note that ``attr_mod`` and ``genetic_code`` are optional if the annotation uses standard attribute field names and the standard genetic code, which the example below does):

>>> my_minorg.reference = 'TAIR10' ## incorrect
AttributeError: can't set attribute
>>> my_minorg.add_reference('TAIR10', '/path/to/TAIR10.fasta', '/path/to/TARI10.gff3', genetic_code = 1, atr_mod = {}) ## correct


