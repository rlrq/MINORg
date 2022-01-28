Examples
========

In all the following examples, the current directory/working directory is presumed to contain all files in https://github.com/rlrq/MINORg/tree/master/examples.

Example 1
---------

Generate gRNA in the coding regions of a reference gene.

* We will use the gene 'AT1G33560'
* As the default feature type for generating gRNA in is "CDS" (coding region), we do not need to explicitly specify it

**Python code**

>>> from minorg.MINORg import MINORg
>>> my_minorg = MINORg(prefix = "example1")
>>> my_minorg.add_reference("TAIR10", "subset_ref_TAIR10.fasta", "subset_ref_TAIR10.gff")
>>> my_minorg.genes = "AT1G33560"
>>> my_minorg.query_reference = True
>>> my_minorg.seq() ## generate target sequence(s), in this case that of AT1G33560
>>> my_minorg.grna() ## generate all possible gRNA in target sequence(s)
>>> my_minorg.filter_feature() ## filter gRNA by feature type
>>> my_minorg.write_pass_grna_fasta() ## write gRNA that pass all filters (in this case, only 1 filter was used)
>>> my_minorg.pass_fasta
'/path/to/current/directory/example1/example1_gRNA_pass.fasta'


**CLI code (bash)**

.. code-block:: bash
   
   $ minorg --prefix example 1 \
            --assembly ./subset_ref_TAIR10.fasta --annotation ./subset_ref_TAIR10.gff \ ## reference genome
            --indv ref --gene AT1G33560 \ ## query genome and gene
            --gc-min 0 --gc-max 1 \ ## accept all gRNA regardless of GC content
            --auto ## auto-approve minimum gRNA sets

Note that the CLI code invokes the full programme, so minimum gRNA set(s) will also be generated. Nevertheless, you can still find a FASTA file of all passing gRNA (that are within the desired feature AND have desired GC content) at '/path/to/current/directory/example1/example1_gRNA_pass.fasta'.


**CLI code (bash)**: using reference genome aliases

.. code-block:: bash
   
   $ minorg --prefix example 1 \
            --reference-set ./arabidopsis_genomes.txt --ref TAIR10 \ ## reference genome
            --indv ref --gene AT1G33560 \ ## query genome and gene
            --gc-min 0 --gc-max 1 \ ## accept all gRNA regardless of GC content
            --auto ## auto-approve minimum gRNA sets

Here, we use ``--reference-set`` to give MINORg a mapping of reference genome FASTA and GFF files to more easy to write aliases. In this case, the alias for the combination of './subset_ref_TAIR10.fasta' and './subset_ref_TAIR10.gff' is 'TAIR10', which we pass to ``--ref``.


**CLI code (bash)**: using config.ini

Open ./sample_config.ini and replace '/path/to/arabidopsis_genomes.txt' with the absolute path to ./araidopsis_genomes.txt. Save it.

.. code-block:: bash
   
   $ export MINORG_CONFIG=./sample_config.ini
   $ minorg --prefix example 1 \
            --indv ref --gene AT1G33560 \ ## query genome and gene
            --gc-min 0 --gc-max 1 \ ## accept all gRNA regardless of GC content
            --auto ## auto-approve minimum gRNA sets

If we open './sample_config.ini', we will find the following lines::

  reference sets = athaliana:/path/to/athaliana_genomes.txt
                   arabidopsis:/<absolute path>/arabidopsis_genomes.txt

These lines tell MINORg to use 'arabidopsis' as an alias for arabidopsis_genomes.txt. ::

  reference set = arabidopsis

This line tells MINORg to pass arabidopsis_genomes.txt (mapped from alias 'arabidopsis') to ``--reference-set``. ::

  reference = TAIR10

And finally this line tells MINORg to use the FASTA and GFF combination mapped to the alias 'TAIR10' in arabidopsis_genomes.txt as the reference genome.


Example 2
---------

Description of desired outcome

Python
CLI
