Tutorial
========

In all the following tutorial, the current directory/working directory is presumed to contain all files in https://github.com/rlrq/MINORg/tree/master/examples. If you have not downloaded the files, please do so and navigate to the directory that contains them.


Command line
------------

Note that all command line code snippets in the following tutorial are for **bash terminal**. You may have to adapt them according to your operating system.

Getting started
~~~~~~~~~~~~~~~

Let us begin with the simplest MINORg execution:

.. code-block:: bash
   
   $ minorg --target ./sample_CDS.fasta --directory ./example00
   Final gRNA sequence(s) have been written to /tmp/tmptrr0sca_/minorg/minorg_gRNA_final.fasta
   Final gRNA sequence ID(s), gRNA sequence(s), and target(s) have been written to /tmp/tmptrr0sca_/minorg/minorg_gRNA_final.map
   
   1 mutually exclusive gRNA set(s) requested. 1 set(s) found.
   Output files have been generated in /path/to/current/directory/example00

The above combination of arguments tells MINORg to generate gRNA from targets in a user-provided FASTA file (``--target ./sample_CDS.fasta``) and to output files into directory ``--directory ./example00``. By default, MINORg generates 20 bp gRNA using NGG PAM.



Python package
--------------

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


