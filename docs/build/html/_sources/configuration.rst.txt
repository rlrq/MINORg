Configuration
=============

Note that all command line code snippets in the following tutorial are for **bash terminal**. You may have to adapt them according to your operating system.


config file
-----------

The config.ini file is used for two main purposes: setting default value and assigning short aliases to long values. All entries in the config.ini file that end with 'alias' or 'sets' are used to assign aliases, while the rest are for setting default values. Use the following to export the config file location as an environment variable:

.. code-block:: bash
   
   $ export MINORG_CONFIG=/path/to/config.ini

If you'd like the config file to apply to all users, consider setting it globally.


Alias lookups
-------------

To fully understand how aliases are assigned, we shall use the sample_config.ini file (provided at https://github.com/rlrq/MINORg/blob/master/examples/sample_config.ini) and zoom in on reference genomes as an example. First, let use look at ``reference sets`` under the section ``[lookup]``. Here, we assign aliases for lookup files (which at the command line we can specify using ``--reference-set``). ::
  
  reference sets = athaliana:/path/to/athaliana_genomes.txt
                   arabidopsis:/path/to/arabidopsis_genomes.txt

Each lookup file contains mapping of aliases to reference genome files as well as meta data. Let us look at the file athaliana_genomes.txt (provided at https://github.com/rlrq/MINORg/blob/master/examples/athaliana_genomes.txt). ::

  tair10;TAIR10	/path/to/subset_ref_TAIR10.fasta	/path/to/subset_ref_TAIR10.gff	1	
  araport11	/path/to/subset_ref_Araport11.fasta	/path/to/subset_ref_Araport11.gff	1	

This is a tab-separated file where:

* column 1: semicolon-separated alias(es)
* column 2: path to genome FASTA file
* column 3: path to genome GFF3 file
* column 4: NCBI genetic code number or name (optional if standard genetic code)
* column 5: mapping of nonstandard GFF3 attribute field names to standard field names (optional if standard)

By using ``--reference-set athaliana`` in the command line execution, we can use ``--reference tair10`` (or ``--reference TAIR10``) to tell MINORg to use '/path/to/subset_ref_TAIR10.fasta' and '/path/to/subset_ref_TAIR10.gff' as reference assembly and annotation files respectively without having to type their paths out. In fact, we can even specify multiple reference genomes using ``--reference tair10,araport11`` and making sure to use a comma to separate reference genome aliases.

But there's more! We can set 'athaliana' as the default reference set AND 'tair10' as the default reference in the ``[data]`` section of the config file::

  reference = TAIR10
  reference set = athaliana

This way, unless you wish to use a different reference genome from the default, you won't have to type ``--reference-set athaliana --reference tair10`` either! This is particularly useful if you primarily design gRNA for only a single species. Of course, if you wish, you can combine all reference set files into a single massive file containing the mapping information for all possible reference genomes instead of having multiple files. However, I personally find it easier to maintain smaller files with descriptive names.

You can view the reference genomes available using:

.. code-block:: bash
   
   $ minorg --references ## prints contents of default reference set file
   
   Valid genome aliases (defined in /path/to/athaliana_genomes.txt):
   
   <semicolon-separated genome alias(es)>	<FASTA file>	<GFF3 file>	<NCBI genetic code>	<attribute name mapping>
   tair10;TAIR10	/path/to/subset_ref_Araport11.fasta	/path/to/subset_ref_Araport11.gff
   araport11	/path/to/subset_ref_TAIR10.fasta	/path/to/subset_ref_TAIR10.gff
   
   $ minorg --reference-set arabidopsis --references ## prints contents of arabidopsis_genomes.txt instead
   
   Valid genome aliases (defined in /path/to/arabidopsis_genomes.txt):
   
   <semicolon-separated genome alias(es)>	<FASTA file>	<GFF3 file>	<NCBI genetic code>	<attribute name mapping>
   tair10;TAIR10	/path/to/subset_ref_Araport11.fasta	/path/to/subset_ref_Araport11.gff
   araport11	/path/to/subset_ref_TAIR10.fasta	/path/to/subset_ref_TAIR10.gff
   araly2;alyrata2	/path/to/subset_ref_Araly2.fasta	/path/to/subset_ref_Araly2.gff
   araha1;ahalleri1	/path/to/subset_ref_Araha1.fasta	/path/to/subset_ref_Araha1.gff


The same logic applies as well to ``cluster sets``\ -``cluster set`` (``--cluster-set``\ -``--cluster``\ -``--clusters``) and ``genome sets``\ -``genome set`` (``--genome-set``\ -``--indv``\ -``--genomes``), with the caveat that there is no option to set default clusters or query genomes.

+----------------+--------------+-------------+---------------------------------------------+
|**reference**   |**cluster**   |**genome**   |**description**                              |
+----------------+--------------+-------------+---------------------------------------------+
|reference sets  |cluster sets  |genome sets  |[config] assign alias to lookup files        |
+----------------+--------------+-------------+---------------------------------------------+
|reference set   |cluster set   |genome set   |[config] set default lookup file             |
+----------------+--------------+-------------+---------------------------------------------+
|reference       |              |             |[config] set default value                   |
+----------------+--------------+-------------+---------------------------------------------+
|-\-reference-set|-\-cluster-set|-\-genome-set|[CLI] specify lookup file (alias or path)    |
+----------------+--------------+-------------+---------------------------------------------+
|-\-reference    |-\-cluster    |-\-indv      |[CLI] specify reference/cluster/indv (alias) |
+----------------+--------------+-------------+---------------------------------------------+
|-\-references   |-\-clusters   |-\-genomes   |[CLI] print contents of lookup file to screen|
+----------------+--------------+-------------+---------------------------------------------+



Lookup file formats
~~~~~~~~~~~~~~~~~~~

cluster
+++++++

The lookup files specified in ``cluster sets`` should look like this::

  RPS4;TTR1;RPS4_TTR1	AT5G45050,AT5G45060,AT5G45200,AT5G45210,AT5G45220,AT5G45230,AT5G45240,AT5G45250
  RPS6	AT5G46260,AT5G46270,AT5G46450,AT5G46470,AT5G46490,AT5G46510,AT5G46520

This is a tab-separated file where:

* column 1: semicolon-separated alias(es)
* column 2: comma-separated gene IDs of genes in the cluster


genome
++++++

The lookup files specified in ``genome sets`` should look like this::

  9654	/path/to/subset_9654.fasta
  9655	/path/to/subset_9655.fasta
  9944	/path/to/subset_9944.fasta
  9947	/path/to/subset_9947.fasta

This is a tab-separated file where:

* column 1: semicolon-separated alias(es)
* column 2: path to FASTA file


reference
+++++++++

For ease of reference, the format for reference lookup files specified in ``reference sets`` is reproduced here::
  
  tair10;TAIR10	/path/to/subset_ref_TAIR10.fasta	/path/to/subset_ref_TAIR10.gff	1	
  araport11	/path/to/subset_ref_Araport11.fasta	/path/to/subset_ref_Araport11.gff	1	

This is a tab-separated file where:

* column 1: semicolon-separated alias(es)
* column 2: path to genome FASTA file
* column 3: path to genome GFF3 file
* column 4: NCBI genetic code number or name (optional if standard genetic code)
* column 5: mapping of nonstandard GFF3 attribute field names to standard field names (optional if standard)

..
   However, having multiple files where each file is maintained by a separate user is one way users can update their list of reference genomes without needing a common file accessible by everyone.
