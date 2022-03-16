Configuration
=============

Note that all command line code snippets in the following tutorial are for **bash terminal**. You may have to adapt them according to your operating system.


config file
-----------

The config.ini file is used for two main purposes: setting default values and assigning short aliases to long values. All entries in the config.ini file that end with 'alias' or 'sets' are used to assign aliases, while the rest are for setting default values. Use the following to export the config file location as an environment variable:

.. code-block:: bash
   
   $ export MINORG_CONFIG=/path/to/config.ini

If you'd like the config file to apply to all users, consider setting it globally.


Alias lookups
-------------

There are two types of alias lookups in MINORg: 1-level lookups and 2-level lookups.

1-level lookup
~~~~~~~~~~~~~~

1-level lookups are defined directly in the config file. Parameters that use 1-level lookups are:

* ``--assembly``
* ``--annotation``
* ``--db``
* ``--attr-mod``
* ``--domain``

To understand how aliases are assigned using the config file, let us take a look at ``assembly alias`` under section ``[lookup]``. Here, we assign aliases for FASTA files of reference assemblies in the format <semicolon-separated alias>:<value>. ::

  assembly alias = tair10;TAIR10:/path/to/subset_ref_TAIR10.fasta
                   araport11:/path/to/subset_ref_Araport11.fasta
                   araly2:/path/to/subset_ref_Araly2.fasta
                   araha1:/path/to/subset_ref_Araha1.fasta

With this setup, we can use ``--assembly TAIR10`` instead of ``--assembly /path/to/subset_ref_TAIR10.fasta`` when building our MINORg command (although either way is acceptable). This applies to ``annotation alias`` (``--annotation``), ``rps database alias`` (``--db``), and ``gff attribute modification presets`` (``--attr-mod``) as well.

+----------------+--------------+----------------+----------------+---------------------------------------------+
|**assembly**    |**annotation**|**rps database**|**attr mod**    |**description**                              |
+----------------+--------------+----------------+----------------+---------------------------------------------+
|assembly alias  |annotation    |rps database    |gff attribute   |[config] assign alias to values              |
|                |alias         |alias           |modification    |                                             |
|                |              |                |presets         |                                             |
+----------------+--------------+----------------+----------------+---------------------------------------------+
|-\-assembly     |-\-annotation |-\-db           |-\-attr-mod     |[CLI] specify value (alias or raw)           |
+----------------+--------------+----------------+----------------+---------------------------------------------+

Unlike the above parameters, ``domain alias`` is specified in the reverse format using commas instead of semi-colons (<value>:<comma-separated aliases) for readability reasons, as most PSSM-Ids are comprised of the same number of digits. Therefore, the ``domain alias`` in section ``[lookup]`` looks like this instead::
  
  domain alias = 366714:TIR
                 366375:NB-ARC,NBS
                 375519:Rx_N

As with the above parameters, you can use ``--domain TIR`` instead of ``--domain 366714`` (although either way is acceptable). Do note, however, that if you use ``--domain TIR``, newly generate file and sequence names will use ``TIR`` instead of ``366714`` if the domain is to be included in the name.


2-level lookup
~~~~~~~~~~~~~~

2-level lookups are defined using a file that maps aliases to values. This means that there is a first layer of lookup that maps aliases to those files, and then a second layer of lookup in those files that maps aliases to values. Parameters that use 2-level lookups are:

* ``--reference``
* ``--cluster``
* ``--indv``

To understand how aliases are assigned, we shall use the sample_config.ini file (provided at https://github.com/rlrq/MINORg/blob/master/examples/sample_config.ini) and zoom in on reference genomes as an example. First, let use look at ``reference sets`` under the section ``[lookup]``. Here, we assign aliases for lookup files (which at the command line we can specify using ``--reference-set``). ::
  
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

By using ``--reference-set athaliana`` in the command line execution, we can use ``--reference TAIR10`` (or ``--reference tair10``) to tell MINORg to use '/path/to/subset_ref_TAIR10.fasta' and '/path/to/subset_ref_TAIR10.gff' as reference assembly and annotation files respectively without having to type their paths out. In fact, we can even specify multiple reference genomes using ``--reference TAIR10,araport11`` and making sure to use a comma to separate reference genome aliases.

But there's more! We can set 'athaliana' as the default reference set AND 'TAIR10' as the default reference in the ``[data]`` section of the config file::

  reference = TAIR10
  reference set = athaliana

This way, unless you wish to use a different reference genome from the default, you won't have to type ``--reference-set athaliana --reference TAIR10`` either! This is particularly useful if you primarily design gRNA for only a single species. Of course, if you wish, you can combine all reference set files into a single massive file containing the mapping information for all possible reference genomes instead of having multiple files. However, I personally find it easier to maintain smaller files with descriptive names.

You can view the reference genomes in the default reference set using:

.. code-block:: bash
   
   $ minorg --references ## prints contents of default reference set file
   
   Valid genome aliases (defined in /path/to/athaliana_genomes.txt):
   
   <semicolon-separated genome alias(es)>	<FASTA file>	<GFF3 file>	<NCBI genetic code>	<attribute name mapping>
   tair10;TAIR10	/path/to/subset_ref_Araport11.fasta	/path/to/subset_ref_Araport11.gff
   araport11	/path/to/subset_ref_TAIR10.fasta	/path/to/subset_ref_TAIR10.gff

To view the reference genomes in a non-default reference set, use:

.. code-block:: bash
   
   $ minorg --reference-set arabidopsis --references ## prints contents of arabidopsis_genomes.txt instead
   
   Valid genome aliases (defined in /path/to/arabidopsis_genomes.txt):
   
   <semicolon-separated genome alias(es)>	<FASTA file>	<GFF3 file>	<NCBI genetic code>	<attribute name mapping>
   tair10;TAIR10	/path/to/subset_ref_Araport11.fasta	/path/to/subset_ref_Araport11.gff
   araport11	/path/to/subset_ref_TAIR10.fasta	/path/to/subset_ref_TAIR10.gff
   araly2;alyrata2	/path/to/subset_ref_Araly2.fasta	/path/to/subset_ref_Araly2.gff
   araha1;ahalleri1	/path/to/subset_ref_Araha1.fasta	/path/to/subset_ref_Araha1.gff

Note that you can also provide an alias mapping file that is not in the config file to MINORg by specifying the path to the file instead of using a non-existent alias (e.g. ``--reference-set /path/to/arabidopsis_genomes.txt``).

Parameters
++++++++++

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
|-\-reference    |-\-cluster    |-\-indv      |[CLI] specify reference/cluster/indv         |
|                |              |             |(comma-separated alias(es))                  |
+----------------+--------------+-------------+---------------------------------------------+
|-\-references   |-\-clusters   |-\-genomes   |[CLI] print contents of lookup file to screen|
+----------------+--------------+-------------+---------------------------------------------+

Alternative parameters
++++++++++++++++++++++

Do note that, because of the nature of these lookups, you cannot simply provide the value(s) mapped to an alias to ``--reference``, ``--cluster``, or ``--indv``. If the desired files/genes are not specified in any mapping file, you will have to use the following alternatives:

+-------------+--------------------+---------------------------------------------------------+
|             |**using aliases**   |**not using aliases**                                    |
+-------------+--------------------+---------------------------------------------------------+
|**reference**|-\-reference        |-\-assembly <path to FASTA> -\-annotation <path to GFF3> |
|             |<alias(es)>         |-\-genetic-code <number or name> -\-attr-mod <attrinute  |
|             |                    |modiication>*                                            |
+-------------+--------------------+---------------------------------------------------------+
|**cluster**  |-\-cluster          |-\-gene <comma-separated gene IDs>**                     |
|             |<alias(es)>         |                                                         |
+-------------+--------------------+---------------------------------------------------------+
|**indv**     |-\-indv <alias(es)> |-\-query <path to FASTA>***                              |
+-------------+--------------------+---------------------------------------------------------+

\* ``--genetic-code`` and ``--attr-mod`` are optional if the reference genome uses the standard genetic code and standard GFF attribute field names respectively. Do note that you CANNOT SPECIFY MULTIPLE reference genomes if not using aliases.

\*\* If not using aliases, each cluster must be processed separately (i.e. a different MINORg execution for each cluster), as MINORg has no way of knowing which gene belongs to which cluster if you use ``--gene``.

\*\*\* You may specify multiple query FASTA files by using ``--query <FASTA>`` as many times as needed (e.g. ``--query <FASTA 1> --query <FASTA 2> --query <FASTA 3>``).


2-level lookup file formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
* column 5: mapping of nonstandard GFF3 attribute field names to standard field names (optional if standard) (see :ref:`Parameters:Attribute modification format (reference lookup file)` for format)

..
   However, having multiple files where each file is maintained by a separate user is one way users can update their list of reference genomes without needing a common file accessible by everyone.
