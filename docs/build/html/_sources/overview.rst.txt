Overview
========

MINORg (**MI**\ nimum **NO**\ n-**R**\ eference **g**\ RNA) is a 4-part programme created to design a minimum number of gRNA to cover multiple non-reference targets. Nevertheless, MINORg is also capable of designing gRNA for one target as well as for reference genes. It is available as both a command line programme as well as a Python package.

The 4 broad steps of MINORg are as follows:
* seq: generation of target sequences
* grna: generation of all potential gRNA from target sequences
* filter: filtering of potential gRNA by GC, off-target, and/or within-feature
* minimumset: generation of minimum set(s) of gRNA that cover all target sequences


Command line
------------

Each of the subcommands (``seq``, ``grna``, ``filter``, and ``minimumset``) can be separately executed at the command line using:

.. code-block:: bash
   
   $ minorg <subcommmand> <arguments>

If no subcommand is specified, MINORg will default to the full programme.

.. code-block:: bash
   
   $ minorg <arguments>

To view the help page for each subcommand (which describes the subcommand and its parameters), use:

.. code-block:: bash
   
   $ minorg <subcommand> --help

To view the help page for the full programme, use:

.. code-block:: bash
   
   $ minorg full --help

Using the following will print a help page that lists common parameters and valid subcommand names:

.. code-block:: bash
   
   $ minorg --help

Using a config.ini file, the command line version of MINORg allows users to supply short aliases in place of file names and/or combinations of parameters, as well as set default values for some parameters (such as the reference genome). For example, with the appropriate config.ini setup and lookup files, ``--reference TAIR10`` can be used in place of ``--assembly /path/to/TAIR10/genome.fasta --annotation /path/to/TAIR10/genome.gff3``, and ``--indv ler1`` can be used in place of ``--query /path/to/ler1.fasta``. For details on how to set up a config.ini file, see :ref:`Tutorial:Command line` in the :ref:`Tutorial:Tutorial` section.

The command line interface is a wrapper for the Python package and is built using Typer (https://github.com/tiangolo/typer).


Python package
--------------

Unlike the command line version, the Python package does not support aliases or preset parameter combinations beyond optionally reading default values from a config.ini file. Nevertheless, most arguments at the command line have equivalents as attributes in the class :class:`~minorg.MINORg.MINORg`. For a more comprehensive list of the similarities and differences, please refer to the :ref:`Overview:Parameters` section below.


Parameters
----------

Several CLI arguments have no equivalents in the Python module as they were intended to simplify the building of commands for users who have little to no experience with coding. Users of the Python package are assumed to be comfortable with generating their own preset parameter combinations.

The table below lists the major similarities and differences between CLI arguments and the Python package's MINORg class attributes (note that some attriutes are in fact properties, but they setting them should be no different from setting attributes).

+---------------+----------------------+----------------------+-------------------------+
|Category       |CLI arguments         |Python attributes     |Description              |
+---------------+----------------------+----------------------+-------------------------+
|General        |directory             |directory             |output directory         |
|               +----------------------+----------------------+-------------------------+
|               |prefix                |prefix                |output file/directory    |
|               |                      |                      |prefix                   |
|               +----------------------+----------------------+-------------------------+
|               |thread                |thread                |threads                  |
+---------------+----------------------+----------------------+-------------------------+
|Executable     |blastn                |blastn                |local BLAST's blastn     |
|               +----------------------+----------------------+-------------------------+
|               |rpsblast              |rpsblast              |local BLAST's            |
|               |                      |                      |rpsblast/rpsblast+       |
|               +----------------------+----------------------+-------------------------+
|               |mafft                 |mafft                 |MAFFT                    |
+---------------+----------------------+----------------------+-------------------------+
|Reference      |reference <alias>     |reference <Reference  |reference genome         |
|genomes        |                      |object>               |                         |
|               +----------------------+----------------------+-------------------------+
|(CLI: seq,     |assembly              |                      |reference genome FASTA   |
|full;          +----------------------+----------------------+-------------------------+
|               |annotation            |                      |reference genome GFF     |
|Python: seq,   +----------------------+----------------------+-------------------------+
|filter)        |attr_mod              |                      |mapping for non-standard |
|               |                      |                      |GFF attribute field names|
|               +----------------------+----------------------+-------------------------+
|               |genetic_code          |                      |NCBI genetic code number |
|               |                      |                      |or name                  |
|               +----------------------+----------------------+-------------------------+
|               |ext_gene              |                      |FASTA file of genes to   |
|               |                      |                      |add to reference genome  |
|               +----------------------+----------------------+-------------------------+
|               |ext_cds               |                      |FASTA file of CDS of     |
|               |                      |                      |genes to add to reference|
|               |                      |                      |genome                   |
+---------------+----------------------+----------------------+-------------------------+
|[seq]          |gene                  |gene\ **s**           |gene IDs                 |
|               +----------------------+----------------------+-------------------------+
|target         |cluster               |                      |cluster aliases          |
|definition     +----------------------+----------------------+-------------------------+
|               |indv                  |                      |individuals to discover  |
|               |                      |                      |targets in               |
|               +----------------------+----------------------+-------------------------+
|               |target                |target                |FASTA file of sequences  |
|               |                      |                      |to find gRNA in          |
|               +----------------------+----------------------+-------------------------+
|               |query                 |query                 |FASTA file(s) to discover|
|               |                      |                      |targets in               |
|               +----------------------+----------------------+-------------------------+
|               |domain <alias>        |                      |aliases of domains to    |
|               |                      |                      |find gRNA in             |
|               +----------------------+----------------------+-------------------------+
|               |domain <Pssm-Id>      |pssm_ids              |Pssm-Id(s) of domains to |
|               |                      |                      |find gRNA in             |
|               +----------------------+----------------------+-------------------------+
|               |                      |domain_name           |human-readable domain    |
|               |                      |                      |name used in sequence and|
|               |                      |                      |file names in place of   |
|               |                      |                      |Pssm-Ids                 |
+---------------+----------------------+----------------------+-------------------------+
|[seq]          |minid                 |minid                 |minimum hit % identity   |
|               +----------------------+----------------------+-------------------------+
|inferring      |minlen                |minlen                |minimum merged hits      |
|homologues     |                      |                      |length                   |
|               +----------------------+----------------------+-------------------------+
|from BLASTN    |mincdslen             |mincdslen             |minimum CDS length of    |
|hits           |                      |                      |merged hits              |
|               +----------------------+----------------------+-------------------------+
|               |check_recip           |check_recip           |execute reciprocal check |
|               +----------------------+----------------------+-------------------------+
|               |relax_recip           |relax_recip           |execute relaxed          |
|               |                      |                      |reciprocal check         |
|               +----------------------+----------------------+-------------------------+
|               |merge_within          |merge_within          |maximum distance between |
|               |                      |                      |hits for merging         |
|               +----------------------+----------------------+-------------------------+
|               |check_id_before_merge |check_id_before_merge |filter hits by % identity|
|               |                      |                      |before merging           |
+---------------+----------------------+----------------------+-------------------------+
|[seq]          |db                    |db                    |path to local RPS-BLAST  |
|               |                      |                      |database                 |
|RPS-BLAST      +----------------------+----------------------+-------------------------+
|options        |remote_rps            |remote_rps            |use remote RPS-BLAST     |
|               |                      |                      |database (currently      |
|               |                      |                      |non-functional)          |
+---------------+----------------------+----------------------+-------------------------+
|[grna]         |pam                   |pam                   |PAM pattern              |
|               +----------------------+----------------------+-------------------------+
|               |length                |length                |gRNA length              |
+---------------+----------------------+----------------------+-------------------------+
|[filter]       |gc_min                |gc_min                |minimum GC content       |
|               +----------------------+----------------------+-------------------------+
|GC             |gc_max                |gc_max                |maximum GC content       |
+---------------+----------------------+----------------------+-------------------------+
|[filter]       |feature               |feature               |GFF3 feature type        |
|               +----------------------+----------------------+-------------------------+
|feature        |max_insertion         |max_insertion         |maximum allowable        |
|               |                      |                      |insertion in feature     |
|               |                      |                      |                         |
|               +----------------------+----------------------+-------------------------+
|               |min_within_n          |min_within_n          |minimum number of        |
|               |                      |                      |reference genes which    |
|               |                      |                      |features overlap with    |
|               |                      |                      |gRNA range in alignment  |
|               +----------------------+----------------------+-------------------------+
|               |min_within_fraction   |min_within_fraction   |minimum fraction of      |
|               |                      |                      |reference genes which    |
|               |                      |                      |features overlap with    |
|               |                      |                      |gRNA range in alignment  |
+---------------+----------------------+----------------------+-------------------------+
|[filter]       |background            |background            |FASTA files in which to  |
|               |                      |                      |search for potential     |
|background     |                      |                      |off-targets              |
|               +----------------------+----------------------+-------------------------+
|               |screen_reference      |screen_reference      |include reference genomes|
|               |                      |                      |in search for potential  |
|               |                      |                      |off-targets              |
|               +----------------------+----------------------+-------------------------+
|               |                      |mask                  |FASTA files of additional|
|               |                      |                      |sequences to mask        |
|               +----------------------+----------------------+-------------------------+
|               |unmask_ref            |                      |unmask reference genes   |
|               +----------------------+----------------------+-------------------------+
|               |mask_gene             |                      |additional genes to mask |
|               +----------------------+----------------------+-------------------------+
|               |unmask_gene           |                      |genes to unmask          |
|               +----------------------+----------------------+-------------------------+
|               |mask_cluster          |                      |additional clusters to   |
|               |                      |                      |mask                     |
|               +----------------------+----------------------+-------------------------+
|               |unmask_cluster        |                      |clusters to unmask       |
|               +----------------------+----------------------+-------------------------+
|               |ot_pamless            |ot_pamless            |ignore absense of PAM for|
|               |                      |                      |potential off-targets    |
|               +----------------------+----------------------+-------------------------+
|               |ot_mismatch           |ot_mismatch           |minimum acceptable       |
|               |                      |                      |mismatches for           |
|               |                      |                      |off-targets              |
|               +----------------------+----------------------+-------------------------+
|               |ot_gap                |ot_gap                |minimum acceptable gaps  |
|               |                      |                      |for off-targets          |
|               +----------------------+----------------------+-------------------------+
|               |skip_bg_check         |                      |skip off-target check    |
+---------------+----------------------+----------------------+-------------------------+
|[filter]       |exclude               |exclude               |FASTA file of gRNA       |
|exclude        |                      |                      |sequences to exclude     |
+---------------+----------------------+----------------------+-------------------------+
|[minimumset]   |accept_invalid        |accept_invalid        |score 'NA' as 'pass'     |
|               +----------------------+----------------------+-------------------------+
|               |accept_feature_unknown|accept_feature_unknown|score 'NA' as 'pass' for |
|               |                      |                      |feature check            |
|               +----------------------+----------------------+-------------------------+
|               |                      |accept_invalid_field  |score 'NA' as 'pass' if  |
|               |                      |                      |all entries for a check  |
|               |                      |                      |are 'NA'                 |
|               +----------------------+----------------------+-------------------------+
|               |sets                  |sets                  |number of gRNA sets      |
|               |                      |                      |                         |
|               +----------------------+----------------------+-------------------------+
|               |auto                  |auto                  |generate sets without    |
|               |                      |                      |require manual user      |
|               |                      |                      |confirmation for each set|
+---------------+----------------------+----------------------+-------------------------+
