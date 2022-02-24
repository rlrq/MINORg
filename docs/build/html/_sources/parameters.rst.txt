Parameters
==========

In addition to a brief overview of the equivalent parameters between CLI and Python versions of MINORg, this section provides additional information about the input format and utility of some of the more cryptic parameters.

.. contents:: Contents
   :local:
   :depth: 3


CLI vs Python
-------------

Several CLI arguments have no equivalents in the Python module as they were intended to simplify the building of commands for users who have little to no experience with coding. Users of the Python package are assumed to be comfortable with generating their own preset parameter combinations.

The table below lists the major similarities and differences between CLI arguments and the Python package's MINORg class attributes (note that some attributes are in fact properties, but with the exception of ``reference``, setting them should be no different from setting attributes).

Note: Some parameters only apply to specific subcommands (in addition to, of course, the full programme). The relevant subcommands will be indicated in square brackets in the 'Category' column.

Python attributes in the table below indicated with an asterisk (*) should be set using a dedicated method.

+---------------+---------------------------+---------------------------+-------------------------+
|**Category**   |**CLI arguments**          |**Python attributes**      |**Description**          |
+---------------+---------------------------+---------------------------+-------------------------+
|General        |-\-directory               |directory                  |output directory         |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-prefix                  |prefix                     |output file/directory    |
|               |                           |                           |prefix                   |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-thread                  |thread                     |threads                  |
+---------------+---------------------------+---------------------------+-------------------------+
|Executable     |-\-blastn                  |blastn                     |local BLAST's blastn     |
|               +---------------------------+---------------------------+-------------------------+
|(path to       |-\-rpsblast                |rpsblast                   |local BLAST's            |
|executable     |                           |                           |rpsblast/rpsblast+       |
|               +---------------------------+---------------------------+-------------------------+
|if not in      |-\-mafft                   |mafft                      |MAFFT                    |
|command-search +---------------------------+---------------------------+-------------------------+
|               |-\-bedtools                |bedtools                   |EXCEPTION: path to       |
|path)          |                           |                           |directory containing     |
|               |                           |                           |BEDTools executables     |
+---------------+---------------------------+---------------------------+-------------------------+
|Reference      |-\-:ref:`reference         |:ref:`reference            |reference genome         |
|genomes        |<parameters:Reference>`    |<parameters:Reference>`\ * |                         |
|               +---------------------------+---------------------------+-------------------------+
|(CLI: seq,     |-\-assembly                |                           |reference genome FASTA   |
|full;          +---------------------------+---------------------------+-------------------------+
|               |-\-annotation              |                           |reference genome GFF     |
|Python: seq,   +---------------------------+---------------------------+-------------------------+
|filter)        |-\-:ref:`attr_mod          |                           |mapping for non-standard |
|               |<parameters:Attribute      |                           |GFF attribute field names|
|               |modification>`             |                           |                         |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-genetic-code            |                           |NCBI genetic code number |
|               |                           |                           |or name                  |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-extend-gene             |                           |FASTA file of genes to   |
|               |                           |                           |add to reference genome  |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-extend-cds              |                           |FASTA file of CDS of     |
|               |                           |                           |genes to add to reference|
|               |                           |                           |genome                   |
+---------------+---------------------------+---------------------------+-------------------------+
|[seq]          |-\-gene                    |gene\ **s**                |gene IDs                 |
|               +---------------------------+---------------------------+-------------------------+
|target         |-\-cluster                 |                           |cluster aliases          |
|definition     +---------------------------+---------------------------+-------------------------+
|               |-\-indv                    |                           |individuals to discover  |
|               |                           |                           |targets in               |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-target                  |target                     |FASTA file of sequences  |
|               |                           |                           |to find gRNA in          |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-query                   |query                      |FASTA file(s) to discover|
|               |                           |                           |targets in               |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-domain <alias>          |                           |aliases of domains to    |
|               |                           |                           |find gRNA in             |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-domain <Pssm-Id>        |pssm_ids                   |Pssm-Id(s) of domains to |
|               |                           |                           |find gRNA in             |
|               +---------------------------+---------------------------+-------------------------+
|               |                           |domain_name                |human-readable domain    |
|               |                           |                           |name used in sequence and|
|               |                           |                           |file names in place of   |
|               |                           |                           |Pssm-Ids                 |
+---------------+---------------------------+---------------------------+-------------------------+
|[seq]          |-\-minid                   |minid                      |minimum hit % identity   |
|               +---------------------------+---------------------------+-------------------------+
|inferring      |-\-minlen                  |minlen                     |minimum merged hits      |
|homologues     |                           |                           |length                   |
|               +---------------------------+---------------------------+-------------------------+
|from BLASTN    |-\-mincdslen               |mincdslen                  |minimum CDS length of    |
|hits           |                           |                           |merged hits              |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-check-recip             |check_recip                |execute reciprocal check |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-relax-recip             |relax_recip                |execute relaxed          |
|               |                           |                           |reciprocal check         |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-merge-within            |merge_within               |maximum distance between |
|               |                           |                           |hits for merging         |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-check-id-before-merge   |check_id_before_merge      |filter hits by % identity|
|               |                           |                           |before merging           |
+---------------+---------------------------+---------------------------+-------------------------+
|[seq]          |-\-db                      |db                         |path to local RPS-BLAST  |
|               |                           |                           |database                 |
|RPS-BLAST      +---------------------------+---------------------------+-------------------------+
|options        |-\-remote-rps              |remote_rps                 |use remote RPS-BLAST     |
|               |                           |                           |database (currently      |
|               |                           |                           |non-functional)          |
+---------------+---------------------------+---------------------------+-------------------------+
|[grna]         |-\-:ref:`pam               |:ref:`pam <parameters:pam>`|PAM pattern              |
|               |<parameters:pam>`          |                           |                         |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-length                  |length                     |gRNA length              |
+---------------+---------------------------+---------------------------+-------------------------+
|[filter]       |-\-gc-min                  |gc_min                     |minimum GC content       |
|               +---------------------------+---------------------------+-------------------------+
|GC             |-\-gc-max                  |gc_max                     |maximum GC content       |
+---------------+---------------------------+---------------------------+-------------------------+
|[filter]       |-\-feature                 |feature                    |GFF3 feature type        |
|               +---------------------------+---------------------------+-------------------------+
|feature        |-\-max-insertion           |max_insertion              |maximum allowable        |
|               |                           |                           |insertion in feature     |
|               |                           |                           |                         |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-min-within-n            |min_within_n               |minimum number of        |
|               |                           |                           |reference genes which    |
|               |                           |                           |features overlap with    |
|               |                           |                           |gRNA range in alignment  |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-min-within-fraction     |min_within_fraction        |minimum fraction of      |
|               |                           |                           |reference genes which    |
|               |                           |                           |features overlap with    |
|               |                           |                           |gRNA range in alignment  |
+---------------+---------------------------+---------------------------+-------------------------+
|[filter]       |-\-background              |background                 |FASTA files in which to  |
|               |                           |                           |search for potential     |
|background     |                           |                           |off-targets              |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-screen-reference        |screen_reference           |include reference genomes|
|               |                           |                           |in search for potential  |
|               |                           |                           |off-targets              |
|               +---------------------------+---------------------------+-------------------------+
|               |                           |mask                       |FASTA files of additional|
|               |                           |                           |sequences to mask        |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-unmask-ref              |                           |unmask reference genes   |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-mask-gene               |                           |additional genes to mask |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-unmask-gene             |                           |genes to unmask          |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-mask-cluster            |                           |additional clusters to   |
|               |                           |                           |mask                     |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-unmask-cluster          |                           |clusters to unmask       |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-ot-pamless              |ot_pamless                 |ignore absense of PAM for|
|               |                           |                           |potential off-targets    |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-ot-mismatch             |ot_mismatch                |minimum acceptable       |
|               |                           |                           |mismatches for           |
|               |                           |                           |off-targets              |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-ot-gap                  |ot_gap                     |minimum acceptable gaps  |
|               |                           |                           |for off-targets          |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-skip-bg-check           |                           |skip off-target check    |
+---------------+---------------------------+---------------------------+-------------------------+
|[filter]       |-\-exclude                 |exclude                    |FASTA file of gRNA       |
|exclude        |                           |                           |sequences to exclude     |
+---------------+---------------------------+---------------------------+-------------------------+
|[minimumset]   |-\-accept-invalid          |accept_invalid             |score 'NA' as 'pass'     |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-accept-feature-unknown  |accept_feature_unknown     |score 'NA' as 'pass' for |
|               |                           |                           |feature check            |
|               +---------------------------+---------------------------+-------------------------+
|               |                           |accept_invalid_field       |score 'NA' as 'pass' if  |
|               |                           |                           |all entries for a check  |
|               |                           |                           |are 'NA'                 |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-sets                    |sets                       |number of gRNA sets      |
|               |                           |                           |                         |
|               +---------------------------+---------------------------+-------------------------+
|               |-\-auto                    |auto                       |generate sets without    |
|               |                           |                           |require manual user      |
|               |                           |                           |confirmation for each set|
+---------------+---------------------------+---------------------------+-------------------------+


Parameter types
---------------

Flags and arguments
~~~~~~~~~~~~~~~~~~~

Flag
++++

Flags are parameters that do not take values.

**CLI**: ``--auto``

For example:

.. code-block:: bash

   $ minorg <other arguments> --auto

Simply using ``--auto`` tells MINORg to automate set generation.


**Python**: ``auto``

In Python, flags are raised by setting the value of their attributes to ``True`` or ``False``. For example:

>>> from minorg.MINORg import MINORg
>>> my_minorg = MINORg()
>>> my_minorg.auto = True ## raise flag for parameter 'auto'


Argument
++++++++

These parameters take values.

**CLI**: all parameters that are not flags

.. code-block:: bash
   
   $ minorg <other arguments> --prefix my_minorg

``--prefix my_minorg`` tells MINORg to use 'my_minorg' as a prefix for output files and directories.

**Python**: all parameters that are not flags

>>> from minorg.MINORg import MINORg
>>> my_minorg = MINORg()
>>> my_minorg.prefix = 'my_minorg' ## tells MINORg to use 'my_minorg' as prefix for output files and directories


Paths
~~~~~

| **CLI**: As all paths will be resolved to absolute paths, relative paths are acceptable. Nevertheless, do be careful with relative paths and NEVER use them in the config file or in lookup files.
| **Python**: Paths are NOT resolved (except directory and config file). Absolute paths are STRONGLY RECOMMENDED. Be careful with relative paths.

Executables
+++++++++++

Default values for executables may be specified in the config file (see :ref:`Configuration:Configuration` for more on the config file).


blastn, rpsblast/rpsblast+, MAFFT
_________________________________

**CLI**: ``--blastn``, ``--rpsblast``, ``--maff``

.. code-block:: bash

   $ minorg <other arguments> --blastn /usr/bin/blastn

**Python**: ``blastn``, ``rpsblast``, ``mafft``

>>> from minorg.MINORg import MINORg
>>> my_minorg = MINORg()
>>> my_minorg.blastn = '/usr/bin/blastn' ## tells MINORg where the blastn executable is

If an executable is in the command-search path, specifying these parameters is optional, although you may, if you desire, specify the command itself (e.g. 'blastn' instead of '/usr/bin/blastn'). If not, the **path to the executable** is required.

To determine if blastn and rpsblast (or rpsblast+ depending on your BLAST+ version) in the command-search path, execute::

  blastn -version

If it prints something like ::

  blastn: 2.6.0+
   Package: blast 2.6.0, build Jan 15 2017 17:12:27

then 'blastn' IS in your command-search path. Repeat this with 'rpsblast' and/or 'rpsblast+'.

BEDTools
________

**CLI**: ``bedtools``

.. code-block:: bash

   $ minorg <other arguments> --bedtools /path/to/bedtools2/bin/

**Python**: ``bedtools``

>>> from minorg.MINORg import MINORg
>>> my_minorg = MINORg()
>>> my_minorg.bedtools = '/path/to/bedtools2/bin/' ## tells MINORg where the BEDTools executables are

If bedtools is in the command-search path, you should NOT use this parameter. If not, the path to the **directory containing the BEDTools executables** is required.

To determine if the BEDTools executables are in your command-search path, execute::

.. code-block::
   
   bedtools --version

If it prints something like ::

  bedtools v2.26.0

then 'bedtools' is in your command-search path.


Alias lookup
~~~~~~~~~~~~

Note that aliases are **case-sensitive**.

1-level lookup
++++++++++++++

See also: :ref:`configuration:1-level lookup`

1-level lookup parameters have preset values mapped to aliases defined in a configuration file. Users may use either the alias(es) or provide raw values.

| **CLI**: ``--assembly``, ``--annotation``, ``--db``, ``attr-mod``, ``--domain``
| **Python**: Does not support aliases. Raw values only.


2-level lookup
++++++++++++++

See also: :ref:`configuration:2-level lookup`

2-level lookup parameters use a combination of 2 parameters. The first parameter (suffixed with `set`) specifies a file containing alias mapping information for the second parameter (not suffixed). Aliases for the first parameter are defined in a configuration file, and functions effectively the same way a 1-level lookup parameter does. The second parameter reads alias mapping information from the file specified by the first parameter. Unlike the first parameters, users may only use alias(es)--raw values are not allowed. To specify raw values, different parameters must be used (see :ref:`configuration:Alternative parameters` for which).

| **CLI**: ``--reference-set``\ -``reference``, ``--cluster-set``\ -``--cluster``, ``--genome-set``\ -``--indv``
| **Python**: Does not support aliases. Raw values only.


Predefined lookup
+++++++++++++++++

Predefined lookup parameters are built into the programme. Users may use either the alias(es) or raw values.

| **CLI**: ``--pam``
| **Python**: ``pam``


Raw values
++++++++++

All other parameters are raw values only.


Multiple values
~~~~~~~~~~~~~~~

Comma-separated (CLI)
+++++++++++++++++++++

**CLI**: ``--reference``, ``--cluster``, ``--gene``, ``--indv``

Comma-separated multiple value arguments accept multiple values for a single argument so long as the values are comma-separated. For example, multiple genes can be specified using ``--gene 'geneA,geneB,geneC'``.


Multi-argument (CLI)
++++++++++++++++++++

**CLI**: ``--reference``, ``--cluster``, ``--gene``, ``--indv``, ``--query``, ``--feature``, ``--ext-gene``, ``--ext-cds``, ``--mask-gene``, ``--unmask-gene``, ``--mask-cluster``, ``--unmask-cluster``
.. , ``--ot-indv`` (not implemented)

Multi-argument parameters accept multiple values by re-using a parameter. For example, multiple genes can be specified using ``--gene geneA --gene geneB --gene geneC``.

(Note that some parameters can be both comma-separated AND multi-argument, and that these features can be combined. For example, ``--gene geneA --gene geneB,geneC`` is also valid.)


Multi-value list (Python)
+++++++++++++++++++++++++

**Python**: ``genes``

Multiple values for a single parameter may be provided to MINORg in a list. For example:

>>> from minorg.MINORg import MINORg
>>> my_minorg = MINORg()
>>> my_minorg.genes = ['geneA'] ## specify a single value
>>> my_minorg.genes = ['geneA', 'geneB', 'geneC'] ## specify multiple values


Multi-value dictionary (Python)
+++++++++++++++++++++++++++++++

**Python**: ``query``, ``background``




Parameter descriptions
----------------------

Reference
~~~~~~~~~

**Type**: :ref:`Parameters:Argument`, :ref:`Parameters:2-level lookup`

| **CLI**: ``--reference``
| **Python**: ``reference``
| **Config file**:

  | set default: ``reference`` (section ``[data]``)
  | set default set: ``reference set`` (section ``[data]``)
  | assign aliases to sets: ``reference sets`` (section ``[lookup]``)




TODO: also link to attribute modification section below when describing setting reference for Python


Attribute modification
~~~~~~~~~~~~~~~~~~~~~~

**Type**: :ref:`Parameters:Argument`, :ref:`Parameters:1-level lookup`

| **CLI**: ``--attr-mod``
| **Python**: NA (see argument ``attr_mod`` of :meth:`~minorg.MINORg.MINORg.add_reference` instead)
| **Config file**:

  | set default: ``gff attribute modification`` (section ``[data]``)
  | assign aliases: ``gff attribute modification presets`` (section ``[lookup]``)


This parameter tells MINORg how to map non-standard GFF3 field names to standard GFF3 field names. This feature was originally developed when I tried to retrieve sequences using the IRGSP-1.0 annotation for rice (*Oryza sativa* subsp. Nipponbare) and discovered that it uses 'Locus_id' instead of 'Parent' for mRNA annotations.

The input given to ``--attr-mod`` should follow this format (with quotes):

  ‘<feature>:<standard>=<nonstandard>,<standard>=<nonstandard>;<feature>:<standard>=<nonstandard>’

Examples:

  ``--attr-mod 'mRNA:Parent=Locus_id,ID=transcript_id;CDS:Parent=transcript_id'``
    'Locus_id' and 'transcript_id' are non-standard field names for
    fields 'Parent' and 'ID' respectively for the feature type 'mRNA',
    and 'transcript_id' is the non-standard name for the field 'Parent' for the feature type 'CDS'.

  ``--attr-mod 'all:ID=id'``
    'id' is the non-standard field name for the field 'ID' for all feature types.

See http://gmod.org/wiki/GFF3 for standard attribute field names (see section titled ‘Column 9: “attributes”’).


Extended genome
~~~~~~~~~~~~~~~

**Type**: :ref:`Parameters:Argument`, :ref:`Parameters:Raw values`, :ref:`Parameters:Multi-argument (CLI)`, :ref:`Parameters:Multi-value list (Python)`

| **CLI**: ``--extend-gene``, ``--extend-cds``
| **Python**: ``ext_gene``, ``ext_cds``

These parameters accept FASTA files and allow MINORg to infer coding regions (CDS) from genomic (``--extend-gene``, ``ext_gene``) and CDS-only (``--extend-cds``, ``ext_cds``) sequences. They should be used when you do not have a GFF3 annotation file for your desired genes, but DO have the above mentioned sequences. MINORg will align gene and CDS-only sequences using MAFFT to generate a GFF3 annotation file with inferred intron-exon boundaries. These genes will then be added to the reference genome **and you can use their gene IDs as you would reference gene IDs**. You may provide multiple files to each parameter--MINORg will process them all simultaneously.

For MINORg to map the CDS-only sequences to the correct gene sequences, CDS-only sequences should be named according to the the format: '<gene ID>.<CDS ID>'

For example, given the following CDS sequences::

  >geneA.1
  ATGATGATGATGATGATGATGATGTAA
  >geneA.two
  ATGATGATGATGATGATGATGTAA
  >geneA.foo.bar
  ATGATGATGATGATGATGTAA
  >geneB.1
  ATGAAAAAAAAAAAAAAAAAATAA

And the following gene sequences::

  >geneA
  ATGATGATGATGATGATGATGATGTAA
  >geneA.foo
  ATGATGATGATGATGATGATGATGTAA
  >geneB
  ATGAAAAAAAAAAAAAAAAAAAAAAAATAA

CDS sequences ``geneA.1`` and ``geneA.two`` will be mapped to gene sequence ``geneA``, ``geneA.foo.bar`` will be mapped to ``geneA.foo``, and ``geneB.1`` will be mapped to ``geneB``. Note that ``geneA.1`` and ``geneA.two`` will be treated as different isoforms of the gene ``geneA``. 


PAM
~~~

**Type**: :ref:`Parameters:Argument`, :ref:`Parameters:Predefined lookup`

| **CLI**: ``--pam``
| **Python**: ``pam``
| **Config file**:

  | set default: ``pam`` (section ``[grna]``)
  | assign aliases: ``pam alias`` (section ``[lookup]``) (not yet implemented)


By default, MINORg designs gRNA for SpCas9 systems (i.e. 3' NGG PAM). You may specify other PAM patterns for non-SpCas9 systems using ``--pam``. It is recommended that any PAM pattern that uses special characters be enclosed in quotes, as it may lead to unexpected behaviour otherwise at the terminal.

Under the hood, MINORg uses regex to match PAM sites. Therefore, it is in theory possible to utilise the full suite of Python regex syntax to customise your PAM pattern. However, do take care to avoid using  ``.`` as a wildcard, as MINORg uses this character to determine where gRNA is relative to a PAM pattern.


Ambiguous bases and repeats
+++++++++++++++++++++++++++

Unlike many gRNA designers, MINORg accepts ambiguous bases (see: https://genome.ucsc.edu/goldenPath/help/iupac.html for IUPAC codes) as well as variable number of repeats.

  Example: The pattern 'R{1,2}T' (where 'R' means 'A' or 'G', and {1,2} means either 1 to 2 repetitions
  of the character right before it) will match 'AT', 'GT', 'AAT', 'AGT', 'GAT', and 'GGT'.


Spacers and 3' or 5' PAM
++++++++++++++++++++++++

In the absence of 'N' in the PAM pattern, MINORg will assume 3' PAM with 1 spacer base (such as in the 3' 'NGG' of SpCas9). If a pattern includes an 'N' at either end, MINORg will assume that the gRNA is directly adjacent to the 'N' base of the pattern. To specify a 5' PAM in the absence of 'N' in the PAM pattern, '.' should be inserted where the gRNA is.

  Example 1: ``--pam .NGG`` and ``--pam NGG`` and ``--pam GG`` are functionally identical.
  The latter two will be expanded to the most explicit pattern: ``.NGG``.

  Example 2: If a CRISPR system uses 'GG' PAM with NO spacer 'N' base, the PAM pattern has to be
  specified to MINORg as ``--pam .GG``. Otherwise, MINORg will insert a spacer 'N' base, giving rise
  to the incorrect explicit pattern of ``.NGG`` instead.

  Example 3: AacCas12b uses a 5' PAM with the pattern 'TTN', which can be specified to MINORg as
  ``--pam TTN`` or ``--pam TTN.``, where ``.`` indicates where the gRNA is.
  ``.`` is optional as this PAM pattern (TTN) includes 'N' at the end.
  Therefore, MINORg will infer a 5' PAM.

  Example 4: Cas12a uses a 5' PAM with the pattern 'TTTV', which can be specified to MINORg as
  ``--pam TTTV.`` or ``--pam 'T{3}V.'``, where ``.`` indicates where the gRNA is.
  As the PAM pattern does not include 'N', the gRNA position MUST be explicitly indicated using ``.``.
  If ``--pam TTTV`` is (incorrectly) used, MINORg will default to a 3' PAM AND add a spacer base,
  expanding it to the undesired explicit pattern ``.NTTTV`` .

For a PAM-less search, use: ``--pam .`` or ``--pam '.'``.


Preset PAM patterns
+++++++++++++++++++

MINORg comes with several preset PAM patterns for different CRISPR systems.

  For example: ``--pam SpCas9`` and ``--pam .NGG`` are functionally identical.

+-------------+----------------+--------------------------------------+
|**alias(es)**|**PAM sequence**|**Notes**                             |
|             |(explicit)      |                                      |
+-------------+----------------+--------------------------------------+
|SpCas9 OR    |.NGG            |default                               |
|spcas9       |                |                                      |
+-------------+----------------+--------------------------------------+
|SaCas9T OR   |.NGRRT          |                                      |
|sacas9t      |                |                                      |
+-------------+----------------+--------------------------------------+
|SaCas9N OR   |.NGRRN          |                                      |
|sacas9n      |                |                                      |
+-------------+----------------+--------------------------------------+
|NmeCas9 OR   |.NNNNGATT       |                                      |
|nmecas9      |                |                                      |
+-------------+----------------+--------------------------------------+
|CjCas9 OR    |.NNNNRYAC       |                                      |
|cjcas9       |                |                                      |
+-------------+----------------+--------------------------------------+
|StCas9 OR    |.NNAGAAW        |                                      |
|stcas9       |                |                                      |
+-------------+----------------+--------------------------------------+
|Cas12a OR    |TTTV.           |5' PAM                                |
|cas12a       |                |                                      |
+-------------+----------------+--------------------------------------+
|AacCas12b OR |TTN.            |5' PAM                                |
|aaccas12b    |                |                                      |
+-------------+----------------+--------------------------------------+
|BhCas12b OR  |DTTN.           |5' PAM                                |
|bhcas12b     |                |                                      |
+-------------+----------------+--------------------------------------+

..
   |Cas14ds OR   |.TTTA           |T-rich PAM for dsDNA cleavage (no PAM |
   |cas14ds      |                |required for ssDNA)                   |
   +-------------+----------------+--------------------------------------+



RPS-BLAST local database
~~~~~~~~~~~~~~~~~~~~~~~~
**Type**: :ref:`Parameters:Argument`, :ref:`Parameters:1-level lookup`

| **CLI**: ``--db``
| **Python**: ``db``

  | set default: ``rps database`` (section ``[data]``)
  | assign aliases: ``rps database alias`` (section ``[lookup]``)

The latest CDD database may be downloaded at ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.targ.gz. As the CDD database is regularly updated, the PSSM-Id for a domain shown at the CDD website is subject to change. Thus, I also recommend downloading ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz, which contains information that maps PSSM-Ids to domain accession IDs as well as domain names of the database version at the point of downloading.

Note: As the local database itself consists of multiple files with different extensions, the path provided to this parameter is not to any single file. For example, given the following file structure::

  /
  +-- root/
      |-- other_files/
      |-- rps_db/
          |-- Cdd.aux
          |-- Cdd.freq
          |-- Cdd.loo
          |-- Cdd.phr
          |-- Cdd.pin
          |-- Cdd.psd
          |-- Cdd.psi
          |-- Cdd.psq
          +-- Cdd.rps

where the database is contained in the directory ``/root/rsp_db/``, the appropriate path to pass to this parameter is: ``/root/rps_db/Cdd``, where the trailing 'Cdd' is the prefix of all of the database's files


RPS-BLAST remote database
~~~~~~~~~~~~~~~~~~~~~~~~~
**Type**: :ref:`Parameters:Flag`

| **CLI**: ``--remote-rps``
| **Python**: ``remote_rps``

  | set default: ``remote rps`` (section ``[data]``)

While it is in theory possible to use the remote CDD database & servers instead of local ones, the ``--remote`` option for the 'rpsblast'/'rpsblast+' command from the BLAST+ package has never worked for me. In any case, if your version of local rpsblast is able to access the remote database, you can use ``--remote-rps`` instead of ``--db /path/to/rpsblast/db``.
