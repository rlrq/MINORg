Overview
========

Parameters
----------

Several CLI arguments have no equivalents in the Python module as they were intended to simplify the building of commands for users who have little to no experience with coding. Users of the Python package are assumed to comfortable with generating their own preset parameter combinations.


+---------------+---------------------+---------------------+-------------------------+
|Category       |CLI arguments        |Python attributes    |Description              |
+---------------+---------------------+---------------------+-------------------------+
|General        |directory            |directory            |output directory         |
|               +---------------------+---------------------+-------------------------+
|               |prefix               |prefix               |output file/directory    |
|               |                     |                     |prefix                   |
|               +---------------------+---------------------+-------------------------+
|               |thread               |thread               |threads                  |
+---------------+---------------------+---------------------+-------------------------+
|Executable     |blastn               |blastn               |local BLAST's blastn     |
|               +---------------------+---------------------+-------------------------+
|               |rpsblast             |rpsblast             |local BLAST's            |
|               |                     |                     |rpsblast/rpsblast+       |
|               +---------------------+---------------------+-------------------------+
|               |mafft                |mafft                |MAFFT                    |
+---------------+---------------------+---------------------+-------------------------+
|Reference      |reference <alias>    |reference <Reference |reference genome         |
|genomes        |                     |object>              |                         |
|(CLI:          +---------------------+---------------------+-------------------------+
|seq, full;     |assembly             |                     |reference genome FASTA   |
|Python: seq,   +---------------------+---------------------+-------------------------+
|filter)        |annotation           |                     |reference genome GFF     |
+---------------+---------------------+---------------------+-------------------------+
|[seq]          |gene                 |genes                |gene IDs                 |
|target         +---------------------+---------------------+-------------------------+
|definition     |cluster              |                     |cluster aliases          |
|               +---------------------+---------------------+-------------------------+
|               |indv                 |                     |individuals to discover  |
|               |                     |                     |targets in               |
|               +---------------------+---------------------+-------------------------+
|               |target               |target               |FASTA file of sequences  |
|               |                     |                     |to find gRNA in          |
|               +---------------------+---------------------+-------------------------+
|               |query                |query                |FASTA file(s) to discover|
|               |                     |                     |targets in               |
|               +---------------------+---------------------+-------------------------+
|               |domain <alias>       |                     |aliases of domains to    |
|               |                     |                     |find gRNA in             |
|               +---------------------+---------------------+-------------------------+
|               |domain <Pssm-Id>     |pssm_ids             |Pssm-Id(s) of domains to |
|               |                     |                     |find gRNA in             |
+---------------+---------------------+---------------------+-------------------------+
|[seq]          |minid                |minid                |minimum hit % identity   |
|inferring      |                     |                     |                         |
|homologues from+---------------------+---------------------+-------------------------+
|BLASTN hits    |minlen               |minlen               |minimum merged hits      |
|               |                     |                     |length                   |
|               +---------------------+---------------------+-------------------------+
|               |mincdslen            |mincdslen            |minimum CDS length of    |
|               |                     |                     |merged hits              |
|               +---------------------+---------------------+-------------------------+
|               |check_recip          |check_recip          |execute reciprocal check |
|               |                     |                     |                         |
|               +---------------------+---------------------+-------------------------+
|               |relax_recip          |relax_recip          |execute relaxed          |
|               |                     |                     |reciprocal check         |
|               +---------------------+---------------------+-------------------------+
|               |merge_within         |merge_within         |maximum distance between |
|               |                     |                     |hits for merging         |
|               +---------------------+---------------------+-------------------------+
|               |check_id_before_merge|check_id_before_merge|filter hits by % identity|
|               |                     |                     |before merging           |
+---------------+---------------------+---------------------+-------------------------+
|               |                     |                     |                         |
|               |                     |                     |                         |
+---------------+---------------------+---------------------+-------------------------+
|               |                     |                     |                         |
|               |                     |                     |                         |
+---------------+---------------------+---------------------+-------------------------+
|               |                     |                     |                         |
|               |                     |                     |                         |
+---------------+---------------------+---------------------+-------------------------+
|               |                     |                     |                         |
|               |                     |                     |                         |
+---------------+---------------------+---------------------+-------------------------+
|               |                     |                     |                         |
+---------------+---------------------+---------------------+-------------------------+
