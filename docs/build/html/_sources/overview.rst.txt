Overview
========

(Note that all command line code snippets in the following tutorial are for **bash terminal**. You may have to adapt them according to your operating system.)

MINORg (**MI**\ nimum **NO**\ n-**R**\ eference **g**\ RNA) is a 4-part programme created to design a minimum number of gRNA to cover multiple non-reference targets. Nevertheless, MINORg is also capable of designing gRNA for one target as well as for reference genes or from user-provided targets. It is available as both a command line programme as well as a Python package.

The 4 broad steps of MINORg are as follows:

* seq: generation of target sequences
* grna: generation of all potential gRNA from target sequences
* filter: filtering of potential gRNA by GC, off-target, and/or within-feature
* minimumset: generation of minimum set(s) of gRNA that cover all target sequences


Command line
------------

The command line interface is a wrapper for the Python package and is built using Typer (https://github.com/tiangolo/typer).


Subcommands
~~~~~~~~~~~

Each of the subcommands (``seq``, ``grna``, ``filter``, and ``minimumset``) can be separately executed at the command line using:

.. code-block:: bash
   
   $ minorg <subcommmand> <arguments>

If no subcommand is specified, MINORg will default to the full programme.

.. code-block:: bash
   
   $ minorg <arguments>

Help page
~~~~~~~~~

To view the help page for each subcommand (which describes the subcommand and its parameters), use:

.. code-block:: bash
   
   $ minorg <subcommand> --help

To view the help page for the full programme, use:

.. code-block:: bash
   
   $ minorg full --help

Using the following will print a help page that lists common parameters and valid subcommand names:

.. code-block:: bash
   
   $ minorg --help

Aliases
~~~~~~~

Using a config.ini file, the command line version of MINORg allows users to supply short aliases in place of file names and/or combinations of parameters, as well as set default values for some parameters (such as the reference genome). For example, with the appropriate config.ini setup and lookup files, ``--reference TAIR10`` can be used in place of ``--assembly /path/to/TAIR10.fasta --annotation /path/to/TAIR10.gff3``, and ``--indv ler1`` can be used in place of ``--query /path/to/ler1.fasta``. For details on how to set up a config.ini file, see :ref:`Configuration:Configuration` in the :ref:`Tutorial:Tutorial` section.



Python package
--------------

Unlike the command line version, the Python package does not support aliases or preset parameter combinations beyond optionally reading default values from a config.ini file. Nevertheless, most arguments at the command line have equivalents as attributes in the class :class:`~minorg.MINORg.MINORg`. For a more comprehensive list of the similarities and differences, please refer to :ref:`Parameters:CLI vs Python`.
