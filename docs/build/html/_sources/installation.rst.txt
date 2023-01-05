Installation
============

Running in Docker
-----------------

A Docker image with MINORg (along with all of its dependencies) pre-installed on an linux system with Python 3.9 is available at: https://hub.docker.com/r/rlrq/minorg

You may pull the latest version using:

.. code-block::
   
   docker pull rlrq/minorg

This image comes bundled with the Cdd v3.18 database for domain search, which when unzipped will take up 4.8G. A Cdd database is required for domain search, but domain search is fully optional. You can pull a Docker image without the Cdd database using:

.. code-block::
   
   docker pull rlrq/minorg-lite

Note that some multithreading processes are incompatible with Docker and so are disabled in the Docker image, meaning that MINORg may take longer to run. If do not wish to use Docker or if you wish to install MINORg on your system instead, please refer to the rest of the Installation section for how to install MINORg.

Note for Windows users
----------------------

As one of the dependencies is not available for Windows, it is strongly recommended that Windows users use a Unix emulator such as "Ubuntu on Windows" (available for Windows 10 and up). You may follow the instructions at https://mafft.cbrc.jp/alignment/software/ubuntu_on_windows.html under 'Installation of Ubuntu' here for a succint guide OR https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#1-overview for a more thorough walkthrough. After, you may install the following dependencies using guides for 'Ubuntu' or 'Linux'.


Dependencies
------------

Python 3.6
++++++++++

MINORg is written and tested using Python 3.6.9, and minimally requires Python 3.6.

You may download and install Python here: https://www.python.org/downloads/


MAFFT
+++++

MAFFT is a multiple sequence alignment software that is capable of handling large datasets.

Follow one of the following links most applicable to your setup:

* Mac OS X: https://mafft.cbrc.jp/alignment/software/macosx.html
* Windows using Ubuntu emulator: https://mafft.cbrc.jp/alignment/software/ubuntu_on_windows.html
* Linux: https://mafft.cbrc.jp/alignment/software/linux.html
* From source: https://mafft.cbrc.jp/alignment/software/source.html

In all cases, do take note of the full path name to the MAFFT executable if the programme is not in your command-search path.

To determine if **mafft** is in your command-search path, execute

.. code-block::
   
   mafft --version

If it prints something like ::

  v7.427 (2019/Mar/24)

then it IS in your command-search path. If not, do take note of the full path name to the mafft executable as you may need to pass it to MINORg using ``--mafft <path to mafft>`` (at command line) or ``my_minorg.mafft = '<path to mafft>'`` (in Python).

BLAST+
++++++

BLAST+ is a suite of command-line tools that allows users to perform BLAST searches locally on their own machines. MINORg requires two of the tools: blastn and rpsblast (or rpsblast+)

Follow one of the following links most applicable to your setup:

* Mac OS X: https://www.ncbi.nlm.nih.gov/books/NBK569861/#_intro_Installation_MacOSX_
* Linux/Unix: https://www.ncbi.nlm.nih.gov/books/NBK52640/
* RedHat Linux: https://www.ncbi.nlm.nih.gov/books/NBK569861/#_intro_Installation_RedHat_Linux_
* From source: https://www.ncbi.nlm.nih.gov/books/NBK569861/#_intro_Installation_Source_tarball_


To determine if **blastn** is in your command-search path, execute at the command line::

  blastn -version

If it prints something like ::

  blastn: 2.6.0+
   Package: blast 2.6.0, build Jan 15 2017 17:12:27

then it IS in your command-search path. If not, do take note of the full path name to the blastn executable as you may need to pass it to MINORg using ``--blastn <path to blastn>`` (at command line) or ``my_minorg.blastn = '<path to blastn>'`` (in Python).

To determine if **rpsblast** is in your command-search path, execute at the command line::

  rpsblast -version

If it prints something like ::

  rpsblast: 2.11.0+
   Package: blast 2.11.0, build Oct 6 2020 03:24:05

then it IS in your command-search path. If not, execute the following to determine if your version of BLAST+ comes with 'rpsblast+' instead of 'rpsblast'::

  rpsblast+ -version

If it prints something like ::

  rpsblast+: 2.6.0+
   Package: blast 2.6.0, build Jan 15 2017 17:12:27

then it IS in your command-search path. If not, do take note of the full path name to the rpsblast (or rpsblast+ depending on version) executables as you may need to pass it to MINORg using ``--rpsblast <path to rpsblast or rpsblast+>`` (at command line) or ``my_minorg.rpsblast = '<path to rpsblast or rpsblast+>'`` (in Python).



BEDTools
++++++++

BEDTools is a "swiss-army knife of tools for a wide-range of genomic analysis tasks." It is not available for Windows.

You may follow the instructions to install it here: https://bedtools.readthedocs.io/en/latest/content/installation.html

To determine if BEDTools is in your command-search path, execute at the command line::

  bedtools --version

If it prints something like ::

  bedtools v2.26.0

then it IS in your command-search path. If not, do take note of the full path name to the directory containing BEDTools executables as you may need to pass it to MINORg using ``--bedtools <path>`` (at command line) or ``my_minorg.bedtools = '<path>'`` (in Python).


pysam dependencies
++++++++++++++++++

'pysam' is part of the 'pybedtools' package that MINORg uses. You DO NOT need to install pysam separately as it will be installed together with MINORg, but you MAY need to install some of its dependencies, as there are a handful that are not automatically installed with it. You may install them using your OS's package manager.

Some dependencies include (but may not be limited to):

* curses
  
  * Debian/Ubuntu: libncurses5-dev
  * RPM-based linux distributions: ncurses-devel
    
* zlib
  
  * Debian/Ubuntu: zlib1g-dev
  * RPM-based linux distributions or Cygwin: zlib-devel
    
* libbzip2
  
  * Debian/Ubuntu: libbz2-dev
  * RPM-based linux distributions or Cygwin: bzip2-devel
    
* liblzma
  
  * Debian/Ubuntu: liblzma-dev
  * RPM-based linux distributions or Cygwin: xz-devel
  * MacOS via Homebrew: xz

Do note that you may fail at installing pysam as part of MINORg's dependencies when installing MINORg according to :ref:`Installation:Install MINORg` if you are still missing some dependencies. Simply install the missing package described in the error message and try installing MINORg again.


Install MINORg
--------------

Test version can be installed from testpypi using:

.. code-block::

   python3 -m pip install --upgrade --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ minorg

