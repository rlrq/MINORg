## THINGS TO NOTE WHEN MIGRATING
## name changes:
## - config.reference_ext --> MINORg.assembly
## - config.annotation_ext --> MINORg.annotation
## removed:
## - config.bed_red, config.reference_indexed, bed_ext

## TODO: figure out which args to keep in MINORgCLI.args and which to incorporate directly into MINORg

## remove this for release
import sys
sys.path.append("/mnt/chaelab/rachelle/scripts/minorgpy")

import os
import copy
import shutil
import logging
import warnings
import itertools
import tempfile
import regex as re
from Bio import SearchIO, Seq

from minorg.log import MINORgLogger

from minorg.mafftcommandline_add import MafftCommandline
from minorg.filter_grna import make_target_feature_ranges_function, within_feature
from minorg.minimum_set import get_minimum_set

from minorg.index import IndexedFasta
from minorg.annotation import GFF
from minorg.parse_config import (
    Param,
    Params,
    get_val_none,
    get_val_default,
    parse_multiline_multikey_sdict,
    mv_dir_overwrite
)
# from minorg.filter_grna import make_exclude_function as mk_excl_fn
from minorg.constants import (
    INDV_GENOMES_ALL,
    INDV_GENOMES_NONE,
    INDV_GENOMES_REF,
    REFERENCED_ALL,
    REFERENCED_NONE
)
from minorg.extract_homologue import find_homologue_multiindv
from minorg.extend_reference import extend_reference
from minorg.reference import AnnotatedFasta
from minorg.functions import (
    split_none,
    cat_files,
    non_string_iter,
    assign_alias,
    blast,
    adjusted_feature_ranges,
    ranges_union,
    ranges_intersect,
    ranges_subtract
)
from minorg.generate_grna import find_multi_gRNA
from minorg.fasta import dict_to_fasta, fasta_to_dict
from minorg.filter_grna import mask_identical
from minorg.grna import gRNAHits
from minorg.pam import PAM

from typing import Optional, Type, Dict

LOGGING_LEVEL = logging.DEBUG

#######################
##  PARSE FUNCTIONS  ##
#######################

def parse_lookup(iterable, lookup, return_first = False):
    '''
    commas not allowed in values
    '''
    iterable = split_none(iterable)
    if (not isinstance(lookup, dict)) or iterable is None: return None
    mapped = [lookup.get(val, val) for val in iterable]
    if return_first: return mapped[0]
    else: return mapped

def valid_readable_file(pathname):
    pathname = os.path.abspath(pathname)
    if not os.path.exists(pathname):
        raise InvalidPath(pathname)
    if not os.path.isfile(pathname):
        raise InvalidFile(pathname)
    if not os.access(pathname, os.R_OK):
        raise UnreadableFile(pathname)
    return pathname

def valid_aliases(aliases, lookup: dict, raise_error: bool = True, message: str = None,
                  param: Param = None, none_value = None, all_value = None, clear_value = None,
                  display_cmd = None, additional_message = None, return_mapping = False):
    '''
    Generates appropriate error + message if alias(es) is/are invalid.
    Throws the error generated if requested (i.e. raise_error = True).
    Requires alias(es) (str of single alias or iterable of multiple aliases) + set lookup dictionary
    '''
    if isinstance(aliases, str): aliases = [aliases]
    ## return None without raising error if none_value is only value provided
    if none_value is not None and len(aliases) == 1 and aliases[0] == none_value: return None
    ## identify unknown aliases
    unknown = []
    mapping = {}
    for alias in aliases:
        if alias not in lookup and alias != all_value and alias != clear_value: unknown.append(alias)
        else: mapping[alias] = lookup[alias]
    if unknown: ## if there are unmappable values (i.e. unknown aliases, generate error message)
        if (message is None and param is not None and param.description is not None):
            message1 = lambda x: ( f"Unrecognised {param.description}" + \
                                   f" {'/'.join(param.names)}:" + \
                                   f" {', '.join(sorted(x))}" )
            if display_cmd is not None:
                message2 = lambda x: ( message1(x) +
                                       (f"\nFor a list of valid {param.description}"
                                        f"{'' if param.alias_value_description is None else ' and ' + param.alias_value_description}"
                                        f", execute 'minorg {display_cmd}'." ) )
            else: message2 = lambda x: message1(x)
            message3 = lambda x: ( message2(x) +
                                   (f"\nYou may update the list of valid {param.description}"
                                    f"{'' if not param.alias_value_description else ' and ' + param.alias_value_description} by modifying the config file located at {params.config_file}"))
            if none_value is not None:
                message4 = lambda x: ( message3(x) +
                                       (f"\nNote that the default empty value '{none_value}' should not"
                                        " be used with other values.") )
            else: message4 = lambda x: message3(x)
            if additional_message is not None: message = lambda x: message4(x) + '\n' + additional_message
            else: message = message4
        error = typer.BadParameter( f"Unrecognised alias(es): {sorted(unknown)}\n"
                                    if message is None else (message(sorted(unknown)) + '\n') )
        if raise_error: raise error
        elif return_mapping: return mapping
        else: return None
    elif return_mapping: return mapping
    else: return None


###############
##  CLASSES  ##
###############

class PathHandler:
    """
    Tracks new files/directories and temporary files.
    
    Attributes:
        tmp (bool): whether temporary directory is used
        keep_tmp (bool): whether to retain temporary files/directories
        out_dir (str): absolute path to final output directory
        new_dirs (list of str): paths of newly created directories
        tmp_files (list of str): paths of temporary files/directories
        tmp_dir (str): path of temporary directory
        directory (str): path of directory being written into.
            If tmp=True, this will point to a temporary directory.
            Else, this will be the output directory.
    """
    def __init__(self, tmp = False, keep_tmp = False, directory = None):
        """
        Create a PathHandler object.
        
        Arguments:
            tmp (bool): write files to temporary directory
                (to be deleted or moved to final output directory using resolve)
            keep_tmp (bool): retain temporary files/directories
            directory (str): final output directory (will be directly written into if tmp=False)
        """
        self.tracebacklimit = (None if LOGGING_LEVEL <= logging.DEBUG else 0)
        self.out_dir = None if directory is None else os.path.abspath(directory) ## final output directory
        self.tmp = (tmp if self.out_dir is not None else True)
        self.keep_tmp = keep_tmp
        self.new_dirs = []
        self.tmp_files = set()
        ## create temporary directory
        if self.tmp:
            ## contents to be moved to self.out_dir upon fulfilment of command using self.resolve
            self.tmpdir = tempfile.mkdtemp()
            ## self.directory is used outside PathHandler objects.
            ## PathHandler object will handle cleanup if self.directory points to self.tmpdir
            self.directory = self.tmpdir
            print(f"tmpdir: {self.tmpdir}")
        elif self.out_dir is not None:
            self.tmpdir = None
            self.directory = self.out_dir
            self.mkdir(self.directory)
            print(f"Files will be written directly in: {self.directory}")
        else:
            print("Directory (directory = <path>) or temporary directory (tmp = True) required.")
            exit
    @property
    def tracebacklimit(self):
        return self._tracebacklimit
    @tracebacklimit.setter
    def tracebacklimit(self, limit):
        self._tracebacklimit = limit
        import sys
        sys.tracebacklimit = self._tracebacklimit
    ## update output directory
    def mv(self, directory):
        """
        Update output directory.
        
        Contents of the current output directory will also be moved to the new directory.
        
        Arguments:
            directory (str): required, path to new output directory
        """
        if directory is None: return
        directory = os.path.abspath(directory)
        if directory == self.out_dir: return
        ## if files are being written directly to output directory, move to new destination
        if self.directory == self.out_dir:
            mv_dir_overwrite(self.directory, directory)
            self.tmp_files = set(re.sub('^' + re.escape(self.directory), directory, fname)
                                 for fname in self.tmp_files)
            self.new_dirs = [re.sub('^' + re.escape(self.directory), directory, dirname)
                             for dirname in self.new_dirs]
            self.directory = directory
        self.out_dir = directory
    ## make directory (or directories)
    def mkdir(self, *directories, tmp = False):
        """
        Create directory/directories.
        
        Arguments:
            *directories (str): path
            tmp (bool): mark directory as temporary (for deletion when self.rm_tmpfiles is called)
        """
        output = []
        for directory in directories:
            if not os.path.isabs(directory):
                directory = os.path.abspath(os.path.join(self.directory, directory))
            output.append(directory)
            if not os.path.exists(directory) and not directory in self.new_dirs:
                os.mkdir(directory)
                self.new_dirs.append(directory)
                if tmp: self.tmp_files.add(directory)
        return output[0] if len(directories) == 1 else output
    ## generate fname (usage: self.mkfname('ref', 'tmp.txt'))
    def mkfname(self, *path, tmp = False):
        """
        Generate path.
        
        If path provided is not absolute, self.directory is assumed to be the root directory.
        
        Arguments:
            *path (str): required, path to output file
                (e.g. self.mkfname('tmp', 'tmp.fasta') --> <self.directory>/tmp/tmp.fasta)
            tmp (bool): mark file as temporary (for deletion when self.rm_tmpfiles is called)
        
        Returns
        -------
        str
            path
        """
        if os.path.isabs(os.path.join(*path)): path = os.path.join(*path)
        else: path = os.path.join(self.directory, *path)
        if tmp: self.tmp_files.add(path)
        return path
    ## mkfname, except it also creates path to directory + empty file if they don't already exist
    def reserve_fname(self, *path, tmp = False, newfile = False):
        """
        Generate new file.
        
        Operates exactly as :meth:`PathHandler.mkfname`, with the additional options of creating an empty file or clearing an existing file.
        
        Arguments:
            *path (str): path to output file
            tmp (bool): mark file as temporary (for deletion when self.rm_tmpfiles is called)
            newfile (bool): clear an existing file at destination if it already exists
        
        Returns
        -------
        str
            path
        """
        fname = self.mkfname(*path)
        fdir = os.path.dirname(fname)
        os.makedirs(fdir, exist_ok = True)
        if newfile:
            open(fname, 'w').close()
        else:
            open(fname, 'a').close()
        if tmp: self.tmp_files.add(fname)
        return fname
    def rm_tmpfiles(self, *fnames):
        """
        Delete all files and directories marked as temporary (in self.tmp_files).
        
        Directories will only be deleted if they are empty.
        
        Arguments:
            *fnames (str): optional, path. If provided, only the specified files/directories will be deleted if they are also marked as temporary. Otherwise, all temporary files and directories are deleted.
        """
        if not self.keep_tmp:
            to_remove = set(fnames if fnames else self.tmp_files)
            files_to_remove = [fname for fname in to_remove if
                               (os.path.exists(fname) and os.path.isfile(fname))]
            dir_to_remove = [fname for fname in to_remove if
                             (os.path.exists(fname) and os.path.isdir(fname))]
            for fname in files_to_remove:
                os.remove(fname)
            for dirname in dir_to_remove:
                if not os.listdir(dirname):
                    os.rmdir(dirname)
            ## update tmp_files
            self.tmp_files -= to_remove
        return
    ## move output files to final output directory
    def resolve(self):
        """
        If self.tmp, move output files to final output directory.
        """
        ## move files to final destination (lol)
        if self.tmp:
            if self.out_dir:
                ## create final output directory if doesn't exist
                if not os.path.exists(self.out_dir):
                    os.makedirs(self.out_dir, exist_ok = True)
                ## remove tmp files
                self.rm_tmpfiles()
                ## copy items
                mv_dir_overwrite(self.tmpdir, self.out_dir)
                print(f"Output files have been generated in {self.out_dir}")
            ## remove tmpdir
            shutil.rmtree(self.tmpdir)

# _verbose = True
# _directory = "/mnt/chaelab/rachelle/scripts/minorgpy/test/minorg"
# _config = "/mnt/chaelab/rachelle/scripts/minorgpy/config.ini"
# _prefix = "test"
# _tmp = False
# _keep_tmp = True
# _keep_on_crash = True
# _thread = 1
# _ext_cds = "/mnt/chaelab/rachelle/data/NLR/RPW8/merged-KZ10_Ms0_TW_onlyRPW8_curated-Ms0_onlyRPW8_CDSonly.fasta"
# _ext_genome = "/mnt/chaelab/rachelle/data/NLR/RPW8/merged-KZ10_Ms0_TW_onlyRPW8_curated-Ms0_onlyRPW8_CDScomplete.fasta"
# _query = {10015: "/mnt/chaelab/shared/anna_lena/contigs/10015.contigs.fasta"}
# _query_reference = True

# x = MINORg(verbose = _verbose, directory = _directory, config = _config, prefix = _prefix,
#            tmp = _tmp, keep_tmp = _keep_tmp, thread = _thread)
# x.genes = ["AT3G44400", "AT3G44480", "AT3G44670", "AT3G44630"]
# x.subset_annotation()
# # [y.get_attr("ID") for y in x.reference["TAIR10"].annotation.get_subfeatures("AT3G44400")]
# x.ext_cds = _ext_cds
# x.ext_gene = _ext_genome
# x.extend_reference()
# # print("extended genes:", [y.get_attr("ID") for y in x.reference["Extended"].annotation.get_features("gene")])
# # print("mRNA pep:", x.get_reference_seq("CDS", adj_dir = True, translate = True, isoform_lvl = "mRNA"))
# # x.get_reference_seq("CDS", adj_dir = True, translate = True, seqid_template = "$source|$isoform|$feature|$n", apply_template_to_dict = True)
# x.query = _query
# x.query_reference = _query_reference
# x.target()
# x.generate_grna()
# x.filter_background()
# x.grna_hits.filter_seqs("background")
# x.filter_gc()
# x.grna_hits.filter_seqs("GC")
# x.grna_hits.filter_seqs("background", "GC")
# x.align_reference_and_targets()
# x.filter_feature()
# x.grna_hits.filter_hits("feature")
# x.minimumset()



# _verbose = True
# _directory = "/mnt/chaelab/rachelle/scripts/minorgpy/test/minorg_nrg"
# _config = "/mnt/chaelab/rachelle/scripts/minorgpy/config.ini"
# _prefix = "test"
# _tmp = False
# _keep_tmp = True
# _keep_on_crash = True
# _thread = 1

# # x = MINORg(verbose = _verbose, directory = _directory, config = _config, prefix = _prefix,
# #            tmp = _tmp, keep_tmp = _keep_tmp, thread = _thread)

# from MINORg import MINORg
# x = MINORg(directory = "/mnt/chaelab/rachelle/scripts/minorgpy/test/minorg_nrg", prefix = "test",
#            tmp = False, keep_tmp = True, thread = 1)
# x.genes = ["AT5G66900", "AL8G44500.v2.1"]
# x.add_reference("TAIR10", "/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta", "/mnt/chaelab/rachelle/scripts/minorgpy/test/full/test_full_NRG1_1/reduced_ann/reduced.TAIR10.gff", replace = True)
# x.add_reference("araly2", "/mnt/chaelab/rachelle/data/Alyrata/v2.1/Alyrata_384_v1.fa", "/mnt/chaelab/rachelle/scripts/minorgpy/test/full/test_full_NRG1_1/reduced_ann/reduced.araly2.gff")
# x.subset_annotation()
# # [y.get_attr("ID") for y in x.reference["TAIR10"].annotation.get_subfeatures(*x.genes)]
# x.query_reference = True
# x.pssm_ids = "366375" # nb-arc
# x.rpsblast = "rpsblast+"
# x.db = "/mnt/chaelab/shared/blastdb/Cdd.v3.18/Cdd"
# x.seq()
# x.grna()

# x.screen_reference = True
# x.filter_background()
# x.grna_hits.filter_seqs("background")
# ## 352 without domain, 118 in NB-ARC 366375
# x.filter_gc()
# x.grna_hits.filter_seqs("GC")
# ## 370 without domain, 114 in NB-ARC 366375
# x.grna_hits.filter_seqs("background", "GC")
# ## 321 without domain, 106 in NB-ARC 366375
# x.align_reference_and_targets()
# x.filter_feature()
# x.grna_hits.filter_hits("feature")
# ## 344 without domain, 104 in NB-ARC 366375
# x.passed_grna
# ## 278 without domain, 89 in NB-aRC 366375
# x.minimumset()

class MINORg (PathHandler):
    """
    Tracks parameters, intermediate files, and output file for reuse.
    
    >>> from minorg.MINORg import MINORg
    >>> my_minorg = MINORg(directory = "/my/output/directory", 
                           prefix = "test", tmp = False, keep_tmp = True, thread = 1)
    >>> my_minorg.add_reference("TAIR10", "/path/to/TAIR10_Chr.all.fasta", "/path/to/TAIR10_GFF3.genes.gff", replace = True)
    >>> my_minorg.add_reference("araly2", "/path/to/Alyrata_384_v1.fa", "/path/to/Alyrata_384_v2.1.gene.gff3")
    >>> my_minorg.genes = ["AT5G66900", "AL8G44500.v2.1"]
    >>> my_minorg.subset_annotation()
    >>> [y.get_attr("ID") for y in my_minorg.reference["TAIR10"].annotation.get_subfeatures(*x.genes)]
    [['AT5G66900.1']]
    >>> [y.get_attr("ID") for y in my_minorg.reference["TAIR10"].annotation.get_subfeatures(*x.genes)]
    [['AL8G44500.t1.v2.1']]
    >>> my_minorg.query_reference = True
    >>> my_minorg.seq()
    >>> my_minorg.grna()
    >>> my_minorg.screen_reference = True
    >>> my_minorg.filter_background()
    >>> my_minorg.grna_hits.filter_seqs("background")
    gRNAHits(len = 352)
    >>> my_minorg.filter_gc()
    >>> my_minorg.grna_hits.filter_seqs("GC")
    gRNAHits(len = 370)
    >>> my_minorg.grna_hits.filter_seqs("background", "GC")
    gRNAHits(len = 321)
    >>> my_minorg.filter_feature() ## by default, MINORg only retains gRNA in CDS
    Max acceptable insertion length: 15
    >>> my_minorg.grna_hits.filter_hits("feature")
    gRNAHits(len = 344)
    >>> my_minorg.passed_grna
    gRNAHits(len = 278)
    >>> my_minorg.minimumset()
    >>> my_minorg.resolve()
    
    Attributes:
        prefix (str): [general] prefix for output directories and files
        thread (int): [general] maximum number of threads for parallel processing

        blastn (str): [executable] path to blastn executable
        rpsblast (str): [executable] path to rpsblast or rpsblast+ executable
        mafft (str): [executable] path to mafft executable

        db (str): [RPS-BLAST option] path to local RPS-BLAST database
        remote_rps (bool): [RPS-BLAST option] use remote database instead of local database for RPS-BLAST

        pssm_ids (list of str):
            [seq] list of Pssm-Ids of domain(s) for homologue search. 
            If multiple Pssm-Ids are provided, overlapping domains will be merged.
        domain_name (str): human-readable domain name used in sequence and file names in place of Pssm-Ids
        genes (list or str): [seq] list of target gene IDs
        query_reference (bool): [seq] include reference genes as targets

        minlen (int): [seq: homologue] minimum homologue length (bp)
        minid (float): [seq: homologue] minimum hit % identity
        mincdslen (int):
            [seq: homologue] minimum number of bases in homologue aligned with reference gene CDS
        merge_within (int): [seq: homologue] maximum distance (bp) between hits for merging
        check_recip (bool): [seq: homologue] execute reciprocal check
        relax_recip (bool): [seq: homologue] execute relaxed reciprocal check
        check_id_before_merge (bool):
            [seq: homologue] filter out hits by % identity before merging into potential homologues

        length (int): [grna] gRNA length (bp)
        pam (str): [grna] PAM pattern

        screen_reference (bool): [filter: background] include reference genome for screening
        ot_pamless (bool): [filter: background] ignore absence of PAM when assessing off-target gRNA hits
        offtarget (func):
            [filter: background] function that accepts Biopython's HSP and QueryResult objects and determines whether an off-target gRNA hit is problematic.
            If not set by user, ot_mismatch and ot_gap will be used
            ## TODO!!! integrate this somehow?
        ot_mismatch (int):
            [filter: background] minimum number of mismatches allowed for off-target gRNA hits
        ot_gap (int): [filter: background] minimum number of gaps allowed for off-target gRNA hits

        gc_min (float):
            [filter: GC] minimum GC content (between 0 and 1, where 0 is no GC content and 1 is all GC)
        gc_max (float):
            [filter: GC] maximum GC content (betweew 0 and 1, where 0 is no GC content and 1 is all GC)

        feature (str): [filter: feature] GFF3 feature within which gRNA are to be designed (default="CDS")
        max_insertion (int): [filter: feature] maximum allowable insertion size in feature (bp)
        min_within_n (int):
            [filter: feature] minimum number of reference genes which feature a gRNA must align within
        min_within_fraction (float): 
            [filter: feature] minimum fraction of reference genes which feature a gRNA must align within
            (between 0 and 1, where 0 is none and 1 is all; if 0, min_within_n will be set to 1)

        sets (int): [minimumset] number of sets to generate
        auto (bool): [minimumset] generate sets without requiring manual user confirmation for each set
        accept_invalid (bool): [minimumset] score 'NA' as 'pass'
        accept_feature_unknown (bool): [minimumset] score 'NA' as 'pass' for feature check
        accept_invalid_field (bool): [minimumset] score 'NA' as 'pass' if all entries for a check are 'NA'
        pass_map (str):
            [minimumset] path to output .map file for gRNA that pass all valid checks
            (autogenerated by MINORg if not provided)
        pass_fasta (str):
            [minimumset] path to output .fasta file for gRNA that pass all valid checks
            (autogenerated by MINORg if not provided)
        final_map (str):
            [minimumset] path to output .map file for final gRNA set(s)
            (autogenerated by MINORg if not provided)
        final_fasta (str):
            [minimumset] path to output .fasta file for final gRNA set(s)
            (autogenerated by MINORg if not provided)
    """
    def __init__(self, directory = None, config = None, prefix = None,
                 thread = None, keep_tmp = False, **kwargs):
        """Create a MINORg object.
        
        Arguments:
            directory (str): path to output directory
            config (str): path to config.ini file
            prefix (str): prefix for output directories and files
            thread (int): maximum number of threads for parallel processing
            keep_tmp (bool): retain temporary files
            **kwargs: other arguments supplied to parent class :class:`PathHandler`
        
        Returns
        -------
        :class:`~minorg.MINORg.MINORg`
            a MINORg object
        """
        self.verbose = False
        super().__init__(directory = directory, keep_tmp = keep_tmp, **kwargs)
        
        ## parse config.ini file if provided (if None, default values will be stored in Params obj)
        self.params = Params(config)
        
        ## handle log file
        self.logfile = MINORgLogger(level = LOGGING_LEVEL)
        
        ## parameter things
        ## status tracking
        self.subset_gid = []
        self.subset_ref = []
        ## general output params
        self.prefix = get_val_default(prefix, default = self.params.prefix.default)
        ## general execution params
        self.thread = get_val_default(thread, default = self.params.thread.default)
        ## binaries/executables
        self.rpsblast = self.params.rpsblast.default
        self.blastn = self.params.blastn.default
        self.mafft = self.params.mafft.default
        ## RPS-BLAST database
        self.db = self.params.db.default
        self.remote_rps = self.params.remote_rps.default
        ## data/annotation
        ## (sets and aliases not supported for MINORg object beyond retrieval of default reference)
        self.reference = {} ## use MINORg.clear_reference(), MINORg.add_reference(), and MINORg.remove_reference() to update this attribute
        self._parse_reference()
        ## annotation format options
        ## seqid template stuff
        self.ref_seqid_template = "Reference|${source}|${domain}|${feature}|${complete}|${gene}|${isoform}"
        self.ref_seqid_gene_pattern = "(?<=\\|stitched\\||\\|complete\\|).+(?=\\|\d+-\d+)"
        self.ref_seqid_feature_pattern = "(?<=\\|)[^|]+(?=\\|stitched\\||\\|complete\\|)"
        self.ref_seqid_source_pattern = "(?<=^Reference\\|)[^|]+(?=\\|)"
        ## homologue params
        self.pssm_ids = [] ## if multiple pssm_ids are provided, overlapping domains will be combined and output will not distinguish between one pssm_id or another
        self._domain_name = None ## human-readable name for self.pssm_ids
        self.gff_domain = None ## gff file for domains, generated by MINORg
        ## previously self.query_map storing [(<alias>, <path>)]. Now: {<alias>: <path>}
        ## self.query_map is now @property
        self.query = {}
        self.genes = [] ## [<geneA>, <geneB>]
        self.query_reference = False ## if True, reference genes will also be queried for gRNA
        self.ref_gene = None ## fasta file, generated by MINORg
        self.ref_cds = None ## fasta file, generated by MINORg
        self.minlen = self.params.minlen.default
        self.minid = self.params.minid.default
        self.mincdslen = self.params.mincdslen.default
        self._check_recip = self.params.check_recip.default
        self.relax_recip = self.params.relax_recip.default
        self.check_id_before_merge = self.params.check_id_before_merge.default
        self.merge_within = self.params.merge_within.default
        ## generate gRNA params
        self.target = None ## output by self.seq()
        self._length = self.params.length.default
        self._pam = PAM(pam = self.params.pam.default, gRNA_length = self.params.length.default)
        ## filter params
        self.mask = []
        self.masked = {}
        self.grna_hits = None ## gRNAHits object, generated by MINORg.generate_grna()
        self.grna_fasta = None ## fasta file, generated by MINORg.write_all_grna_fasta() if not specified by user
        self.grna_map = None ## .map file, generated by MINORg.write_all_grna_map() if not specified by user
        self.screen_reference = True
        self._offtarget = None ## user-provided function that takes in pam, QueryResult obj, and HSP object and returns whether the off-target hit is acceptable. If not set, defaults to using self.ot_mismatch and self.ot_gap (combined with OR) to set lower limits for acceptable off-target
        self.background = {} ## {<alias>: <fname>}b
        self.ot_mismatch = self.params.ot_mismatch.default
        self.ot_gap = self.params.ot_gap.default
        self.ot_pamless = self.params.ot_pamless.default
        self.gc_min = self.params.gc_min.default
        self.gc_max = self.params.gc_max.default
        self.alignment = None ## fasta file, generated by MINORg
        self.max_insertion = self.params.max_insertion.default
        self.min_within_n = self.params.min_within_n.default
        self.min_within_fraction = self.params.min_within_fraction.default
        self._feature = self.params.feature.default
        ## minimum set
        self.sets = 1
        self.auto = True
        self.pass_map = None ## .map file, generated by MINORg.write_pass_grna_map() if not specified
        self.pass_fasta = None ## fasta file, generated by MINORg.write_pass_grna_fasta() if not specified
        self.final_map = None ## .map file, generated by MINORg.write_all_grna_map() if not specified
        self.final_fasta = None ## fasta file, generated by MINORg.write_all_grna_fasta() if not specified
        self.accept_invalid = False
        self.accept_feature_unknown = False
        self.accept_invalid_field = True

        ## move logfile
        self.logfile.update_filename(self.mkfname(f".log"))
    
    ## getters
    @property
    def query_map(self):
        """
        Mapping of alias to FASTA for queries.
        
        :getter: Returns list of (<alias>, </path/to/FASTA>)
        :type: list
        """
        return [(str(k), v) for k, v in self.query.items()]
    @property
    def domain_name(self):
        """
        Domain name for use in sequence IDs and file names
        
        :getter: Returns domain name if set, else 'gene' if not :attr:`MINORg.pssm_ids`, else dash-separated :attr:`MINORg.pssm_ids`
        :setter: Sets domain name
        :type: str
        """
        if self._domain_name is not None: return self._domain_name
        elif not self.pssm_ids: return "gene"
        elif non_string_iter(self.pssm_ids): return '-'.join(self.pssm_ids)
        else: return self.pssm_ids
    @property
    def check_recip(self):
        """
        Whether to execute reciprocal check
        
        :getter: Returns True if MINORg.check_recip=True OR MINORg.relax_recip=True
        :setter: Sets this check
        :type: bool
        """
        ## return True if self.relax_recip is raised
        return self.relax_recip or self._check_recip
    @property
    def assemblies(self):
        """
        Reference assemblies
        
        :getter: Returns {<alias>: </path/to/FASTA>}
        :type: dict
        """
        return {alias: ref.fasta for alias, ref in self.reference.items()}
    @property
    def annotations(self):
        """
        Reference annotations
        
        :getter: Returns {<alias>: </path/to/GFF3>}
        :type: dict
        """
        return {alias: ref._gff for alias, ref in self.reference.items()}
    @property
    def attr_mods(self):
        """
        Attribute modifications
        
        :getter: Returns {<reference alias>: <dict of attribute modifications>}
        :type: dict
        """
        return {alias: ref.attr_mod for alias, ref in self.reference.items()}
    @property
    def offtarget(self): ## TODO!!! integrate this somehow?
        """
        Function to assess off-target goodness of alignment
        
        :getter: Returns function that accepts HSP and QueryResult
        :type: func
        """
        if self._offtarget:
            return self._offtarget
        else:
            return mk_excl_fn(pam = self.pam, outside_targets = None, indexed_fa = self.assembly,
                              max_mismatch = max(1, self.ot_mismatch)-1,
                              max_gap = max(1, self.ot_gap)-1)
    @property
    def pam(self):
        """
        PAM pattern
        
        :getter: Returns expanded PAM pattern (includes explicit position of gRNA; e.g. '.NGG')
        :setter: Sets PAM pattern and parses it into expanded PAM pattern
        :type: str
        """
        return self._pam
    @property
    def length(self):
        """
        gRNA length
        
        :getter: Returns gRNA length (bp)
        :setter: Sets gRNA length and updates MINORg.PAM
        :type: int
        """
        return self._length
    @property
    def feature(self): ## returns first feature if multiple are present
        """
        Target feature
        
        :getter: Returns first feature if multiple have been provided
        :setter: Sets target feature
        :type: str
        """
        if non_string_iter(self._feature):
            if len(self._feature) == 1: return tuple(self._feature)[0]
            elif len(self._feature) == 0: return None
            else: raise MessageError("Too many features")
        else: return self._feature
    @property
    def features(self):
        """
        Target features
        
        :getter: Returns list of target feature(s)
        :type: list
        """
        if non_string_iter(self._feature): return self._feature
        else: return [self._feature]
    @property
    def passed_grna(self):
        """
        gRNA that have passed background, GC, and feature checks
        
        Relevant attributes: accept_invalid, accept_invalid_field
        
        :getter: Returns gRNA that have passed checks
        :type: gRNAHits
        """
        filtered = self.grna_hits.filter_seqs("background", "GC",
                                              accept_invalid = self.accept_invalid,
                                              accept_invalid_field = self.accept_invalid_field)
        filtered = filtered.filter_hits("CDS", "feature",
                                        accept_invalid = self.accept_invalid,
                                        accept_invalid_field = self.accept_invalid_field,
                                        exclude_empty_seqs=True)
        return filtered
    
    ## setters
    @domain_name.setter
    def domain_name(self, val):
        self._domain_name = val
    @offtarget.setter
    def offtarget(self, f):
        self._offtarget = f
    @pam.setter
    def pam(self, pam):
        self._pam = PAM(pam = pam, gRNA_length = self.length)
    @length.setter
    def length(self, length):
        self._length = length
        self.pam = self.pam
    @feature.setter
    def feature(self, val):
        if non_string_iter(val): self._feature = val
        else: self._feature = [val]
    @check_recip.setter
    def check_recip(self, val):
        self._check_recip = val
    
    ## update 'mv' so that it also moves logfile
    def mv(self, directory):
        """
        Update output directory.
        
        Contents of the current output directory will also be moved to the new directory.
        
        Arguments:
            directory (str): required, path to new output directory
        """
        super().mv(directory)
        ## move logfile
        self.logfile.move(self.directory, self.prefix)
    
    ## update mkfname so it adds a prefix to each file
    def mkfname(self, *path, tmp = False, prefix = True):
        """
        Generate new file name.
        
        If path provided is not absolute, self.directory is assumed to be the root directory.
        
        Arguments:
            *path (str): required, path to output file (e.g. self.mkfname('tmp', 'tmp.fasta'))
            tmp (bool): mark file as temporary (for deletion when self.rm_tmpfiles is called)
            prefix (bool): prefix self.prefix to basename
        """
        if path and prefix:
            path_last = path[-1]
            fname = self.prefix + '_' + os.path.basename(path_last)
            path = path[:-1] + type(path)([os.path.join(os.path.dirname(path_last), fname)])
        return super().mkfname(*path, tmp = tmp)
    
    def results_fname(self, *path, reserve = False, newfile = False, **kwargs):
        """
        Generate new file name in results directory (<self.directory>/<self.prefix>).
        
        Operates exactly as :meth:`MINORg.mkfname`, with the additional options of creating an empty file or clearing an existing file.
        
        Arguments:
            *path (str): path to output file
            reserve (bool): create empty file at destination if it does not exist
            newfile (bool): clear an existing file at destination if it already exists
            **kwargs: other arguments supplied to :meth:`MINORg.mkfname`
        """
        fpath = os.path.dirname(self.mkfname(self.prefix, *path))
        ## create directory if directory path does not currently exist
        if not os.path.exists(fpath): self.mkdir(fpath)
        ## return relevant path
        if reserve or newfile: return self.reserve_fname(self.prefix, *path, **kwargs, newfile = newfile)
        else: return self.mkfname(self.prefix, *path, **kwargs)
    
    def get_ref_seqid(self, seqid, attr):
        """
        Extract attribute of sequence from sequence name (for sequences generated by self.target)
        
        Arguments:
            seqid (str): required, sequence name
            attr (str): required, attribute to return.
                Valid attributes: 'gene', 'feature', 'source'.
        
        Returns
        -------
        str
        """
        if attr == "gene": return re.search(self.ref_seqid_gene_pattern, seqid).group(0)
        elif attr == "feature": return re.search(self.ref_seqid_feature_pattern, seqid).group(0)
        elif attr == "source": return re.search(self.ref_seqid_source_pattern, seqid).group(0)
        else: raise MessageError(f"Unknown attribute: {attr}")
    
    ####################
    ##  ADD FUNCIONS  ##
    ####################
    
    def add_background(self, fname, alias = None):
        """
        Add file for background filter.
        
        Arguments:
            fname (str): required, path to file
            alias (str): optional, alias for background file.
                Used when writing mask report and for removing background.
        """
        if not alias:
            alias = f"bg_{str(len(self.background) + 1).zfill(3)}"
        self.background[alias] = IndexedFasta(fname)
        return
    
    def remove_background(self, alias):
        """
        Remove file for background filter.
        
        Arguments:
            alias (str): required, alias of background file to remove
        """
        if alias not in self.background:
            raise MessageError(f"'{alias}' is not a valid background alias.")
        del self.background[alias]
        return
    
    def parse_grna_map_from_file(self, fname):
        """
        Read candidate gRNA from .map file output by MINORg.
        
        Arguments:
            fname (str): required, path to MINORg .map file
        """
        grna_hits = gRNAHits()
        grna_hits.parse_from_mapping(fname, version = None)
        self.grna_hits = grna_hits
    
    def rename_grna(self, fasta):
        """
        Rename candidate gRNA according to FASTA file.
        
        Arguments:
            fasta (str): required, path to FASTA file
        """
        self.grna_hits.rename_seqs(fasta)
    
    def write_all_grna_fasta(self, fout = None):
        """
        Write FASTA file of all candidate gRNA.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.directory if not provided.
        """
        if fout is None: fout = get_val_default(self.grna_fasta, self.results_fname("gRNA_all.fasta"))
        self.grna_fasta = fout
        self.grna_hits.write_fasta(fout, write_all = True)
    
    def write_all_grna_map(self, fout = None, write_checks = True):
        """
        Write .map file of all candidate gRNA.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.directory if not provided.
            write_checks (bool): write all check statuses
        """
        if fout is None: fout = get_val_default(self.grna_map, self.results_fname("gRNA_all.map"))
        self.grna_map = fout
        self.grna_hits.write_mapping(fout, version = 2, write_all = True, write_checks = write_checks)
    
    def write_pass_grna_fasta(self, fout = None):
        """
        Write FASTA file of gRNA that pass all valid checks.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.directory if not provided.
        """
        if fout is None: fout = get_val_default(self.pass_fasta, self.results_fname("gRNA_pass.fasta"))
        self.passed_grna.write_fasta(fout, write_all = True)
    
    def write_pass_grna_map(self, fout = None):
        """
        Write .map file of gRNA that pass all valid checks.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.directory if not provided.
        """
        if fout is None: fout = get_val_default(self.pass_map, self.results_fname("gRNA_pass.map"))
        self.passed_grna.write_mapping(fout, version = 2, write_all = True, write_checks = False)
    
    ##########################
    ##  PARSERS FOR CONFIG  ##
    ##########################
    
    def _parse_reference(self, reference = None, reference_set = None):
        """
        Parse reference assembly and annotation read from config.ini file into self.params.
        
        Arguments:
            reference (str): alias of reference genome
            reference_set (str): alias of or path to reference set file
        """
        ## parse reference lookup
        val_ref_set = parse_lookup(get_val_default(reference_set, self.params.reference_set.default),
                                   self.params.reference_sets, return_first = True)
        if val_ref_set is not None:
            if os.path.exists(os.path.abspath(val_ref_set)):
                try:
                    with open(val_ref_set, 'r') as f:
                        reference_aliases = {k: v.split('\t') for k, v in
                                             parse_multiline_multikey_sdict(f.read(),
                                                                            kv_sep = '\t').items()}
                        reference_aliases = {k: v + (4 - len(v))*['']
                                             for k, v in reference_aliases.items()}
                except:
                    raise InputFormatError(error_src = "file", hint = f"File: {val_ref_set}")
            else:
                warnings.warn( ( "Unrecognised reference set lookup file alias or non-existent file:"
                                 f" {val_ref_set}. You may manually provide reference genomes using"
                                 " MINORg.add_reference(<alias>, <path to FASTA>, <path to GFF3>") )
        ## get reference
        val_ref = get_val_default(reference, self.params.reference.default)
        if val_ref:
            none_val = INDV_GENOMES_NONE
            mapped = valid_aliases(aliases = val_ref, lookup = reference_aliases,
                                   none_value = INDV_GENOMES_NONE, all_value = INDV_GENOMES_ALL,
                                   param = self.params.reference, return_mapping = True)
            for alias, ref in mapped.items():
                fasta, ann, genetic_code, attr_mod = reference_aliases[alias]
                self.add_reference(alias, fasta, ann, genetic_code = genetic_code,
                                   attr_mod = self.params.parse_attr_mod(attr_mod))
        return
    
    ############################
    ##  REFERENCE/ANNOTATION  ##
    ##       FUNCTIONS        ##
    ############################
    
    ## takes path to FASTA and GFF/BED file and indexes them before storing.
    def add_reference(self, alias, fasta, annotation, genetic_code = 1, attr_mod = {},
                      memsave = True, replace = False):
        """
        Add reference genome.
        
        Arguments:
            alias (str): required, name of fasta-annotation combo
            fasta (str): required, path to fasta file
            annotation (str): required path to GFF3 file
            genetic_code (str or int): NCBI genetic code (default=1)
            attr_mod (dict): mapping for non-standard attribute field names in GFF3 file (default={})
            replace (bool): replace any existing reference with same alias (default=False)
        """
        if alias in self.reference and not replace:
            raise Exception(f"ID '{alias}' is already in use.")
        self.reference[alias] = AnnotatedFasta(fasta, annotation, attr_mod = attr_mod,
                                               genetic_code = genetic_code, name = alias, memsave = memsave)
        return
    def remove_reference(self, alias):
        """
        Remove reference genome.
        
        Arguments:
            alias (str): required, alias of reference genome
        """
        if alias in self.reference:
            del self.reference[alias]
        return
    def clear_reference(self):
        """
        Remove all reference genomes.
        """
        self.reference = {}
        return
    def _subset_annotation(self, *ids, quiet = True, preserve_order = True):
        if ids:
            ref_to_subset = []
            ## if reference has not been subset before, queue it
            if set(self.reference.keys()) != self.subset_ref:
                ref_to_subset = [ref for ref in self.reference if ref not in self.subset_ref]
            ## if genes are different from previous subset_annotation execution, queue ALL references
            if set(ids) != self.subset_gid:
                ref_to_subset = list(self.reference.keys())
            ## start subsetting
            self.mkdir("subset_ann", tmp = True)
            mk_fname = lambda x: self.reserve_fname("subset_ann", f"subset.{x}.gff",
                                                    tmp = True, newfile = True)
            for alias in ref_to_subset:
                ref = self.reference[alias]
                ref.subset_annotation(*ids, fout = mk_fname(alias), preserve_order = preserve_order)
            self.subset_ref = set(self.reference.keys())
            self.subset_gid = set(ids)
        elif not quiet:
            typer.echo("MINORg._subset_annotation requires ids.")
        return
    def subset_annotation(self, quiet = True, preserve_order = True):
        """
        Subset annotations of all reference genomes according to self.genes.
        
        Reduces annotation lookup time.
        
        Arguments:
            quiet (bool): silence printing of non-essential messages
            sort (bool): sort subset data
        """
        if self.genes:
            self._subset_annotation(*self.genes, quiet = quiet, preserve_order = preserve_order)
        elif not quiet:
            typer.echo("MINORg.subset_annotation requires MINORg.genes.")
        return
    def extend_reference(self, ext_gene, ext_cds):
        """
        Extend reference genomes from FASTA files of gene and CDS sequences.
        
        Sequence IDs of gene sequences in ext_gene file(s) should not be repeated and should not contain '|'. Sequences in ext_cds file(s) should be named <sequence ID in ext_gene file(s)>.<int> to indicate their parent gene. Corresponding gene and CDS sequences will be aligned with mafft for inference of CDS range in genes.
        
        Arguments:
            ext_gene (str or list): required, path to file or list of paths to files of gene sequences
            ext_cds (str or list): required, path to file or list of paths to files of CDS sequences
        """
        # if not (self.ext_gene and self.ext_cds):
        #     return
        # ext_gene = self.ext_gene
        # ext_cds = self.ext_cds
        if not non_string_iter(ext_gene): ext_gene = [ext_gene]
        if not non_string_iter(ext_cds): ext_cds = [ext_cds]
        tmp_dir = self.mkdir("tmp_extend_aln")
        ext_fasta = self.mkfname("ext.fasta", prefix = False)
        ext_gff = self.mkfname("ext.gff", prefix = False)
        extend_reference(ext_gene, ext_cds, ext_fasta, ext_gff, mafft = self.mafft,
                         # feature_type = "gene", subfeature_type = "CDS",
                         thread = self.thread, directory = tmp_dir, tmp = True, logger = self.logfile)
        self.add_reference("Extended", ext_fasta, ext_gff)
        return
    def _get_reference_seq(self, ids, *features, adj_dir=False, by_gene=False, mktmp=None, ref=None,
                           seqid_template="Reference|$source|$gene|$isoform|$feature|$n",
                           translate=False, isoform_lvl=None, gff_domain=None, fout=None,
                           # apply_template_to_dict=False,
                           **fouts):
        '''
        Separately processes each feature (*features)
        Returns (<dict of seqs of feature1>, <dict of seqs of feature2>...)
        Path to output files for each feature can optionally be provided (**fouts)
            (e.g. _self.get_reference_seq("CDS", CDS = <path to file>) --> returns dict AND writes to file)
        Use fout to get sequence of features with the the specified ids regardless of feature type
        '''
        output = {}
        mktmp = (lambda x: self.mkfname(x)) if mktmp is None else mktmp
        ref_aliases = (list(self.reference.keys()) if ref is None
                       else ref if non_string_iter(ref) else [ref])
        if fout is not None:
            fouts[None] = fout
            features = features + (None,)
        import string
        for feature in features:
            fout = fouts.get(feature, None)
            if fout is None: continue
            feature_seqs = {}
            for ref_alias in ref_aliases:
                ref = self.reference[ref_alias]
                ## process each gene separately
                for feature_id in ids:
                    ## fill in gene and source
                    ft_seqid_template = string.Template(seqid_template).safe_substitute(gene = feature_id,
                                                                                        source = ref_alias)
                    ## get features at isoform_lvl if specified
                    if isoform_lvl:
                        feature_ids = [x.get_attr("ID", fmt=list)[0] for x in
                                       ref.annotation.subset(feature_id, features = isoform_lvl,
                                                             subfeatures = True, sort = False)]
                    else:
                        feature_ids = [feature_id]
                    seqs = ref.get_feature_seq(*feature_ids, fout = None, feature_type = feature,
                                               fout_gff = gff_domain, adj_dir = adj_dir,
                                               translate = translate, by_gene = by_gene, mktmp = mktmp,
                                               db = self.db, rpsblast = self.rpsblast, pssm_id = self.pssm_ids,
                                               seqid_template = ft_seqid_template,
                                               apply_template_to_dict = True)
                    feature_seqs = {**feature_seqs,
                                    **dict(itertools.chain(*map(lambda x: x.items(), seqs.values())))}
            if fout:
                dict_to_fasta(feature_seqs, fout)
            output = {**output, **feature_seqs}
        return output
    def get_reference_seq(self, *features, adj_dir=False, by_gene=False, ref=None,
                          seqid_template="Reference|$source|$gene|$isoform|$feature|$n",
                          translate=False, isoform_lvl=None, fout=None,
                          # apply_template_to_dict=False,
                          **fouts):
        """
        Get sequence(s) of reference genes or isoforms.
        
        If self.pssm_ids is provided, sequences are restricted to the relevant domain(s).
        
        Arguments:
            *features (str): optional, GFF3 feature type(s) to retrieve
                If not specified, sequence of the gene or isoform will be retrieved directly.
            ref (str): optional, alias of reference genome from which to extract sequences.
                If not specified, all reference genomes will be searched for genes.
            seqid_template (str): optional, template for output sequence name.
                Template will be parsed by strings.Template.
                The default template is "Reference|$source|$gene|$isoform|$feature|$n".
                
                    - $source: reference genome alias
                    - $gene: gene/isoform ID
                    - $isoform: isoform ID if by_gene = False, else same as $gene
                    - $feature: GFF3 feature type
                    - $n: if multiple domains are present, they will be numbered according to 
                      proximity to 5' of sense strand
            
            adj_dir (bool): output sense strand
            by_gene (bool): merge sequences from all isoforms of a gene
            translate (bool): translate sequence. Should be used with adj_dir.
            fout (str): optional, path to output file if features is not specified.
            **fouts (str): optional, path to output files for each feature if features is specified.
                Example: self.get_reference_seq("CDS", CDS = <path to file>) returns dict AND writes to file
        
        Returns
        -------
        dict
            sequence(s) of reference genes or isoforms (specified in self.genes), grouped by feature.
            Format: {<feature>: {<seqid>: <seq>}}
        """
        return self._get_reference_seq(self.genes, *features, adj_dir = adj_dir,
                                       by_gene = by_gene, ref = ref,
                                       seqid_template = seqid_template, fout = fout,
                                       translate = translate, isoform_lvl = isoform_lvl,
                                       # apply_template_to_dict=False,
                                       **fouts)
    
    #########################
    ##  SUBCMD: HOMOLOGUE  ##
    #########################
    
    def _homologue(self, ids, fout, ref_dir = "ref", quiet = True, domain_name = None,
                   minlen = None, minid = None, mincdslen = None, merge_within = None,
                   check_recip = None, relax_recip = None, check_id_before_merge = None):
        """
        Discovery homologues based on self's attributes
        """
        if ids is None:
            raise Exception("Requires ids (list of gene IDs)")
        if not all(map(lambda gid: gid in self.subset_gid, ids)):
            self._subset_annotation(*ids)
        ## instantiate variables
        printi = ((lambda msg: None) if quiet else (print))
        seqid_template = "Reference|${source}|${domain}|${feature}|${complete}|${gene}|${isoform}"
        domain_name = get_val_default(domain_name, self.domain_name)
        make_refname = lambda suf, **kwargs: self.results_fname(ref_dir, f"ref_{domain_name}_{suf}",
                                                                reserve = True, **kwargs)
        ## generate fname
        fasta_gene = make_refname(f"gene.fasta")
        fasta_cds = make_refname(f"CDS.fasta")
        gff_domain = (None if not self.pssm_ids else make_refname(f"domain.gff"))
        ## get reference sequences
        printi("Extracting reference sequences")
        self.get_reference_seq("gene", "CDS", adj_dir = True, by_gene = True, mktmp = make_refname,
                               seqid_template = seqid_template, gene = fasta_gene, CDS = fasta_cds,
                               gff_domain = gff_domain)
        ## store reference data location
        self.ref_gene = fasta_gene
        self.ref_cds = fasta_cds
        self.gff_domain = gff_domain
        ## get homologues
        printi("Finding homologues")
        if self.query:
            find_homologue_multiindv(fasta_queries = self.query_map, fout = fout,
                                     ## execution options
                                     # directory = self.results_fname(ref_dir, prefix = False),
                                     directory = self.results_fname(prefix = False),
                                     threads = self.thread, keep_tmp = self.keep_tmp,
                                     ## blastn options
                                     blastn = self.blastn,
                                     fasta_complete = self.ref_gene, fasta_cds = self.ref_cds,
                                     ## reciprocal blast options
                                     genes = ids, 
                                     check_reciprocal = get_val_default(check_recip, self.check_recip),
                                     relax = get_val_default(relax_recip, self.relax_recip),
                                     ## reference options (for reciprocal blast)
                                     # reference = self.reference,
                                     gff = self.annotations, fasta_ref = self.assemblies,
                                     attribute_mod = self.attr_mods,
                                     ## filtering options
                                     min_len = get_val_default(minlen, self.minlen),
                                     min_id = get_val_default(minid, self.minid),
                                     min_cds_len = get_val_default(mincdslen, self.mincdslen),
                                     check_id_before_merge = get_val_default(check_id_before_merge,
                                                                             self.check_id_before_merge),
                                     merge_within_range = get_val_default(merge_within, self.merge_within))
        return
        
    def seq(self, minlen = None, minid = None, mincdslen = None, quiet = True,
            check_recip = None, relax_recip = None, check_id_before_merge = None):
        """
        Identify targets and generate self.target file based on self's attributes.
        
        All arguments are optional. If not provided, the corresponding value stored in self's attributes will be used.
        
        If self.pssm_ids has been provided, either 1) self.rps_hits
        OR 2) self.db AND self.rpsblast are required.
        
        >>> my_minorg = MINORg(directory = "/my/output/directory", 
                               prefix = "test", tmp = False, keep_tmp = True, thread = 1)
        >>> my_minorg.add_reference("TAIR10", "/path/to/TAIR10_Chr.all.fasta", "/path/to/TAIR10_GFF3.genes.gff", replace = True)
        >>> my_minorg.genes = ["AT5G66900"]
        >>> my_minorg.subset_annotation()
        >>> my_minorg.pssm_ids = "366375" # NB-ARC domain
        >>> my_minorg.query_reference = True
        >>> my_minorg.seq()
        Exception: If self.pssm_ids has been provided, self.db (path to RPS-BLAST database) AND self.rpsblast (path to rpsblast or rpsblast+ executable OR command name if available at command line) are required
        >>> my_minorg.db = "/path/to/rpsblast/database"
        >>> my_minorg.rpslast = "rpsblast+"
        >>> my_minorg.seq()
        >>> my_minorg.target
        '/my/output/directory/test/test_3663775_targets.fasta'
        
        Arguments:
            minlen (int): optional. See attributes.
            minid (float): optional. See attributes.
            mincdslen (in): optional. See attributes.
            check_recip (bool): optional. See attributes.
            relax_recip (bool): optional. See attributes.
            check_id_before_merge (bool): optional. See attributes.
        """
        if not self.genes:
            raise Exception("Requires self.genes (list of gene IDs)")
        if self.pssm_ids and not (self.db and self.rpsblast):
            raise Exception(("If self.pssm_ids has been provided, self.db (path to RPS-BLAST database)"
                             " AND self.rpsblast (path to rpsblast or rpsblast+ executable"
                             " OR command name if available at command line)"
                             " are required"))
        ## get kwargs
        kwargs = locals()
        del kwargs["self"]
        ## create output file
        fout = self.results_fname(f"{self.domain_name}_targets.fasta", newfile = True)
        ## get homologues
        self._homologue(self.genes, fout, ref_dir = "ref", **kwargs)
        ## append reference genes to targets file
        if self.query_reference:
            tmp = self.results_fname(f"tmp.fasta", tmp = True)
            cat_files([fout, self.ref_gene], tmp, remove = False)
            os.rename(tmp, fout)
        self.target = fout
        return
    
    # def homologue(self, **kwargs):
    #     return self.seq(**kwargs)
    
    # def homolog(self, **kwargs):
    #     return self.seq(**kwargs)
    
    ####################
    ##  SUBCMD: gRNA  ##
    ####################
    
    def _generate_grna(self, fname, pam = None, length = None):
        if pam is None:
            pam = self.pam
        if length is None:
            length = self.length
        return find_multi_gRNA(fname, pam = pam, gRNA_len = length)
    
    def grna(self):
        """
        Generate all possible gRNA from self.target file based on self's attributes: pam, length.
        """
        self.grna_hits = self._generate_grna(self.target, pam = self.pam, length = self.length)
        self.write_all_grna_fasta()
        # self.grna_hits.write_fasta(self.grna_fasta, write_all = True)
        return
    
    ######################
    ##  SUBCMD: FILTER  ##
    ######################
    
    def _mask_ontargets(self, *mask_fnames):
        tmp_f = self.reserve_fname("tmp.tsv", newfile = True, tmp = True)
        self.masked = {}
        ## mask function
        def _mask_identical(alias_subject):
            fnames = [IndexedFasta(subject).filename for subject in assign_alias(alias_subject).values()]
            return {subject: tuple(itertools.chain(*[mask_identical(str(mask_fname), subject, tmp_f,
                                                                    blastn = self.blastn)
                                                     for mask_fname in mask_fnames if mask_fname]))
                    for subject in fnames}
        ## mask in reference
        ## (we don't use ranges of self.genes directly because some genes may be identical
        ##  but NOT in self.genes and we want to mask those too)
        self.masked = {**self.masked, **_mask_identical(self.assemblies)}
        ## mask in query
        self.masked = {**self.masked, **_mask_identical(self.query)}
        ## mask in backgrounds
        self.masked = {**self.masked, **_mask_identical(self.background)}
        return
    
    def write_mask_report(self, fout):
        """
        Write mask report to file.
        
        Details masked regions as well as which sequence each region is identical to.
        
        Arguments:
            fout (str): required, path to output file
        """
        with open(fout, 'w') as f:
            inv_fnames = {fname: alias for alias, fname in
                          {**{alias: ref.fasta for alias, ref in self.reference.items()},
                           **self.query,
                           **assign_alias(self.background)}.items()}
            f.write("#alias\tfname\n")
            for fname, alias in inv_fnames.items():
                f.write(f"#{alias}\t{fname}\n")
            f.write('\t'.join(["source", "molecule", "start", "end", "masked"]) + '\n')
            for fname, masked in self.masked.items():
                alias = inv_fnames[fname]
                for entry in masked:
                    f.write('\t'.join([str(alias),
                                       entry.molecule, str(entry.start), str(entry.end), entry.masked]) + '\n')
        return
    
    def _is_offtarget_pos(self, hsp, ifasta):
        if ifasta.filename in self.masked:
            for masked_target in self.masked[ifasta.filename]:
                same_molecule = hsp.hit_id == masked_target.molecule
                if same_molecule:
                    if masked_target.within((hsp.hit_start, hsp.hit_end)): return False
        return True
    
    def is_offtarget_pos(self, hsp, ifasta):
        """
        Check whether gRNA hit is outside masked target regions.
        
        Arguments:
            hsp (Biopython HSP): required
            ifasta (IndexedFasta): required, subject of BLAST search
        
        Returns
        -------
        bool
            whether a HSP hit is outside of a masked region
        """
        return self._is_offtarget_pos(hsp, ifasta)
    
    def _is_offtarget_aln(self, hsp, query_result,
                          mismatch_check = True, gap_check = True, fully_aligned_check = True):
        qlen = query_result.seq_len
        mismatch_excess = (max(1, self.ot_mismatch) - 1) - (qlen - hsp.ident_num)
        gap_excess = (max(1, self.ot_gap) - 1) - hsp.gap_num
        unaligned_excess = qlen - hsp.query_end + hsp.query_start
        mismatch_pass = (mismatch_excess >= 0) if mismatch_check else True
        gap_pass = (gap_excess >= 0) if gap_check else True
        fully_aligned = unaligned_excess == 0 if fully_aligned_check else True
        return fully_aligned and mismatch_pass and gap_pass
    
    def is_offtarget_aln(self, hsp, query_result, **kwargs):
        """
        Check whether (potentially off-target) gRNA hit aligns too well.
        
        Used in combination with is_offtarget_pos to determine if an off-target gRNA hit could be problematic.
        
        Arguments:
            hsp (Biopython HSP): required
            query_result (Biopython QueryResult): required
            **kwargs: other arguments supplied passed to :meth:`MINORg._is_offtarget_aln`
        
        Returns
        -------
        bool
            whether a HSP hit meets the threshold for problematic off-target gRNA hit
        """
        return self._is_offtarget_aln(hsp, query_result, **kwargs)
    
    def _has_pam(self, hsp, query_result, ifasta):
        pam_pre = self.pam.five_prime()
        pam_post = self.pam.three_prime()
        pam_pre_max = self.pam.five_prime_maxlen()
        pam_post_max = self.pam.three_prime_maxlen()
        seq = ifasta[hsp.hit_id]
        qlen = query_result.seq_len
        ## get maximum PAM distance (unused gap allowance + uanligned gRNA length)
        gap_excess = (max(1, self.ot_gap) - 1) - hsp.gap_num
        len_excess = qlen - (hsp.query_end - hsp.query_start)
        pre_excess = hsp.query_start
        post_excess = qlen - hsp.query_end
        pre_len_max = pam_pre_max + gap_excess + pre_excess
        post_len_max = pam_post_max + gap_excess + post_excess
        ## get the sequences flanking the gRNA hit
        if hsp.hit_strand == 1:
            pre_seq = seq[hsp.hit_start - pre_len_max : hsp.hit_start]
            post_seq = seq[hsp.hit_end : hsp.hit_end + post_len_max]
        else:
            post_seq = seq[hsp.hit_start - post_len_max: hsp.hit_start].reverse_complement()
            pre_seq = seq[hsp.hit_end : hsp.hit_end + pre_len_max].reverse_complement()
        ## assess whether pre-gRNA and post-gRNA regions match PAM pattern
        has_pam_pre = False if not pam_pre else bool(re.search(pam_pre, str(pre_seq)))
        has_pam_post = False if not pam_post else bool(re.search(pam_post, str(post_seq)))
        return has_pam_pre or has_pam_post
    
    def is_offtarget_pam(self, hsp, query_result, ifasta):
        """
        Check for presence of PAM in close proximity and correct orientation to gRNA hit.
        
        Arguments:
            hsp (Biopython HSP): required
            query_result (Biopython QueryResult): required
            ifasta (IndexedFasta): required
        
        Returns
        -------
        bool
            whether a HSP hit is in close proximity and correct orientation to a PAM site
        """
        if self.ot_pamless: return True
        return self._has_pam(hsp, query_result, ifasta)
    
    ## function to judge if off-target crosses threshold
    def is_offtarget_hit(self, hsp, query_result, ifasta):
        """
        Check whether a gRNA hit aligns too well outside of masked target regions.
        
        Determined by is_offtarget_pos AND is_offtarget_aln AND is_offtarget_pam for a given HSP, QueryResult, and FASTA file combination.
        
        Arguments:
            hsp (Biopython HSP): required
            query_result (Biopython QueryResult): required
            ifasta (IndexedFasta): required
        
        Returns
        -------
        bool
            whether a gRNA hit is off-target and problematic
        """
        ifasta = IndexedFasta(ifasta)
        return self.is_offtarget_pos(hsp, ifasta) and self.is_offtarget_aln(hsp, query_result) and \
            self.is_offtarget_pam(hsp, query_result, ifasta)
    
    def _offtarget_hits(self, fasta_subject, keep_output = False, fout = None):
        if not fout:
            fout = self.mkfname("tmp.xml")
            keep_output = False
        if self.grna_fasta is None:
            self.write_all_grna_fasta()
        from Bio.Blast.Applications import NcbiblastnCommandline
        blast(NcbiblastnCommandline, header = '', fout = fout, outfmt = 5,
              query = self.grna_fasta, subject = IndexedFasta(fasta_subject).filename,
              cmd = self.blastn, task = "blastn-short",
              # mt_mode = (0 if self.thread == 1 else 1),
              num_threads = self.thread)
        offtarget = set(hsp.query_id
                        for query_result in SearchIO.parse(fout, "blast-xml")
                        for hit in query_result for hsp in hit
                        if self.is_offtarget_hit(hsp, query_result, fasta_subject))
        if not keep_output: os.remove(fout)
        return offtarget
    
    def filter_background(self, keep_blast_output = False, mask_reference = True,
                          *other_mask_fnames):
        """
        Set background filter check for candidate gRNAs.
        
        Masks target sequences in all FASTA files to be screened for off-target, BLASTs all candidate gRNAs to those FASTA files, and assesses each BLAST hit individually for whether they could potentially be problematic.
        
        Relevant attributes: screen_reference, ot_pamless, ot_mismatch, ot_gap

        Relevant methods: add_background, remove_background
        
        Arguments:
            keep_blast_output (bool): retain BLAST output file. Default behaviour deletes it.
            mask_reference (bool): mask reference genes (default=True)
            *other_mask_fnames (str): optional, paths to other FASTA files not in self.background that are also to be screened for off-target
        """
        if not self.grna_hits:
            raise MessageError( ("MINORg.filter_background requires self.grna_hits"
                                 " (generated with self.generate_grna)") )
        if self.grna_fasta is None:
            self.grna_fasta = self.results_fname("gRNA_all.fasta") ## write gRNA to file so we can BLAST it
            self.grna_hits.write_fasta(self.grna_fasta, write_all = True)
        if mask_reference and not self.ref_gene and self.genes:
            ref_to_mask = self.mkfname("tmp_ref_to_mask.fasta", tmp = True)
            self.get_reference_seq(fout = ref_to_mask)
            self.ref_gene = ref_to_mask
        self._mask_ontargets(self.target,
                             *([self.ref_gene] if mask_reference and self.genes else []),
                             *self.mask,
                             *other_mask_fnames)
        excl_seqid = set()
        def get_excl_seqids(alias_ifasta, descr = None):
            output = set()
            for alias, ifasta in alias_ifasta.items():
                ot_hits = self.mkfname(f"hit_{descr}_{alias}.xml" if descr else f"hit_{alias}.xml")
                output |= self._offtarget_hits(ifasta, fout = ot_hits,
                                               keep_output = keep_blast_output)
            return output
        if self.screen_reference and self.reference:
            excl_seqid |= get_excl_seqids(self.reference, descr = "ref")
        if self.background:
            excl_seqid |= get_excl_seqids(self.background, descr = "bg")
        if self.query:
            excl_seqid |= get_excl_seqids(self.query, descr = "query")
        ## parse bg check status
        grna_seqs = fasta_to_dict(self.grna_fasta)
        grna_screened = tuple(grna_seqs.values())
        grna_failed = set(str(grna_seqs[seqid]).upper() for seqid in excl_seqid)
        grna_passed = set(str(seq).upper() for seq in grna_screened if str(seq).upper() not in grna_failed)
        ## record bg check status
        self.grna_hits.set_seqs_check("background", True, grna_passed)
        self.grna_hits.set_seqs_check("background", False, grna_failed)
        return
    
    def filter_gc(self, gc_min = None, gc_max = None):
        """
        Set GC check for candidate gRNA.
        
        Arguments:
            gc_min (float): optional. See attributes.
            gc_max (float): optional. See attributes.
        """
        gc_min = get_val_default(gc_min, self.gc_min)
        gc_max = get_val_default(gc_max, self.gc_max)
        for grna_seq in self.grna_hits.gRNAseqs.values():
            grna_seq.set_gc_check(gc_min = gc_min, gc_max = gc_max)
        return
    
    def _adjust_feature_range(self, seq, gene, ref_alias, feature, domain = False, **kwargs):
        gene_ann = self.reference[ref_alias].annotation.get_id(gene, output_list = False)
        feature_ranges = ranges_union([[(x.start, x.end)] for x in
                                        self.reference[ref_alias].annotation.get_subfeatures_full(gene,
                                                                                                  feature_types = feature)])
        if domain and self.gff_domain:
            domain_anns = GFF(fname = self.gff_domain, quiet = True, fmt = "GFF")
            domain_ann = [entry for entry in domain_anns.get_id(gene, output_list = True)
                          if entry.source == ref_alias]
            seq_range = ranges_intersect([(gene_ann.start, gene_ann.end)],
                                         (ranges_union([[(x.start, x.end)] for x in domain_ann])))[0]
            feature_ranges = ranges_intersect([seq_range], feature_ranges)
        else:
            seq_range = (gene_ann.start, gene_ann.end)
        return adjusted_feature_ranges(seq, seq_range, feature_ranges, strand = gene_ann.strand, **kwargs)
    
    def align_reference(self, fout):
        """
        Align reference genes.
        
        CDS are first aligned. Full genes are then added to the alignment.
        
        Arguments:
            fout (str): required, path to output FASTA file
        """
        if not self.ref_cds and self.ref_gene:
            raise MessageError( ("MINORg.align_reference_and_targets requires"
                                 " MINORg.ref_cds and MINORg.ref_gene.") )
        tmp_f = self.results_fname("tmp_align_reference.fasta", tmp = True)
        ## align CDS
        with open(fout, 'w') as f:
            stdout, stderr = MafftCommandline(self.mafft, input = self.ref_cds, # quiet = True,
                                              thread = self.thread)()
            f.write(stdout)
        ## align gene
        with open(tmp_f, 'w') as f:
            stdout, stderr = MafftCommandline(self.mafft, add = self.ref_gene, # quiet = True,
                                              thread = self.thread, input = fout)()
            f.write(stdout)
        ## adjust CDS' alignments relative to reference based on reference gene annotations
        aln = fasta_to_dict(tmp_f)
        seqids_cds = list(fasta_to_dict(self.ref_cds).keys())
        seqids_gene = list(fasta_to_dict(self.ref_gene).keys())
        def same_source_and_gene(seqid1, seqid2):
            return ((self.get_ref_seqid(seqid1, "source") == self.get_ref_seqid(seqid2, "source"))
                    and (self.get_ref_seqid(seqid1, "gene") == self.get_ref_seqid(seqid2, "gene")))
        for seqid_cds in seqids_cds:
            seqid_gene = [seqid for seqid in seqids_gene if same_source_and_gene(seqid_cds, seqid)]
            if not seqid_gene: continue
            seqid_gene = seqid_gene[0]
            seq_gene = aln[seqid_gene]
            source = self.get_ref_seqid(seqid_cds, "source")
            gene = self.get_ref_seqid(seqid_cds, "gene")
            adj_ranges = self._adjust_feature_range(seq_gene, gene, source, "CDS",
                                                    domain = bool(self.gff_domain), subtract_gaps = False)
            new_cds_seq = Seq.Seq('-'*(len(seq_gene)))
            for r in adj_ranges:
                start, end = r
                new_cds_seq = new_cds_seq[:start] + seq_gene[start:end] + new_cds_seq[end:]
            aln[seqid_cds] = new_cds_seq
        dict_to_fasta(aln, fout)
        os.remove(tmp_f)
        return
    
    ## untested
    def align_reference_and_targets(self, domain_name = None):
        """
        Align reference genes and targets. Path to FASTA file generated will be stored in self.alignment.
        
        Arguments:
            domain_name (str): optional, domain name to be used when naming output file
        """
        if not self.ref_cds and self.ref_gene and self.target:
            raise MessageError( ("MINORg.align_reference_and_targets requires"
                                 " MINORg.ref_cds, MINORg.ref_gene, and MINORg.targets.") )
        fasta_aln = self.results_fname(f"{get_val_default(domain_name, self.domain_name)}_mafft.fasta")
        tmp_f = self.results_fname("tmp_align_reference_and_targets.fasta", tmp = True)
        ## align reference
        self.align_reference(tmp_f)
        ## add non-reference seqs to alignment
        ref_seqids = set(fasta_to_dict(self.ref_gene).keys())
        nonref_seqs = {seqid: seq for seqid, seq in fasta_to_dict(self.target).items()
                       if seqid not in ref_seqids}
        if nonref_seqs:
            tmp_targets_nr = self.results_fname("tmp2.fasta", tmp = True)
            dict_to_fasta(nonref_seqs, tmp_targets_nr)
            with open(fasta_aln, "w+") as f:
                stdout, stderr = MafftCommandline(self.mafft, add = tmp_targets_nr, # quiet = True,
                                                  thread = self.thread, input = tmp_f,
                                                  adjustdirectionaccurately = True)()
                f.write(stdout)
            self.rm_tmpfiles(tmp_targets_nr, tmp_f)
        else:
            os.rename(tmp_f, fasta_aln)
        self.alignment = fasta_aln
        return
    
    def filter_feature(self, max_insertion = None, min_within_n = None, min_within_fraction = None,
                       alignment_rvs_pattern = "^_R_$seqid$$"):
        """
        Set feature check for candidate gRNAs.
        
        Range(s) for desired feature of inferred homologue targets discoverd by MINORg will be inferred based on alignment with reference genes based on self's attributes: max_insertion
        
        Arguments:
            max_insertion (int): optional. See attributes.
            min_within_n (int): optional. See attributes.
            min_within_fraction (float): optional. See attributes.
        """
        max_insertion = get_val_default(max_insertion, self.max_insertion)
        min_within_n = get_val_default(min_within_n, self.min_within_n)
        min_within_fraction = get_val_default(min_within_fraction, self.min_within_fraction)
        if not self.alignment:
            self.align_reference_and_targets()
        alignment = fasta_to_dict(self.alignment)
        genes = {'|'.join(seqid.split('|')[5:-1]):
                 seqid for seqid in alignment
                 if (seqid.split('|')[0] == "Reference" and seqid.split('|')[3] == "gene")}
        gene_anns = {alias: {gene: ref.annotation.get_id(gene, output_list = False) for gene in self.genes}
                     for alias, ref in self.reference.items()}
        if self.pssm_ids and self.gff_domain:
            ## All fields in domain GFF file will be identical to that of the gene with exception
            ##   of feature ('domain'), ranges (whatever the start and end of the domain are),
            ##   and source ('minorg').
            domain_anns = GFF(fname = self.gff_domain, quiet = True, fmt = "GFF")
        if self.feature is None:
            feature_ranges = {source: {gene: [(feature.start, feature.end)]
                                       for gene, feature in gene_ann.items()}
                              for source, gene_ann in gene_anns.items()}
        else:
            feature_ranges = {source: {gene:
                                       ranges_union([[(x.start, x.end) for x in
                                                      ref.annotation.get_subfeatures_full(gene, feature_types = self.features)]])
                                       for gene in genes}
                              for source, ref in self.reference.items()}
        def adjust_feature_ranges(gene, seqid, **kwargs):
            source = self.get_ref_seqid(seqid, attr = "source")
            gene_ann = gene_anns[source]
            gene_feature_ranges = feature_ranges[source]
            if self.gff_domain:
                domain_ann = [entry for entry in domain_anns.get_id(gene, output_list = True)
                              if entry.source == source]
                return adjusted_feature_ranges(alignment[seqid],
                                               ranges_intersect([(gene_ann[gene].start, gene_ann[gene].end)],
                                                                (ranges_union([[(x.start, x.end)] for x in
                                                                               domain_ann])))[0],
                                               gene_feature_ranges[gene],
                                               strand = gene_ann[gene].strand,
                                               **kwargs)
            else:
                return adjusted_feature_ranges(alignment[seqid],
                                               (gene_ann[gene].start, gene_ann[gene].end),
                                               gene_feature_ranges[gene],
                                               strand = gene_ann[gene].strand,
                                               **kwargs)
        from string import Template
        def is_rvs(aln_seqid, seqid):
            return re.search(Template(alignment_rvs_pattern).substitute(seqid = re.escape(seqid)), aln_seqid)
        feature_only_ranges = {seqid: adjust_feature_ranges(gene, seqid, subtract_gaps = True)
                               for gene, seqid in genes.items()}
        feature_gaps_ranges = {seqid: ranges_subtract(adjust_feature_ranges(gene, seqid, subtract_gaps = False),
                                                      feature_only_ranges[seqid])
                               for gene, seqid in genes.items()}
        ## define acceptable ranges in targets
        get_target_feature_ranges = make_target_feature_ranges_function(feature_only_ranges,
                                                                        feature_gaps_ranges,
                                                                        max_insertion = max_insertion)
        targets_feature_ranges = {seqid: get_target_feature_ranges(alignment[seqid], seqid)
                                  for seqid in alignment}
        ## iterate through all gRNAs
        for gRNA_seq, coverage in self.grna_hits.hits.items():
            ## assess each hit
            for gRNA_hit in coverage:
                seq_id_seq = [(seq_id, seq) for seq_id, seq in alignment.items()
                              if (seq_id == gRNA_hit.target_id or is_rvs(seq_id, gRNA_hit.target_id))]
                ## if target not among user-specified targets or not in alignment, skip to next gRNA hit
                if len(seq_id_seq) == 0: continue
                seq_id, seq = [(seq_id, seq) for seq_id, seq in alignment.items()
                               if (seq_id == gRNA_hit.target_id or is_rvs(seq_id, gRNA_hit.target_id))][0]
                seq_plus = seq_id == gRNA_hit.target_id
                gRNA_hit.set_parent_sense('+' if seq_plus else '-') ## used to tie break set cover
                gRNA_range = gRNA_hit.range if (not is_rvs(seq_id, gRNA_hit.target_id)) \
                             else gRNA_hit.reverse_range
                ## append gRNAHit object if within at least 1 feature
                if (within_feature(targets_feature_ranges[seq_id], seq, gRNA_range,
                                   min_within_n = min_within_n, min_within_fraction = min_within_fraction)):
                    # new_coverage.append(gRNA_hit)
                    gRNA_hit.set_feature_check(True)
                else:
                    gRNA_hit.set_feature_check(False)
        return
    
    ##########################
    ##  SUBCMD: MINIMUMSET  ##
    ##########################
    
    def _check_valid_status(self, check_name, type, descr):
        """
        Checks if check contains only NAs, if some entries are NAs and some are set, or if all are set.
        
        Arguments:
            check_name (str): required, check name
            type (str): required, type of feature for which check is relevant. Should be 'hit' or 'seq'.
            descr (str): required, description of check
        """
        if type == "hit":
            valid_check = self.grna_hits.valid_hit_check(check_name)
        else:
            valid_check = self.grna_hits.valid_seq_check(check_name)
        if not valid_check:
            self.logfile.warning(f"Ignoring unset check: {descr}")
            # warnings.warn(f"Ignoring unset check: {descr}")
            return
        if type == "hit":
            check_values = set(x.check(check_name) for x in self.grna_hits.flatten_hits())
        else:
            check_values = set(x.check(check_name) for x in self.grna_hits.flatten_gRNAseqs())
        if None in check_values:
            # warings.warn(f"The {descr} status of at least one gRNA {type} is not known.")
            self.logfile.warning(f"The {descr} status of at least one gRNA {type} is not known.")
            if self.accept_invalid:
                print( (f"The programme will assume that these {type}s are VALID"
                        " and treat then as viable candidates.\n") )
            else:
                print( (f"The programme will assume that these {type}s are INVALID"
                        " and exclude them as candidates.\n") )
        return
    
    def minimumset(self, sets = None, manual = None, fasta = None, fout_fasta = None, fout_map = None,
                   report_full_path = True):
        """
        Generate minimum set(s) of gRNA required to cover all targets.
        
        Arguments:
            sets (int): optional. See attributes.
            manual (bool): manually approve all gRNA sets. Defaults to NOT self.auto if not provided.
            fasta (str): optional, path to FASTA file. Used for renaming gRNA.
            fout_fasta (str): optional, path to output FASTA file.
                Autogenerated using self.prefix and self.directory if not provided.
            fout_map (str): optional, path to output .map file.
                Autogenerated using self.prefix and self.directory if not provided.
            report_full_path (bool): print full path to output files upon successful writing
        """
        if sets is None: sets = self.sets
        ## assume all targets in self.grna_hits are to be, well, targeted
        targets = set(hit.target_id for hit in self.grna_hits.flatten_hits())
        ## check if statuses has been set. If not, warn user.
        self._check_valid_status("GC", "seq", "%GC")
        self._check_valid_status("background", "seq", "off-target")
        self._check_valid_status("feature", "hit", "within-feature")
        ## filter grna by checks
        grna_hits = copy.deepcopy(self.passed_grna)
        ## remove sequences not included in fasta file from grna_hits if fasta file provided
        if fasta:
            grna_hits.remove_seqs(*[seq for seq in grna_hits.seqs() if
                                    str(seq).upper() not in map(lambda x: str(x).upper(),
                                                                fasta_to_dict(fasta).values())])
            grna_hits.rename_seqs(fasta)
        ## warn user if desired number of sets cannot be returned
        if len(grna_hits) < sets:
            # warnings.warn(f"The gRNA sequences cannot cover all target sequences the desired number of times ({len(grna_hits)} valid gRNA, {sets} set(s) requested).\n")
            self.logfile.warning(f"The gRNA sequences cannot cover all target sequences the desired number of times ({len(grna_hits)} valid gRNA, {sets} set(s) requested).\n")
            return
        ## start generating sets
        grna_sets = []
        while len(grna_sets) < sets:
            ## get a (minimum) set of gRNA sequences
            seq_set = get_minimum_set(grna_hits, set_num = len(grna_sets) + 1,
                                      manual_check = (manual if manual is not None else not self.auto),
                                      targets = targets)
            ## if valid set returned
            if seq_set:
                grna_sets.append(seq_set) ## add to existing list of sets
                grna_hits.remove_seqs(seq_set) ## remove seqs in seq_set so they're not repeated
            else:
                # warnings.warn(f"The gRNA sequences cannot cover all target sequences the desired number of times ({sets}). (Failed at set {len(grna_sets) + 1} of {sets})\n")
                self.logfile.warning(f"The gRNA sequences cannot cover all target sequences the desired number of times ({sets}). (Failed at set {len(grna_sets) + 1} of {sets})\n")
                break
        ## write
        if not fout_fasta: fout_fasta = get_val_default(self.final_fasta,
                                                        self.results_fname("gRNA_final.fasta"))
        if not fout_map: fout_map = get_val_default(self.final_map,
                                                    self.results_fname("gRNA_final.map"))
        ## write gRNA fasta file and gRNA mapping
        if grna_sets:
            self.grna_hits.write_fasta(fout_fasta, seqs = itertools.chain(*grna_sets), fasta = fasta)
            print("Final gRNA sequence(s) have been written to"
                  f" {fout_fasta if report_full_path else os.path.basename(fout_fasta)}")
            self.grna_hits.write_mapping(fout_map, sets = grna_sets, fasta = fasta, version = 2)
            print("Final gRNA sequence ID(s), gRNA sequence(s), and target(s) have been written to"
                  f" {fout_map if report_full_path else os.path.basename(fout_map)}")
        ## print summary
        print(f"\n{sets} mutually exclusive gRNA set(s) requested. {len(grna_sets)} set(s) found.")
        return




# :param directory: path to output directory
# :type directory: str
# :param config: path to config.ini file
# :type config: str
# :param prefix: prefix for output files and directories
# :type prefix: str
# :param thread: maximum number of threads for parallel processing
# :type thread: int
# :param keep_tmp: retain temporary files
# :type keep_tmp: bool

# :ivar prefix: [general] (str) prefix for output directories and files
# :ivar thread: [general] (int) maximum number of threads for parallel processing

# >>> from minorg.MINORg import MINORg
# >>> my_minorg = MINORg()
# '{"oh no"}'


# """
# Create a MINORg object.

# Arguments:
# - directory - path to output directory (string)
# - config - path to config.ini file (string)
# - prefix - prefix for output directories and files (string)
# - thread - maximum number of threads for parallel processing (integer)
# - keep_tmp - retain temporary files (bool)

# Attributes:
# - prefix - prefix for output directories and files (string)
# - thread - maximum number of threads for parallel processing (integer)
# - rpsblast - path to rpsblast or rpsblast+ executable (string)
# - blastn - path to blastn executable (string)
# - mafft - path to mafft executable (string)
# - db - path to local RPS-BLAST database (string)
# - remote_rps - use remote database instead of local database for RPS-BLAST (bool)
# - pssm_ids - list of Pssm-Ids of domain(s) for homologue search. If multiple Pssm-Ids are provided, overlapping domain hits will be merged. (list of string)
# - domain_name - human-readable domain name used in sequence and file names in place of Pssm-Ids (string)
# - genes - list of target gene IDs (list of string)
# - query_reference - include reference genes as targets for gRNA discovery (bool)

# Attibutes (homologue discovery):
# - minlen - minimum homologue length (bp) (int)
# - minid - minimum hit % identity (float)
# - mincdslen - minimum number of bases in homologue aligned with reference gene CDS (int)
# - check_recip - execute reciprocal check (bool)
# - relax_recip - execute relaxed reciprocal check (bool)
# - check_id_before_merge - filter out hits by % identity before merging into potential homologues (bool)
# - merge_within - maximum distance (bp) between hits for merging (int)

# Attributes (gRNA):
# - length - gRNA length (bp) (int)
# - pam - PAM pattern (string)

# Attributes (background filter):
# - screen_reference - include reference genome for screening (bool)
# - ot_pamless - ignore absence of PAM when assessing off-target gRNA hits (bool)
# - offtarget - function that accepts Biopython's HSP and QueryResult objects and determines whether an off-target gRNA hit is problematic; if not set by user, ot_mismatch and ot_gap will be used (func) ## TODO!!! integrate this somehow?
# - ot_mismatch - minimum number of mismatches allowed for off-target gRNA hits (int)
# - ot_gap - minimum number of gaps allowed for off-target gRNA hits (int)

# Attributes (GC filter):
# - gc_min - minimum GC content (between 0 and 1, where 0 is no GC content and 1 is all GC) (float)
# - gc_max - maximum GC content (betweew 0 and 1, where 0 is no GC content and 1 is all GC) (float)

# Attributes (feature filter):
# - feature - GFF3 feature within which gRNA are to be designed (string)
# - max_insertion - maximum allowable insertion size in feature (bp) (int)
# - min_within_n - minimum number of reference genes which feature a gRNA must align within (bool)
# - min_within_fraction - minimum fraction of reference genes which feature a gRNA must align within (between 0 and 1, where 0 is none and 1 is all; if 0, min_within_n will be set to 1) (float)

# Attributes (minimum set):
# - sets - number of sets to generate (int)
# - auto - generate sets without requiring manual user confirmation for each set (bool)
# - accept_invalid - score 'NA' as 'pass' (bool)
# - accept_feature_unknown - score 'NA' as 'pass' for feature check (bool)
# - accept_invalid_field - score 'NA' as 'pass' if all entries for a check are 'NA' (bool)
# - pass_map - path to output .map file for gRNA that pass all valid checks (autogenerated by MINORg if not provided) (string)
# - pass_fasta - path to output .fasta file for gRNA that pass all valid checks (autogenerated by MINORg if not provided) (string)
# - final_map - path to output .map file for final gRNA set(s) (autogenerated by MINORg if not provided) (string)
# - final_fasta - path to output .fasta file for final gRNA set(s) (autogenerated by MINORg if not provided) (string)
# """
