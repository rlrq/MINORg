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
from Bio.Seq import Seq

from minorg import (
    _logging_level,
    _warning,
    MINORgWarning,
    MINORgError
)

from minorg.log import MINORgLogger

from minorg.mafftcommandline_add import MafftCommandline
from minorg.filter_grna import (
    mask_identical,
    make_target_feature_ranges_function,
    within_feature
)
from minorg.minimum_set import (
    all_best_nr,
    all_best_pos,
    make_get_minimum_set
)

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
from minorg.extract_homologue import (
    merge_hits_and_filter,
    recip_blast_multiref
)
from minorg.extend_reference import extend_reference
from minorg.reference import AnnotatedFasta
from minorg.functions import (
    split_none,
    is_empty_file,
    cat_files,
    non_string_iter,
    assign_alias,
    adjusted_feature_ranges,
    ranges_union,
    ranges_intersect,
    ranges_subtract,
    imap_progress
)
from minorg.blast import blast, BlastHSP
from minorg.searchio import searchio_parse
from minorg.generate_grna import find_multi_gRNA
from minorg.fasta import dict_to_fasta, fasta_to_dict
from minorg.grna import gRNAHits
from minorg.pam import PAM
from minorg.ot_regex import OffTargetExpression

from typing import Optional, Type, Dict

warnings.showwarning = _warning
warnings.filterwarnings("always", category = MINORgWarning)

LOGGING_LEVEL = _logging_level


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
        self.directory = None if directory is None else os.path.abspath(directory) ## final output directory
        self.tmp = (tmp if self.directory is not None else True)
        self.keep_tmp = keep_tmp
        self.new_dirs = []
        self.tmp_files = set()
        ## create temporary directory
        if self.tmp:
            ## contents to be moved to self.directory upon fulfilment of command using self.resolve
            self.tmpdir = tempfile.mkdtemp()
            ## self.active_directory is used outside PathHandler objects.
            ## PathHandler object will handle cleanup if self.active_directory points to self.tmpdir
            self.active_directory = self.tmpdir
            print(f"tmpdir: {self.tmpdir}")
        elif self.directory is not None:
            self.tmpdir = None
            self.active_directory = self.directory
            self.mkdir(self.active_directory)
            print(f"Files will be written directly in: {self.active_directory}")
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
    # ## update output directory
    # def mv(self, directory):
    #     """
    #     Update output directory.
        
    #     Contents of the current output directory will also be moved to the new directory if files are
    #     already being written to output directory (i.e. not currently using a tmp dir)
        
    #     Arguments:
    #         directory (str): required, path to new output directory
    #     """
    #     if directory is None: return
    #     directory = os.path.abspath(directory)
    #     if directory == self.directory: return
    #     ## if files are being written directly to output directory, move to new destination
    #     if self.active_directory == self.directory:
    #         mv_dir_overwrite(self.active_directory, directory)
    #         self.tmp_files = set(re.sub('^' + re.escape(self.active_directory), directory, fname)
    #                              for fname in self.tmp_files)
    #         self.new_dirs = [re.sub('^' + re.escape(self.active_directory), directory, dirname)
    #                          for dirname in self.new_dirs]
    #         self.active_directory = directory
    #     self.directory = directory
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
                directory = os.path.abspath(os.path.join(self.active_directory, directory))
            output.append(directory)
            if not os.path.exists(directory):
                def make_parent(d):
                    if not os.path.exists(os.path.dirname(d)):
                        make_parent(os.path.dirname(d))
                    else:
                        os.mkdir(d)
                make_parent(directory)
                self.new_dirs.append(directory)
                if tmp: self.tmp_files.add(directory)
        return output[0] if len(directories) == 1 else output
    ## generate fname (usage: self.mkfname('ref', 'tmp.txt'))
    def mkfname(self, *path, tmp = False):
        """
        Generate path.
        
        If path provided is not absolute, self.active_directory is assumed to be the root directory.
        
        Arguments:
            *path (str): required, path to output file
                (e.g. self.mkfname('tmp', 'tmp.fasta') --> <self.active_directory>/tmp/tmp.fasta)
            tmp (bool): mark file as temporary (for deletion when self.rm_tmpfiles is called)
        
        Returns:
            str: path
        """
        if os.path.isabs(os.path.join(*path)): path = os.path.join(*path)
        else: path = os.path.join(self.active_directory, *path)
        if tmp: self.tmp_files.add(path)
        return path
    ## mkfname, except it also creates path to directory + empty file if they don't already exist
    def reserve_fname(self, *path, tmp = False, newfile = False, **kwargs):
        """
        Generate new file.
        
        Operates exactly as :meth:`PathHandler.mkfname`, with the additional options of creating an empty file or clearing an existing file.
        
        Arguments:
            *path (str): path to output file
            tmp (bool): mark file as temporary (for deletion when self.rm_tmpfiles is called)
            newfile (bool): clear an existing file at destination if it already exists
            **kwargs: other arguments to pass to self.mkfname
        
        Returns: 
            str: path
        """
        fname = self.mkfname(*path, **kwargs)
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
        ## remove tmp files
        self.rm_tmpfiles()
        ## move files to final destination (lol)
        if self.tmp:
            if self.directory:
                ## create final output directory if doesn't exist
                if not os.path.exists(self.directory):
                    os.makedirs(self.directory, exist_ok = True)
                ## move logfile
                self.logfile.move(self.directory, self.prefix)
                ## copy items (but not the whole directory, just directory minorg created)
                mv_dir_overwrite(self.mkfname(self.tmpdir, self.prefix),
                                 self.mkfname(self.directory, self.prefix),
                                 ignore = [self.directory])
                ## update file locations
                self.new_dirs = [re.sub('^' + re.escape(self.active_directory), self.directory, dirname)
                                 for dirname in self.new_dirs]
                print(f"Output files have been generated in {self.directory}")
            else:
                while True:
                    inpt = input( ("No output directory specified. Project will be deleted."
                                   " Proceed? (Y/N) ") )
                    if inpt in {'N', 'n'}:
                        return
                    elif inpt not in {'Y', 'y'}:
                        print(f"Invalid input: '{inpt}'")
                    else:
                        break
            ## remove tmpdir
            shutil.rmtree(self.tmpdir)
            self.active_directory = self.directory

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
# x.add_reference("/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta", "/mnt/chaelab/rachelle/scripts/minorgpy/test/full/test_full_NRG1_1/reduced_ann/reduced.TAIR10.gff", alias = "TAIR10", replace = True)
# x.add_reference("/mnt/chaelab/rachelle/data/Alyrata/v2.1/Alyrata_384_v1.fa", "/mnt/chaelab/rachelle/scripts/minorgpy/test/full/test_full_NRG1_1/reduced_ann/reduced.araly2.gff", alias = "araly2")
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
# x.valid_grna()
# ## 278 without domain, 89 in NB-aRC 366375
# x.minimumset()

# import sys
# sys.path.append("/mnt/chaelab/rachelle/scripts/minorgpy")
# from minorg.MINORg import MINORg
# x = MINORg()
# x.parse_grna_map_from_file("/mnt/chaelab/rachelle/scripts/minorgpy/test/minimumset/test_gRNA_all.map")
# x.minimumset()
# x.final_map
# x.resolve() ## input 'n'
# x.directory = "/mnt/chaelab/rachelle/scripts/minorgpy/test/minimumset/test_newcheck"
# x.resolve()
# x.final_map ## unchanged :(

class MINORg (PathHandler):
    """
    Tracks parameters, intermediate files, and output file for reuse.
    
    >>> from minorg.MINORg import MINORg
    >>> my_minorg = MINORg(directory = "/my/output/directory", 
                           prefix = "test", tmp = False, keep_tmp = True, thread = 1)
    >>> my_minorg.add_reference("/path/to/TAIR10_Chr.all.fasta", "/path/to/TAIR10_GFF3.genes.gff", alias = "TAIR10", replace = True)
    >>> my_minorg.add_reference("/path/to/Alyrata_384_v1.fa", "/path/to/Alyrata_384_v2.1.gene.gff3", alias = "araly2")
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
    gRNAHits(gRNA = 352)
    >>> my_minorg.filter_gc()
    >>> my_minorg.grna_hits.filter_seqs("GC")
    gRNAHits(gRNA = 370)
    >>> my_minorg.grna_hits.filter_seqs("background", "GC")
    gRNAHits(gRNA = 321)
    >>> my_minorg.filter_feature() ## by default, MINORg only retains gRNA in CDS
    Max acceptable insertion length: 15
    >>> my_minorg.grna_hits.filter_hits("feature")
    gRNAHits(gRNA = 344)
    >>> my_minorg.valid_grna()
    gRNAHits(gRNA = 278)
    >>> my_minorg.minimumset()
    >>> my_minorg.resolve()
    
    Attributes:
        directory (str): [general] final output directory
        auto_update_files (bool): [general] whether to update gRNA .map and FASTA files automatically
            after a subcommand is called
        prefix (str): [general] prefix for output directories and files
        thread (int): [general] maximum number of threads for parallel processing

        blastn (str): [executable] path to blastn executable
        rpsblast (str): [executable] path to rpsblast or rpsblast+ executable
        mafft (str): [executable] path to mafft executable
        bedtools (str): [executable] path to directory containing BEDTools executables.
            Use ONLY IF BEDTools is not in your command-search path.

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
        pam (:class:`~minorg.pam.PAM`): [grna] PAM pattern

        screen_reference (bool): [filter: background] include reference genome for screening
        ot_pamless (bool): [filter: background] ignore absence of PAM when assessing off-target gRNA hits
        offtarget (func):
            [filter: background] function that accepts Biopython's HSP and QueryResult objects and determines whether an off-target gRNA hit is problematic.
            If not set by user, ot_pattern will be used.
            If both offtarget and ot_pattern are not set, ot_mismatch and ot_gap will be used.
            ## TODO!!! integrate this somehow?
        ot_pattern (:class:`~minorg.ot_regex.OffTargetExpression`):
            [filter: background] pattern specifying maximum number of gap(s) and/or mismatch(es)
            within given range(s) of a gRNA in off-target regions to disqualify the gRNA
        ot_unaligned_as_mismatch (bool):
            [filter: background] treat unaligned positions as mismatches (used with ot_pattern)
        ot_unaligned_as_gap (bool):
            [filter: background] treat unaligned positions as gaps
            (specifically as insertions; used with ot_pattern)
        ot_mismatch (int):
            [filter: background] minimum number of mismatches allowed for off-target gRNA hits
        ot_gap (int): [filter: background] minimum number of gaps allowed for off-target gRNA hits
        mask (list): [filter: background] FASTA file of additional sequences to mask in background

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
        
        exclude (str): [filter: exclude] path to FASTA file containing gRNA sequences to exclude

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
    def __init__(self, directory = None, config = None, prefix = "minorg",
                 thread = None, keep_tmp = False, cli = False,
                 auto_update_files = True, **kwargs):
        """Create a MINORg object.
        
        Arguments:
            directory (str): path to output directory
            config (str): path to config.ini file
            prefix (str): prefix for output directories and files (default='minorg')
            thread (int): maximum number of threads for parallel processing
            keep_tmp (bool): retain temporary files
            cli (bool): whether MINORg is called at command line.
                Not currently used, but may be used in future to determine what to print.
            auto_update_files (bool): whether to update gRNA .map and FASTA files automatically
                after a subcommand is called
            **kwargs: other arguments supplied to parent class :class:`PathHandler`
        
        Returns: 
            :class:`~minorg.MINORg.MINORg`: a MINORg object
        """
        self.verbose = False
        if directory is None:
            print(f"Setting current directory ({os.getcwd()}) as output directory")
            directory = os.getcwd()
        super().__init__(directory = directory, keep_tmp = keep_tmp, **kwargs)
        
        ## track whether MINORg is being used at command-line.
        ## this may be used to determine which statements to print
        self._cli = cli
        
        ## parse config.ini file if provided (if None, default values will be stored in Params obj)
        self.params = Params(config)
        
        ## handle log file
        self.logfile = MINORgLogger(level = LOGGING_LEVEL)
        
        ## parameter things
        ## file update behaviour
        self.auto_update_files = auto_update_files
        ## status tracking
        self.subset_gid = []
        self.subset_ref = []
        self._genes_updated_since_alignment = False
        ## general output params
        self.prefix = get_val_default(prefix, default = self.params.prefix.default)
        ## general execution params
        self.thread = get_val_default(thread, default = self.params.thread.default)
        ## binaries/executables
        self.rpsblast = self.params.rpsblast.default
        self.blastn = self.params.blastn.default
        self.mafft = self.params.mafft.default
        self.bedtools = self.params.bedtools.default
        ## RPS-BLAST database
        self.db = self.params.db.default
        # self.db_versions = self.params.db_versions.default ## TODO to allow users to specify domain name and let MINORg automatically fetch pssm-ids
        self.remote_rps = self.params.remote_rps.default
        ## data/annotation
        ## (sets and aliases not supported for MINORg object beyond retrieval of default reference)
        # self.reference = {} ## use MINORg.clear_reference(), MINORg.add_reference(), and MINORg.remove_reference() to update this attribute
        self._reference = {}
        self._parse_reference()
        ## annotation format options
        ## seqid template stuff
        self.ref_seqid_template = "Reference|${source}|${domain}|${n}|${feature}|${complete}|${gene}|${range}"
        self.ref_seqid_gene_pattern = "(?<=\\|stitched\\||\\|complete\\|).+(?=\\|\d+-\d+)"
        self.ref_seqid_feature_pattern = "(?<=\\|)[^|]+(?=\\|stitched\\||\\|complete\\|)"
        self.ref_seqid_source_pattern = "(?<=^Reference\\|)[^|]+(?=\\|)"
        ## homologue params
        # self.domains = {} ## stores domain: {pssmids} mapping. TODO to allow users to specify domain name and let MINORg automatically fetch pssm-ids
        self.pssm_ids = [] ## if multiple pssm_ids are provided, overlapping domains will be combined and output will not distinguish between one pssm_id or another
        self._domain_name = None ## human-readable name for self.pssm_ids
        self.gff_domain = None ## gff file for domains, generated by MINORg
        ## previously self.query_map storing [(<alias>, <path>)]. Now: {<alias>: <path>}
        ## self.query_map is now @property
        self.query = {}
        self._genes = [] ## [<geneA>, <geneB>]
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
        self.grna_eqv = None ## .eqv file, generated by MINORg.write_all_grna_eqv() if not specified by user
        self.screen_reference = True
        self._offtarget = None ## user-provided function that takes in pam, QueryResult obj, and HSP object and returns whether the off-target hit is acceptable. If not set, defaults to using self.ot_mismatch and self.ot_gap (combined with OR) to set lower limits for acceptable off-target
        self.background = {} ## {<alias>: <fname>}
        self._ot_function = self._is_offtarget_aln_nopattern
        self._ot_pattern = None
        self._ot_unaligned_as_mismatch = True
        self._ot_unaligned_as_gap = False
        self.ot_mismatch = self.params.ot_mismatch.default ## ignored if self.ot_pattern != None
        self.ot_gap = self.params.ot_gap.default ## ignored if self.ot_pattern != None
        self.ot_pamless = self.params.ot_pamless.default
        self.gc_min = self.params.gc_min.default
        self.gc_max = self.params.gc_max.default
        self.alignment = None ## fasta file, generated by MINORg
        self.max_insertion = self.params.max_insertion.default
        self.min_within_n = self.params.min_within_n.default
        self.min_within_fraction = self.params.min_within_fraction.default
        self._feature = self.params.feature.default
        self.exclude = None
        ## minimum set
        self.sets = self.params.sets.default
        self.auto = self.params.auto.default
        self._priority = "pos"
        self.prioritise_nr = self.params.prioritise_nr.default
        self.pass_map = None ## .map file, generated by MINORg.write_pass_grna_map() if not specified
        self.pass_fasta = None ## fasta file, generated by MINORg.write_pass_grna_fasta() if not specified
        self.pass_eqv = None ## .eqv file, generated by MINORg.write_pass_grna_eqv() if not specified
        self.final_map = None ## .map file, generated by MINORg.write_all_grna_map() if not specified
        self.final_fasta = None ## fasta file, generated by MINORg.write_all_grna_fasta() if not specified
        self.accept_invalid = False
        self.accept_feature_unknown = False
        self.accept_invalid_field = True

        ## move logfile
        self.logfile.update_filename(self.mkfname(f".log"))
    
    ## getters
    @property
    def genes(self):
        """
        Target gene IDs
        
        :getter: Returns list of target gene IDs (str)
        :type: list
        """
        return self._genes
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
    def reference(self):
        return self._reference
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
            else: raise MINORgError("Too many features")
        else: return self._feature
    @property
    def features(self):
        """
        Target features
        
        :getter: Returns list of target feature(s)
        :type: list
        """
        if non_string_iter(self._feature): return list(self._feature)
        else: return [self._feature]
    @property
    def ot_pam(self):
        """
        Remove gRNA from candidates if off-target hit has PAM site nearby
        
        :getter: Returns whether to consider presence of PAM site for off-target filtering
        :type: bool
        """
        return not self.ot_pamless
    @property
    def ot_pattern(self):
        """
        Off-target mismatch/gap pattern
        
        :getter: Returns OffTargetExpression object
        :type: `:class:~minorg.ot_regex.OffTargetExpression`
        """
        return self._ot_pattern
    @property
    def ot_unaligned_as_mismatch(self):
        """
        Treat unaligned positions as mismatch (used with ot_pattern)
        
        :getter: Returns whether to count unaligned positions as mismatch
        :type: bool
        """
        return self._ot_unaligned_as_mismatch
    @property
    def ot_unaligned_as_gap(self):
        """
        Treat unaligned positions as gaps (specifically insertions; used with ot_pattern)
        
        :getter: Returns whether to count unaligned positions as gaps
        :type: bool
        """
        return self._ot_unaligned_as_gap
    @property
    def passed_grna(self):
        """
        gRNA that have passed background, GC, and feature checks
        
        Relevant attributes: accept_invalid, accept_invalid_field
        
        :getter: Returns gRNA that have passed standard checks
        :type: gRNAHits
        """
        filtered = self.grna_hits.filter_seqs("background", "GC", "exclude",
                                              accept_invalid = self.accept_invalid,
                                              accept_invalid_field = self.accept_invalid_field,
                                              quiet = True)
        filtered = filtered.filter_hits("CDS", "feature",
                                        accept_invalid = (self.accept_invalid or self.accept_feature_unknown),
                                        accept_invalid_field = self.accept_invalid_field,
                                        exclude_empty_seqs=True, quiet = True)
        return filtered
    @property
    def prioritise_nr(self):
        """
        Prioritise non-redundancy for minimum set generation
        
        :getter: Returns whether non-redundancy is prioritised
        :type: bool
        """
        return self._priority == "nr"
    @property
    def prioritize_nr(self):
        """
        Prioritise non-redundancy for minimum set generation
        
        :getter: Returns whether non-redundancy is prioritised
        :type: bool
        """
        return self.prioritise_nr
    @property
    def prioritise_pos(self):
        """
        Prioritise position (i.e. 5' gRNA) for minimum set generation
        
        :getter: Returns whether position (i.e. 5' gRNA) is prioritised
        :type: bool
        """
        return self._priority == "pos"
    @property
    def prioritize_pos(self):
        """
        Prioritise position (i.e. 5' gRNA) for minimum set generation
        
        :getter: Returns whether position (i.e. 5' gRNA) is prioritised
        :type: bool
        """
        return self.prioritise_pos
    
    ## setters
    @genes.setter
    def genes(self, val):
        if val != self.genes:
            self._genes_updated_since_alignment = True
        self.genes = val
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
    def check_recip(self, val: bool):
        self._check_recip = val
    @ot_pam.setter
    def ot_pam(self, val: bool):
        self.ot_pamless = (not val)
    @ot_pattern.setter
    def ot_pattern(self, val: str):
        ## sets self._ot_function to self._is_offtarget_aln_pattern if is valid pattern,
        ## else to self._is_offtarget_aln_nopattern if not valid
        if val is None: self._ot_function = self._is_offtarget_aln_nopattern
        elif isinstance(val, OffTargetExpression):
            self._ot_pattern = val
            self._ot_function = self._is_offtarget_aln_pattern
        elif isinstance(val, str):
            self._ot_pattern = OffTargetExpression(val,
                                                   unaligned_as_mismatch = self.ot_unaligned_as_mismatch,
                                                   unaligned_as_gap = self.ot_unaligned_as_gap)
            self._ot_function = self._is_offtarget_aln_pattern
    @ot_unaligned_as_mismatch.setter
    def ot_unaligned_as_mismatch(self, val: bool):
        self._ot_unaligned_as_mismatch = val
        if self.ot_pattern is not None:
            self._ot_pattern = OffTargetExpression(self.ot_patern.pattern,
                                                   unaligned_as_mismatch = val,
                                                   unaligned_as_gap = self.ot_unaligned_as_gap)
    @ot_unaligned_as_gap.setter
    def ot_unaligned_as_gap(self, val: bool):
        self._ot_unaligned_as_gap = val
        if self.ot_pattern is not None:
            self._ot_pattern = OffTargetExpression(self.ot_patern.pattern,
                                                   unaligned_as_mismatch = self.ot_unaligned_as_mismatch,
                                                   unaligned_as_gap = val)
    @prioritise_nr.setter
    def prioritise_nr(self, val):
        if not isinstance(val, bool):
            raise MINORgError("True/False value only")
        self._priority = "nr" if val else "pos"
    @prioritize_nr.setter
    def prioritize_nr(self, val):
        self.prioritise_nr = val
    @prioritise_pos.setter
    def prioritise_pos(self, val):
        if not isinstance(val, bool):
            raise MINORgError("True/False value only")
        self._priority = "pos" if val else "nr"
    @prioritize_pos.setter
    def prioritize_pos(self, val):
        self.prioritise_pos = val
    
    ## update 'mv' so that it also moves logfile
    def _mv_frm_tmp(self, directory):
        """
        Update output directory.
        
        Contents of the current output directory will also be moved to the new directory if files are
        already being written to output directory (i.e. not currently using a tmp dir)
        
        Arguments:
            directory (str): required, path to new output directory
        """
        src_dir = self.active_directory
        ## move MINORg contents
        if directory is None: return
        directory = os.path.abspath(directory)
        if directory == self.directory: return
        ## if files are being written directly to output directory, move to new destination
        if self.active_directory == self.directory:
            ## make output directory
            if not os.path.exists(directory):
                os.makedirs(directory, exist_ok = True)
            ## copy items (but not the whole directory, just directory minorg created)
            minorg_dir = self.mkfname(self.active_directory, self.prefix, prefix = False)
            if os.path.exists(minorg_dir):
                out_dir = self.mkfname(directory, self.prefix, prefix = False)
                os.makedirs(out_dir, exist_ok = True)
                mv_dir_overwrite(minorg_dir, out_dir, ignore = [out_dir])
            ## update some file paths
            self.tmp_files = set(re.sub('^' + re.escape(self.active_directory), directory, fname)
                                 for fname in self.tmp_files)
            self.new_dirs = [re.sub('^' + re.escape(self.active_directory), directory, dirname)
                             for dirname in self.new_dirs]
            self.active_directory = directory
            ## move logfile
            self.logfile.move(directory, self.prefix)
        self.directory = directory
        ## update stored filenames
        self._update_stored_fnames(src_dir, self.active_directory)
    
    ## update mkfname so it adds a prefix to each file
    def mkfname(self, *path, tmp = False, prefix = True):
        """
        Generate new file name.
        
        If path provided is not absolute, self.active_directory is assumed to be the root directory.
        
        Arguments:
            *path (str): required, path to output file (e.g. self.mkfname('tmp', 'tmp.fasta'))
            tmp (bool): mark file as temporary (for deletion when self.rm_tmpfiles is called)
            prefix (bool): prefix self.prefix to basename
        """
        if path and prefix:
            path_last = path[-1]
            path_last_base = os.path.basename(path_last)
            ## add underscore between prefix and basename if basename doesn't start with '.'
            if path_last_base[0] != '.':
                fname = self.prefix + '_' + path_last_base
            else:
                fname = self.prefix + path_last_base
            path = path[:-1] + type(path)([os.path.join(os.path.dirname(path_last), fname)])
        return super().mkfname(*path, tmp = tmp)
    
    def results_fname(self, *path, reserve = False, newfile = False, **kwargs):
        """
        Generate new file name in results directory (<self.active_directory>/<self.prefix>).
        
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
    
    def resolve(self) -> None:
        """
        Wrapper for super().resolve() that also updates stored filenames.
        Removes logfile if not in CLI mode.
        """
        src_dir = self.active_directory
        dst_dir = self.directory
        ## resolve
        print("resolving")
        super().resolve()
        ## update file paths if dst is not None (i.e. if resolve doesn't delete everything)
        if dst_dir:
            ## update stored filenames
            self._update_stored_fnames(src_dir, dst_dir)
        if not self._cli:
            os.remove(self.logfile.filename)
    
    def _update_stored_fnames(self, src_dir, dst_dir) -> None:
        """
        Updates paths of files (potentially) generated by MINORg if they are in src_dir to dst_dir.
        Attributes updated are: gff_domain, ref_gene, ref_cds, target, grna_fasta,
            grna_map, alignment, pass_map, pass_fasta, final_map, final_fasta
        
        Arguments:
            src_dir (str): path to source directory
            dst_dir (str) path to destination directory
        """
        src_dir_split = os.path.normpath(src_dir).split(os.path.sep)
        dst_dir_split = os.path.normpath(dst_dir).split(os.path.sep)
        def update_fname(fname):
            if not fname: return fname
            fname_split = os.path.normpath(fname).split(os.path.sep)
            if fname_split[:len(src_dir_split)] == src_dir_split:
                return os.path.join(*(dst_dir_split + fname_split[len(src_dir_split):]))
            return fname
        ## update
        for attr in ("gff_domain", "ref_gene", "ref_cds", "target", "grna_fasta", "grna_map",
                     "alignment", "pass_map", "pass_fasta", "final_map", "final_fasta"):
            try:
                setattr(self, attr, update_fname(getattr(self, attr)))
            except Exception as e:
                print(attr, getattr(self, attr))
                raise e
        return

    ###############
    ##  GETTERS  ##
    ###############
    
    def get_ref_seqid(self, seqid, attr):
        """
        Extract attribute of sequence from sequence name (for sequences generated by self.target)
        
        Arguments:
            seqid (str): required, sequence name
            attr (str): required, attribute to return.
                Valid attributes: 'gene', 'feature', 'source'.
        
        Returns: 
            str
        """
        if attr == "gene": return re.search(self.ref_seqid_gene_pattern, seqid).group(0)
        elif attr == "feature": return re.search(self.ref_seqid_feature_pattern, seqid).group(0)
        elif attr == "source": return re.search(self.ref_seqid_source_pattern, seqid).group(0)
        else: raise MINORgError(f"Unknown attribute: {attr}")
    
    def valid_grna(self, *check_names, accept_invalid = None, accept_invalid_field = None):
        """
        gRNA that have passed background, GC, and feature checks, as well as any user-set checks
        
        Relevant attributes: accept_invalid, accept_invalid_field
        
        Arguments:
            *check_names (str): optional, checks to include.

                If not specified: all checks, both standard and non-standard, will be used.
        
                If 'standard': only background, GC, and feature checks will be used.
                Can be used in combination with non-standard check names.
        
                If 'nonstandard': only non-standard checks will be used.
                Can be used in combination with standard check names.
        
        Returns:
            :class:`~minorg.grna.gRNAHits`
        """
        accept_invalid = (accept_invalid if accept_invalid is not None
                          else self.accept_invalid)
        accept_invalid_field = (accept_invalid_field if accept_invalid_field is not None
                                else self.accept_invalid_field)
        ## expand check names
        extant_check_names = self.check_names(standard = True, nonstandard = True, hit = True, seq = True)
        if not check_names:
            check_names = copy.copy(extant_check_names)
        else:
            if "standard" in check_names:
                check_names.remove("standard")
                check_names.extend(self.check_names(standard = True, nonstandard = False))
            if "nonstandard" in check_names:
                check_names.remove("nonstandard")
                check_names.extend(self.check_names(standard = False, nonstandard = True))
        ## check for unknown check names
        unknown_check_names = [check_name for check_name in check_names if check_name not in extant_check_names]
        if unknown_check_names:
            warnings.warn( f"Unknown check name(s): {','.join(unknown_check_names)}" , MINORgWarning)
        check_names = [check_name for check_name in check_names if check_name in extant_check_names]
        ## filter
        seqs_checks = set(check_names).intersection({"background", "GC", "exclude"})
        hits_checks = set(check_names) - seqs_checks
        filtered = self.grna_hits
        if seqs_checks:
            filtered = filtered.filter_seqs(*seqs_checks,
                                            accept_invalid = self.accept_invalid,
                                            accept_invalid_field = self.accept_invalid_field,
                                            quiet = True, report_invalid_field = True)
        if hits_checks:
            filtered = filtered.filter_hits(*hits_checks,
                                            accept_invalid = (self.accept_invalid
                                                              or self.accept_feature_unknown),
                                            accept_invalid_field = self.accept_invalid_field,
                                            exclude_empty_seqs = True, quiet = True,
                                            report_invalid_field = True)
        return filtered
    def check_names(self, seq = True, hit = True, standard = True, nonstandard = True) -> list:
        """
        Get gRNA check names
        
        Arguments:
            seq (bool): get check names for gRNA sequences (gRNASeq objects)
            hit (bool): get check names for gRNA hits (gRNAHit objects)
            standard (bool): get standard check names
            nonstandard (bool): get nonstandard check names
        
        Returns: 
            list: List of check names (str)
        """
        if not (seq or hit):
            raise MINORgError("Either one or both of 'seq=True' or 'hit=True' is required.")
        if not (standard or nonstandard):
            raise MINORgError("Either one or both of 'standard=True' or 'nonstandard=True' is required.")
        check_names = (self.grna_hits.check_names if (seq and hit)
                       else self.grna_hits.check_names_seqs if seq
                       else self.grna_hits.check_names_hits)
        if standard and nonstandard:
            return check_names
        else:
            standard_check_names = {"background", "GC", "exclude", "CDS", "feature"}
            if standard:
                return list(set(check_names).intersection(standard_check_names))
            elif nonstandard:
                return list(set(check_names) - standard_check_names)
    
    ####################
    ##  ADD FUNCIONS  ##
    ####################
    
    def add_query(self, fname, alias = None):
        """
        Add file to be queried.
        
        Arguments:
            fname (str): required, path to file
            alias (str): optional, alias for query file
        """
        if not alias:
            max_id = max([int(re.search("\d+$", y).group(0)) for y in
                          ["query_0"] + [x for x in self.query if re.search("^query_\d+$", x)]])
            alias = f"query_{str(max_id + 1).zfill(3)}"
        elif alias in self.query:
            raise MINORgError(f"ID '{alias}' is already in use.")
        self.query[alias] = IndexedFasta(str(fname))
        return
    
    def remove_query(self, alias):
        """
        Remove file for querying.
        
        Arguments:
            alias (str): required, alias of query file to remove
        """
        if alias not in self.query:
            raise MINORgError(f"'{alias}' is not a valid query alias.")
        del self.query[alias]
        return

    def add_background(self, fname, alias = None):
        """
        Add file for background filter.
        
        Arguments:
            fname (str): required, path to file
            alias (str): optional, alias for background file.
                Used when writing mask report and for removing background.
        """
        if not alias:
            max_id = max([int(re.search("\d+$", y).group(0)) for y in
                          ["bg_0"] + [x for x in self.background if re.search("^bg_\d+$", x)]])
            alias = f"bg_{str(max_id + 1).zfill(3)}"
        elif alias in self.background:
            raise MINORgError(f"ID '{alias}' is already in use.")
        self.background[alias] = IndexedFasta(str(fname))
        return
    
    def remove_background(self, alias):
        """
        Remove file for background filter.
        
        Arguments:
            alias (str): required, alias of background file to remove
        """
        if alias not in self.background:
            raise MINORgError(f"'{alias}' is not a valid background alias.")
        del self.background[alias]
        return
    
    def parse_grna_map_from_file(self, fname):
        """
        Read candidate gRNA from .map file output by MINORg.
        
        Arguments:
            fname (str): required, path to MINORg .map file
        """
        grna_hits = gRNAHits()
        # grna_hits.parse_from_mapping(fname, version = None)
        grna_hits.parse_from_mapping(fname, targets = self.target)
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
                Autogenerated using self.prefix and self.active_directory if not provided.
        """
        if fout is None: fout = get_val_default(self.grna_fasta, self.results_fname("gRNA_all.fasta"))
        self.reserve_fname(fout, prefix = False)
        self.grna_fasta = fout
        self.grna_hits.write_fasta(fout, write_all = True)
    
    def write_all_grna_map(self, fout = None, write_checks = True):
        """
        Write .map file of all candidate gRNA.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.active_directory if not provided.
            write_checks (bool): write all check statuses
        """
        if fout is None: fout = get_val_default(self.grna_map, self.results_fname("gRNA_all.map"))
        self.reserve_fname(fout, prefix = False)
        self.grna_map = fout
        self.grna_hits.write_mapping(fout, version = 5, write_all = True,
                                     write_checks = write_checks)
    
    def write_all_grna_eqv(self, fout = None):
        """
        Write .eqv file (grouping gRNA by equivalent coverage) of all candidate gRNA.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.active_directory if not provided.
        """
        if fout is None: fout = self.results_fname("gRNA_all.eqv")
        self.reserve_fname(fout, prefix = False)
        self.grna_eqv = fout
        self.grna_hits.write_equivalents(fout, write_all = True)
    
    def write_pass_grna_fasta(self, fout = None):
        """
        Write FASTA file of gRNA that pass all valid checks.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.active_directory if not provided.
        """
        if fout is None: fout = get_val_default(self.pass_fasta, self.results_fname("gRNA_pass.fasta"))
        self.reserve_fname(fout, prefix = False)
        self.pass_fasta = fout
        self.valid_grna().write_fasta(fout, write_all = True)
    
    def write_pass_grna_map(self, fout = None):
        """
        Write .map file of gRNA that pass all valid checks.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.active_directory if not provided.
        """
        if fout is None: fout = get_val_default(self.pass_map, self.results_fname("gRNA_pass.map"))
        self.reserve_fname(fout, prefix = False)
        self.pass_map = fout
        self.valid_grna().write_mapping(fout, version = 5, write_all = True, write_checks = False)
    
    def write_pass_grna_eqv(self, fout = None):
        """
        Write .eqv file (grouping gRNA by equivalent coverage) of gRNA that pass all valid checks.
        
        Arguments:
            fout (str): optional, absolute path to output file.
                Autogenerated using self.prefix and self.active_directory if not provided.
        """
        if fout is None: fout = self.results_fname("gRNA_pass.eqv")
        self.reserve_fname(fout, prefix = False)
        self.pass_eqv = fout
        self.valid_grna().write_equivalents(fout, write_all = True)
    
    def write_pass_grna_files(self, fasta = None, map = None, eqv = None):
        self.write_pass_grna_map(fout = map)
        self.write_pass_grna_fasta(fout = fasta)
        self.write_pass_grna_eqv(fout = eqv)
    
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
                                 " MINORg.add_reference(<path to FASTA>, <path to GFF3>, alias = '<alias>')"),
                               MINORgWarning)
        ## get reference
        val_ref = get_val_default(reference, self.params.reference.default)
        if val_ref:
            none_val = INDV_GENOMES_NONE
            mapped = valid_aliases(aliases = val_ref, lookup = reference_aliases,
                                   none_value = INDV_GENOMES_NONE, all_value = INDV_GENOMES_ALL,
                                   param = self.params.reference, return_mapping = True)
            for alias, ref in mapped.items():
                fasta, ann, genetic_code, attr_mod = reference_aliases[alias]
                self.add_reference(fasta, ann, alias = alias, genetic_code = genetic_code,
                                   attr_mod = self.params.parse_attr_mod(attr_mod))
        return
    
    ############################
    ##  REFERENCE/ANNOTATION  ##
    ##       FUNCTIONS        ##
    ############################
    
    ## takes path to FASTA and GFF/BED file and indexes them before storing.
    def add_reference(self, fasta, annotation, alias = None, genetic_code = 1, attr_mod = {},
                      memsave = True, replace = False):
        """
        Add reference genome.
        
        Arguments:
            fasta (str): required, path to fasta file
            annotation (str): required path to GFF3 file
            alias (str): optional, name of fasta-annotation combo
            genetic_code (str or int): NCBI genetic code (default=1)
            attr_mod (dict): mapping for non-standard attribute field names in GFF3 file (default={})
            replace (bool): replace any existing reference with same alias (default=False)
        """
        if alias in self.reference and not replace:
            raise MINORgError(f"ID '{alias}' is already in use.")
        if not alias:
            max_id = max([int(re.search("\d+$", y).group(0)) for y in
                          ["Reference_0"] + [x for x in self.reference if re.search("^Reference_\d+$", x)]])
            alias = f"Reference_{str(max_id + 1).zfill(3)}"
        self._reference[alias] = AnnotatedFasta(fasta, annotation, attr_mod = attr_mod,
                                                genetic_code = genetic_code, name = alias, memsave = memsave)
        return
    def remove_reference(self, alias):
        """
        Remove reference genome.
        
        Arguments:
            alias (str): required, alias of reference genome
        """
        if alias in self.reference:
            del self._reference[alias]
        return
    def clear_reference(self):
        """
        Remove all reference genomes.
        """
        self._reference = {}
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
        self.add_reference(ext_fasta, ext_gff, alias = "Extended")
        return
    def _get_reference_seq(self, ids, *features, adj_dir=True, by_gene=False, mktmp=None, ref=None,
                           seqid_template="Reference|$source|$gene|$isoform|$feature|$n|$complete",
                           translate=False, isoform_lvl=None, gff_domain=None, fout=None, complete=False,
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
            feature_seqs = {}
            for ref_alias in ref_aliases:
                ref = self.reference[ref_alias]
                ## process each gene separately
                for feature_id in ids:
                    ## fill in gene and source
                    ft_seqid_template = string.Template(seqid_template).safe_substitute(
                        gene = feature_id,
                        source = ref_alias,
                        domain = (self.domain_name if self.domain_name else
                                  ','.join(map(str, self.pssm_ids)) if self.pssm_ids
                                  else "gene"))
                    ## get features at isoform_lvl if specified
                    if isoform_lvl:
                        feature_ids = [x.get_attr("ID", fmt=list)[0] for x in
                                       ref.annotation.subset(feature_id, features = isoform_lvl,
                                                             subfeatures = True, sort = False)]
                    else:
                        feature_ids = [feature_id]
                    seqs = ref.get_feature_seq(*feature_ids, fout = None, feature_type = feature,
                                               fout_gff = gff_domain, adj_dir = adj_dir, complete = complete,
                                               translate = translate, by_gene = by_gene, mktmp = mktmp,
                                               db = self.db, rpsblast = self.rpsblast, pssm_id = self.pssm_ids,
                                               remote_rps = self.remote_rps,
                                               seqid_template = ft_seqid_template,
                                               apply_template_to_dict = True)
                    feature_seqs = {**feature_seqs,
                                    **dict(itertools.chain(*map(lambda x: x.items(), seqs.values())))}
            if fout:
                dict_to_fasta(feature_seqs, fout)
            output = {**output, **feature_seqs}
        return output
    def get_reference_seq(self, *features, adj_dir=True, by_gene=False, ref=None,
                          seqid_template="Reference|$source|$gene|$isoform|$feature|$n|$complete",
                          translate=False, isoform_lvl=None, fout=None, complete=False,
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
                The default template is "Reference|$source|$gene|$isoform|$feature|$n|$complete".
                
                    - $source: reference genome alias
                    - $gene: gene/isoform ID
                    - $isoform: isoform ID if by_gene = False, else same as $gene
                    - $feature: GFF3 feature type
                    - $n: if multiple domains are present, they will be numbered according to 
                      proximity to 5' of sense strand
                    - $complete: 'complete' if ``complete=True`` else 'stitched'
            
            adj_dir (bool): output sense strand
            by_gene (bool): merge sequences from all isoforms of a gene
            isoform_lvl (bool): if GFF features 'mRNA' or 'protein' specified, 
                sequences will be separated by mRNA/protein isoforms
            complete (bool): merge output range(s) together into single range.
                Output will stll be a list of tuple. (e.g. [(<smallest start>, <largest end>)])
            translate (bool): translate sequence. Should be used with adj_dir.
            fout (str): optional, path to output file if features is not specified.
            **fouts (str): optional, path to output files for each feature if features is specified.
                Example: self.get_reference_seq("CDS", CDS = <path to file>) returns dict AND writes to file
        
        Returns:
            dict: sequence(s) of reference genes or isoforms (specified in self.genes), grouped by feature.
                  Format: {<feature>: {<seqid>: <seq>}}
        """
        return self._get_reference_seq(self.genes, *features, adj_dir = adj_dir,
                                       by_gene = by_gene, ref = ref, complete = complete,
                                       seqid_template = seqid_template, fout = fout,
                                       translate = translate, isoform_lvl = isoform_lvl,
                                       # apply_template_to_dict=False,
                                       **fouts)
    def generate_ref_gene_cds(self, ref_dir = "ref", quiet = True, domain_name = None):
        """
        Get and store self.genes' gene and CDS sequences, and GFF data for domains.
        
        The default sequence ID template is:
        "Reference|$source|$domain|$n|$feature|$complete|$gene|$range".
        
            - $source: reference genome alias
            - $domain: PSSM ID or domain name ('gene' if not specified)
            - $n: if multiple domains are present, they will be numbered according to 
              proximity to 5' of sense strand
            - $feature: GFF3 feature type
            - $complete: 'complete' if sequence includes intervening sequences not of feature type $feature.
              'stitched' if sequence is concatenated from features of feature type $feature.
            - $gene: gene/isoform ID
            - $range: range of sequence in gene
            
        """
        if not all(map(lambda gid: gid in self.subset_gid, self.genes)):
            self.subset_annotation()
        ## instantiate variables
        printi = ((lambda msg: None) if quiet else (print))
        # seqid_template = "Reference|${source}|${domain}|${n}|${feature}|${complete}|${gene}|${isoform}"
        seqid_template = self.ref_seqid_template
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
        self.logfile.devsplain(f"{self.genes}, {self.ref_gene}, {self.ref_cds}, {self.gff_domain}")
        return
    
    #########################
    ##  SUBCMD: HOMOLOGUE  ##
    #########################


    def _homologue_indv(self, fout, fasta_query, quiet = True, lvl = 0,
                        check_recip = None, relax_recip = None,
                        **for_merge_hits_and_filter):
        """
        Find homologue in single FASTA file.

        1. Execute BLASTN of reference sequences (``fasta_complete`` and ``fasta_cds``) against ``fasta_query``.
        2. :func:`~minorg.extract_homologue.merge_hits_and_filter`: Merge hits based on proximity to each other, filtering by length and % identity, to generate candidate homologues.
        3. :func:`~minorg.extract_homologue.recip_blast_multiref`: If check_reciprocal=True, execute BLASTN of candidate homologues to reference genome.
            - :func:`~minorg.extract_homologue.recip_blast_multiref`: Remove candidate homolougues if the hit with the best bitscore is NOT to a gene in ``gene``.
            - :func:`~minorg.extract_homologue.recip_blast_multiref`: If relax=False, candidate homologues which best bitscore hit overlaps with gene in ``gene`` AND ALSO a gene NOT IN ``gene`` will be removed.

        Arguments:
            fout (str): required, path to output FASTA file in which to write homologues
            fasta_query (str): required, path to FASTA file in which to search for homologues
            quiet (bool): print only essential messages
            lvl (int): optional, indentation level of printed messages
            **for_merge_hits_and_filter: additional arguments for 
                :func:`~minorg.extract_homologue.for_merge_hits_and_filter`
        """
        tmp_pref = os.path.basename(fout)
        ## execute BLASTN of reference genes against query
        ## note: directory must be valid path
        tmp_dir = self.results_fname(f"homologue_{tmp_pref}", prefix = False, tmp = True)
        tsv_blast_ref = self.reserve_fname(self.mkfname(tmp_dir,
                                                        f"{tmp_pref}_tmp_blastn_ref.tsv"))
        tsv_blast_cds = self.reserve_fname(self.mkfname(tmp_dir,
                                                        f"{tmp_pref}_tmp_blastn_cds.tsv"))
        from Bio.Blast.Applications import NcbiblastnCommandline
        blast(blastf = NcbiblastnCommandline, cmd = self.blastn,
              header = None, ## default fields
              fout = tsv_blast_ref, query = self.ref_gene, subject = fasta_query)
        blast(blastf = NcbiblastnCommandline, cmd = self.blastn,
              header = None, ## default fields
              # header = "qseqid,sseqid,pident,length,sstart,send",
              fout = tsv_blast_cds, query = self.ref_cds, subject = fasta_query)
        ## check for hits
        # if not parse_get_data(tsv_blast_ref)[1]:
        if is_empty_file(tsv_blast_ref):
            # raise MINORgError("No blast hits during homologue search, exiting programme.")
            ## write empty file
            with open(fout, "w+") as f:
                pass
            return
        ## extract homologue
        merge_hits_and_filter(blast6_fname = tsv_blast_ref, fout = fout, fasta = fasta_query,
                              blast6cds_fname = tsv_blast_cds, quiet = quiet, **for_merge_hits_and_filter)
        ## check reciprocal
        check_recip = get_val_default(check_recip, self.check_recip)
        if check_recip and self.genes:
            recip_blast_multiref(fasta_target = fout,
                                 genes = self.genes, blastn = self.blastn, bedtools = self.bedtools,
                                 gff = self.annotations, fasta_ref = self.assemblies,
                                 attribute_mod = self.attr_mods,
                                 relax = get_val_default(check_recip, self.relax_recip),
                                 directory = tmp_dir,
                                 keep_tmp = self.keep_tmp, lvl = lvl, quiet = quiet)
        ## remove tmp files
        if not self.keep_tmp:
            self.rm_tmpfiles(tsv_blast_ref, tsv_blast_cds, tmp_dir)
        return

    def _homologue(self, fout, ref_dir = "ref", quiet = True, domain_name = None,
                   minlen = None, minid = None, mincdslen = None, merge_within = None,
                   check_recip = None, relax_recip = None, check_id_before_merge = None):
        """
        Discovery homologues based on self's attributes
        """
        printi = ((lambda msg: None) if quiet else (print))
        ## get reference data
        self.generate_ref_gene_cds(ref_dir = ref_dir, quiet = quiet, domain_name = domain_name)
        ## define args
        args = [[i, self.results_fname(f"{i}_targets.fasta"),
                 (fasta_query if not isinstance(fasta_query, IndexedFasta)
                  else fasta_query.filename)]
                for i, fasta_query in self.query_map]
        ## define function for multithreading
        def f(args):
            indv_i, fout, fasta_query = args
            self._homologue_indv(fout = fout, fasta_query = fasta_query, indv_i = indv_i,
                                 check_recip = get_val_default(check_recip, self.check_recip),
                                 relax_recip = get_val_default(check_recip, self.relax_recip),
                                 ## merge_hits_and_filter options
                                 min_len = get_val_default(minlen, self.minlen),
                                 min_id = get_val_default(minid, self.minid),
                                 min_cds_len = get_val_default(mincdslen, self.mincdslen),
                                 check_id_before_merge = get_val_default(check_id_before_merge,
                                                                         self.check_id_before_merge),
                                 merge_within_range = get_val_default(merge_within, self.merge_within))
            return fout
        ## get homologues
        printi("Finding homologues")
        if self.query:
            self.logfile.devsplain(f"chkpt1: {[x[0] for x in self.query_map]}")
            tmp_fouts = imap_progress(f, args, threads = self.thread)
            cat_files(tmp_fouts, fout, remove = True)
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
        >>> my_minorg.add_reference("/path/to/TAIR10_Chr.all.fasta", "/path/to/TAIR10_GFF3.genes.gff", alias = "TAIR10", replace = True)
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
        if self.pssm_ids and not (self.rpsblast and self.db):
            raise Exception(("If self.pssm_ids has been provided,"
                             " self.rpsblast (path to rpsblast or rpsblast+ executable"
                             " OR command name if available at command line)"
                             " AND self.db (path to RPS-BLAST database)"
                             " are required"))
        ## get kwargs
        kwargs = locals()
        del kwargs["self"]
        ## create output file
        fout = self.results_fname(f"{self.domain_name}_targets.fasta", newfile = True)
        ## get homologues
        self._homologue(fout, ref_dir = "ref", **kwargs)
        ## append reference genes to targets file
        if self.query_reference:
            tmp = self.results_fname(f"tmp.fasta", tmp = True)
            cat_files([fout, self.ref_gene], tmp, remove = False)
            os.rename(tmp, fout)
        self.target = fout
        ## check if any target sequences found
        if not(fasta_to_dict(fout)):
            warnings.warn("No target sequences found.", MINORgWarning)
            if not self._cli and not self.query_reference:
                print( ("Did you mean to query reference genes?"
                        " If so, don't forget to use '<MINORg object>.query_reference = True'.") )
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
        if self.auto_update_files:
            self.write_all_grna_map()
            self.write_all_grna_eqv()
        # self.grna_hits.write_fasta(self.grna_fasta, write_all = True)
        return
    
    ######################
    ##  SUBCMD: FILTER  ##
    ######################
    
    def _mask_ontargets(self, *mask_fnames):
        tmp_f = self.reserve_fname("masking_tmp.blast.tsv", newfile = True, tmp = True)
        ## reset
        self.masked = {}
        ## mask function
        def _mask_identical(alias_subject, descr = None):
            fnames = [IndexedFasta(subject).filename for subject in assign_alias(alias_subject).values()]
            fnames = [fname for fname in fnames if fname not in self.masked] ## avoid repeats
            msg = ((lambda curr, last: f"{descr}: {curr}/{last} done.") if descr is not None
                   else (lambda curr, last: f"{curr}/{last} done."))
            args = [(mask_fname, subject, self.results_fname("masking", f"tmp_{i}_{j}.tsv", reserve = True))
                    for i, mask_fname in enumerate(mask_fnames) if mask_fname
                    for j, subject in enumerate(fnames)]
            ## define function for multithreading
            def f(args):
                mask_fname, subject, tmp_f = args
                output = mask_identical(mask_fname, subject, tmp_f, blastn = self.blastn)
                os.remove(tmp_f)
                return [subject, output]
            outputs = imap_progress(f, args, threads = self.thread, msg = msg)
            output = {subject:
                      tuple(set(itertools.chain(*[masked for fname, masked in outputs if fname == subject])))
                      for subject in fnames}
            return output
        ## mask in reference
        ## (we don't use ranges of self.genes directly because some genes may be identical
        ##  but NOT in self.genes and we want to mask those too)
        self.masked = {**self.masked, **_mask_identical(self.assemblies, descr = "reference")}
        ## mask in query
        self.masked = {**self.masked, **_mask_identical(self.query, descr = "query")}
        ## mask in backgrounds
        self.masked = {**self.masked, **_mask_identical(self.background, descr = "background")}
        ## slate 'masking' directory for removal
        self.results_fname("masking", prefix = False, tmp = True)
        return
    
    def write_mask_report(self, fout):
        """
        Write mask report to file.
        
        Details masked regions as well as which sequence each region is identical to.
        
        Arguments:
            fout (str): required, path to output file
        """
        with open(fout, 'w') as f:
            inv_fnames = {str(fname): alias for alias, fname in
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
        
        Returns: 
            bool: whether a HSP hit is outside of a masked region
        """
        return self._is_offtarget_pos(hsp, ifasta)
    
    def _is_offtarget_aln_pattern(self, hsp, query_result, **kwargs):
        return self.ot_pattern(BlastHSP(hsp, query_result, None))
    
    def _is_offtarget_aln_nopattern(self, hsp, query_result,
                                    mismatch_check = True, gap_check = True, fully_aligned_check = True):
        # ## note: XX_pass --> XX exceeds threshold for non-problematic off-target
        # qlen = query_result.seq_len
        # mismatch_excess = (max(1, self.ot_mismatch) - 1) - (qlen - hsp.ident_num)
        # gap_excess = (max(1, self.ot_gap) - 1) - hsp.gap_num
        # unaligned = qlen - hsp.query_end + hsp.query_start
        # mismatch_pass = (mismatch_excess >= 0) if mismatch_check else True
        # gap_pass = (gap_excess >= 0) if gap_check else True
        # fully_aligned = unaligned == 0 if fully_aligned_check else True
        # return fully_aligned and mismatch_pass and gap_pass
        ## note XX_problematic --> XX meets problematic off-target threshold
        qlen = query_result.seq_len
        aligned = max(hsp.query_start, hsp.query_end) - min(hsp.query_start, hsp.query_end)
        unaligned = qlen - aligned
        mismatch_excess_problematic_quota = (max(1, self.ot_mismatch) - 1) - (hsp.mismatch_num)
        gap_excess_problematic_quota = (max(1, self.ot_gap) - 1) - hsp.gap_num
        mismatch_problematic = (mismatch_excess_problematic_quota >= 0)
        gap_problematic = (gap_excess_problematic_quota >= 0)
        align_problematic = (((mismatch_excess_problematic_quota + gap_excess_problematic_quota)
                              - unaligned) >= 0)
        ## combine all checks
        is_problematic = True
        if gap_check: is_problematic = is_problematic and gap_problematic
        if mismatch_check: is_problematic = is_problematic and mismatch_problematic
        if fully_aligned_check: is_problematic = is_problematic and align_problematic
        return is_problematic
    
    def is_offtarget_aln(self, hsp, query_result, **kwargs):
        """
        Check whether (potentially off-target) gRNA hit aligns too well.
        
        Used in combination with is_offtarget_pos to determine if an off-target gRNA hit could be problematic.
        
        Arguments:
            hsp (Biopython HSP): required
            query_result (Biopython QueryResult): required
            **kwargs: other arguments supplied passed to :meth:`MINORg._is_offtarget_aln`
        
        Returns: 
            bool: whether a HSP hit meets the threshold for problematic off-target gRNA hit
        """
        # return self._is_offtarget_aln(hsp, query_result, **kwargs)
        return self._ot_function(hsp, query_result, **kwargs)
    
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
        
        Returns: 
            bool: whether a HSP hit is in close proximity and correct orientation to a PAM site
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
        
        Returns: 
            bool: whether a gRNA hit is off-target and problematic
        """
        ifasta = IndexedFasta(ifasta)
        return self.is_offtarget_pos(hsp, ifasta) and self.is_offtarget_aln(hsp, query_result) and \
            self.is_offtarget_pam(hsp, query_result, ifasta)
    
    def _offtarget_hits(self, fasta_subject, keep_output = False, fout = None, thread = None):
        thread = thread if thread is not None else self.thread
        if not fout:
            fout = self.mkfname("tmp.tsv")
            keep_output = False
        if self.grna_fasta is None:
            self.write_all_grna_fasta()
        from Bio.Blast.Applications import NcbiblastnCommandline
        fields = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", ## standard fields part 1
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", ## standard fields part 2
                  "gaps", "qlen", "btop", "nident"] ## custom fields
        blast(NcbiblastnCommandline, header = fields, fout = fout, outfmt = 6,
              query = self.grna_fasta, subject = IndexedFasta(fasta_subject).filename,
              cmd = self.blastn, task = "blastn-short",
              # mt_mode = (0 if self.thread == 1 else 1),
              num_threads = thread)
        offtarget = set(hsp.query_id
                        for query_result in searchio_parse(fout, "blast-tab", fields = fields)
                        for hit in query_result for hsp in hit
                        if self.is_offtarget_hit(hsp, query_result, fasta_subject))
        if not keep_output: os.remove(fout)
        return offtarget
    
    def filter_background(self, *other_mask_fnames, keep_blast_output = None, mask_reference = True):
        """
        Set background filter check for candidate gRNAs.
        
        Masks target sequences in all FASTA files to be screened for off-target, BLASTs all candidate gRNAs to those FASTA files, and assesses each BLAST hit individually for whether they could potentially be problematic.
        
        Relevant attributes: screen_reference, ot_pamless, ot_mismatch, ot_gap

        Relevant methods: add_background, remove_background
        
        Arguments:
            *other_mask_fnames (str): optional, paths to other FASTA files not in self.background that are also to be screened for off-target
            keep_blast_output (bool): retain BLAST output file. 
                If not provided, defaults to self.keep_tmp. ``False`` deletes it.
            mask_reference (bool): mask reference genes (default=True)
        """
        if not self.grna_hits:
            raise MINORgError( ("MINORg.filter_background requires gRNA"
                                " (generated with self.grna())") )
        if keep_blast_output is None: keep_blast_output = self.keep_tmp
        if self.grna_fasta is None:
            self.write_all_grna_fasta() ## write gRNA to file so we can BLAST it
        if mask_reference and not self.ref_gene and self.genes:
            ref_to_mask = self.mkfname("tmp_ref_to_mask.fasta", tmp = True)
            self.get_reference_seq(fout = ref_to_mask)
            self.ref_gene = ref_to_mask
        print("Masking on-targets")
        self._mask_ontargets(self.target,
                             *([self.ref_gene] if mask_reference and self.genes else []),
                             *self.mask,
                             *other_mask_fnames)
        self.write_mask_report(self.results_fname("masked.txt"))
        print("Finding off-targets")
        excl_seqid = set()
        def get_excl_seqids(alias_ifasta, descr = None, screened = set()):
            fname = lambda x: x if not isinstance(x, IndexedFasta) else x.filename
            blast_thread = max(1, int(self.thread/len(alias_ifasta)))
            msg = ((lambda curr, last: f"{descr}: {curr}/{last} done.") if descr is not None
                   else (lambda curr, last: f"{curr}/{last} done."))
            args = [(ifasta,
                     self.mkfname(f"hit_{descr}_{alias}.tsv" if descr else f"hit_{alias}.tsv"))
                    for alias, ifasta in alias_ifasta.items()
                    if fname(ifasta) not in screened]
            ## define function for multithreading
            def f(args):
                ifasta, fout = args
                output = self._offtarget_hits(ifasta, fout = fout, keep_output = keep_blast_output,
                                              thread = blast_thread)
                return output
            outputs = imap_progress(f, args, threads = self.thread, msg = msg)
            screened |= set(fname(ifasta) for ifasta in alias_ifasta.values())
            return set(itertools.chain(*outputs)), screened
        screened = set()
        if (self.screen_reference or self.query_reference) and self.reference:
            to_excl, screened = get_excl_seqids(self.reference, descr = "reference", screened = screened)
            excl_seqid |= to_excl
        if self.query:
            to_excl, screened = get_excl_seqids(self.query, descr = "query", screened = screened)
            excl_seqid |= to_excl
        if self.background:
            to_excl, screened = get_excl_seqids(self.background, descr = "background", screened = screened)
            excl_seqid |= to_excl
        ## parse bg check status
        grna_seqs = fasta_to_dict(self.grna_fasta)
        grna_screened = tuple(grna_seqs.values())
        grna_failed = set(str(grna_seqs[seqid]).upper() for seqid in excl_seqid)
        grna_passed = set(str(seq).upper() for seq in grna_screened if str(seq).upper() not in grna_failed)
        ## record bg check status
        self.grna_hits.set_seqs_check("background", True, grna_passed)
        self.grna_hits.set_seqs_check("background", False, grna_failed)
        ## update .map and FASTA files
        if self.auto_update_files:
            self.write_all_grna_map()
            self.write_pass_grna_files()
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
        ## update .map and FASTA files
        if self.auto_update_files:
            self.write_all_grna_map()
            self.write_pass_grna_files()
        return
    
    def _adjust_feature_range(self, seq, gene, ref_alias, feature, domain = False, **kwargs):
        gene_ann = self.reference[ref_alias].annotation.get_id(gene, output_list = False)
        feature_ranges = ranges_union([[(x.start, x.end)] for x in
                                        self.reference[ref_alias].annotation.get_subfeatures_full(gene, feature_types = feature)])
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
        if not (self.ref_cds and self.ref_gene):
            print( ("MINORg.align_reference_and_targets requires"
                    " MINORg.ref_cds and MINORg.ref_gene."
                    " Generating required files.") )
            self.generate_ref_gene_cds()
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
            new_cds_seq = Seq('-'*(len(seq_gene)))
            for r in adj_ranges:
                start, end = r
                new_cds_seq = new_cds_seq[:start] + seq_gene[start:end] + new_cds_seq[end:]
            aln[seqid_cds] = new_cds_seq
        dict_to_fasta(aln, fout)
        os.remove(tmp_f)
        self._genes_updated_since_alignment = False
        return
    
    ## untested
    def align_reference_and_targets(self, domain_name = None, realign_only_if_updated = True):
        """
        Align reference genes and targets. Path to FASTA file generated will be stored in self.alignment.
        
        Arguments:
            domain_name (str): optional, domain name to be used when naming output file
        """
        # if not self.ref_cds and self.ref_gene and self.target:
        #     raise MINORgError( ("MINORg.align_reference_and_targets requires"
        #                         " MINORg.ref_cds, MINORg.ref_gene, and MINORg.targets.") )
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
        if not self.alignment or self._genes_updated_since_alignment:
            self.align_reference_and_targets()
        alignment = fasta_to_dict(self.alignment)
        genes = {'|'.join(seqid.split('|')[6:-1]):
                 seqid for seqid in alignment
                 if (seqid.split('|')[0] == "Reference" and seqid.split('|')[4] == "gene")}
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
        # print("domain anns:", '\n'.join(x.generate_str() for x in domain_anns))
        # print("feature_ranges:", feature_ranges)
        def adjust_feature_ranges(gene, seqid, **kwargs):
            source = self.get_ref_seqid(seqid, attr = "source")
            gene_ann = gene_anns[source]
            gene_feature_ranges = feature_ranges[source]
            if self.pssm_ids and self.gff_domain:
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
        # check = [0]
        def is_rvs(aln_seqid, seqid):
            filled_rvs_template = Template(alignment_rvs_pattern).substitute(seqid = re.escape(seqid))
            # if check[0] == 0 and '7273.tig' in aln_seqid:
            #     print(aln_seqid, seqid, filled_rvs_template, re.escape(aln_seqid), re.search(filled_rvs_template, re.escape(aln_seqid)), re.search(filled_rvs_template, aln_seqid))
            #     check[0] = 1
            return re.search(filled_rvs_template, aln_seqid)
        feature_only_ranges = {seqid: adjust_feature_ranges(gene, seqid, subtract_gaps = True)
                               for gene, seqid in genes.items()}
        # print("feature_only_ranges:", feature_only_ranges)
        feature_gaps_ranges = {seqid: ranges_subtract(adjust_feature_ranges(gene, seqid, subtract_gaps = False),
                                                      feature_only_ranges[seqid])
                               for gene, seqid in genes.items()}
        # print("feature_gaps_ranges:", feature_gaps_ranges)
        ## define acceptable ranges in targets
        get_target_feature_ranges = make_target_feature_ranges_function(feature_only_ranges,
                                                                        feature_gaps_ranges,
                                                                        max_insertion = max_insertion)
        targets_feature_ranges = {seqid: get_target_feature_ranges(alignment[seqid], seqid)
                                  for seqid in alignment}
        # print("targets_feature_ranges:", targets_feature_ranges)
        # self.logfile.devsplain(f"{targets_feature_ranges}, {len(self.grna_hits.hits.items())}")
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
        ## update .map and FASTA files
        if self.auto_update_files:
            self.write_all_grna_map()
            self.write_pass_grna_files()
        return
    
    def filter_exclude(self):
        """
        Set exclude check for candidate gRNAs. Overwrites existing exclude check using self.exclude.
        
        If self.exclude is set, gRNA which sequences appear in the file at self.exclude will fail this check.
        """
        if self.exclude is not None:
            to_exclude = set(str(x).upper() for x in fasta_to_dict(self.exclude).values())
            self.grna_hits.set_all_seqs_check_by_function("exclude",
                                                          lambda grna_seq: str(grna_seq) not in to_exclude)
            ## update .map and FASTA files
            if self.auto_update_files:
                self.write_all_grna_map()
                self.write_pass_grna_files()
        return
    
    def filter(self, background_check = True, feature_check = True, gc_check = True):
        """
        Execute self.filter_background, self.filter_feature, and self.filter_gc based on self's attributes.
        
        Arguments:
            background_check (bool): filter gRNA by off-target (default=True)
            feature_check (bool): filter gRNA by within-feature (default=True)
            gc_check (bool): filter gRNA by GC content (default=True)
        """
        ## lower auto update flag so we don't re-write the files three times
        auto_update = self.auto_update_files
        self.auto_update_files = False
        if background_check:
            self.logfile.devsplain("Filtering background")
            self.filter_background()
        if feature_check and self.genes:
            self.logfile.devsplain("Filtering feature")
            self.filter_feature()
        if gc_check:
            self.logfile.devsplain("Filtering GC")
            self.filter_gc()
        ## return flag to original status and update files if flag is raised
        self.auto_update_files = auto_update
        ## update .map and FASTA files
        if self.auto_update_files:
            self.write_all_grna_map()
            self.write_pass_grna_files()
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
    
    def minimumset(self, sets = None, manual = None, fasta = None,
                   fout_fasta = None, fout_map = None, fout_eqv = None,
                   report_full_path = True, all_checks = True, nonstandard_checks = [],
                   exclude_check = True, gc_check = True, background_check = True, feature_check = True):
        """
        Generate minimum set(s) of gRNA required to cover all targets.
        
        Arguments:
            sets (int): optional. See attributes.
            manual (bool): manually approve all gRNA sets. Defaults to NOT self.auto if not provided.
            fasta (str): optional, path to FASTA file. Used for renaming gRNA.
            fout_fasta (str): optional, path to output FASTA file.
                Autogenerated using self.prefix and self.active_directory if not provided.
            fout_map (str): optional, path to output .map file.
                Autogenerated using self.prefix and self.active_directory if not provided.
            fout_eqv (str): optional, path to output .eqv file of passing gRNA.
                Autogenerated using self.prefix and self.active_directory if not provided.
            report_full_path (bool): print full path to output files upon successful writing
            all_checks (bool): include all checks for filtering (including any in user-added columns)
            nonstandard_checks (list): non-standard checks to include for filtering (subject to self.accept_invalid_field)
            exclude_check (bool): include 'exclude' field for filtering (subject to self.accept_invalid_field)
            gc_check (bool): include 'GC' field for filtering (subject to self.accept_invalid_field)
            background_check (bool): include 'background' field for filtering
                (subject to self.accept_invalid_field)
            feature_check (bool): include 'feature' field for filtering (subject to self.accept_invalid_field)
        """
        if not self.grna_hits:
            raise MINORgError( ("MINORg.minimumset requires gRNA."
                                " You may read a mapping file using"
                                " self.parse_grna_map_from_file('<path to file>')"
                                " or use self.grna() to generate gRNA.") )
        if sets is None: sets = self.sets
        ## assume all targets in self.grna_hits are to be, well, targeted
        targets = set(hit.target_id for hit in self.grna_hits.flatten_hits())
        ## set exclude check
        self.filter_exclude()
        ## check if statuses has been set. If not, warn user.
        # if all_checks or exclude_check:
        #     self._check_valid_status("exclude", "seq", "exclude")
        if all_checks or gc_check:
            self._check_valid_status("GC", "seq", "%GC")
        if all_checks or background_check:
            self._check_valid_status("background", "seq", "off-target")
        if all_checks or feature_check:
            self._check_valid_status("feature", "hit", "within-feature")
        if all_checks:
            for check_name in self.check_names(hit = False, seq = True,
                                               standard = False, nonstandard = True):
                self._check_valid_status(check_name, "seq", check_name)
            for check_name in self.check_names(hit = True, seq = False,
                                               standard = False, nonstandard = True):
                self._check_valid_status(check_name, "hit", check_name)
        ## filter grna by checks
        grna_hits = copy.deepcopy(self.valid_grna(*([] if all_checks else
                                                    ((["exclude"] if exclude_check else []) +
                                                     (["GC"] if gc_check else []) +
                                                     (["background"] if background_check else []) +
                                                     (["feature"] if feature_check else []) +
                                                     nonstandard_checks))))
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
        ## start generating sets
        grna_sets = []
        grna_hits_for_setcover = copy.deepcopy(grna_hits)
        ## for tie breaker functions (only actually used when !self.prioritise_nr):
        ## - cov: dictionary of {'<gRNA seq>': [<gRNAHit items>]} for unselected gRNA
        ## - all_cov: dictionary of {'<gRNA seq>': [<gRNAHit items>]} for all gRNA
        ## - covered: set of {<IDs of targets covered by already chosen gRNA>}
        if not self.prioritise_nr: ## tie break with pos (favour 5')
            tie_breaker = lambda cov, all_cov, covered: tuple(all_best_nr(all_best_pos(cov, all_cov, covered),
                                                                          all_cov, covered).items())[0]
        else: ## tie break with non-redundancy in coverage
            # tie_breaker = lambda cov, all_cov, covered: tuple(all_best_pos(all_best_nr(cov, all_cov, covered),
            #                                                                all_cov, covered).items())[0]
            tie_breaker = lambda *args: tuple()
        get_minimum_set = make_get_minimum_set(
            grna_hits_for_setcover, prioritise_nr = self.prioritise_nr,
            manual_check = (manual if manual is not None else not self.auto),
            targets = targets, tie_breaker = tie_breaker, suppress_warning = True)
        while len(grna_sets) < sets:
            ## get a (minimum) set of gRNA sequences
            seq_set = get_minimum_set()
            ## if valid set returned
            if seq_set:
                grna_sets.append(seq_set) ## add to existing list of sets
            else:
                # warnings.warn(f"The gRNA sequences cannot cover all target sequences the desired number of times ({sets}). (Failed at set {len(grna_sets) + 1} of {sets})\n")
                self.logfile.warning(f"The gRNA sequences cannot cover all target sequences the desired number of times ({sets}). (Failed at set {len(grna_sets) + 1} of {sets})\n")
                break
        ## write
        if not fout_fasta: fout_fasta = get_val_default(self.final_fasta,
                                                        self.results_fname("gRNA_final.fasta"))
        if not fout_map: fout_map = get_val_default(self.final_map,
                                                    self.results_fname("gRNA_final.map"))
        if not fout_eqv: fout_eqv = get_val_default(self.pass_eqv,
                                                    self.results_fname("gRNA_pass.eqv"))
        ## write gRNA fasta file and gRNA mapping
        if grna_sets:
            seqs_to_write = itertools.chain(*grna_sets)
            grna_hits.write_fasta(fout_fasta, seqs = seqs_to_write, fasta = fasta)
            print("\nFinal gRNA sequence(s) have been written to"
                  f" {fout_fasta if report_full_path else os.path.basename(fout_fasta)}")
            self.final_fasta = fout_fasta
            grna_hits.write_mapping(fout_map, sets = grna_sets, fasta = fasta, version = 5)
            print("Final gRNA sequence ID(s), gRNA sequence(s), and target(s) have been written to"
                  f" {fout_map if report_full_path else os.path.basename(fout_map)}")
            self.final_map = fout_map
            grna_hits_for_setcover.write_equivalents(fout_eqv, fasta = fasta, write_all = True)
            print("Groupings of passing gRNA with equivalent coverage have been written to"
                  f" {fout_eqv if report_full_path else os.path.basename(fout_eqv)}")
        ## print summary
        print(f"\n{sets} mutually exclusive gRNA set(s) requested. {len(grna_sets)} set(s) found.")
        return
    
    ####################
    ##  SUBCMD: FULL  ##
    ####################
    
    def full(self, manual = None, background_check = True, feature_check = True, gc_check = True):
        """
        Execute full MINORg programme using self's attributes as parameter values.
        
        Arguments:
            manual (bool): manually approve all gRNA sets. Defaults to NOT self.auto if not provided.
            background_check (bool): filter gRNA by off-target (default=True)
            feature_check (bool): filter gRNA by within-feature (default=True)
            gc_check (bool): filter gRNA by GC content (default=True)
        """
        with warnings.catch_warnings():
            warnings.filterwarnings("error", message = "No target sequences found.",
                                    category = MINORgWarning)
            print("query reference:", self.query_reference)
            if not self.target or self._genes_updated_since_alignment:
                try:
                    self.seq(quiet = False)
                except MINORgWarning:
                    warnings.warn("No targets found.", MINORgWarning)
                    return
        self.grna()
        if len(self.grna_hits) == 0:
            warnings.warn("No gRNA found.", MINORgWarning)
            return
        ## filter
        self.filter(background_check = background_check,
                    feature_check = feature_check,
                    gc_check = gc_check)
        ## write file only if not already automatically updated
        if not self.auto_update_files:
            self.write_all_grna_map()
            self.write_pass_grna_files()
        ## abort if no gRNA pass
        if len(self.valid_grna()) == 0:
            warnings.warn("No valid gRNA after filtering.", MINORgWarning)
            return
        ## minimumset
        self.minimumset(report_full_path = False, manual = manual)
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
