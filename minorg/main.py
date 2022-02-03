#!/usr/bin/env python3

## TODO: allow users to provide GFF instead of GFF -> BED (done)
## TODO: 2021/06/10: implement --valid-indv/--valid-genome/--valid-acc, --valid-cluster to allow user to check alias validity (print location of FASTA file for --valid-genome, do member_callback for --valid-cluster) before running
## TODO: 2021/06/10: disable Enum for --indv. Since users are now allowed to submit custom genome lookup files using --genome-lookup, autofill shouldn't be used for --indv unless it can dynamically update to accommodate user's --genome-lookup input even before command submission (which is impossible lol). Basically, restructure the whole --indv thing to use the same kind of code as --cluster
## TODO: 2021/09/07: don't merge ext_gene files w/ reference (or ext_bed files w/ bed). Find some way to allow both to be processed whenever blastn or whatever needs to be done to reference
## TODO: 2021/10/21: update all functions to use config.reference_ext and config.annotation_ext (dictionary of {source_id: path_to_fasta/bed}) instead of args.reference in order to accommodate multiple reference files (since we don't want to merge ext_gene with the reference). I've modified some stuff and test with homologue but not the other subcommands.
## TODO: 2021/10/21: log user-provided indvs that cannot be mapped. Actually, log user ANY ALIAS that cannot be mapped.
## TODO: 2021/10/21: update everything to use gff3 instead of bed
## TODO: 2021/10/21: enable remote rpsblast (not that i've been able to get it to find the remote database but hey we should at least let people have the option)
## TODO: 2021/10/22: DONE. for the generate_grna subcommand, either integrate within_cds (rename that to within_feature) with execute_grna or update execute_filter so it can filter by within feature and not everything else then use execute_filter after execute_grna.
## TODO: 2021/10/26: enable other gRNA filters. Also, implement within feature filter for generate_gRNA.
## TODO: 2021/10/26: clean up execute_grna so it pipes into a more generic function that doesn't require config/args
## TODO: 2021/10/31: use external txt file to track reference sequences ("^Reference|...") so changing naming format or having '|' in seqid is less likely to break the code. Easier maintenance. (columns can be seqid, gene, isoform, source, molecule, comma-separated ranges)
## TODO: 2021/11/08: enable all homologue discovery options for filtering in background as well
## TODO: 2021/11/09: ensure all arguments required for --indv are also enabled for --ot-indv (extract the query mapping part of check_homologue_args into separate function reusable for --ot-indv)
## TODO: 2021/11/18: create Minorg object to facilitate passing of args (and storing parsed args (e.g. Minrog.query_map stores query filepaths mapped to aliases so we don't have to attach it to config)).
## TODO: 2021/11/19: for filter, accept masked_report.tsv so we don't have to re-mask everything?


import os
import sys
import click
import typer
import inspect
import logging
import itertools

from typing import List, Optional
from pathlib import Path
from argparse import Namespace
# from datetime import datetime

## just cheat it a little so the minorg.X imports work without this actually being installed
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from minorg.functions import (
    split_list as split_callback_list#,
    # get_count_dict, cat_files, assign_alias
)
# from minorg.grna import gRNAHits
# from minorg.annotation import reduce_ann
# from minorg.pam import PAM
from minorg.log import MINORgLogger
from minorg.MINORgCLI import MINORgCLI

# from minorg.extend_reference import extend_reference_cli

from minorg.parse_config import (
    # Config,
    Params, OptionParams,
    # SetCoverAlgo, IndvGenomesAll, IndvGenomesAllClear, Lookup,
    # LogFile,
    # parse_multiline_multikey_sdict
)

from minorg.exceptions import (
    MINORgError,
    InputFormatError,
    InvalidPath,
    InvalidFile,
    UnreadableFile
)

# from minorg.display import (
#     print_indent as printi,
#     print_overwrite_multi as printom
# )

default_sub_cmd = "full"
app_main = typer.Typer()
app_sub = typer.Typer()

logging_level = logging.DEBUG

## parse config
try:
    CONFIG = os.environ["MINORG_CONFIG"]
except KeyError:
    print("key error!")
    # ## TODO: set env variable so we can remove this KeyError exception handling thing
    # CONFIG = "/mnt/chaelab/rachelle/scripts/minorgpy/config.ini"
    CONFIG = None

## generate Params, OParmas, and Config objects
oparams = OptionParams() ## namespace for some pre-defined sets of parameter options
# params = Params(CONFIG)
# config = Config(params, keep_on_crash = True)
minor_g = MINORgCLI(config = CONFIG, keep_on_crash = False)
params = minor_g.params

## argument parsing functions
def positive_callback(val):
    """
    Callback that checks if float or integer value is greater than zero.
    
    Raises
    ------
    typer.BadParameter
        If value is not positive
    """
    if val > 0:
        return val
    else:
        raise typer.BadParameter( f"Invalid value: {val}. Positive value required." )

def non_negative_callback(val):
    """
    Callback that checks if float or integer value is greater than or equal to zero.
    
    Raises
    ------
    typer.BadParameter
        If value is negative
    """
    if val >= 0:
        return val
    else:
        raise typer.BadParameter( f"Invalid value: {val}. Non-negative value required." )

def zero_to_one_callback(val):
    """
    Callback that checks if float or integer value is between 0 and 1 (inclusive).
    
    Raises
    ------
    typer.BadParameter
        If value is less than zero or greater than 1
    """
    val = non_negative_callback(val)
    if val <= 1:
        return val
    else:
        raise typer.BadParameter( f"Invalid value: {val}. Must be between 0 and 1 (inclusive)." )

# ## TODO: make sure all tmp files are removed if not using --keep-on-crash (2021/05/12)
# @app_sub.command("homologue")
# @app_sub.command("homolog")
@app_sub.command("target")
@app_sub.command("seq")
def seq(
        ## general options
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        blastn: str = typer.Option(*params.blastn(), **params.blastn.options),
        rpsblast: str = typer.Option(*params.rpsblast(), **params.rpsblast.options),
        mafft: str = typer.Option(*params.mafft(), **params.mafft.options),
        thread: int = typer.Option(*params.thread(), **params.thread.options),
        remote_rps: bool = typer.Option(*params.remote_rps(), **params.remote_rps.options),
        
        ## target definition options
        gene: Optional[List[str]] = typer.Option(*params.gene(), **params.gene.options,
                                                 callback = split_callback_list),
        cluster: Optional[List[str]] = typer.Option(*params.cluster(), **params.cluster.options,
                                                    callback = split_callback_list),
        indv: Optional[List[str]] = typer.Option(*params.indv(), **params.indv.options,
                                                 callback = split_callback_list),
        target: Path = typer.Option(*params.target(), **params.target.options,
                                    **oparams.file_valid),
        query: Optional[List[Path]] = typer.Option(*params.query(), **params.query.options,
                                                   **oparams.file_valid),
        domain: str = typer.Option(*params.domain(), **params.domain.options,
                                   callback = minor_g.domain_callback),
        feature: Optional[List[str]] = typer.Option(*params.feature(), **params.feature.options,
                                                    callback = split_callback_list),
        minid: float = typer.Option(*params.minid(), **params.minid.options, callback = non_negative_callback),
        minlen: int = typer.Option(*params.minlen(), **params.minlen.options, callback = positive_callback),
        mincdslen: int = typer.Option(*params.mincdslen(), **params.mincdslen.options,
                                      callback = positive_callback),
        check_recip: bool = typer.Option(*params.check_recip(), **params.check_recip.options),
        relax_recip: bool = typer.Option(*params.relax_recip(), **params.relax_recip.options),
        merge_within: int = typer.Option(*params.merge_within(), **params.merge_within.options,
                                         callback = non_negative_callback),
        check_id_before_merge: bool = typer.Option(*params.check_id_before_merge(),
                                                   **params.check_id_before_merge.options),
        cluster_set: str = typer.Option(*params.cluster_set(), **params.cluster_set.options,
                                        **oparams.file_valid, is_eager = True,
                                        callback = minor_g.cluster_set_callback),
        genome_set: str = typer.Option(*params.genome_set(), **params.genome_set.options,
                                       **oparams.file_valid, is_eager = True,
                                       callback = minor_g.genome_set_callback),
        reference_set: str = typer.Option(*params.reference_set(), **params.reference_set.options,
                                          **oparams.file_valid, is_eager = True,
                                          callback = minor_g.reference_set_callback),
        reference: Optional[List[str]] = typer.Option(*params.reference(), **params.reference.options,
                                                      callback = split_callback_list),
        assembly: Path = typer.Option(*params.assembly(), **params.assembly.options,
                                      **oparams.file_valid),
        annotation: Path = typer.Option(*params.annotation(), **params.annotation.options,
                                        **oparams.file_valid),
        db: str = typer.Option(*params.db(), **params.db.options, callback = minor_g.db_callback),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = minor_g.attr_mod_callback),
        ext_gene: Optional[List[Path]] = typer.Option(*params.ext_gene(), **params.ext_gene.options,
                                                      **oparams.file_valid),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options,
                                                     **oparams.file_valid),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        genetic_code: str = typer.Option(*params.genetic_code(), **params.genetic_code.options,
                                         callback = minor_g.genetic_code_callback),
        
        ## user lookups (this is right at the end because we might need the other args)
        ## genomes is not "eager" because we have to wait for --genome-lookup to be processed
        genomes: bool = typer.Option(*params.genomes(), **params.genomes.options,
                                     callback = minor_g.make_genomes_callback("independent")),
        references: bool = typer.Option(*params.references(), **params.reference.options,
                                        callback = minor_g.references_callback),
        ## clusters & members are not "eager" because we have to wait for --cluster-lookup to be processed
        clusters: bool = typer.Option(*params.clusters(), **params.clusters.options,
                                      callback = minor_g.clusters_callback),
        members: str = typer.Option(*params.members(), **params.members.options,
                                    callback = minor_g.members_callback) ):
    """
    Subcommand seq.
    
    Generate target sequences. 
    If query/queries is/are non-reference genomes, conducts homologue discovery and outputs homologues.
    Else if querying reference genome, reference genes will be output.
    """
    
    ## check validity of args
    args = Namespace(**locals())
    minor_g.parse_args(args, "seq")
    
    ## call execute_homologue.
    minor_g.subcmd_seq()
    minor_g.resolve()
    return

@app_sub.command("grna")
def generate_grna(
        ## general options
        rpsblast: str = typer.Option(*params.rpsblast(), **params.rpsblast.options),
        remote_rps: bool = typer.Option(*params.remote_rps(), **params.remote_rps.options),
        mafft: str = typer.Option(*params.mafft(), **params.mafft.options),
        thread: int = typer.Option(*params.thread(), **params.thread.options),
        
        ## output files: option A
        out_map: Path = typer.Option(*params.out_map(), resolve_path = True),
        out_fasta: Path = typer.Option(*params.out_fasta(), resolve_path = True),
        
        ## output files: option B
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        
        ## gRNA options
        pam: str = typer.Option(*params.pam(), **params.pam.options),
        length: str = typer.Option(*params.length(), **params.length.options),
        span_junction: bool = typer.Option(*params.span_junction(), **params.span_junction.options),
        
        ## target definition options
        target: Path = typer.Option(*params.target(), **params.target.options,
                                    **oparams.file_valid),
        
        ## target definition options if finding gRNA in reference genes
        gene: Optional[List[str]] = typer.Option(*params.gene(), **params.gene.options,
                                                 callback = split_callback_list),
        cluster: Optional[List[str]] = typer.Option(*params.cluster(), **params.cluster.options,
                                                    callback = split_callback_list),
        domain: str = typer.Option(*params.domain(), **params.domain.options,
                                   callback = minor_g.domain_callback),
        reference: Optional[List[str]] = typer.Option(*params.reference(), **params.reference.options,
                                                      callback = split_callback_list),
        assembly: Path = typer.Option(*params.assembly(), **params.assembly.options,
                                      **oparams.file_valid),
        annotation: Path = typer.Option(*params.annotation(), **params.annotation.options,
                                        **oparams.file_valid),
        db: str = typer.Option(*params.db(), **params.db.options, callback = minor_g.db_callback),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = minor_g.attr_mod_callback),
        ext_gene: Optional[List[Path]] = typer.Option(*params.ext_gene(), **params.ext_gene.options,
                                                        **oparams.file_valid),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options,
                                                     **oparams.file_valid),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        genetic_code: str = typer.Option(*params.genetic_code(), **params.genetic_code.options,
                                         callback = minor_g.genetic_code_callback),
        
        ## basic filtering
        gc_min: float = typer.Option(*params.gc_min(), **params.gc_min.options,
                                     callback = zero_to_one_callback),
        gc_max: float = typer.Option(*params.gc_max(), **params.gc_max.options,
                                     callback = zero_to_one_callback),
        feature: Optional[List[str]] = typer.Option(*params.feature(), **params.feature.options,
                                                    callback = split_callback_list),
        
        ## user lookups (this is right at the end because we might need the other args)
        ## genomes is not "eager" because we have to wait for --genome-lookup to be processed
        references: bool = typer.Option(*params.references(), **params.reference.options,
                                        callback = minor_g.references_callback),
        ## clusters & members are not "eager" because we have to wait for --cluster-lookup to be processed
        clusters: bool = typer.Option(*params.clusters(), **params.clusters.options,
                                      callback = minor_g.clusters_callback),
        members: str = typer.Option(*params.members(), **params.members.options,
                                    callback = minor_g.members_callback),
        ## more user lookups, this time for the lookups themselves haha
        cluster_set: str = typer.Option(*params.cluster_set(), **params.cluster_set.options,
                                        **oparams.file_valid, is_eager = True,
                                        callback = minor_g.cluster_set_callback),
        reference_set: str = typer.Option(*params.reference_set(), **params.reference_set.options,
                                          **oparams.file_valid, is_eager = True,
                                          callback = minor_g.reference_set_callback)):
    
    """
    Subcommand grna.
    
    Generate gRNA from either user-provided FASTA file or reference genes
    """
    
    typer.echo("generating grna")
    
    ## check validity of args
    args = Namespace(**locals())
    minor_g.parse_args(args, "grna")
    minor_g.logfile.devsplain(f"post check: {str(vars(args))}")
    
    ## generate gRNA
    minor_g.subcmd_grna()
    minor_g.resolve()
    return

@app_sub.command("filter")
@app_sub.command("check")
def filter_grna(
        ## which filters
        check_all: bool = typer.Option(*params.check_all(), **params.check_all.options),
        gc_check: bool = typer.Option(*params.gc_check(), **params.gc_check.options),
        feature_check: bool = typer.Option(*params.feature_check(), **params.feature_check.options),
        background_check: bool = typer.Option(*params.background_check(), **params.background_check.options),
        reset_checks: bool = typer.Option(*params.reset_checks(), **params.reset_checks.options),
        
        ## general options
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        mafft: str = typer.Option(*params.mafft(), **params.mafft.options),
        thread: int = typer.Option(*params.thread(), **params.thread.options),
        blastn: str = typer.Option(*params.blastn(), **params.blastn.options),
        rpsblast: str = typer.Option(*params.rpsblast(), **params.rpsblast.options),
        remote_rps: bool = typer.Option(*params.remote_rps(), **params.remote_rps.options),
        db: str = typer.Option(*params.db(), **params.db.options, callback = minor_g.db_callback),
        genetic_code: str = typer.Option(*params.genetic_code(), **params.genetic_code.options,
                                         callback = minor_g.genetic_code_callback),
        
        ## input files
        map: Path = typer.Option(*params.map(), **params.map.options, **oparams.file_valid),
        grna: Path = typer.Option(*params.grna(), **params.grna.options, **oparams.file_valid),
        alignment: Path = typer.Option(*params.alignment(), **params.alignment.options, **oparams.file_valid),
        target: Path = typer.Option(*params.target(), **params.target.options, **oparams.file_valid),
        reference: Optional[List[str]] = typer.Option(*params.reference(), **params.reference.options,
                                                      callback = split_callback_list),
        assembly: Path = typer.Option(*params.assembly(), **params.assembly.options,
                                      **oparams.file_valid),
        annotation: Path = typer.Option(*params.annotation(), **params.annotation.options,
                                        **oparams.file_valid),
        ext_gene: Optional[List[Path]] = typer.Option(*params.ext_gene(), **params.ext_gene.options,
                                                      **oparams.file_valid),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options,
                                                     **oparams.file_valid),
        background: Optional[List[Path]] = typer.Option(*params.background(), **params.background.options,
                                                        **oparams.file_valid),
        mask: Optional[List[Path]] = typer.Option(*params.mask(), **params.mask.options, **oparams.file_valid),
        
        ## format options
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = minor_g.attr_mod_callback),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        
        ## GC filter options
        gc_min: float = typer.Option(*params.gc_min(), **params.gc_min.options,
                                     callback = zero_to_one_callback),
        gc_max: float = typer.Option(*params.gc_max(), **params.gc_max.options,
                                     callback = zero_to_one_callback),
        
        ## within feature filter options
        feature: Optional[List[str]] = typer.Option(*params.feature(), **params.feature.options,
                                                    callback = split_callback_list),
        max_insertion: int = typer.Option(*params.max_insertion(), **params.max_insertion.options,
                                          callback = non_negative_callback),
        min_within_n: int = typer.Option(*params.min_within_n(), **params.min_within_n.options,
                                         callback = positive_callback),
        min_within_fraction: float = typer.Option(*params.min_within_fraction(),
                                                  **params.min_within_fraction.options,
                                                  callback = zero_to_one_callback),
        
        ## output options
        out_map: Path = typer.Option(*params.out_map(), resolve_path = True),
        out_fasta: Path = typer.Option(*params.out_fasta(), resolve_path = True),
        in_place: bool = typer.Option(*params.in_place(), **params.in_place.options),
        
        ## background filter options
        screen_reference: bool = typer.Option(*params.screen_reference(), **params.screen_reference.options),
        unmask_ref: bool = typer.Option(*params.unmask_ref(), **params.unmask_ref.options),
        ot_indv: Optional[List[str]] = typer.Option(*params.ot_indv(), **params.ot_indv.options,
                                                    callback = split_callback_list),
        # by_indv: bool = typer.Option(*params.by_indv(), **params.by_indv.options),
        ## masking options
        mask_gene: Optional[List[str]] = typer.Option(*params.mask_gene(), **params.mask_gene.options,
                                                      callback = split_callback_list),
        unmask_gene: Optional[List[str]] = typer.Option(*params.unmask_gene(), **params.unmask_gene.options,
                                                        callback = split_callback_list),
        mask_cluster: Optional[List[str]] = typer.Option(*params.mask_cluster(),
                                                         **params.mask_cluster.options,
                                                         callback = split_callback_list),
        unmask_cluster: Optional[List[str]] = typer.Option(*params.unmask_cluster(),
                                                           **params.unmask_cluster.options,
                                                           callback = split_callback_list),
        domain: str = typer.Option(*params.domain(), **params.domain.options,
                                   callback = minor_g.domain_callback),
        # mask_homologue: bool = typer.Option(*params.mask_homologue(), **params.mask_homologue.options),
        # minid: float = typer.Option(*params.minid(), **params.minid.options,
        #                             callback = non_negative_callback),
        # minlen: int = typer.Option(*params.minlen(), **params.minlen.options, callback = positive_callback),
        # mincdslen: int = typer.Option(*params.mincdslen(), **params.mincdslen.options,
        #                               callback = positive_callback),
        # check_recip: bool = typer.Option(*params.check_recip(), **params.check_recip.options),
        # relax_recip: bool = typer.Option(*params.relax_recip(), **params.relax_recip.options),
        # merge_within: int = typer.Option(*params.merge_within(), **params.merge_within.options,
        #                                  callback = non_negative_callback),
        # check_id_before_merge: bool = typer.Option(*params.check_id_before_merge(),
        #                                        **params.check_id_before_merge.options),
        ## off-target threshold options
        pam: str = typer.Option(*params.pam(), **params.pam.options),
        ot_pamless: bool = typer.Option(*params.ot_pamless(), **params.ot_pamless.options),
        ot_mismatch: int = typer.Option(*params.ot_mismatch(), **params.ot_mismatch.options,
                                        callback = non_negative_callback),
        ot_gap: int = typer.Option(*params.ot_gap(), **params.ot_gap.options, callback = non_negative_callback),
        
        ## exclude filter
        exclude: Path = typer.Option(*params.exclude(), **params.exclude.options, **oparams.file_valid),
        
        ## user lookups (this is right at the end because we might need the other args)
        ## genomes is not "eager" because we have to wait for --genome-lookup to be processed
        genomes: bool = typer.Option(*params.genomes(), **params.genomes.options,
                                     callback = minor_g.make_genomes_callback("independent")),
        references: bool = typer.Option(*params.references(), **params.reference.options,
                                        callback = minor_g.references_callback),
        ## clusters & members are not "eager" because we have to wait for --cluster-lookup to be processed
        clusters: bool = typer.Option(*params.clusters(), **params.clusters.options,
                                      callback = minor_g.clusters_callback),
        members: str = typer.Option(*params.members(), **params.members.options,
                                    callback = minor_g.members_callback),
        ## more user lookups, this time for the lookups themselves haha
        cluster_set: str = typer.Option(*params.cluster_set(), **params.cluster_set.options,
                                        **oparams.file_valid, is_eager = True,
                                        callback = minor_g.cluster_set_callback),
        genome_set: str = typer.Option(*params.genome_set(), **params.genome_set.options,
                                       **oparams.file_valid, is_eager = True,
                                       callback = minor_g.genome_set_callback),
        reference_set: str = typer.Option(*params.reference_set(), **params.reference_set.options,
                                          **oparams.file_valid, is_eager = True,
                                          callback = minor_g.reference_set_callback)):
    
    """
    Subcommand filter.
    
    Filter gRNA by checks (background, feature, and GC)
    """
    
    typer.echo("filtering")
    
    ## check validity of args
    args = Namespace(**locals())
    minor_g.parse_args(args, "filter")
    minor_g.logfile.devsplain(f"post check: {str(vars(args))}")
    minor_g.subcmd_filter()
    minor_g.resolve()
    # ## check arguments are in order
    # check_reference_args(args)
    # check_filter_args(args, standalone = True)
    # ## parse [mask|unmask]_cluster -> genes
    # gene_sets = parse_genes_for_filter(args)
    # ## move log file
    # config.logfile.move(config.directory, args.prefix)
    # config.logfile.args_expanded([args, params])
    
    # # print("post check:", vars(args))
    # config.logfile.devsplain(f"post check: {str(vars(args))}")
    
    # if prefix: config.prefix = prefix
    # if directory: config.out_dir = directory
    
    # ## extend genome if ext_gene and ext_cds are provided
    # extend_reference_cli(args, config)
    
    # ## filter BED/GFF for relevant entries (reduces get_seq search time)
    # if tuple(itertools.chain(*gene_sets.values())):
    #     make_reduced_ann(args, config, gene_sets)
    
    # ## get sequences to be masked
    # to_mask = assign_alias(args.mask, lambda i: f"Usermask_{str(i+1).zfill(2)}")
    # for prefix, genes in gene_sets.items():
    #     if not genes: continue
    #     output_homologue = execute_homologue(args, config, params, prefix, genes,
    #                                          indv = {"ref"}, for_masking = True)
    #     to_mask[prefix] = output_homologue[2]
    #     # if screen_ref and not unmask_ref:
    #     #     to_mask[f"{prefix}_ref"] = output_homologue[4]
    #     # if mask_homologue:
    #     #     to_mask[f"{prefix}_nonref"] = output_homologue[2]
    
    # ## parse mapping file into gRNAHits object
    # grna_hits = gRNAHits()
    # grna_hits.parse_from_mapping(args.mapping, targets = args.target, version = None)
    # ## rename gRNA according to FASTA file if provided
    # if args.grna is not None: grna_hits.assign_gRNAseq_id(args.grna)
    # ## reset checks if requested
    # if args.reset_checks: grna_hits.clear_checks()
    
    # ## call execute_filter
    # # val_or_empty = lambda b, v: [v] if b else []
    # # checks = list(itertools.chain(*[val_or_empty(b, v)
    # #                                 for b, v in [(args.gc_check, "GC"),
    # #                                              (args.feature_check, "feature"),
    # #                                              (args.background_check, "background")]]))
    # grna_pass, map_pass = execute_filter(args, config, directory, prefix, grna_hits,
    #                                      args.mapping, args.grna, args.alignment, args.target,
    #                                      config.annotation_red, fasta_exclude = args.exclude,
    #                                      # checks = checks,
    #                                      domain_gff_bed = None, to_mask = to_mask)
    # config.resolve()
    return

@app_sub.command("minimumset")
def minimumset(
        
        ## data input
        map: Path = typer.Option(*params.map(), **params.map.options, **oparams.file_valid),
        grna: Path = typer.Option(*params.grna(), *params.fasta.names, **oparams.file_valid),
        exclude: Path = typer.Option(*params.exclude(), **oparams.file_valid),
        
        ## output files: option A
        out_map: Path = typer.Option(*params.out_map(), resolve_path = True),
        out_fasta: Path = typer.Option(*params.out_fasta(), resolve_path = True),
        
        ## output files: option B
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        
        ## process options
        sets: int = typer.Option(*params.sets(), callback = positive_callback),
        # sc_algorithm: SetCoverAlgo = typer.Option(*params.sc_algorithm()),
        
        ## flags
        auto: bool = typer.Option(*params.auto(), **params.auto.options),
        accept_invalid: bool = typer.Option(*params.accept_invalid(), **params.accept_invalid.options),
        accept_feature_unknown: bool = typer.Option(*params.accept_feature_unknown(),
                                                    **params.accept_feature_unknown.options) ):
    
    """
    Subcommand minimumset.
    
    Generate minimum set(s) of gRNA required to cover all targets from mapping file and FASTA file of gRNA.
    Requires mapping file (generated by minorg's 'full' subcommand) and a FASTA file of gRNA sequences.
    gRNA sequences not present in the mapping file will be ignored.
    """
    args = Namespace(**locals())
    minor_g.parse_args(args, "minimumset")
    minor_g.subcmd_minimumset()
    minor_g.resolve()
    return

## TODO: wrapper for the minorg.py 'full' function that parses --acc in a way that works by subsetting VdW's FASTA (as is currently the implementation for find_gRNA.sh), then passing that into this 'full' function?? (But then how to prevent this from writing a log that says '-q <file generated by the wrapping function>'?
@app_sub.command("full")
def full(
        
        ## general options
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        blastn: str = typer.Option(*params.blastn(), **params.blastn.options),
        rpsblast: str = typer.Option(*params.rpsblast(), **params.rpsblast.options),
        # rpsblast: Path = typer.Option(*params.rpsblast(), **params.rpsblast.options, **oparams.file_valid),
        mafft: str = typer.Option(*params.mafft(), **params.mafft.options),
        thread: int = typer.Option(*params.thread(), **params.thread.options),
        remote_rps: bool = typer.Option(*params.remote_rps(), **params.remote_rps.options),
        
        ## target definition options
        gene: Optional[List[str]] = typer.Option(*params.gene(), **params.gene.options,
                                                 callback = split_callback_list),
        cluster: Optional[List[str]] = typer.Option(*params.cluster(), **params.cluster.options,
                                                    callback = split_callback_list),
        indv: Optional[List[str]] = typer.Option(*params.indv(), **params.indv.options,
                                                 callback = split_callback_list),
        target: Path = typer.Option(*params.target(), **params.target.options),
        query: Optional[List[Path]] = typer.Option(*params.query(), **params.query.options,
                                                   **oparams.file_valid),
        domain: str = typer.Option(*params.domain(), **params.domain.options,
                                   callback = minor_g.domain_callback),
        minid: float = typer.Option(*params.minid(), **params.minid.options, callback = non_negative_callback),
        minlen: int = typer.Option(*params.minlen(), **params.minlen.options, callback = positive_callback),
        mincdslen: int = typer.Option(*params.mincdslen(), **params.mincdslen.options,
                                      callback = positive_callback),
        check_recip: bool = typer.Option(*params.check_recip(), **params.check_recip.options),
        relax_recip: bool = typer.Option(*params.relax_recip(), **params.relax_recip.options),
        merge_within: int = typer.Option(*params.merge_within(), **params.merge_within.options,
                                         callback = non_negative_callback),
        check_id_before_merge: bool = typer.Option(*params.check_id_before_merge(),
                                                   **params.check_id_before_merge.options),
        reference: Optional[List[str]] = typer.Option(*params.reference(), **params.reference.options,
                                                      callback = split_callback_list),
        assembly: Path = typer.Option(*params.assembly(), **params.assembly.options,
                                      **oparams.file_valid),
        annotation: Path = typer.Option(*params.annotation(), **params.annotation.options,
                                        **oparams.file_valid),
        db: str = typer.Option(*params.db(), **params.db.options, callback = minor_g.db_callback),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = minor_g.attr_mod_callback),
        ext_gene: Optional[List[Path]] = typer.Option(*params.ext_gene(), **params.ext_gene.options,
                                                        **oparams.file_valid),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options,
                                                     **oparams.file_valid),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        genetic_code: str = typer.Option(*params.genetic_code(), **params.genetic_code.options,
                                         callback = minor_g.genetic_code_callback),
        
        ## gRNA generation options
        pam: str = typer.Option(*params.pam(), **params.pam.options),
        length: str = typer.Option(*params.length(), **params.length.options),
        
        ## filter gRNA options (GC)
        gc_min: float = typer.Option(*params.gc_min(), **params.gc_min.options,
                                     callback = non_negative_callback),
        gc_max: float = typer.Option(*params.gc_max(), **params.gc_max.options,
                                     callback = non_negative_callback),
        
        ## filter gRNA options (feature)
        feature: Optional[List[str]] = typer.Option(*params.feature(), **params.feature.options,
                                                    callback = split_callback_list),
        max_insertion: int = typer.Option(*params.max_insertion(), **params.max_insertion.options,
                                          callback = non_negative_callback),
        min_within_n: int = typer.Option(*params.min_within_n(), **params.min_within_n.options,
                                         callback = positive_callback),
        min_within_fraction: float = typer.Option(*params.min_within_fraction(),
                                                  **params.min_within_fraction.options,
                                                  callback = zero_to_one_callback),
        
        ## filter gRNA options (off-target)
        background: Optional[List[Path]] = typer.Option(*params.background(), **params.background.options,
                                                        **oparams.file_valid),
        screen_reference: bool = typer.Option(*params.screen_reference(), **params.screen_reference.options),
        unmask_ref: bool = typer.Option(*params.unmask_ref(), **params.unmask_ref.options),
        mask_gene: Optional[List[str]] = typer.Option(*params.mask_gene(), **params.mask_gene.options,
                                                      callback = split_callback_list),
        unmask_gene: Optional[List[str]] = typer.Option(*params.unmask_gene(), **params.unmask_gene.options,
                                                        callback = split_callback_list),
        mask_cluster: Optional[List[str]] = typer.Option(*params.mask_cluster(),
                                                         **params.mask_cluster.options,
                                                         callback = split_callback_list),
        unmask_cluster: Optional[List[str]] = typer.Option(*params.unmask_cluster(),
                                                           **params.unmask_cluster.options,
                                                           callback = split_callback_list),
        # mask_homologue: bool = typer.Option(*params.mask_homologue(), **params.mask_homologue.options),
        ot_pamless: bool = typer.Option(*params.ot_pamless(), **params.ot_pamless.options),
        ot_mismatch: int = typer.Option(*params.ot_mismatch(), **params.ot_mismatch.options,
                                        callback = non_negative_callback),
        ot_gap: int = typer.Option(*params.ot_gap(), **params.ot_gap.options, callback = non_negative_callback),
        ot_indv: Optional[List[str]] = typer.Option(*params.ot_indv(), **params.ot_indv.options,
                                                    callback = split_callback_list),
        skip_bg_check: bool = typer.Option(*params.skip_bg_check(), **params.skip_bg_check.options),

        ## filter gRNA (other)
        exclude: Path = typer.Option(*params.exclude(), **params.exclude.options, **oparams.file_valid),
        
        ## gRNA minimum set options
        accept_invalid: bool = typer.Option(*params.accept_invalid(), **params.accept_invalid.options),
        accept_feature_unknown: bool = typer.Option(*params.accept_feature_unknown(),
                                                    **params.accept_feature_unknown.options),
        sets: int = typer.Option(*params.sets(), **params.sets.options, callback = positive_callback),
        # sc_algorithm: SetCoverAlgo = typer.Option(*params.sc_algorithm(), **params.sc_algorithm.options),
        auto: bool = typer.Option(*params.auto(), **params.auto.options),
        output_ver: int = typer.Option(*params.output_ver(), **params.output_ver.options),
        
        ## user lookups (this is right at the end because we might need the other args)
        ## genomes is not "eager" because we have to wait for --genome-lookup to be processed
        genomes: bool = typer.Option(*params.genomes(), **params.genomes.options,
                                     callback = minor_g.make_genomes_callback("independent")),
        references: bool = typer.Option(*params.references(), **params.reference.options,
                                        callback = minor_g.references_callback),
        ## clusters & members are not "eager" because we have to wait for --cluster-lookup to be processed
        clusters: bool = typer.Option(*params.clusters(), **params.clusters.options,
                                      callback = minor_g.clusters_callback),
        members: str = typer.Option(*params.members(), **params.members.options,
                                    callback = minor_g.members_callback),
        ## more user lookups, this time for the lookups themselves haha
        cluster_set: str = typer.Option(*params.cluster_set(), **params.cluster_set.options,
                                        **oparams.file_valid, is_eager = True,
                                        callback = minor_g.cluster_set_callback),
        genome_set: str = typer.Option(*params.genome_set(), **params.genome_set.options,
                                       **oparams.file_valid, is_eager = True,
                                       callback = minor_g.genome_set_callback),
        reference_set: str = typer.Option(*params.reference_set(), **params.reference_set.options,
                                          **oparams.file_valid, is_eager = True,
                                          callback = minor_g.reference_set_callback)):
    """
    Subcommand full.
    
    Executes commands homologue, grna, filter, and minimumset in sequence to
    generate minimum set(s) of gRNA required to cover all targets.
    """
    
    ## check validity of args
    args = Namespace(**locals())
    minor_g.parse_args(args, "full")
    minor_g.subcmd_full()
    minor_g.resolve()
    minor_g.logfile.devsplain("heh")
    return
    

#######################
##    COMMON ARGS    ##
##  + PARSE SUB-CMD  ##
#######################

## note: any updates to arguments for app_sub's sub_main must be reflected in arguments for app_main's main func
## requirements for successful parsing of defualt sub cmd using app_main's main function
## - varname must be one of its CLI keywords for both app_sub's sub_main and app_main's main
@app_sub.callback()
def sub_main(
        quiet: bool = typer.Option(*params.quiet()),
        # directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        # prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        version: bool = typer.Option(*params.version(), callback = minor_g.version_callback, is_eager = True)):
    
    """
    Sub-application for redirecting to correct subcommand.
    
    All subcommands are grouped under this application.
    """
    minor_g.logfile.devsplain("in sub main")
    
    ## config
    if quiet: minor_g.verbose = False
    # if prefix: config.prefix = prefix
    # if directory: config.out_dir = directory
    
    return

## used to redirect to default subcommand if subcommand not provided
## requirements for successful parsing of defualt sub cmd using app_main's main function
## - varname must be one of its CLI keywords for both app_sub's sub_main and app_main's main
## - all args in app_main's main must have empty defaults so error will only be thrown in app_sub's sub_main
##  - exception made for 'prefix', which needs to be defined early to generate log file name
## ## parameter checking will only be done using app_sub's sub_main so that main & sub cmds have the same manual
@app_main.command(context_settings = {"allow_extra_args": True, "ignore_unknown_options": True})
def main(ctx: typer.Context,
         quiet: bool = typer.Option(*params.quiet(None)),
         # directory: str = typer.Option(*params.directory(None)),
         # prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
         version: bool = typer.Option(*params.version(None)),
         keep_on_crash: bool = typer.Option(default = False), ## remove this option once done?
         keep_all: bool = typer.Option(default = False),
         _help: bool = typer.Option(*params.help(None))): ## catch help
    """
    Wrapper main application.
    
    An additional application over app_sub that enables default subcommand by
    redirecting an execution that does not specify a subcommand to subcommand full.
    
    Also handles some execution-level shared arguments
    (``--version``, ``--keep-on-crash``, ``--keep-all``).
    """
    minor_g.keep_tmp = keep_all
    
    ## regenerate/reformat args for sub_main
    def regenerate_main_args(local_args):
        
        sub_main_args = inspect.getfullargspec(sub_main).annotations
        format_key = lambda x: f"-{x}" if len(x) == 1 else f"--{x}"
        
        output = []
        for k, v in local_args.items():
            ## ignore local vars not used in sub_main
            if (k not in sub_main_args): continue
            ## format flags if raised else ignore
            elif (sub_main_args[k] is bool) and (v is not None): output.append(format_key(k))
            ## format keyword args if value provided else ignore
            elif (v is not None): output.extend([format_key(k), v])
        return output
    
    ## this condition should always be fulfilled since app_main doesn't use sub cmds
    if ctx.invoked_subcommand is None:
        
        sub_cmd_args = ctx.args
        
        ## if '--help' raised
        if _help:
            ## append --help to end of args (position relative to subcommand is impt!)
            ## - if the first args in sub_cmd_args starts with '-' (i.e. is argument. checking '-' is only possible because this programme doesn't take positional arguments), clear sub_cmd_args so this function doesn't raise 'no such option' error if the keyword argument belongs to a subcommand. if not, we assume that the first element in sub_cmd_args is the subcommand name and we clear all args from second argument onward
            ##  - if no subcommand provided, it defaults to sub_main's help page
            ##  - if subcommand provided, it uses subcommand's help page
            if sub_cmd_args and sub_cmd_args[0][0] == '-': sub_cmd_args = ["--help"]
            else: sub_cmd_args = sub_cmd_args[:1] + ["--help"]
            ## remove tmpdir since programme exists after printing help manual
            minor_g.cleanup()
        
        ## if subcmd and/or subcmd args provided
        ## - handle unexpected subcommands (i.e. no subcommand or incorrect subcommand name)
        ## elif no subcmd args provided at all, just pass as it is to sub_main and raise sub-main's help page
        elif (sub_cmd_args and sub_cmd_args[0] not in {cmd.name for cmd in app_sub.registered_commands}):
            ## if subcommand not provided, use default subcommand
            ## checking '-' is only possible because this programme doesn't take positional arguments
            if sys.argv[1][0] == '-':
                sub_cmd_args = [default_sub_cmd] + sub_cmd_args
            ## else if subcommand is incorrect, retain sub_cmd_args to trigger the appropriate error message
        
        ## combine <main cmd name>, <sub_main's args>, <sub_cmd's name + args>
        sys.argv = [sys.argv[0]] + regenerate_main_args(locals()) + sub_cmd_args
        ## assign default subcmd if no subcmd or args provided at all, else set subcmd to user-provided subcmd
        minor_g.subcmd = 'full' if not sub_cmd_args else sub_cmd_args[0]
        if keep_on_crash: minor_g.keep_on_crash = True
        minor_g.logfile.devsplain(sys.argv)
        
        ## execute subcommand
        app_sub()
    
    return

    
if __name__ == "__main__":
    
    try:
        minor_g.logfile.update_filename(minor_g.mkfname(os.path.basename(minor_g.logfile.filename)))
        minor_g.logfile.devsplain(minor_g.logfile.filename)
        minor_g.logfile.args(["raw", sys.argv])
    except Exception as e:
        print("Unable to instantiate logger.")
        minor_g.cleanup()
        raise e
    
    ## try-except to handle cleanup of tmpdir upon crash etc.
    try:
        minor_g.raw_args = sys.argv
        app_main()
    except SystemExit as e:
        minor_g.cleanup()
        if e.code != 0: ## click returns 0 upon successful execution; catch and ignore
            raise e
    except MINORgError as e:
        print(e.message)
        minor_g.cleanup()
    except Exception as e:
        typer.echo("problem!")
        minor_g.cleanup()
        raise e

## cleanup for sphinx
minor_g.keep_on_crash = False
minor_g.cleanup()
