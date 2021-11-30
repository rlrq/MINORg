#!/usr/bin/env python3

## TODO: allow users to provide GFF instead of GFF -> BED
## TODO: 2021/06/10: implement --valid-indv/--valid-genome/--valid-acc, --valid-cluster to allow user to check alias validity (print location of FASTA file for --valid-genome, do member_callback for --valid-cluster) before running
## TODO: 2021/06/10: disable Enum for --indv. Since users are now allowed to submit custom genome lookup files using --genome-lookup, autofill shouldn't be used for --indv unless it can dynamically update to accommodate user's --genome-lookup input even before command submission (which is impossible lol). Basically, restructure the whole --indv thing to use the same kind of code as --cluster
## TODO: 2021/09/07: don't merge ext_genome files w/ reference (or ext_bed files w/ bed). Find some way to allow both to be processed whenever blastn or whatever needs to be done to reference
## TODO: 2021/10/21: update all functions to use config.reference_ext and config.annotation_ext (dictionary of {source_id: path_to_fasta/bed}) instead of args.reference in order to accommodate multiple reference files (since we don't want to merge ext_genome with the reference). I've modified some stuff and test with homologue but not the other subcommands.
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

from minorg.functions import get_count_dict, cat_files, reduce_ann, assign_alias, gRNAHits
from minorg.log import MINORgLogger

## import subcommand functions
from minorg.subcmd_homologue import execute_homologue
from minorg.subcmd_grna import execute_grna
from minorg.subcmd_filter import execute_filter
from minorg.subcmd_minimumset import execute_minimumset

from minorg.extend_genome import extend_genome
from minorg.extract_homologue import get_ref_by_genes_resolve

from minorg.parse_config import (
    Config, Params, OptionParams,
    SetCoverAlgo, IndvGenomesAll, IndvGenomesAllClear, Lookup,
    LogFile,
    parse_multiline_multikey_sdict
)

from exceptions import (
    MessageError,
    InputFormatError,
    InvalidPath,
    InvalidFile,
    UnreadableFile
)

from display import (
    print_indent as printi,
    print_overwrite_multi as printom
)

__version_main__ = "3.1"
__version_full__ = "2.1"
__version_homologue__ = "1.0"
__version_grna__ = "1.0"
__version_filter__ = "1.0"
__version_minimumset__ = "1.1"

default_sub_cmd = "full"
app_main = typer.Typer()
app_sub = typer.Typer()

logging_level = logging.DEBUG

## parse config
try:
    CONFIG = os.environ["MINORG_CONFIG"]
except KeyError:
    # ## TODO: set env variable so we can remove this KeyError exception handling thing
    # CONFIG = "/mnt/chaelab/rachelle/scripts/minorgpy/config.ini"
    CONFIG = None

## generate Params, OParmas, and Config objects
params = Params(CONFIG)
oparams = OptionParams() ## namespace for some pre-defined sets of parameter options
config = Config(params, keep_on_crash = True)

## argument parsing functions
def is_true(value):
    return value is True

def valid_readable_file(pathname):
    pathname = os.path.abspath(pathname)
    if not os.path.exists(pathname):
        raise InvalidPath(pathname)
    if not os.path.isfile(pathname):
        raise InvalidFile(pathname)
    if not os.access(pathname, os.R_OK):
        raise UnreadableFile(pathname)
    return pathname

def version_callback(value: bool):
    if value:
        typer.echo(
            f"minorg version: {__version_main__}\n"
            f"minorg full: {__version_full__}\n"
            f"minorg homologue: {__version_homologue__}\n"
            f"minorg grna: {__version_grna__}\n"
            f"minorg filter: {__version_filter__}\n"
            f"minorg minimumset: {__version_minimumset__}"
        )
        config.cleanup()
        raise typer.Exit()

def split_none(iterable):
    return None if ( iterable is None or len(iterable) == 0 ) else \
        (iterable if isinstance(iterable, str) else ','.join(iterable)).split(',')

def parse_lookup(iterable, lookup, return_first = False):
    '''
    commas not allowed in values
    '''
    iterable = split_none(iterable)
    if (not isinstance(lookup, dict)) or iterable is None: return None
    mapped = [lookup.get(val, val) for val in iterable]
    if return_first: return mapped[0]
    else: return mapped

# def parse_cluster_set(val):
#     return params.cluster_mapping.get(val, val)

def parse_cluster_set(val):
    return parse_lookup(val, params.cluster_mapping, return_first = True)

def get_cluster_genes(val):
    cluster_genes = config.cluster_aliases.get(val, None)
    if cluster_genes is None:
        raise click.UsageError( (f"'{val}' is not a valid cluster name in"
                                 f" lookup file '{config.cluster_set}'") )
    else:
        return tuple(cluster_genes.split(','))

def parse_genes_from_args(args):
    '''
    To be called after 'gene' or 'cluster' have passed check_homologue_args checks
    Returns: {<prefix>: (gene1, gene2, gene3)}
    '''
    if args.target is not None:
        ## set None as val; execute_homologue will exit if it detects that val is None
        return {args.prefix: None}
    elif args.cluster is None:
        return {args.prefix: tuple(args.gene)}
    else:
        output = {}
        for cluster in args.cluster:
            output[f"{args.prefix}_{cluster}"] = tuple(get_cluster_genes(cluster))
        return output
    return

def parse_genes_for_filter(args, priority = None):
    '''
    To be called after 'mask', 'unmask', 'mask_cluster', 'unmask_cluster', and 'gene' 
        have passed check_filter_args checks
    Returns {<mask/unmask/gene>: (gene1, gene2, gene3)}
    '''
    output = {}
    genes = (set() if args.gene is None else set(args.gene))
    mask_genes = set((tuple() if args.mask_gene is None else tuple(args.mask_gene)) + \
                     (tuple() if args.mask_cluster is None else
                      tuple(itertools.join(*[get_cluster_genes(cluster)
                                             for cluster in args.mask_cluster]))))
    unmask_genes = set((tuple() if args.unmask_gene is None else tuple(args.unmask_gene)) + \
                       (tuple() if args.unmask_cluster is None else
                        tuple(itertools.join(*[get_cluster_genes(cluster)
                                               for cluster in args.unmask_cluster]))))
    ## parse '.' for 'all genes passed to --gene'
    def expand_genes(genes_to_expand):
        if '.' in genes_to_expand and '-' in genes_to_expand:
            raise MessageError("ERROR: '-' and '.' are mutually exclusive for --mask-gene and --unmask-gene")
        if '.' in genes_to_expand:
            genes_to_expand -= {'.'}
            genes_to_expand |= genes
        elif '-' in genes_to_expand:
            genes_to_expand -= {'-'}
            genes_to_expand -= genes
        return genes_to_expand
    mask_genes = expand_genes(mask_genes)
    unmask_genes = expand_genes(unmask_genes)
    ## check for overlaps between mask and unmask if priority is not set
    if priority is None:
        overlap_genes = mask_genes.intersection(unmask_genes)
        if overlap_genes:
            raise MessageError( ("ERROR: The following genes were provided to both"
                                 f" --mask-gene and --unmask-gene: {','.join(sorted(overlap_genes))}") )
    ## generate output dictionary
    output["gene"] = tuple(genes)
    if priority == "mask":
        output["mask"] = tuple(mask_genes)
        output["unmask"] = tuple(unmask_genes - mask_genes)
    elif priority == "unmask":
        output["mask"] = tuple(mask_genes - unmask_genes)
        output["unmask"] = tuple(unmask_genes)
    else:
        output["mask"] = tuple(mask_genes)
        output["unmask"] = tuple(unmask_genes)
    return output

def split_callback_str(val):
    return ','.join(split_none(val))

def split_callback_list(val):
    return split_none(val)

def attr_mod_callback(val):
    if isinstance(val, str):
        try:
            from re import search
            return params.parse_attr_mod(search(r"""[^'"]+""", val).group(0))
        except:
            if val in params.attr_mod_presets:
                raise InputFormatError(message = ( "Unable to parse the GFF attribute modification preset"
                                                   f" ({params.attr_mod_presets[val]}) retrieved"
                                                   " from the config file ({params.config_file})"
                                                   " via alias '{val}'." ),
                                       hint = params.attr_mod.help())
            else: raise InputFormatError(error_src = "for parameter --attr-mod", hint = params.attr_mod.help())
    else: return params.attr_mod.default

def reference_required(args, msg):
    if args.reference is None and (args.assembly is None or args.annotation is None):
        raise click.UsageError( ("'-r <reference genome alias>' OR"
                                 " '--assembly <reference assembly alias or path to FASTA file>"
                                 " --annotation <reference annotation alias or GFF3 file or"
                                 " BED file converted from GFF3 using gff2bed>'is required if using"
                                 f" {msg}.") )
    
def domain_callback(val: str):
    config.raw_domain = val
    if val is not None:
        parsed_domain = params.parse_domain(val)
        if str(parsed_domain) == val and val != "gene":
            typer.echo(f"Parsing input to --domain ({val}) as PSSM-Id")
        return parsed_domain
    
def alias_or_file_callback(val: str, lookup: dict, lookup_fail_message = None):
    mapped_val = parse_lookup(val, lookup, return_first = True)
    # if mapped_val == val and lookup_fail_message: typer.echo(lookup_fail_message(val))
    if mapped_val is not None: return valid_readable_file(mapped_val)

def make_alias_or_file_callback(lookup: dict, param = None,
                                lookup_fail_message = lambda x, y: (f"Unrecognised alias passed to {x}."
                                                                    f" Parsing '{y}' as path.")):
    lookup_fail_message_updated = lambda x: lookup_fail_message('/'.join(param.names), x)
    return (lambda val: alias_or_file_callback(val, lookup = lookup,
                                               lookup_fail_message = lookup_fail_message_updated))


def check_reference_args(args):
    if args.assembly is not None and args.annotation is not None:
        assembly_alias_or_file_callback = make_alias_or_file_callback(params.assembly_aliases, params.assembly)
        assembly_mapped = assembly_alias_or_file_callback(args.assembly)
        if assembly_mapped is not None:
            assembly_mapped = os.path.abspath(assembly_mapped)
            params.indv_genomes["ref"] = assembly_mapped
        ann_alias_or_file_callback = make_alias_or_file_callback(params.annotation_aliases, params.annotation)
        ann_mapped = ann_alias_or_file_callback(args.annotation)
        if ann_mapped is not None:
            ann_mapped = os.path.abspath(ann_mapped)
        config.set_reference(assembly_mapped, ann_mapped)
        args.assembly = assembly_mapped
        args.annotation = ann_mapped
        config.reference_aliases["Reference"] = (args.assembly, args.annotation)
    if args.reference:
        none_val = '-'
        valid_aliases(aliases = args.reference, lookup = config.reference_aliases,
                      none_value = '-', all_value = "all",
                      param = params.reference, display_cmd = "--references",
                      additional_message = ("Alternatively, provide a FASTA file of a genome assembly using"
                                            " '--assembly <path to FASTA file>' and a GFF3 file of genome"
                                            " annotations using '--annotation <path to GFF3 file>'."))
        if set(args.reference) == {none_val}:
            args.reference = None
        else:
            for alias in args.reference:
                fasta, ann = config.reference_aliases[str(alias)]
                config.extend_reference(alias, fasta, ann)
    if ((args.reference and args.assembly and args.annotation) or
        (not args.reference and not (args.assembly and args.annotation))):
        raise click.UsageError((f"Either '--reference <alias>' OR"
                                " '--assembly <FASTA> --annotation <GFF3>' is required"))

## set ref: <reference genome> in params.indv_genomes
def assembly_callback(val):
    return val
    # alias_or_file_callback = make_alias_or_file_callback(params.assembly_aliases, params.assembly)
    # mapped_val = alias_or_file_callback(val)
    # if mapped_val is not None:
    #     params.indv_genomes["ref"] = mapped_val
    # config.set_reference(mapped_val, None)
    # return mapped_val

def annotation_callback(val):
    return val
    # alias_or_file_callback = make_alias_or_file_callback(params.annotation_aliases, params.annotation)
    # mapped_val = alias_or_file_callback(val)
    # config.set_reference(None, mapped_val)
    # return mapped_val

def reference_annotation_callback(val):
    return val

def positive_callback(val):
    if val > 0:
        return val
    else:
        raise typer.BadParameter( f"Invalid value: {val}. Positive value required." )

def non_negative_callback(val):
    if val >= 0:
        return val
    else:
        raise typer.BadParameter( f"Invalid value: {val}. Non-negative value required." )

def zero_to_one_callback(val):
    val = non_negative_callback(val)
    if val <= 1:
        return val
    else:
        raise typer.BadParameter( f"Invalid value: {val}. Must be between 0 and 1 (inclusive)." )

def reference_set_callback(val: str):
    print("reference_set_callback")
    val = parse_lookup(val, params.reference_sets, return_first = True)
    if val is not None:
        if os.path.exists(os.path.abspath(val)):
            config.reference_set = os.path.abspath(val)
            try:
                with open(config.reference_set, 'r') as f:
                    config.reference_aliases = {**config.reference_aliases,
                                                **{k: v.split('\t') for k, v in
                                                   parse_multiline_multikey_sdict(f.read(),
                                                                                  kv_sep = '\t').items()}}
            except:
                raise InputFormatError(error_src = "file", hint = f"File: {config.reference_set}")
        else:
            raise typer.BadParameter( ( "Unrecognised reference set lookup file alias or non-existent file:"
                                        f" {val}"
                                        "\nFor a list of valid reference set lookup file aliases"
                                        " and their locations, execute 'minorg --references'."
                                        "\nAlternatively, please provide a valid file of",
                                        " reference alias-FASTA-GFF3 mapping."
                                        "\nHelp message for --reference-set: "
                                        f"{params.reference_set.help}") )
    return val

def cluster_set_callback(val: str):
    val = parse_lookup(val, params.cluster_sets, return_first = True)
    if val is None:
        config.cluster_aliases = {}
    elif os.path.exists(os.path.abspath(val)):
        config.cluster_set = os.path.abspath(val)
        try:
            with open(config.cluster_set, 'r') as f:
                config.cluster_aliases = parse_multiline_multikey_sdict(f.read(), kv_sep = '\t')
        except:
            raise InputFormatError(error_src = "file", hint = f"File: {config.cluster_set}")
    else:
        raise typer.BadParameter( ( f"Unrecognised cluster set lookup file alias or non-existent file: {val}"
                                    "\nFor a list of valid cluster set lookup file aliases and their locations,"
                                    " execute 'minorg --clusters'."
                                    "\nAlternatively, please provide a valid file of",
                                    " cluster alias-member mapping."
                                    "\nHelp message for --cluster-set: "
                                    f"{params.cluster_set.help}") )
    return val

def genome_set_callback(val: str):
    val = parse_lookup(val, params.genome_sets, return_first = True)
    if val is None:
        config.genome_aliases = params.indv_genomes
    elif os.path.exists(os.path.abspath(val)):
        config.genome_set = os.path.abspath(val)
        try:
            with open(config.genome_set, 'r') as f:
                config.genome_aliases = parse_multiline_multikey_sdict(f.read(), kv_sep = '\t')
        except:
            raise InputFormatError(error_src = "file", hint = f"File: {config.genome_set}")
    else:
        raise typer.BadParameter( ( f"Unrecognised genome set lookup file alias or non-existent file: {val}"
                                    "\nFor a list of valid genome set lookup file aliases and their locations,"
                                    " execute 'minorg --genomes'."
                                    "\nAlternatively, please provide a valid file of",
                                    " genome alias-filename mapping."
                                    "\nHelp message for --genome-set: "
                                    f"{params.genome_set.help}") )
    return val

def references_callback(value: bool):
    if value:
        if config.reference_set is None:
            typer.echo( (f"\nNo reference genomes have been defined."
                         " Please use '--reference-set' to specify a file containing"
                         " alias-FASTA-GFF3 mapping OR use '--assembly <FASTA>' and '--annotation <GFF3>'"
                         " to specify genome sequences and annotations respectively." ))
        else:
            typer.echo(f"\nValid genome aliases (defined in {config.reference_set}):\n")
            typer.echo("<semicolon-separated genome alias(es)>\t<FASTA file>\t<GFF3 file>")
            with open(config.reference_set, 'r') as f:
                aliases = [f.read()]
            typer.echo('\n'.join(aliases) + '\n')
        config.cleanup()
        raise typer.Exit()

def make_genomes_callback(mode: str = "independent"):
    def genomes_callback(value: bool):
        if value:
            if config.genome_set is None:
                typer.echo(f"\nValid genome aliases (defined in {params.config_file}):\n")
                aliases = ['\t'.join(kv) for kv in params.indv_genomes.items()]
            else:
                typer.echo(f"\nValid genome aliases (defined in {config.genome_set}):\n")
                with open(config.genome_set, 'r') as f:
                    aliases = [f.read()]
            typer.echo("<semicolon-separated genome alias(es)>\t<FASTA file>")
            if mode == "independent":
                prepend = [".\t<all genomes except reference>",
                           "ref\t<reference genome>"]
            elif mode == "dependent":
                aliases = [".\t<all genomes passed to -g, including reference>",
                           "-\t<exclude all genomes passed to -g, including reference>",
                           "ref\t<reference genome>"]
            typer.echo('\n'.join(prepend + aliases) + '\n')
            config.cleanup()
            raise typer.Exit()
    return genomes_callback

def clusters_callback(value: bool):
    if value:
        if params.cluster_set is not None:
            typer.echo(f"\nThe following information is retrieved from {config.cluster_set}:\n")
            typer.echo("<semicolon-separated cluster alias(es)>\t<comma-separated cluster members>")
            with open(config.cluster_set, 'r') as f:
                typer.echo(f.read())
        else:
            typer.echo( ("\nNo cluster set lookup file has been specified."
                         "\nYou may update the list of cluster set lookup files"
                         " and set a default set lookup file by"
                         f" modifying the relevant fields in the config file located at {params.config_file}"
                         " or provide a set lookup file for a single execution using --cluster-set.") )
        config.cleanup()
        raise typer.Exit()

def members_callback(value: str):
    if value:
        if value not in config.cluster_aliases:
            if params.cluster_set is not None:
                typer.echo( ( f"{value} is not a valid cluster alias in the following cluster set lookup file:"
                              f" {config.cluster_set}" ) )
            else:
                typer.echo( ("\nNo cluster set lookup file has been specified."
                             "\nYou may update the list of cluster set files"
                             " and set a default set lookup file by"
                             f" modifying the relevant fields in the config file located at {params.config_file}"
                             " or provide a set lookup file for a single execution using --cluster-set.") )
        else:
            typer.echo(f"\nThe following information is retrieved from {params.cluster_set}:\n")
            typer.echo(config.cluster_aliases[value] + '\n')
        config.cleanup()
        raise typer.Exit()

def lookup_sets_callback(lookup: Lookup):
    lookup.print_sets_aliases()
    return

# def genome_sets_callback(value: bool):
#     if value:
#         lookup_sets_callback(config.genome_sets)

# def cluster_sets_callback(value: bool):
#     if value:
#         lookup_sets_callback(config.genome_sets)

    
def valid_aliases(aliases, lookup, raise_error = True, message = None, param = None,
                  none_value = None, all_value = None, clear_value = None,
                  display_cmd = None, additional_message = None):
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
    for alias in aliases:
        if alias not in lookup and alias != all_value and alias != clear_value: unknown.append(alias)
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
        else: return None
    else: return None

def make_reduced_ann(args, config, gene_sets):
    config.logfile.splain("Filtering annotation file(s)")
    if gene_sets:
        # ann_red = config.mkfname(f"tmp_reduced.bed", tmp = True)
        # config.annotation_red = ann_red
        # reduce_ann(gff_beds = config.annotation_ext.values(),
        #            ids = tuple(itertools.chain(*tuple(gene_sets.values()))),
        #            fout = ann_red,
        #            mk_tmpf_name = lambda x: config.reserve_fname(f"tmp_reduced.{x}.gff", tmp = True))
        dir_red_ann = config.mkfname("reduced_ann", tmp = True)
        config.annotation_red = reduce_ann(gff_beds = config.annotation_ext,
                                           ids = tuple(itertools.chain(*tuple(gene_sets.values()))),
                                           fout_fmt = "GFF",
                                           mk_tmpf_name = lambda x: config.reserve_fname(dir_red_ann,
                                                                                         f"reduced.{x}.gff",
                                                                                         tmp = True))
        # bed_reds = []
        # for src, bed in config.bed_ext.items():
        #     bed_red = config.mkfname(f"tmp_reduced.{src}.bed", tmp = True)
        #     bed_reds.append(bed_red)
        #     extract_features_and_subfeatures(bed, tuple(itertools.chain(*tuple(gene_sets.values()))),
        #                                      bed_red, quiet = True, fin_fmt = "BED", fout_fmt = "BED")
        # bed_red = config.mkfname(f"tmp_reduced.bed", tmp = True)
        # cat_files(bed_reds, bed_red, remove = True)
        # config.bed_red = bed_red
    return
    
def check_homologue_args(args, standalone = True):
    
    ## check if user has provided indv
    indv_provided = not (len(args.indv) == 1 and args.indv[0] == IndvGenomesAll.none)
    
    ## MUTUALLY EXCLUSIVE ARGS
    ## -q and -t are mutually exclusive
    if args.query and args.target is not None:
        print(args.query, args.target)
        raise click.UsageError("'-q <FASTA file>' and '-t <FASTA file>' are mutually exclusive.")
    ## -g and -c are mutually exclusive
    if args.gene is not None and args.cluster is not None:
        raise click.UsageError("'-g <gene(s)>' and '-c <cluster(s)>' are mutually exclusive.")
    ## if checking args for function 'homologue' and not 'full', then -q/-t are not compatible with -g/-c/-i
    if ( standalone and
         ( ( args.query or args.target is not None) and
           ( args.gene is not None or args.cluster is not None or indv_provided ) ) ):
        raise click.UsageError( ( "'-q <FASTA file>', '-t <FASTA file>', and"
                                  " '[ -g <gene(s)> or -c <cluster(s)> ] -i <individual(s)>'"
                                  " are mutually exclusive.") )
    
    ## REQUIRED ARGS
    ## check that one of -q, -t, -a|-i is provided
    if not args.query and args.target is None and not indv_provided:
        raise click.UsageError( ("One of the following is required:"
                                 " '-q <FASTA file>', '-t <FASTA file>', or '-i <individual(s)>'") )
    ## if -t is not used, check that -g/-c is provided
    if args.target is None:
        if args.gene is None and args.cluster is None:
            raise click.UsageError("'-g <gene(s)>' or '-c <cluster(s)>' is required if not using -t.'")
    ## if -g/-c is provided, check that --ref and --annotation are also provided
    if args.gene is not None or args.cluster is not None:
        reference_required(args, "'-g <gene(s)>' or '-c <cluster(s)>'")
        # if args.reference is None and (args.assembly is None or args.annotation is None):
        #     raise click.UsageError( ("'-r <reference genome alias>' OR"
        #                              " '--assembly <reference assembly alias or path to FASTA file>"
        #                              " --annotation <reference annotation alias or GFF3 file or"
        #                              " BED file converted from GFF3 using gff2bed>'is required if using"
        #                              " '-g <gene(s)>' or '-c <cluster(s)>'.") )
        # if args.reference is None:
        #     raise click.UsageError( ("'-r <reference genome FASTA>' is required if using"
        #                              " '-g <gene(s)>' or '-c <cluster(s)>'.") )
        # if args.annotation is None:
        #     raise click.UsageError( ("'--annotation <GFF3 file or BED file converted from GFF3"
        #                              " annotations using gff2bed>' is required if using"
        #                              " '-g <gene(s)>' or '-c <cluster(s)>'.") )
        if args.query is None and not indv_provided:
            raise click.UsageError( ("'-q <FASTA file>' or '-i <individual(s)>'"
                                     " is required if using '-g <gene(s)>' or '-c <cluster(s)>'.") )
    ## if --extend-genome or --extend-cds is provided, check that the other is also provided
    if sum(map(lambda x: x is None, [args.ext_cds, args.ext_genome])) == 1:
        raise click.UsageError( ("'--extend-cds <FASTA file>' and '--extend-genome <FASTA file>'"
                                 " should either be used together or not used at all.") )
    ## check that --mafft is provided
    if not standalone and args.mafft is None:
        raise click.UsageError( "Path to mafft executable is required: --mafft <path>" )
    ## check that --blastn is provided if using -q or '-i <not ref only>'
    if ( (not standalone or args.query or (indv_provided and set(args.indv) != {"ref"}))
         and args.blastn is None):
        raise click.UsageError( "Path to blastn executable is required: --blastn <path>" )
    ## check that --rpsblast is provided if using --domain
    if args.domain is not None and args.domain != "gene":
        if args.rpsblast is None:
            raise click.UsageError( ("'--rpsblast <path to rpsblast or rpsblast+ executable>'"
                                     " is required if using '--domain <domain alias or PSSM-Id>'") )
        if args.db is None:
            raise click.UsageError( ("'--db <alias of or path to local RPS-BLAST database"
                                     " OR name of remote RPS-BLAST database>'"
                                     " is required if using '--domain <domain alias or PSSM-Id>'") )
    
    ## VALID ALIASES (check -r, --annotation, --rps-db, --indv & -c & --attr-mod)
    ## not a callback because it requires config.cluster_aliases generated by cluster_set_callback
    if args.cluster is not None:
        valid_aliases(aliases = args.cluster, lookup = config.cluster_aliases,
                      param = params.cluster, display_cmd = "--clusters",
                      additional_message = ( "Alternatively, manually input the desired genes using"
                                             " '-g <gene(s)>' or provide a different cluster set lookup file"
                                             " using '--cluster-set <path to file>'."))
    # ## not a callback to standardise valid alias checking w/ args.cluster
    # if indv_provided:
    #     # genome_aliases = params.indv_genomes
    #     # if args.genome_set is not None:
    #     #     if args.genome_set in params.genome_sets:
    #     #         ## TODO: warn user if alias in genome_set_aliases is already in params.indv_genomes and handle it so that the alias in genome_set_aliases is given priority
    #     #         genome_set_aliases = {}
    #     #         with open(args.genome_set, 'r') as f:
    #     #             genome_set_aliases = {alias: fname for alias, fname in
    #     #                                   [x.split('\t') for x in f.read().split('\n')]}
    #     #         genome_aliases = {**genome_aliases, **genome_set_aliases}
    #     if "ref" not in config.genome_aliases:
    #         config.genome_aliases["ref"] = args.reference
    #     valid_aliases(aliases = args.indv, lookup = config.genome_aliases,
    #                   none_value = IndvGenomesAll.none.value, all_value = IndvGenomesAll.all.value,
    #                   param = params.indv, display_cmd = "--genomes",
    #                   additional_message = ( "Alternatively, provide a FASTA file of the genome in"
    #                                          " which to query using '-q <path to FASTA file>'."))
    
    ## ARGS DEPENDENT ON VALUE OF OTHER ARGS
    ## raise check_recip if relax_recip is raised
    if args.relax_recip:
        args.check_recip = True
    
    ## generate query mapping
    query_map = []
    if args.query:
        query_map += [[i+1, str(query)] for i, query in enumerate(args.query)]
    if indv_provided:
        indvs_special = {"all", "none", "ref", '.'}
        indvs_ref = set(indv for indv in config.reference_ext) if "ref" in args.indv else set()
        print("indvs_ref:", indvs_ref)
        if "all" in args.indv or '.' in args.indv:
            indvs_genome = set(indv for indv in config.genome_aliases if indv not in indvs_special)
        else:
            indvs_genome = set(args.indv) - indvs_special
        indvs = (set(args.indv) - indvs_special).union(indvs_genome).union(indvs_ref)
        if set(args.indv) != {"ref"}: args.indv = type(args.indv)(indvs)
        valid_aliases(aliases = list(indvs_genome), lookup = config.genome_aliases,
                      none_value = IndvGenomesAll.none.value, all_value = IndvGenomesAll.all.value,
                      param = params.indv, display_cmd = "--genomes",
                      additional_message = ( "Alternatively, provide a FASTA file of the genome in"
                                             " which to query using '-q <path to FASTA file>'."))
        # if IndvGenomesAll.all in args.indv:
        #     indvs_special = {IndvGenomesAll.all, IndvGenomesAll.none}
        #     indvs_all = set(indv for indv in IndvGenomesAll if \
        #                     (indv not in indvs_special.union({IndvGenomesAll.ref})))
        #     indvs = (set(args.indv) - indvs_special).union(indvs_all) ## union in case user also specified 'ref'
        #     args.indv = type(args.indv)(indvs)
        query_map += [[indv, config.genome_aliases[str(indv)]] for indv in indvs_genome] + \
                     [[alias, config.reference_aliases[str(alias)][0]] for alias in indvs_ref]
    ## check that file names are unique
    ## - if multiple queries map to the same file, sort by query id and retain first entry
    ## - if same file provided to args.query and args.indv, priorities args.indv regardless of sort order
    query_count = get_count_dict([x[1] for x in query_map])
    query_multi = set(k for k, v in query_count.items() if v > 1)
    config.query_map = [x for x in query_map if x[1] not in query_multi]
    for fname in query_multi:
        queries = sorted([x[0] for x in query_map if x[1] == fname])
        ## check if IndvGenomesAll object is in queries. If, yes, sort and retain first
        indvgenomesobjs = sorted([x for x in queries if isinstance(x, IndvGenomesAll)])
        if indvgenomesobjs: config.query_map.append([indvgenomesobjs[0], fname])
        else: config.query_map.append([queries[0], fname])
    
    # ## doubles as parameter parsing
    # if args.domain is not None:
    #     parsed_domain = params.parse_domain(args.domain)
    #     if str(parsed_domain) == args.domain: typer.echo(f"Parsing input to --domain ({args.domain}) as PSSM-Id")
    #     args.domain = parsed_domain
    
    return
    # return all(checks)

def check_grna_args(args, standalone = True):
    
    ## MUTUALLY EXCLUSIVE ARGS
    ## -g, -c, and -t are mutually exclusive
    if sum(map(lambda x: x is not None, [args.gene, args.cluster, args.target])) > 1:
        raise click.UsageError( ("'-g <gene(s)>', '-c <cluster(s)>', and '-t <path to FASTA file>'"
                                 " are mutually exclusive.") )
    # if args.target is not None and args.feature is not None:
    #     config.logfile.warning("-t <FASTA> is incompatible with -f <feature>. -f will be ignored.")
    
    ## REQUIRED ARGS
    ## check that --pam is provided
    if args.pam is None:
        raise click.UsageError( "'-p <PAM pattern>' is required." )
    ## check that --length is provided
    if args.length is None:
        raise click.UsageError( "'-l <gRNA length (bp)>' is required." )
    ## if -g/-c is provided, check that --ref and --annotation are also provided
    if (args.gene is not None or args.cluster is not None):
        ## set dummy arguments so execute_homologue doesn't break
        if standalone:
            ## set indv to reference
            args.indv = ["ref"]
            ## set arguments to None or False
            args.minlen = args.minid = args.mincdslen = args.merge_within = 0
            args.max_insertion = args.min_within_fraction = 0
            args.min_within_n = 1
            args.check_recip = args.relax_recip = args.check_id_premerge = False
            args.blastn = args.background = None
        reference_required(args, "'-g <gene(s)>' or '-c <cluster(s)>'")
        # if args.reference is None and (args.assembly is None or args.annotation is None):
        #     raise click.UsageError( ("'-r <reference genome alias>' OR"
        #                              " '--assembly <reference assembly alias or path to FASTA file>"
        #                              " --annotation <reference annotation alias or GFF3 file or"
        #                              " BED file converted from GFF3 using gff2bed>'is required if using"
        #                              " '-g <gene(s)>' or '-c <cluster(s)>'.") )
        # if args.reference is None:
        #     raise click.UsageError( ("'-r <reference genome FASTA>' is required if using"
        #                              " '-g <gene(s)>' or '-c <cluster(s)>'.") )
        # if args.annotation is None:
        #     raise click.UsageError( ("'--annotation <GFF3 file or BED file converted from GFF3"
        #                              " annotations using gff2bed>' is required if using"
        #                              " '-g <gene(s)>' or '-c <cluster(s)>'.") )
    ## check that --rpsblast is provided if using --domain
    if args.domain is not None and args.domain != "gene":
        if args.rpsblast is None:
            raise click.UsageError( ("'--rpsblast <path to rpsblast or rpsblast+ executable>'"
                                     " is required if using '--domain <domain alias or PSSM-Id>'") )
        if args.db is None:
            raise click.UsageError( ("'--db <alias of or path to local RPS-BLAST database"
                                     " OR name of remote RPS-BLAST database>'"
                                     " is required if using '--domain <domain alias or PSSM-Id>'") )
    ## if --extend-genome or --extend-cds is provided, check that the other is also provided
    if sum(map(lambda x: x is None, [args.ext_cds, args.ext_genome])) == 1:
        raise click.UsageError( ("'--extend-cds <FASTA file>' and '--extend-genome <FASTA file>'"
                                 " should either be used together or not used at all.") )
    
    ## VALID ALIASES (check -r, --annotation, --db, -c)
    ## not a callback because it requires config.cluster_aliases generated by cluster_set_callback
    if args.cluster is not None:
        valid_aliases(aliases = args.cluster, lookup = config.cluster_aliases,
                      param = params.cluster, display_cmd = "--clusters",
                      additional_message = ( "Alternatively, manually input the desired genes using"
                                             " '-g <gene(s)>' or provide a different cluster set lookup file"
                                             " using '--cluster-set <path to file>'."))
    
    ## some parsing
    if args.target is not None:
        ## set args.indv to empty list so it doesn't break execute_homologue
        ##  note: execute_homologue is used to retrieve reference genes if requested by user
        args.indv = []
    
def check_filter_args(args, standalone = True):
    
    if standalone:
        if args.check_all:
            args.gc_check = args.feature_check = args.background_check = True
        if not args.gc_check and not args.feature_check and not args.background_check:
            raise click.UsageError( ("At least one of the following is required:"
                                     " '--gc-check', '--feature-check', '--background-check'") )
    else:
        args.gc_check = args.feature_check = args.background_check = True
        if args.target:
            args.feature_check = args.background_check = False
        if args.skip_bg_check:
            args.background_check = False
    
    # ## MUTUALLY EXCLUSIVE ARGS
    # ## -g, -c, and -t are mutually exclusive
    # if sum(map(lambda x: x is not None, [args.gene, args.cluster, args.target])) > 1:
    #     raise click.UsageError( ("'-g <gene(s)>', '-c <cluster(s)>', and '-t <path to FASTA file>'"
    #                              " are mutually exclusive.") )
    # # if args.target is not None and args.feature is not None:
    # #     config.logfile.warning("-t <FASTA> is incompatible with -f <feature>. -f will be ignored.")
    
    ## REQUIRED ARGS
    if standalone:
        ## set dummy arguments so execute_homologue doesn't break
        args.minlen = 1
        args.gene = None
        args.unmask_ref = False
        if args.mapping is None:
            raise click.UsageError("'--mapping <path to minorg .map file>' is required.")
        if args.background_check:
            ## if args.grna is not provided, create FASTA file from mapping file
            if args.grna is None:
                args.grna = config.mkfname("tmp_grna_all.fasta")
                from functions import gRNAHits
                grnahits = gRNAHits()
                grnahits.parse_from_mapping(args.mapping, version = None)
                grnahits.write_fasta(args.grna, write_all = True)
    ## GC filter
    if args.gc_check:
        if args.gc_min is None: args.gc_min = 0
        if args.gc_max is None: args.gc_max = 1
    ## feature filter
    if args.feature_check:
        if not args.feature:
            raise click.UsageError("'--feature <feature>' is required if using '--filter-feature'.")
        if standalone and not args.alignment:
            raise click.UsageError("'--alignment <path to FASTA> is required if using '--filter-feature'.")
        reference_required(args, "'--filter-feature'")
        # if args.reference is None and (args.assembly is None or args.annotation is None):
        #     raise click.UsageError( ("'-r <reference genome alias>' OR"
        #                              " '--assembly <reference assembly alias or path to FASTA file>"
        #                              " --annotation <reference annotation alias or GFF3 file or"
        #                              " BED file converted from GFF3 using gff2bed>'is required if using"
        #                              " '--filter-feature'.") )
        # if not args.annotation:
        #     raise click.UsageError( ("'--annotation <GFF3 file or BED file converted from GFF3"
        #                              " annotations using gff2bed>' is required if using '--filter-feature'.") )
        if args.min_within_n is None: args.min_within_n = 1
        if args.min_within_fraction is None: args.min_within_fraction = 0
    ## background filter
    if args.background_check:
        if args.blastn is None:
            raise click.UsageError( "Path to blastn executable is required: --blastn <path>" )
        if args.unmask_ref:
            if not args.screen_ref:
                raise click.UsageError( "'--unmask-ref' is only used with '--screen-ref'." )
            if args.mask_gene or args.mask_cluster:
                raise click.UsageError( ("'--unmask-ref' should not be used with"
                                         " '--mask-gene' or '--mask-cluster'.") )
            if not standalone and "ref" in args.indv:
                config.logfile.warning( ("'--unmask-ref' will be ignored if 'ref' is provided to '-i'"
                                         " and not excluded using '--ot-indv'. When 'ref' is provided"
                                         " to '-i', the reference genome(s) will be treated the same as"
                                         " any other genomes passed to '-i'.") )
        if args.screen_ref:
            if not (args.reference or (args.assembly and args.annotation)):
                raise click.UsageError("'--reference <reference genome alias>' OR"
                                       " '--assembly <path to FASTA> --annotation <path to GFF3>'"
                                       " is required if using '--screen-ref'.")
        if not args.annotation:
            reference_required(args, "'--unmask-gene', '--mask-cluster', or '--unmask-cluster'")
            # if args.mask_gene or args.unmask_gene or args.mask_cluster or args.unmask_cluster:
            #     raise click.UsageError( ("'--annotation <GFF3 file or BED file converted from GFF3"
            #                              " annotations using gff2bed>' is required if using '--mask-gene',"
            #                              " '--unmask-gene', '--mask-cluster', or '--unmask-cluster'.") )
        if args.ot_mismatch is None: args.ot_mismatch = 1
        if args.ot_gap is None: args.ot_mismatch = 0
        if not args.ot_pamless and args.pam is None:
            raise click.UsageError("'--pam <PAM pattern>' is required if not using '--ot-pamless'.")
    
    # ## check that --pam is provided
    # if args.pam is None:
    #     raise click.UsageError( "'-p <PAM pattern>' is required." )
    # ## check that --length is provided
    # if args.length is None:
    #     raise click.UsageError( "'-l <gRNA length (bp)>' is required." )
    # ## if -g/-c is provided, check that --ref and --annotation are also provided
    # if (args.gene is not None or args.cluster is not None):
    #     ## set dummy arguments so execute_homologue doesn't break
    #     if standalone:
    #         ## set indv to reference
    #         args.indv = ["ref"]
    #         ## set arguments to None or False
    #         args.minlen = args.minid = args.mincdslen = args.merge_within = 0
    #         args.max_insertion = args.min_within_fraction = 0
    #         args.min_within_n = 1
    #         args.check_recip = args.relax_recip = args.check_id_premerge = False
    #         args.blastn = args.background = None
    #     if args.reference is None:
    #         raise click.UsageError( ("'-r <reference genome FASTA>' is required if using"
    #                                  " '-g <gene(s)>' or '-c <cluster(s)>'.") )
    #     if args.annotation is None:
    #         raise click.UsageError( ("'--annotation <GFF3 file or BED file converted from GFF3"
    #                                  " annotations using gff2bed>' is required if using"
    #                                  " '-g <gene(s)>' or '-c <cluster(s)>'.") )
    # ## check that --rpsblast is provided if using --domain
    # if args.domain is not None:
    #     if args.rpsblast is None:
    #         raise click.UsageError( ("'--rpsblast <path to rpsblast or rpsblast+ executable>'"
    #                                  " is required if using '--domain <domain alias or PSSM-Id>'") )
    #     if args.db is None:
    #         raise click.UsageError( ("'--db <alias of or path to local RPS-BLAST database"
    #                                  " OR name of remote RPS-BLAST database>'"
    #                                  " is required if using '--domain <domain alias or PSSM-Id>'") )
    # ## if --extend-genome or --extend-cds is provided, check that the other is also provided
    # if sum(map(lambda x: x is None, [args.ext_cds, args.ext_genome])) == 1:
    #     raise click.UsageError( ("'--extend-cds <FASTA file>' and '--extend-genome <FASTA file>'"
    #                              " should either be used together or not used at all.") )
    
    # ## VALID ALIASES (check -r, --annotation, --db, -c)
    # ## not a callback because it requires config.cluster_aliases generated by cluster_set_callback
    # if args.cluster is not None:
    #     valid_aliases(aliases = args.cluster, lookup = config.cluster_aliases,
    #                   param = params.cluster, display_cmd = "--clusters",
    #                   additional_message = ( "Alternatively, manually input the desired genes using"
    #                                          " '-g <gene(s)>' or provide a different cluster set lookup file"
    #                                          " using '--cluster-set <path to file>'."))
    
    # ## some parsing
    # if args.target is not None:
    #     ## set args.indv to empty list so it doesn't break execute_homologue
    #     ##  note: execute_homologue is used to retrieve reference genes if requested by user
    #     args.indv = []


## TODO: make sure all tmp files are removed if not using --keep-on-crash (2021/05/12)
@app_sub.command("homologue")
@app_sub.command("homolog")
def homologue(
        
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
        # indv: Optional[List[IndvGenomesAll]] = typer.Option(*params.indv(), **params.indv.options,
        #                                                     callback = split_callback_list),
        target: Path = typer.Option(*params.target(), **params.target.options,
                                    **oparams.file_valid),
        query: Optional[List[Path]] = typer.Option(*params.query(), **params.query.options,
                                                   **oparams.file_valid),
        domain: str = typer.Option(*params.domain(), **params.domain.options, callback = domain_callback),
        minid: float = typer.Option(*params.minid(), **params.minid.options, callback = non_negative_callback),
        minlen: int = typer.Option(*params.minlen(), **params.minlen.options, callback = positive_callback),
        mincdslen: int = typer.Option(*params.mincdslen(), **params.mincdslen.options,
                                      callback = positive_callback),
        check_recip: bool = typer.Option(*params.check_recip(), **params.check_recip.options),
        relax_recip: bool = typer.Option(*params.relax_recip(), **params.relax_recip.options),
        merge_within: int = typer.Option(*params.merge_within(), **params.merge_within.options,
                                         callback = non_negative_callback),
        check_id_premerge: bool = typer.Option(*params.check_id_premerge(), **params.check_id_premerge.options),
        cluster_set: str = typer.Option(*params.cluster_set(), **params.cluster_set.options,
                                        **oparams.file_valid, is_eager = True,
                                        callback = cluster_set_callback),
        genome_set: str = typer.Option(*params.genome_set(), **params.genome_set.options,
                                       **oparams.file_valid, is_eager = True,
                                       callback = genome_set_callback),
        # reference_annotation: Optional[List[str]] = typer.Option(*params.reference_annotation(),
        #                                                          **params.reference_annotation.options,
        #                                                          callback = reference_annotation_callback),
        reference_set: str = typer.Option(*params.reference_set(), **params.reference_set.options,
                                          **oparams.file_valid, is_eager = True,
                                          callback = reference_set_callback),
        reference: Optional[List[str]] = typer.Option(*params.reference(), **params.reference.options,
                                                      callback = split_callback_list),
        assembly: str = typer.Option(*params.assembly(), **params.assembly.options,
                                      callback = assembly_callback),
        annotation: str = typer.Option(*params.annotation(), **params.annotation.options,
                                       callback = annotation_callback),
        db: str = typer.Option(*params.db(), **params.db.options),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = attr_mod_callback),
        ext_genome: Optional[List[Path]] = typer.Option(*params.ext_genome(), **params.ext_genome.options,
                                                        **oparams.file_valid),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options,
                                                     **oparams.file_valid),
        ext_reference: Optional[List[Path]] = typer.Option(*params.ext_reference(),
                                                           **params.ext_reference.options,
                                                           **oparams.file_valid),
        ext_annotation: Optional[List[Path]] = typer.Option(*params.ext_annotation(),
                                                            **params.ext_annotation.options,
                                                            **oparams.file_valid),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        
        ## user lookups (this is right at the end because we might need the other args)
        ## genomes is not "eager" because we have to wait for --genome-lookup to be processed
        genomes: bool = typer.Option(*params.genomes(), **params.genomes.options,
                                     callback = make_genomes_callback("independent")),
        references: bool = typer.Option(*params.references(), **params.reference.options,
                                        callback = references_callback),
        ## clusters & members are not "eager" because we have to wait for --cluster-lookup to be processed
        clusters: bool = typer.Option(*params.clusters(), **params.clusters.options,
                                      callback = clusters_callback),
        members: str = typer.Option(*params.members(), **params.members.options,
                                    callback = members_callback) ):
    '''
    Discover homologues in non-reference genomes
    '''
    
    ## check validity of args
    args = Namespace(**locals())
    ## this one can't be parsed using alias_or_file_callback cuz the path isn't actually to a file but a prefix
    args.db = os.path.abspath(parse_lookup(args.db, params.rps_db_aliases, return_first = True))
    ## check arguments are in order
    check_reference_args(args)
    check_homologue_args(args, standalone = True)
    ## parse cluster --> gene
    gene_sets = parse_genes_from_args(args)
    ## move logfile
    config.logfile.move(config, args)
    config.logfile.args_expanded(args, params)
    # ## create LogFile obj
    # config.log_file = LogFile(config, args, params)
    
    ## set config values
    if prefix: config.prefix = prefix
    if directory: config.out_dir = directory
    
    ## extend genome if ext_genome and ext_cds are provided
    extend_genome(args, config)
    
    ## make reduced annotation
    make_reduced_ann(args, config, gene_sets)
    
    ## call execute_homologue.
    for prefix, genes in gene_sets.items():
        ## catch output
        output = execute_homologue(args, config, params, prefix, genes)
    
    config.resolve()
    return
    # pass ## execute find homologue function

@app_sub.command("grna")
def generate_grna(
        ## general options
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        rpsblast: str = typer.Option(*params.rpsblast(), **params.rpsblast.options),
        remote_rps: bool = typer.Option(*params.remote_rps(), **params.remote_rps.options),
        mafft: str = typer.Option(*params.mafft(), **params.mafft.options),
        thread: int = typer.Option(*params.thread(), **params.thread.options),
        # blastn: str = typer.Option(*params.blastn(), **params.blastn.options),
        
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
        cluster_set: str = typer.Option(*params.cluster_set(), **params.cluster_set.options,
                                        **oparams.file_valid, is_eager = True,
                                        callback = cluster_set_callback),
        domain: str = typer.Option(*params.domain(), **params.domain.options, callback = domain_callback),
        # reference_annotation: Optional[List[str]] = typer.Option(*params.reference_annotation(),
        #                                                          **params.reference_annotation.options,
        #                                                          callback = reference_annotation_callback),
        reference_set: str = typer.Option(*params.reference_set(), **params.reference_set.options,
                                          **oparams.file_valid, is_eager = True,
                                          callback = reference_set_callback),
        reference: Optional[List[str]] = typer.Option(*params.reference(), **params.reference.options,
                                                      callback = split_callback_list),
        assembly: str = typer.Option(*params.assembly(), **params.assembly.options,
                                      callback = assembly_callback),
        annotation: str = typer.Option(*params.annotation(), **params.annotation.options,
                                       callback = annotation_callback),
        db: str = typer.Option(*params.db(), **params.db.options),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = attr_mod_callback),
        ext_genome: Optional[List[Path]] = typer.Option(*params.ext_genome(), **params.ext_genome.options,
                                                        **oparams.file_valid),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options,
                                                     **oparams.file_valid),
        ext_reference: Optional[List[Path]] = typer.Option(*params.ext_reference(),
                                                           **params.ext_reference.options,
                                                           **oparams.file_valid),
        ext_annotation: Optional[List[Path]] = typer.Option(*params.ext_annotation(),
                                                            **params.ext_annotation.options,
                                                            **oparams.file_valid),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        
        ## basic filtering
        gc_min: float = typer.Option(*params.gc_min(), **params.gc_min.options,
                                     callback = zero_to_one_callback),
        gc_max: float = typer.Option(*params.gc_max(), **params.gc_max.options,
                                     callback = zero_to_one_callback),
        feature: Optional[List[str]] = typer.Option(*params.feature(), **params.feature.options,
                                                    callback = split_callback_list),
        
        ## ??
        placeholder: str = None):
    
    '''
    Generate gRNA from either user-provided FASTA file or reference genes
    '''
    
    typer.echo("generating grna")
    
    ## check validity of args
    args = Namespace(**locals())
    ## this one can't be parsed using alias_or_file_callback cuz the path isn't actually to a file but a prefix
    args.db = os.path.abspath(parse_lookup(args.db, params.rps_db_aliases, return_first = True))
    ## check arguments are in order
    check_reference_args(args)
    check_grna_args(args, standalone = True)
    ## parse cluster --> gene
    gene_sets = parse_genes_from_args(args)
    ## move log file
    config.logfile.move(config, args)
    config.logfile.args_expanded(args, params)
    
    print("post check:", vars(args))
    
    if prefix: config.prefix = prefix
    if directory: config.out_dir = directory
    
    ## extend genome if ext_genome and ext_cds are provided
    extend_genome(args, config)
    
    ## filter BED/GFF for relevant entries (reduces get_seq search time)
    make_reduced_ann(args, config, gene_sets)
    
    ## generate gRNA
    for prefix, genes in gene_sets.items():
        ## call execute_homologue.
        if args.target is not None:
            set_dir = config.mkfname(prefix)
            set_pref = prefix
            set_target = args.target
            set_ann = config.annotation_red
            set_aln = set_gene = set_cds = set_domain_gff_bed = None
            config.mkdir(set_dir)
        else:
            output_homologue = execute_homologue(args, config, params, prefix, genes)
            ## split into multiple lines because there are just too many variables lol
            set_dir, set_pref = output_homologue[:2]
            set_target, set_aln, set_gene, set_cds, set_ann, set_domain_gff_bed = output_homologue[2:]
        ## call execute_grna.
        output_grna = execute_grna(args, config, set_dir, set_pref, set_target)
        set_grna_all, set_map_all, set_grna_fasta_all = output_grna
        ## call execute_filter.
        set_grna_pass, set_map_pass = execute_filter(args, config, set_dir, set_pref, set_grna_all,
                                                     set_map_all, set_grna_fasta_all, set_aln,
                                                     set_target, set_ann,
                                                     domain_gff_bed = set_domain_gff_bed,
                                                     checks = (["GC"] if target is not None else
                                                               ["GC", "feature"]))
    config.resolve()
    pass ## execute grna generation function

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
        db: str = typer.Option(*params.db(), **params.db.options),
        
        ## input files
        mapping: Path = typer.Option(*params.mapping(), **params.mapping.options, **oparams.file_valid),
        grna: Path = typer.Option(*params.grna(), **params.grna.options, **oparams.file_valid),
        alignment: Path = typer.Option(*params.alignment(), **params.alignment.options, **oparams.file_valid),
        target: Path = typer.Option(*params.target(), **params.target.options, **oparams.file_valid),
        # reference_annotation: Optional[List[str]] = typer.Option(*params.reference_annotation(),
        #                                                          **params.reference_annotation.options,
        #                                                          callback = reference_annotation_callback),
        reference_set: str = typer.Option(*params.reference_set(), **params.reference_set.options,
                                          **oparams.file_valid, is_eager = True,
                                          callback = reference_set_callback),
        reference: Optional[List[str]] = typer.Option(*params.reference(), **params.reference.options,
                                                      callback = split_callback_list),
        assembly: str = typer.Option(*params.assembly(), **params.assembly.options,
                                      callback = assembly_callback),
        annotation: str = typer.Option(*params.annotation(), **params.annotation.options,
                                       callback = annotation_callback),
        ext_genome: Optional[List[Path]] = typer.Option(*params.ext_genome(), **params.ext_genome.options,
                                                        **oparams.file_valid),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options,
                                                     **oparams.file_valid),
        ext_reference: Optional[List[Path]] = typer.Option(*params.ext_reference(),
                                                           **params.ext_reference.options,
                                                           **oparams.file_valid),
        ext_annotation: Optional[List[Path]] = typer.Option(*params.ext_annotation(),
                                                            **params.ext_annotation.options,
                                                            **oparams.file_valid),
        background: Optional[List[Path]] = typer.Option(*params.background(), **params.background.options,
                                                        **oparams.file_valid),
        mask: Optional[List[Path]] = typer.Option(*params.mask(), **params.mask.options, **oparams.file_valid),
        
        ## format options
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = attr_mod_callback),
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
        
        ## background filter options
        screen_ref: bool = typer.Option(*params.screen_ref(), **params.screen_ref.options),
        # unmask_ref: bool = typer.Option(*params.unmask_ref(), **params.unmask_ref.options),
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
        domain: str = typer.Option(*params.domain(), **params.domain.options, callback = domain_callback),
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
        # check_id_premerge: bool = typer.Option(*params.check_id_premerge(),
        #                                        **params.check_id_premerge.options),
        ## off-target threshold options
        pam: str = typer.Option(*params.pam(), **params.pam.options),
        ot_pamless: bool = typer.Option(*params.ot_pamless(), **params.ot_pamless.options),
        ot_mismatch: int = typer.Option(*params.ot_mismatch(), **params.ot_mismatch.options,
                                        callback = non_negative_callback),
        ot_gap: int = typer.Option(*params.ot_gap(), **params.ot_gap.options, callback = non_negative_callback),
        
        ## exclude filter
        exclude: Path = typer.Option(*params.exclude(), **params.exclude.options, **oparams.file_valid),
        
        ## ??
        placeholder: str = None):
    
    '''
    Filter gRNA by checks
    '''
    
    typer.echo("filtering")
    
    ## check validity of args
    args = Namespace(**locals())
    ## check arguments are in order
    check_reference_args(args)
    check_filter_args(args, standalone = True)
    ## parse [mask|unmask]_cluster -> genes
    gene_sets = parse_genes_for_filter(args)
    ## move log file
    config.logfile.move(config, args)
    config.logfile.args_expanded(args, params)
    
    print("post check:", vars(args))
    
    if prefix: config.prefix = prefix
    if directory: config.out_dir = directory
    
    ## extend genome if ext_genome and ext_cds are provided
    extend_genome(args, config)
    
    ## filter BED/GFF for relevant entries (reduces get_seq search time)
    if tuple(itertools.chain(*gene_sets.values())):
        make_reduced_ann(args, config, gene_sets)
    
    ## get sequences to be masked
    to_mask = assign_alias(args.mask, lambda i: f"Usermask_{str(i+1).zfill(2)}")
    for prefix, genes in gene_sets.items():
        if not genes: continue
        output_homologue = execute_homologue(args, config, params, prefix, genes,
                                             indv = {"ref"}, for_masking = True)
        to_mask[prefix] = output_homologue[2]
        # if screen_ref and not unmask_ref:
        #     to_mask[f"{prefix}_ref"] = output_homologue[4]
        # if mask_homologue:
        #     to_mask[f"{prefix}_nonref"] = output_homologue[2]
    
    ## parse mapping file into gRNAHits object
    grna_hits = gRNAHits()
    grna_hits.parse_from_mapping(args.mapping, targets = args.target, version = None)
    ## rename gRNA according to FASTA file if provided
    if args.grna is not None: grna_hits.assign_gRNAseq_id(args.grna)
    ## reset checks if requested
    if args.reset_checks: grna_hits.clear_checks()
    
    ## call execute_filter
    # val_or_empty = lambda b, v: [v] if b else []
    # checks = list(itertools.chain(*[val_or_empty(b, v)
    #                                 for b, v in [(args.gc_check, "GC"),
    #                                              (args.feature_check, "feature"),
    #                                              (args.background_check, "background")]]))
    grna_pass, map_pass = execute_filter(args, config, directory, prefix, grna_hits,
                                         args.mapping, args.grna, args.alignment, args.target,
                                         config.annotation_red, fasta_exclude = args.exclude,
                                         # checks = checks,
                                         domain_gff_bed = None, to_mask = to_mask)
    config.resolve()
    return

@app_sub.command("minimumset")
def minimumset(
        
        ## data input
        mapping: Path = typer.Option(*params.mapping(), **oparams.file_valid),
        grna: Path = typer.Option(*params.grna(), *params.fasta.names, **oparams.file_valid),
        exclude: Path = typer.Option(*params.exclude(), **oparams.file_valid),
        
        ## output files: option A
        out_mapping: Path = typer.Option(*params.out_mapping()),
        out_fasta: Path = typer.Option(*params.out_fasta()),
        
        ## output files: option B
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        
        ## process options
        sets: int = typer.Option(*params.sets(), callback = positive_callback),
        sc_algorithm: SetCoverAlgo = typer.Option(*params.sc_algorithm()),
        
        ## flags
        auto: bool = typer.Option(*params.auto(), **params.auto.options),
        accept_invalid: bool = typer.Option(*params.accept_invalid(), **params.accept_invalid.options),
        accept_feature_unknown: bool = typer.Option(*params.accept_feature_unknown(),
                                                    **params.accept_feature_unknown.options) ):
    
    '''
    Generate minimum set(s) of gRNA required to cover all targets from mapping file and FASTA file of gRNA.
    Requires mapping file (generated by minorg's 'full' subcommand) and a FASTA file of gRNA sequences.
    gRNA sequences not present in the mapping file will be ignored.
    '''
    args = locals()
    ## move log file
    config.logfile.move(config, args)
    config.logfile.args_expanded(args, params)
    
    if prefix: config.prefix = prefix
    if directory: config.out_dir = directory
    
    if ( args.out_mapping is None or args.out_fasta is None ):
        typer.echo(f"Output files will be generated in '{config.out_dir}' with the prefix '{config.prefix}'.")
    
    ## TODO: write log file
    execute_minimum_set(args, config)
    
    # from minimum_set import get_minimum_sets_from_files_and_write
    # ## note that the directory passed to this function is the tmp directory
    # get_minimum_sets_from_files_and_write(mapping = args.mapping, fasta = args.grna, exclude = args.exclude,
    #                                       directory = config.directory, prefix = config.prefix,
    #                                       fout_mapping = args.out_mapping, fout_fasta = args.out_fasta,
    #                                       num_sets = args.sets, sc_algorithm = args.sc_algorithm,
    #                                       output_map_ver = 2, manual_check = (not args.auto),
    #                                       ignore_invalid = args.accept_invalid,
    #                                       accept_unknown_within_feature_status = args.accept_feature_unknown)
    
    # typer.echo(f"boo {sets}")
    # typer.echo(f"verbose: {config.verbose}")
    config.resolve() ## move tmp stuff to output dir
    
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
        # indv: Optional[List[IndvGenomesAll]] = typer.Option(*params.indv(), **params.indv.options,
        #                                                     callback = split_callback_list),
        target: Path = typer.Option(*params.target(), **params.target.options),
        query: Optional[List[Path]] = typer.Option(*params.query(), **params.query.options,
                                                   **oparams.file_valid),
        domain: str = typer.Option(*params.domain(), **params.domain.options, callback = domain_callback),
        minid: float = typer.Option(*params.minid(), **params.minid.options, callback = non_negative_callback),
        minlen: int = typer.Option(*params.minlen(), **params.minlen.options, callback = positive_callback),
        mincdslen: int = typer.Option(*params.mincdslen(), **params.mincdslen.options,
                                      callback = positive_callback),
        check_recip: bool = typer.Option(*params.check_recip(), **params.check_recip.options),
        relax_recip: bool = typer.Option(*params.relax_recip(), **params.relax_recip.options),
        merge_within: int = typer.Option(*params.merge_within(), **params.merge_within.options,
                                         callback = non_negative_callback),
        check_id_premerge: bool = typer.Option(*params.check_id_premerge(), **params.check_id_premerge.options),
        # reference_annotation: Optional[List[str]] = typer.Option(*params.reference_annotation(),
        #                                                          **params.reference_annotation.options,
        #                                                          callback = reference_annotation_callback),
        reference: Optional[List[str]] = typer.Option(*params.reference(), **params.reference.options,
                                                      callback = split_callback_list),
                                      # callback = make_alias_or_file_callback(params.reference_aliases,
                                      #                                        params.reference)),
        assembly: str = typer.Option(*params.assembly(), **params.assembly.options,
                                      callback = assembly_callback),
        annotation: str = typer.Option(*params.annotation(), **params.annotation.options,
                                       callback = annotation_callback),
        db: str = typer.Option(*params.db(), **params.db.options),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = attr_mod_callback),
        ext_genome: Optional[List[Path]] = typer.Option(*params.ext_genome(), **params.ext_genome.options,
                                                        **oparams.file_valid),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options,
                                                     **oparams.file_valid),
        ext_reference: Optional[List[Path]] = typer.Option(*params.ext_reference(),
                                                           **params.ext_reference.options,
                                                           **oparams.file_valid),
        ext_annotation: Optional[List[Path]] = typer.Option(*params.ext_annotation(),
                                                            **params.ext_annotation.options,
                                                            **oparams.file_valid),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        
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
        screen_ref: bool = typer.Option(*params.screen_ref(), **params.screen_ref.options),
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
        sc_algorithm: SetCoverAlgo = typer.Option(*params.sc_algorithm(), **params.sc_algorithm.options),
        auto: bool = typer.Option(*params.auto(), **params.auto.options),
        output_ver: int = typer.Option(*params.output_ver(), **params.output_ver.options),
        
        ## user lookups (this is right at the end because we might need the other args)
        ## genomes is not "eager" because we have to wait for --genome-lookup to be processed
        genomes: bool = typer.Option(*params.genomes(), **params.genomes.options,
                                     callback = make_genomes_callback("independent")),
        references: bool = typer.Option(*params.references(), **params.reference.options,
                                        callback = references_callback),
        ## clusters & members are not "eager" because we have to wait for --cluster-lookup to be processed
        clusters: bool = typer.Option(*params.clusters(), **params.clusters.options,
                                      callback = clusters_callback),
        members: str = typer.Option(*params.members(), **params.members.options,
                                    callback = members_callback),
        ## more user lookups, this time for the lookups themselves haha
        cluster_set: str = typer.Option(*params.cluster_set(), **params.cluster_set.options,
                                        **oparams.file_valid, is_eager = True,
                                        callback = cluster_set_callback),
        genome_set: str = typer.Option(*params.genome_set(), **params.genome_set.options,
                                       **oparams.file_valid, is_eager = True,
                                       callback = genome_set_callback),
        reference_set: str = typer.Option(*params.reference_set(), **params.reference_set.options,
                                          **oparams.file_valid, is_eager = True,
                                          callback = reference_set_callback)):
    
        # genome_set: bool = typer.Option(*params.genome_set(), **params.genome_set.options,
        #                                 callback = genome_sets_callback),
        # cluster_set: bool = typer.Option(*params.cluster_set(), **params.cluster_set.options,
        #                                  callback = cluster_sets_callback)):
    '''
    Executes commands homologue, grna, filter, and minimumset in sequence to
    generate minimum set(s) of gRNA required to cover all targets.
    '''
    
    ## check validity of args
    args = Namespace(**locals())
    ## this one can't be parsed using alias_or_file_callback cuz the path isn't actually to a file but a prefix
    args.db = parse_lookup(args.db, params.rps_db_aliases, return_first = True)
    if args.db is not None: args.db = os.path.abspath(args.db)
    ## check arguments are in order
    check_reference_args(args)
    check_homologue_args(args, standalone = False)
    check_grna_args(args, standalone = False)
    check_filter_args(args, standalone = False)
    ## parse cluster --> gene
    gene_sets = parse_genes_from_args(args)
    gene_sets_for_filter = parse_genes_for_filter(args, priority = "mask")
    ## move log file
    config.logfile.move(config, args)
    config.logfile.args_expanded([args, params])
    
    print("post check:", vars(args))
    
    if prefix: config.prefix = prefix
    if directory: config.out_dir = directory
    
    ## extend genome if ext_genome and ext_cds are provided
    extend_genome(args, config)
    
    ## filter BED/GFF for relevant entries (reduces get_seq search time)
    make_reduced_ann(args, config, gene_sets)
        
    ## iterate through all gene sets
    for prefix, genes in gene_sets.items():
        ## do the thing where we find gRNA.
        ## call execute_homologue.
        if args.target is not None:
            set_dir = config.mkfname(prefix)
            set_pref = prefix
            set_target = args.target
            set_aln = set_gene = set_cds = set_ann = set_domain_gff_bed = None
            config.mkdir(set_dir)
        else:
            ## TODO: i'm not sure what "cds = False" was for. To update execute_homologue if I figure it out
            output_homologue = execute_homologue(args, config, params, prefix, genes)
                                                 # cds = False)
            ## split into multiple lines because there are just too many variables lol
            set_dir, set_pref = output_homologue[:2]
            set_target, set_aln, set_gene, set_cds, set_ann, set_domain_gff_bed = output_homologue[2:]
        ## call execute_grna.
        output_grna = execute_grna(args, config, set_dir, set_pref, set_target)
        set_grna_all, set_map_all, set_grna_fasta_all = output_grna
        ## call execute_filter.
        if args.screen_ref and not args.unmask_ref:
            if set(args.indv) != {"ref"}:
                to_mask = {"targets": set_target, "reference": set_gene}
            else:
                to_mask = {"reference": set_gene}
        else:
            to_mask = {"targets": set_target}
        set_grna_pass, set_map_all = execute_filter(args, config, set_dir, set_pref, set_grna_all,
                                                    set_map_all, set_grna_fasta_all, set_aln,
                                                    to_mask, set_ann,
                                                    domain_gff_bed = set_domain_gff_bed)
        ## call execute_minimumset.
        set_map_final = '_'.join(set_map_all.split('_')[:-1]) + "_final.map"
        set_grna_fasta_final = '_'.join(set_grna_fasta_all.split('_')[:-1]) + "_final.fasta"
        execute_minimumset(args, config, prefix = prefix, directory = set_dir,
                           mapping = set_map_all, grna = set_grna_fasta_all,
                           fout_mapping = set_map_final, fout_fasta = set_grna_fasta_final)
    
    typer.echo("heh")
    config.resolve()
    return

@app_sub.command("test")
def test(suff: int = 1):
    typer.echo("test")
    fname = config.mkfname(f"bruh_{suff}.txt")
    open(fname, 'a').close()
    fname = config.mkfname("hoo boy", f"aiyo_{suff}.txt")
    config.reserve_fname(fname)
    typer.echo("trying to resolve")
    config.resolve()
    typer.echo("resolved")
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
        version: bool = typer.Option(*params.version(), callback = version_callback, is_eager = True)):
    
    '''
    For documentation I guess??
    '''
    print("in sub main")
    
    ## config
    if quiet: config.verbose = False
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
    
    config.keep_tmp = keep_all
    
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
            config.cleanup()
        
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
        config.subcmd = 'full' if not sub_cmd_args else sub_cmd_args[0]
        if keep_on_crash: config.keep_on_crash = True
        print(sys.argv)
        
        ## execute subcommand
        app_sub()
    
    return

# app.set_default_command(full) ## idk how click_default_group is supposed to work :/

# print(main.params)

# print("callback:", dir(app_main.registered_callback))
    
if __name__ == "__main__":
    
    try:
        config.logfile = MINORgLogger(level = logging_level)
        config.logfile.update_filename(config.mkfname(os.path.basename(config.logfile.filename)))
        print(config.logfile.filename)
        config.logfile.args(["raw", sys.argv])
    except Exception as e:
        print("Unable to instantiate logger.")
        config.cleanup()
        raise e
    
    ## try-except to handle cleanup of tmpdir upon crash etc.
    try:
        config.raw_args = sys.argv
        app_main()
    except SystemExit as e:
        if e.code == 0: ## click returns 0 upon successful execution; catch and ignore
            config.cleanup()
        else:
            config.cleanup()
            raise e
    except MessageError as e:
        print(e.message)
        config.cleanup()
    except Exception as e:
        typer.echo("problem!")
        config.cleanup()
        raise e

## testing: ../minorg.py -d . test --quiet
## testing: ../minorg.py --dir ./tx/ty test --suff 2


