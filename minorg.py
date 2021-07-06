#!/usr/bin/env python3

## TODO: allow users to provide GFF instead of GFF -> BED
## TODO: 2021/06/10: implement --valid-indv/--valid-genome/--valid-acc, --valid-cluster to allow user to check alias validity (print location of FASTA file for --valid-genome, do member_callback for --valid-cluster) before running
## TODO: 2021/06/10: disable Enum for --indv. Since users are now allowed to submit custom genome lookup files using --genome-lookup, autofill shouldn't be used for --indv unless it can dynamically update to accommodate user's --genome-lookup input even before command submission (which is impossible lol). Basically, restructure the whole --indv thing to use the same kind of code as --cluster

import os
import sys
import click
import typer
import inspect
import itertools

from typing import List, Optional
from pathlib import Path
from argparse import Namespace
# from datetime import datetime

from scripts.get_ref import get_ref_by_genes_resolve
from functions import extract_features_and_subfeatures, get_count_dict

## import subcommand functions
from subcmd_homologue import execute_homologue
from subcmd_grna import execute_grna
from subcmd_filter import execute_filter
from subcmd_minimumset import execute_minimumset

from parse_config import (
    config, params, oparams,
    SetCoverAlgo, IndvGenomesAll, IndvGenomesAllClear, Lookup,
    LogFile,
    parse_multiline_multikey_sdict
)

# from subcmd_minorg import (
#     execute_homologue,
#     execute_grna,
#     execute_filter,
#     execute_minimumset
# )

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

__version_main__ = "2.1"
__version_full__ = "2.1"
__version_homologue__ = "1.0"
__version_grna__ = "unimplemented"
__version_filter__ = "unimplemented"
__version_minimumset__ = "1.1"

default_sub_cmd = "full"
app_main = typer.Typer()
app_sub = typer.Typer()

# ## import some namespace things
# params = Params(os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini")) ## parse params
# oparams = OptionParams() ## namespace for some pre-defined sets of parameter options
# config = Config(params)

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
            cluster_genes = config.cluster_aliases.get(cluster, None)
            if cluster_genes is None:
                raise click.UsageError( (f"'{args.cluster}' is not a valid cluster name in"
                                         f" lookup file '{config.cluster_set}'") )
            else:
                output[f"{args.prefix}_{cluster}"] = tuple(cluster_genes.split(','))
        return output
    return

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

## set ref: <reference genome> in params.indv_genomes
def reference_callback(val):
    alias_or_file_callback = make_alias_or_file_callback(params.reference_aliases, params.reference)
    mapped_val = alias_or_file_callback(val)
    if mapped_val is not None:
        params.indv_genomes["ref"] = mapped_val
    return mapped_val

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

def cluster_set_callback(val: str):
    val = parse_lookup(val, params.cluster_mapping, return_first = True)
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
    val = parse_lookup(val, params.genome_mapping, return_first = True)
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

def genome_sets_callback(value: bool):
    if value:
        lookup_sets_callback(config.genome_sets)

def cluster_sets_callback(value: bool):
    if value:
        lookup_sets_callback(config.genome_sets)

    
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

def make_reduced_bed(args, gene_sets):
    bed_red = None
    if args.gene or args.cluster:
        bed_red = config.mkfname("tmp_reduced.bed", tmp = True)
        config.bed_red = bed_red
        extract_features_and_subfeatures(args.bed, tuple(itertools.chain(*tuple(gene_sets.values()))),
                                         bed_red, quiet = True, fin_fmt = "BED", fout_fmt = "BED")
    return bed_red
    
def check_homologue_args(args, homologue_discovery_only = True):
    
    ## check if user has provided indv
    indv_provided = not (len(args.indv) == 1 and args.indv[0] is IndvGenomesAll.none)
    
    ## MUTUALLY EXCLUSIVE ARGS
    ## -q and -t are mutually exclusive
    if args.query and args.target is not None:
        print(args.query, args.target)
        raise click.UsageError("'-q <FASTA file>' and '-t <FASTA file>' are mutually exclusive.")
    ## -g and -c are mutually exclusive
    if args.gene is not None and args.cluster is not None:
        raise click.UsageError("'-g <gene(s)>' and '-c <cluster(s)>' are mutually exclusive.")
    ## if checking args for function 'homologue' and not 'full', then -q/-t are not compatible with -g/-c/-i
    if ( homologue_discovery_only and
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
    ## if -g/-c is provided, check that --ref and --bed are also provided
    if args.gene is not None or args.cluster is not None:
        if args.reference is None:
            raise click.UsageError( ("'-r <reference genome FASTA>' is required if using"
                                     " '-g <gene(s)>' or '-c <cluster(s)>'.") )
        if args.bed is None:
            raise click.UsageError( ("'--bed <BED file converted from GFF3 annotations using gff2bed>'"
                                     " is required if using '-g <gene(s)>' or '-c <cluster(s)>'.") )
        if args.query is None and not indv_provided:
            raise click.UsageError( ("'-q <FASTA file>' or '-i <individual(s)>'"
                                     " is required if using '-g <gene(s)>' or '-c <cluster(s)>'.") )
    ## if --extend-genome or --extend-cds is provided, check that the other is also provided
    if sum(map(lambda x: x is None, [args.ext_cds, args.ext_genome])) == 1:
        raise click.UsageError( ("'--extend-cds <FASTA file>' and '--extend-genome <FASTA file>'"
                                 " should either be used together or not used at all.") )
    
    ## VALID ALIASES (check -r, --bed, --rps-db, --indv & -c & --attr-mod)
    ## not a callback because it requires config.cluster_aliases generated by cluster_set_callback
    if args.cluster is not None:
        valid_aliases(aliases = args.cluster, lookup = config.cluster_aliases,
                      param = params.cluster, display_cmd = "--clusters",
                      additional_message = ( "Alternatively, manually input the desired genes using"
                                             " '-g <gene(s)>' or provide a different cluster set lookup file"
                                             " using '--cluster-set <path to file>'."))
    ## not a callback to standardise valid alias checking w/ args.cluster
    if indv_provided:
        valid_aliases(aliases = args.indv, lookup = params.indv_genomes,
                      none_value = IndvGenomesAll.none.value, all_value = IndvGenomesAll.all.value,
                      param = params.indv, display_cmd = "--genomes",
                      additional_message = ( "Alternatively, provide a FASTA file of the genome in"
                                             " which to query using '-q <path to FASTA file>'."))
    
    ## ARGS DEPENDENT ON VALUE OF OTHER ARGS
    ## raise check_recip if relax_recip is raised
    if args.relax_recip:
        args.check_recip = True
    
    ## generate query mapping
    query_map = []
    if args.query:
        query_map += [[i, query] for i, query in enumerate(args.query)]
    if indv_provided:
        if IndvGenomesAll.all in args.indv:
            indvs_special = {IndvGenomesAll.all, IndvGenomesAll.none}
            indvs_all = set(indv for indv in IndvGenomesAll if \
                            (indv not in indvs_special.union({IndvGenomesAll.ref})))
            indvs = (set(args.indv) - indvs_special).union(indvs_all) ## union in case user also specified 'ref'
            args.indv = type(args.indv)(indvs)
        print("chkpt1:", len(config.genome_aliases))
        query_map += [[indv, config.genome_aliases[indv]] for indv in args.indv]
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

## TODO: make sure all tmp files are removed if not using --keep-on-crash (2021/05/12)
@app_sub.command("homologue")
@app_sub.command("homolog")
def homologue(
        
        ## general options
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        
        ## target definition options
        gene: Optional[List[str]] = typer.Option(*params.gene(), **params.gene.options,
                                                 callback = split_callback_list),
        cluster: Optional[List[str]] = typer.Option(*params.cluster(), **params.cluster.options,
                                                    callback = split_callback_list),
        indv: Optional[List[IndvGenomesAll]] = typer.Option(*params.indv(), **params.indv.options,
                                                            callback = split_callback_list),
        target: Path = typer.Option(*params.target(), **params.target.options),
        query: Optional[List[Path]] = typer.Option(*params.query(), **params.query.options),
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
        reference: str = typer.Option(*params.reference(), **params.reference.options,
                                      callback = reference_callback),
        bed: str = typer.Option(*params.bed(), **params.bed.options,
                                callback = make_alias_or_file_callback(params.gff_bed_aliases, params.bed)),
        db: str = typer.Option(*params.db(), **params.db.options),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = attr_mod_callback),
        ext_genome: Optional[List[Path]] = typer.Option(*params.ext_genome(), **params.ext_genome.options),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        rpsblast: str = typer.Option(*params.rpsblast(), **params.rpsblast.options),
        blastn: str = typer.Option(*params.blastn(), **params.blastn.options),
        
        ## user lookups (this is right at the end because we might need the other args)
        ## genomes is not "eager" because we have to wait for --genome-lookup to be processed
        genomes: bool = typer.Option(*params.genomes(), **params.genomes.options,
                                     callback = make_genomes_callback("independent")),
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
    check_homologue_args(args, homologue_discovery_only = True)
    ## parse cluster --> gene
    gene_sets = parse_genes_from_args(args)
    ## create LogFile obj
    config.log_file = LogFile(config, args, params)
    
    ## set config values
    if prefix: config.prefix = prefix
    if directory: config.out_dir = directory
    ## make reduced bed
    make_reduced_bed(args, gene_sets)
    
    ## call execute_homologue.
    for prefix, genes in gene_sets.items():
        execute_homologue(prefix, genes, args, config, params)
    
    config.resolve()
    return
    # pass ## execute find homologue function

@app_sub.command("grna")
def generate_grna(
        ## general options
        directory: Path = typer.Option(*params.directory(), **params.directory.options, **oparams.dir_new),
        prefix: str = typer.Option(*params.prefix(), **params.prefix.options),
        blastn: str = typer.Option(*params.blastn(), **params.blastn.options),
        
        ## gRNA options
        pam: str = typer.Option(*params.pam(), **params.pam.options),
        length: str = typer.Option(*params.length(), **params.length.options),
        span_junction: bool = typer.Option(*params.span_junction(), **params.span_junction.options),
        
        ## target definition options
        target: Path = typer.Option(*params.target(), **params.target.options),
        
        ## target definition options if finding gRNA in reference genes
        gene: Optional[List[str]] = typer.Option(*params.gene(), **params.gene.options,
                                                 callback = split_callback_list),
        cluster: Optional[List[str]] = typer.Option(*params.cluster(), **params.cluster.options,
                                                    callback = split_callback_list),
        cluster_set: str = typer.Option(*params.cluster_set(), **params.cluster_set.options,
                                        **oparams.file_valid, is_eager = True,
                                        callback = cluster_set_callback),
        reference: str = typer.Option(*params.reference(), **params.reference.options,
                                      callback = reference_callback),
        bed: str = typer.Option(*params.bed(), **params.bed.options,
                                callback = make_alias_or_file_callback(params.gff_bed_aliases, params.bed)),
        db: str = typer.Option(*params.db(), **params.db.options),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = attr_mod_callback),
        ext_genome: Optional[List[Path]] = typer.Option(*params.ext_genome(), **params.ext_genome.options),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        rpsblast: Path = typer.Option(*params.rpsblast(), **params.rpsblast.options, **oparams.file_valid),
        
        ## ??
        placeholder: str = None):
    
    typer.echo("generating grna")
    config.resolve()
    pass ## execute grna generation function

@app_sub.command("filter")
def filter_grna():
    typer.echo("filtering")
    config.resolve()
    pass ## exceute filtering function

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
        accept_cds_unknown: bool = typer.Option(*params.accept_cds_unknown(),
                                                **params.accept_cds_unknown.options) ):
    
    '''
    Generate minimum set(s) of gRNA required to cover all targets from mapping file and FASTA file of gRNA.
    Requires mapping file (generated by minorg's 'full' subcommand) and a FASTA file of gRNA sequences.
    gRNA sequences not present in the mapping file will be ignored.
    '''
    
    if prefix: config.prefix = prefix
    if directory: config.out_dir = directory
    
    if ( out_mapping is None or out_fasta is None ):
        typer.echo(f"Output files will be generated in '{config.out_dir}' with the prefix '{config.prefix}'.")
    
    ## TODO: write log file
    
    from scripts.get_minimum_set import get_minimum_sets_from_files_and_write
    ## note that the directory passed to this function is the tmp directory
    get_minimum_sets_from_files_and_write(mapping = mapping, fasta = grna, exclude = exclude,
                                          directory = config.directory, prefix = config.prefix,
                                          fout_mapping = out_mapping, fout_fasta = out_fasta,
                                          num_sets = sets, sc_algorithm = sc_algorithm,
                                          output_map_ver = 2, manual_check = (not auto),
                                          ignore_invalid = accept_invalid,
                                          accept_unknown_within_cds_status = accept_cds_unknown)
    
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
        
        ## target definition options
        gene: Optional[List[str]] = typer.Option(*params.gene(), **params.gene.options,
                                                 callback = split_callback_list),
        cluster: Optional[List[str]] = typer.Option(*params.cluster(), **params.cluster.options,
                                                    callback = split_callback_list),
        # indv: Optional[List[str]] = typer.Option(*params.indv(), **params.indv.options,
        #                                          callback = split_callback_list),
        indv: Optional[List[IndvGenomesAll]] = typer.Option(*params.indv(), **params.indv.options,
                                                            callback = split_callback_list),
        target: Path = typer.Option(*params.target(), **params.target.options),
        query: Optional[List[Path]] = typer.Option(*params.query(), **params.query.options),
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
        reference: str = typer.Option(*params.reference(), **params.reference.options,
                                      callback = reference_callback),
                                      # callback = make_alias_or_file_callback(params.reference_aliases,
                                      #                                        params.reference)),
        bed: str = typer.Option(*params.bed(), **params.bed.options,
                                callback = make_alias_or_file_callback(params.gff_bed_aliases, params.bed)),
        db: str = typer.Option(*params.db(), **params.db.options),
        attr_mod: str = typer.Option(*params.attr_mod(), **params.attr_mod.options,
                                     callback = attr_mod_callback),
        ext_genome: Optional[List[Path]] = typer.Option(*params.ext_genome(), **params.ext_genome.options),
        ext_cds: Optional[List[Path]] = typer.Option(*params.ext_cds(), **params.ext_cds.options),
        sep: str = typer.Option(*params.sep(), **params.sep.options),
        rpsblast: Path = typer.Option(*params.rpsblast(), **params.rpsblast.options, **oparams.file_valid),
        
        ## gRNA generation options
        pam: str = typer.Option(*params.pam(), **params.pam.options),
        length: str = typer.Option(*params.length(), **params.length.options),
        
        ## filter gRNA options
        mismatch: int = typer.Option(*params.mismatch(), **params.mismatch.options,
                                     callback = non_negative_callback),
        gap: int = typer.Option(*params.gap(), **params.gap.options, callback = non_negative_callback),
        exclude: Path = typer.Option(*params.exclude(), **params.exclude.options, **oparams.file_valid),
        accept_invalid: bool = typer.Option(*params.accept_invalid(), **params.accept_invalid.options),
        accept_cds_unknown: bool = typer.Option(*params.accept_cds_unknown(),
                                                **params.accept_cds_unknown.options),
        skip_bg_check: bool = typer.Option(*params.skip_bg_check(), **params.skip_bg_check.options),
        screen_ref: bool = typer.Option(*params.screen_ref(), **params.screen_ref.options),
        unmask_ref: bool = typer.Option(*params.unmask_ref(), **params.unmask_ref.options),
        
        ## gRNA minimum set options
        sets: int = typer.Option(*params.sets(), **params.sets.options, callback = positive_callback),
        sc_algorithm: SetCoverAlgo = typer.Option(*params.sc_algorithm(), **params.sc_algorithm.options),
        auto: bool = typer.Option(*params.auto(), **params.auto.options),
        output_ver: int = typer.Option(*params.output_ver(), **params.output_ver.options),
        
        ## user lookups (this is right at the end because we might need the other args)
        ## genomes is not "eager" because we have to wait for --genome-lookup to be processed
        genomes: bool = typer.Option(*params.genomes(), **params.genomes.options,
                                     callback = make_genomes_callback("independent")),
        ## clusters & members are not "eager" because we have to wait for --cluster-lookup to be processed
        clusters: bool = typer.Option(*params.clusters(), **params.clusters.options,
                                      callback = clusters_callback),
        members: str = typer.Option(*params.members(), **params.members.options,
                                    callback = members_callback),
        ## more user lookups, this time for the lookups themselves haha
        genome_sets: bool = typer.Option(*params.genome_sets(), **params.genome_sets.options,
                                         callback = genome_sets_callback),
        cluster_sets: bool = typer.Option(*params.cluster_sets(), **params.cluster_sets.options,
                                          callback = cluster_sets_callback):
    '''
    Executes commands homologue, grna, filter, and minimumset in sequence to
    generate minimum set(s) of gRNA required to cover all targets.
    '''
    
    ## check validity of args
    args = Namespace(**locals())
    ## this one can't be parsed using alias_or_file_callback cuz the path isn't actually to a file but a prefix
    args.db = os.path.abspath(parse_lookup(args.db, params.rps_db_aliases, return_first = True))
    ## check arguments are in order
    check_homologue_args(args, homologue_discovery_only = False)
    ## parse cluster --> gene
    gene_sets = parse_genes_from_args(args)
    ## TODO: write LOG file
    config.log_file = LogFile(config, args, params)
    
    print("post check:", vars(args))
    
    if prefix: config.prefix = prefix
    if directory: config.out_dir = directory
    
    ## filter BED/GFF for relevant entries (reduces get_seq search time)
    make_reduced_bed(args, gene_sets)
    # if args.gene or args.cluster:
    #     bed_red = config.mkfname("tmp_reduced.bed")
    #     config.bed_red = bed_red
    #     extract_features_and_subfeatures(args.bed, tuple(itertools.chain(*tuple(gene_sets.values()))),
    #                                      bed_red, quiet = True, fin_fmt = "BED", fout_fmt = "BED")
    #     config.tmp_files.append(bed_red)
    
    ## generate wrapper for get_seq for common args
    def get_seq_ref_genes(genes, feature, fasta_out, out_dir, **kwargs):
        get_ref_by_genes_resolve(genes = genes, feature = feature,
                                 out_dir = out_dir, fout = fasta_out, no_bed = True, ## output options
                                 ref_fasta_files = reference, bed = bed_red, ## reference + annotation
                                 attribute_mod = attr_mod, by_gene = True, **kwargs)
    
    ## iterate through all gene sets
    ## TODO: handle situations when args.target is defined
    print(gene_sets)
    for prefix, genes in gene_sets.items():
        ## do the thing where we find gRNA.
        ## call execute_homologue.
        if args.target is not None:
            set_dir = config.mkfname(prefix)
            set_pref = prefix
            set_target = args.target
            set_aln = set_complete = set_cds = None
        else:
            ## TODO: i'm not sure what "cds = False" was for. To update execute_homologue if I figure it out
            set_dir, set_pref, set_target, set_aln, set_complete, set_cds = execute_homologue(prefix, genes,
                                                                                              args, config,
                                                                                              params)
                                                                                              # cds = False)
        ## call execute_grna.
        set_grna_all, set_grna_all = execute_grna(set_dir, set_pref, set_target, args, config)
        # ## call execute_filter.
        # set_gRNA, set_grna_pass, set_map_pass = execute_filter(set_dir, set_pref, set_grna_all, set_map_all,
        #                                                        set_aln, set_complete, set_cds, args, config)
        # ## call execute_minimumset.
        # execute_minimumset(set_grna_pass, set_map_pass, args, config)
        pass
    
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
         _help: bool = typer.Option(*params.help(None))): ## catch help
    
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
            ## - if no subcommand provided, it defaults to sub_main's help page
            ## - if subcommand provided, it uses subcommand's help page
            sub_cmd_args = sub_cmd_args + ["--help"]
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
