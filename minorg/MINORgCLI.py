## imports
import os
import click
import typer
import itertools

from minorg.MINORg import (
    MINORg,
    PathHandler,
    parse_lookup,
)
# from minorg.filter_grna import make_exclude_function as mk_excl_fn
from minorg.constants import (
    INDV_GENOMES_ALL,
    INDV_GENOMES_NONE,
    INDV_GENOMES_REF,
    REFERENCED_ALL,
    REFERENCED_NONE
)
from minorg.exceptions import (
    MessageError,
    InputFormatError,
    InvalidPath,
    InvalidFile,
    UnreadableFile
)
from minorg.parse_config import (
    Param,
    parse_multiline_multikey_sdict
)
from minorg.grna import gRNAHits
from minorg.pam import PAM
from minorg.functions import (
    get_count_dict,
    fasta_to_dict
)

    
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

class MINORgCLI (MINORg):
    
    def __init__(self, config = None, verbose = False, keep_tmp = False, keep_on_crash = False):
        from argparse import Namespace
        super().__init__(tmp = True, keep_tmp = keep_tmp, keep_on_crash = keep_on_crash, config = config)
        self._raw_args = Namespace() ## hidden from user. Stores raw args for posterity/logging
        self._args = Namespace()
        # self.minorg = MINORg(verbose = verbose, config = config, keep_tmp = keep_tmp,
        #                      keep_on_crash = keep_on_crash) ## will not have prefix, directory, and thread
        self.config = config
        self._prefix = None
        ## sets (path to file)
        self.reference_sets = None
        self.cluster_sets = None
        self.genome_sets = None
        ## set aliases
        self.cluster_set = None
        self.reference_set = None
        self.genome_set = None
        ## aliases parsed from sets
        self.reference_aliases = {}
        self.cluster_aliases = {}
        self.genome_aliases = {}
        ## others
        self.gene_sets = {}
        self.mask_gene_sets = {}
        self.background_check = True
        self.feature_check = True
        self.gc_check = True
    
    @property
    def args(self): return self._args
    @args.setter
    def args(self, args):
        self._raw_args = args
        self._args = args
    
    def copy_args(self, *argnames):
        for argname in argnames:
            try:
                setattr(self, argname, getattr(self.args, argname))
            except Exception as e:
                print(argname)
                raise e
        return
    
    ################################
    ##  UPDATE SUBSET_ANNOTATION  ##
    ##   TO USE SELF.GENE_SETS    ##
    ################################
    
    def subset_annotation(self, quiet = True, sort = True):
        if self.gene_sets:
            self._subset_annotation(ids = tuple(itertools.chain(*[x for x in self.gene_sets.values()
                                                                  if x is not None])),
                                    quiet = quiet, sort = sort)
        return
    
    #####################
    ##  ALIAS OR FILE  ##
    ##     PARSERS     ##
    #####################
    
    ## alias_or_file_parser: parse val as alias if valid else parse as path
    ## alias_or_file_callback --> MINORgCLI._alias_or_file_parser
    def _alias_or_file_parser(self, val: str, lookup: dict):
        mapped_val = parse_lookup(val, lookup, return_first = True)
        if mapped_val is not None: return valid_readable_file(mapped_val)
    ## make_alias_or_file_callback --> MINORgCLI._make_alias_or_file_parser
    def _make_alias_or_file_parser(self, lookup: dict, param = None,
                                   lookup_fail_msg = lambda x, y: (f"Unrecognised alias passed to {x}."
                                                                   f" Parsing '{y}' as path.")):
        lookup_fail_msg_updated = lambda x: lookup_fail_msg('/'.join(param.names), x)
        return (lambda val: self._alias_or_file_parser(val, lookup = lookup,
                                                       lookup_fail_msg = lookup_fail_msg_updated))
    
    #####################
    ##  SET CALLBACKS  ##
    #####################
    
    def _parse_file_set(self, val: str, lookup: dict, description: str,
                        cli_lookup_arg: str, lookup_fmt: str, param: Param,
                        num_fields = None,
                        parse_key = lambda k: k, parse_value = lambda v: v):
        val = parse_lookup(val, lookup, return_first = True)
        if val is not None:
            if os.path.exists(os.path.abspath(val)):
                set_fname = os.path.abspath(val)
                try:
                    with open(set_fname, 'r') as f:
                        set_aliases = {parse_key(k): parse_value(v) for k, v in
                                       parse_multiline_multikey_sdict(f.read(),
                                                                      kv_sep = '\t').items()}
                        if num_fields is not None:
                            set_aliases = {k: v + (num_fields - len(v))*['']
                                           for k, v in set_aliases.items()}
                except:
                    raise InputFormatError(error_src = "file", hint = f"File: {set_fname}")
            else:
                raise typer.BadParameter( ( f"Unrecognised {description} set lookup file alias or"
                                            f" non-existent file: {val}"
                                            f"\nFor a list of valid {description} set lookup file aliases"
                                            f" and their locations, execute 'minorg {cli_lookup_arg}'."
                                            "\nAlternatively, please provide a valid file of",
                                            f" {description} {lookup_fmt} mapping."
                                            f"\nHelp message for {param.long}: {param.help}") )
        else:
            set_fname = None
            set_aliases = {}
        return val, set_fname, set_aliases
    
    def reference_set_callback(self, val: str):
        val, set_fname, set_aliases = self._parse_file_set(val, lookup = self.params.reference_sets,
                                                           description = "reference",
                                                           cli_lookup_arg = "--references",
                                                           lookup_fmt = "alias-FASTA-GFF3-genetic code",
                                                           param = self.params.reference_set,
                                                           parse_value = (lambda v: v.split('\t')),
                                                           num_fields = 4)
        self.reference_set = set_fname
        self.reference_aliases = {**self.reference_aliases, **set_aliases}
        # self.args.reference_set = val
        return val
    
    def cluster_set_callback(self, val: str):
        val, set_fname, set_aliases = self._parse_file_set(val, lookup = self.params.cluster_sets,
                                                           description = "cluster",
                                                           cli_lookup_arg = "--clusters",
                                                           lookup_fmt = "alias-members",
                                                           param = self.params.cluster_set)
        self.cluster_set = set_fname
        self.cluster_aliases = {**self.cluster_aliases, **set_aliases}
        return val
    
    def genome_set_callback(self, val: str):
        val, set_fname, set_aliases = self._parse_file_set(val, lookup = self.params.genome_sets,
                                                           description = "genome",
                                                           cli_lookup_arg = "--genomes",
                                                           lookup_fmt = "alias-filename",
                                                           param = self.params.genome_set)
        self.genome_set = set_fname
        self.genome_aliases = {**self.genome_aliases, **set_aliases}
        return val
    
    
    ###########################
    ##  RETRIEVE FROM ALIAS  ##
    ###########################
    
    def get_cluster_genes(self, val: str):
        cluster_genes = self.cluster_aliases.get(val, None)
        if cluster_genes is None:
            raise click.UsageError( (f"'{val}' is not a valid cluster name in"
                                     f" lookup file '{self.cluster_set}'") )
        else:
            return tuple(cluster_genes.split(','))
    
    
    ########################
    ##  LOOKUP CALLBACKS  ##
    ########################
    
    def references_callback(self, value: bool):
        if value:
            if self.reference_set is None:
                typer.echo( (f"\nNo reference genomes have been defined."
                             " Please use '--reference-set' to specify a file containing"
                             " alias-FASTA-GFF3 mapping OR use '--assembly <FASTA>' and '--annotation <GFF3>'"
                             " to specify genome sequences and annotations respectively." ))
            else:
                typer.echo(f"\nValid genome aliases (defined in {self.reference_set}):\n")
                typer.echo("<semicolon-separated genome alias(es)>\t<FASTA file>\t<GFF3 file>\t<NCBI genetic code>")
                with open(self.reference_set, 'r') as f:
                    aliases = [f.read()]
                typer.echo('\n'.join(aliases) + '\n')
            self.cleanup()
            raise typer.Exit()
    
    def make_genomes_callback(self, mode: str = "independent"):
        def genomes_callback(value: bool):
            if value:
                if self.genome_set is None:
                    typer.echo(f"\nValid genome aliases (defined in {self.params.config_file}):\n")
                    aliases = ['\t'.join(kv) for kv in self.params.indv_genomes.items()]
                else:
                    typer.echo(f"\nValid genome aliases (defined in {self.genome_set}):\n")
                    with open(self.genome_set, 'r') as f:
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
                self.cleanup()
                raise typer.Exit()
        return genomes_callback
    
    def clusters_callback(self, value: bool):
        if value:
            if self.params.cluster_set is not None:
                typer.echo(f"\nThe following information is retrieved from {self.cluster_set}:\n")
                typer.echo("<semicolon-separated cluster alias(es)>\t<comma-separated cluster members>")
                with open(self.cluster_set, 'r') as f:
                    typer.echo(f.read())
            else:
                typer.echo( ("\nNo cluster set lookup file has been specified."
                             "\nYou may update the list of cluster set lookup files"
                             " and set a default set lookup file by modifying the relevant fields in"
                             f" the config file located at {self.params.config_file}"
                             " or provide a set lookup file for a single execution using --cluster-set.") )
            self.cleanup()
            raise typer.Exit()

    def members_callback(self, value: str):
        if value:
            if value not in self.cluster_aliases:
                if self.params.cluster_set is not None:
                    typer.echo( ( f"{value} is not a valid cluster alias in the following"
                                  f" cluster set lookup file: {self.cluster_set}" ) )
                else:
                    typer.echo( ("\nNo cluster set lookup file has been specified."
                                 "\nYou may update the list of cluster set files"
                                 " and set a default set lookup file by modifying the relevant fields"
                                 f" in the config file located at {self.params.config_file}"
                                 " or provide a set lookup file for a single execution using --cluster-set.") )
            else:
                typer.echo(f"\nThe following information is retrieved from {self.params.cluster_set}:\n")
                typer.echo(self.cluster_aliases[value] + '\n')
            self.cleanup()
            raise typer.Exit()
    
    
    ######################
    ##  ARGS CALLBACKS  ##
    ######################
        
    def attr_mod_callback(self, val):
        if isinstance(val, str):
            try:
                from regex import search
                return self.params.parse_attr_mod(search(r"""[^'"]+""", val).group(0))
            except:
                if val in self.params.attr_mod_presets:
                    raise InputFormatError(message = ( "Unable to parse the GFF attribute modification preset"
                                                       f" ({self.params.attr_mod_presets[val]}) retrieved"
                                                       " from the config file ({self.params.config_file})"
                                                       " via alias '{val}'." ),
                                           hint = self.params.attr_mod.help())
                else: raise InputFormatError(error_src = "for parameter --attr-mod",
                                             hint = self.params.attr_mod.help())
        else: return self.params.attr_mod.default
    
    def domain_callback(self, val: str):
        if val is not None:
            self.domain_name = val
            parsed_domain = self.params.parse_domain(val)
            if str(parsed_domain) == val and val != "gene":
                typer.echo(f"Parsing input to --domain ({val}) as PSSM-Id")
            self.domain = parsed_domain
            return parsed_domain
    
    def db_callback(self, val: str):
        if val is not None:
            parsed_db = os.path.abspath(parse_lookup(val, self.params.rps_db_aliases, return_first = True))
            self.db = parsed_db
            return parsed_db
    
    def genetic_code_callback(self, val):
        try: parsed_val = int(val)
        except ValueError: parsed_val = val
        self.genetic_code = parsed_val
        return parsed_val
    
    
    #########################
    ##  MULTI-ARGS CHECKS  ##
    #########################
    
    def reference_required(self, msg):
        if self.args.reference is None and (self.args.assembly is None or self.args.annotation is None):
            raise click.UsageError( ("'-r <reference genome alias>' OR"
                                     " '--assembly <reference assembly alias or path to FASTA file>"
                                     " --annotation <reference annotation alias or GFF3 file or"
                                     " BED file converted from GFF3 using gff2bed>'is required if using"
                                     f" {msg}.") )
    
    def check_reference_args(self, require = False):
        if self.args.assembly is not None and self.args.annotation is not None:
            if self.args.reference:
                typer.echo("As --assembly and --annotation are used, --reference will be ignored.")
            self.clear_reference()
            self.add_reference("Reference", os.path.abspath(assembly), os.path.abspath(ann_mapped),
                               genetic_code = self.args.genetic_code, attr_mod = self.args.attr_mod)
            self.args.assembly = assembly_mapped
            self.args.annotation = annotation_mapped
        elif self.args.reference:
            none_val = '-'
            valid_aliases(aliases = self.args.reference, lookup = self.reference_aliases,
                          none_value = none_val, all_value = "all",
                          param = self.params.reference, display_cmd = "--references",
                          additional_message = ("Alternatively, provide a FASTA file of a genome assembly using"
                                                " '--assembly <path to FASTA file>' and a GFF3 file of genome"
                                                " annotations using '--annotation <path to GFF3 file>'."))
            self.clear_reference()
            if set(self.args.reference) == {none_val}:
                self.args.reference = None
            else:
                for alias in self.args.reference:
                    fasta, ann, genetic_code, attr_mod = self.reference_aliases[str(alias)]
                    attr_mod = None if attr_mod == '' else attr_mod
                    self.add_reference(alias, fasta, ann, genetic_code = genetic_code,
                                       attr_mod = self.params.parse_attr_mod(attr_mod))
                    # self.extend_reference(alias, fasta, ann)
        elif require and self.args.target is not None:
            raise click.UsageError((f"Either '--reference <alias>' OR"
                                    " '--assembly <FASTA> --annotation <GFF3>' is required"
                                    " if '--target <FASTA>' is not used."))
        ## if --extend-gene or --extend-cds is provided, check that the other is also provided
        if sum(map(lambda x: x is None, [self.args.ext_gene, self.args.ext_cds])) == 1:
            raise click.UsageError( ("'--extend-gene <FASTA file>' and '--extend-cds <FASTA file>'"
                                     " should either be used together or not used at all.") )
        elif self.args.ext_gene and self.args.ext_cds:
            self.ext_gene = self.args.ext_gene
            self.ext_cds = self.args.ext_cds
            self.extend_reference()
        return
    
    def check_target_args(self, standalone = True):
        ## check if user has provided indv
        indv_provided = not (len(self.args.indv) == 1 and self.args.indv[0] == INDV_GENOMES_NONE)
        ## MUTUALLY EXCLUSIVE ARGS
        ## -q and -t are mutually exclusive
        if self.args.query and self.args.target is not None:
            raise click.UsageError("'-q <FASTA file>' and '-t <FASTA file>' are mutually exclusive.")
        ## -g and -c are mutually exclusive
        if self.args.gene is not None and self.args.cluster is not None:
            raise click.UsageError("'-g <gene(s)>' and '-c <cluster(s)>' are mutually exclusive.")
        ## if checking args for function 'target' and not 'full', then -q/-t are not compatible with -g/-c/-i
        if ( standalone and
             ( ( self.args.query or self.args.target is not None) and
               ( self.args.gene is not None or self.args.cluster is not None or indv_provided ) ) ):
            raise click.UsageError( ( "'-q <FASTA file>', '-t <FASTA file>', and"
                                      " '[ -g <gene(s)> or -c <cluster(s)> ] -i <individual(s)>'"
                                      " are mutually exclusive.") )
        ## REQUIRED ARGS
        ## check that one of -q, -t, -a|-i is provided
        if not self.args.query and self.args.target is None and not indv_provided:
            raise click.UsageError( ("One of the following is required:"
                                     " '-q <FASTA file>', '-t <FASTA file>', or '-i <individual(s)>'") )
        ## if -t is not used, check that -g/-c is provided
        if self.args.target is None:
            if self.args.gene is None and self.args.cluster is None:
                raise click.UsageError("'-g <gene(s)>' or '-c <cluster(s)>' is required if not using -t.'")
        ## if -g/-c is provided, check that --ref and --annotation are also provided
        if self.args.gene is not None or self.args.cluster is not None:
            self.reference_required("'-g <gene(s)>' or '-c <cluster(s)>'")
            if self.args.query is None and not indv_provided:
                raise click.UsageError( ("'-q <FASTA file>' or '-i <individual(s)>'"
                                         " is required if using '-g <gene(s)>' or '-c <cluster(s)>'.") )
        ## check that --mafft is provided
        if not standalone and self.args.mafft is None:
            raise click.UsageError( "Path to mafft executable is required: --mafft <path>" )
        ## check that --blastn is provided if using -q or '-i <not ref only>'
        if ( (not standalone or self.args.query or (indv_provided and set(self.args.indv) != {"ref"}))
             and self.args.blastn is None):
            raise click.UsageError( "Path to blastn executable is required: --blastn <path>" )
        ## check that --rpsblast is provided if using --domain
        if self.args.domain is not None and self.args.domain != "gene":
            if self.args.rpsblast is None:
                raise click.UsageError( ("'--rpsblast <path to rpsblast or rpsblast+ executable>'"
                                         " is required if using '--domain <domain alias or PSSM-Id>'") )
            if self.args.db is None:
                raise click.UsageError( ("'--db <alias of or path to local RPS-BLAST database"
                                         " OR name of remote RPS-BLAST database>'"
                                         " is required if using '--domain <domain alias or PSSM-Id>'") )
        ## VALID ALIASES (check -r, --annotation, --rps-db, --indv & -c & --attr-mod)
        ## not a callback because it requires config.cluster_aliases generated by cluster_set_callback
        if self.args.cluster is not None:
            valid_aliases(aliases = self.args.cluster, lookup = self.cluster_aliases,
                          param = self.params.cluster, display_cmd = "--clusters",
                          additional_message = ( "Alternatively, manually input the desired genes using"
                                                 " '-g <gene(s)>' or provide a different cluster set"
                                                 " lookup file using '--cluster-set <path to file>'."))
        ## ARGS DEPENDENT ON VALUE OF OTHER ARGS
        ## raise check_recip if relax_recip is raised
        if self.args.relax_recip:
            self.args.check_recip = True
        
        ## FORMAT ARGS FOR AND PASS ARGS TO MINORg
        ## self.query and self.query_references handled separately
        ##  generate query mapping (set self.query)
        query_map = []
        if self.args.query:
            query_map += [[i+1, str(query)] for i, query in enumerate(self.args.query)]
        if indv_provided:
            indvs_special = {INDV_GENOMES_NONE, INDV_GENOMES_REF, INDV_GENOMES_ALL}
            if INDV_GENOMES_REF in self.args.indv:
                self.query_reference = True
            if INDV_GENOMES_ALL in self.args.indv:
                indvs_genome = set(indv for indv in self.genome_aliases if indv not in indvs_special)
            else:
                indvs_genome = set(self.args.indv) - indvs_special
            indvs = (set(self.args.indv) - indvs_special).union(indvs_genome)
            if set(self.args.indv) != {INDV_GENOMES_REF}: self.args.indv = type(self.args.indv)(indvs)
            valid_aliases(aliases = list(indvs_genome), lookup = self.genome_aliases,
                          none_value = INDV_GENOMES_NONE, all_value = INDV_GENOMES_ALL,
                          param = self.params.indv, display_cmd = "--genomes",
                          additional_message = ( "Alternatively, provide a FASTA file of the genome in"
                                                 " which to query using '-q <path to FASTA file>'."))
            query_map += [[indv, self.genome_aliases[str(indv)]] for indv in indvs_genome]#  + \
                         # [[alias, self.reference_aliases[str(alias)][0]] for alias in indvs_ref]
        ## check that file names are unique
        ## - if multiple queries map to the same file, sort by query id and retain first entry
        ## - if same file provided to args.query and args.indv, priorities args.indv regardless of sort order
        query_count = get_count_dict([x[1] for x in query_map])
        query_multi = set(k for k, v in query_count.items() if v > 1)
        self.query = {k: v for k, v in query_map if v not in query_multi}
        for fname in query_multi:
            queries = sorted([x[0] for x in query_map if x[1] == fname])
            self.query[queries[0]] = fname
        return
    
    def check_grna_args(self, standalone = True):
        ## MUTUALLY EXCLUSIVE ARGS
        ## -g, -c, and -t are mutually exclusive
        if sum(map(lambda x: x is not None, [self.args.gene, self.args.cluster, self.args.target])) > 1:
            raise click.UsageError( ("'-g <gene(s)>', '-c <cluster(s)>', and '-t <path to FASTA file>'"
                                     " are mutually exclusive.") )
        ## REQUIRED ARGS
        ## check that --pam is provided
        if self.args.pam is None:
            raise click.UsageError( "'-p <PAM pattern>' is required." )
        ## check that --length is provided
        if self.args.length is None:
            raise click.UsageError( "'-l <gRNA length (bp)>' is required." )
        ## if -g/-c is provided, check that --ref and --annotation are also provided
        if (self.args.gene is not None or self.args.cluster is not None):
            self.reference_required("'-g <gene(s)>' or '-c <cluster(s)>'")
        ## check that --rpsblast is provided if using --domain
        if self.args.domain is not None and self.args.domain != "gene":
            if self.args.rpsblast is None:
                raise click.UsageError( ("'--rpsblast <path to rpsblast or rpsblast+ executable>'"
                                         " is required if using '--domain <domain alias or PSSM-Id>'") )
            if self.args.db is None:
                raise click.UsageError( ("'--db <alias of or path to local RPS-BLAST database"
                                         " OR name of remote RPS-BLAST database>'"
                                         " is required if using '--domain <domain alias or PSSM-Id>'") )
        ## VALID ALIASES (check -r, --annotation, --db, -c)
        ## not a callback because it requires config.cluster_aliases generated by cluster_set_callback
        if self.args.cluster is not None:
            valid_aliases(aliases = self.args.cluster, lookup = self.cluster_aliases,
                          param = self.params.cluster, display_cmd = "--clusters",
                          additional_message = ( "Alternatively, manually input the desired genes using"
                                                 " '-g <gene(s)>' or provide a different cluster set lookup file"
                                                 " using '--cluster-set <path to file>'."))
    
    def check_filter_args(self, standalone = True):
        if self.args.screen_reference or self.args.background: self.args.background_check = True
        if standalone:
            if self.args.check_all:
                self.args.gc_check = self.args.feature_check = self.args.background_check = True
            if not self.args.gc_check and not self.args.feature_check and not self.args.background_check:
                raise click.UsageError( ("At least one of the following is required:"
                                         " '--gc-check', '--feature-check', '--background-check'") )
        else:
            self.args.gc_check = self.args.feature_check = self.args.background_check = True
            if self.args.target:
                self.args.feature_check = self.args.background_check = False
            if self.args.skip_bg_check:
                self.args.background_check = False
        self.background_check = self.args.background_check
        self.feature_check = self.args.feature_check
        self.gc_check = self.args.gc_check
        ## REQUIRED ARGS
        if standalone:
            if self.args.map is None:
                raise click.UsageError("'--map <path to minorg .map file>' is required.")
            if self.args.background_check:
                ## if args.grna is not provided, create FASTA file from .map file
                if self.args.grna is None:
                    self.args.grna = self.mkfname("tmp_grna_all.fasta", prefix = False, tmp = True)
                    grnahits = gRNAHits()
                    grnahits.parse_from_mapping(self.args.map, version = None)
                    grnahits.write_fasta(self.args.grna, write_all = True)
            if self.args.in_place:
                self.args.out_map = self.args.map
        ## feature filter
        if self.args.feature_check:
            if not self.args.feature:
                raise click.UsageError("'--feature <feature>' is required if using '--filter-feature'.")
            if standalone and not self.args.alignment:
                raise click.UsageError("'--alignment <path to FASTA> is required if using '--filter-feature'.")
            self.reference_required("'--filter-feature'")
        ## background filter
        if self.args.background_check:
            if self.args.blastn is None:
                raise click.UsageError( "Path to blastn executable is required: --blastn <path>" )
            if self.args.unmask_ref:
                if not self.args.screen_reference:
                    raise click.UsageError( "'--unmask-ref' is only used with '--screen-reference'." )
                if self.args.mask_gene or self.args.mask_cluster:
                    raise click.UsageError( ("'--unmask-ref' should not be used with"
                                             " '--mask-gene' or '--mask-cluster'.") )
                if not standalone and "ref" in self.args.indv:
                    self.logfile.warning( ("'--unmask-ref' will be ignored if 'ref' is provided to '-i'"
                                           " and not excluded using '--ot-indv'. When 'ref' is provided"
                                           " to '-i', the reference genome(s) will be treated the same as"
                                           " any other genomes passed to '-i'.") )
            if self.args.screen_reference:
                if not (self.args.reference or (self.args.assembly and self.args.annotation)):
                    raise click.UsageError("'--reference <reference genome alias>' OR"
                                           " '--assembly <path to FASTA> --annotation <path to GFF3>'"
                                           " is required if using '--screen-reference'.")
            if not self.args.annotation:
                self.reference_required("'--unmask-gene', '--mask-cluster', or '--unmask-cluster'")
            if self.args.ot_mismatch is None: self.args.ot_mismatch = 1
            if self.args.ot_gap is None: self.args.ot_mismatch = 0
            if not self.args.ot_pamless and self.args.pam is None:
                raise click.UsageError("'--pam <PAM pattern>' is required if not using '--ot-pamless'.")
            if self.args.pam is not None:
                if standalone:
                    ## set minimum length for least restrictive regex
                    self.args.pam = PAM(pam = self.args.pam, gRNA_length = 1)
                else:
                    self.args.pam = PAM(pam = self.args.pam)
    
    
    ##########################
    ##  MULTI-ARGS PARSERS  ##
    ##########################
    
    def parse_cluster(self):
        '''
        Maps cluster aliases to genes and returns {<prefix>: [<genes>]}
        '''
        output = {}
        for cluster in self.args.cluster:
            output[f"{self.master_prefix}_{cluster}"] = tuple(self.get_cluster_genes(cluster))
        return output
    
    def parse_genes(self):
        '''
        To be called after 'gene' or 'cluster' have passed check_target_args checks
        Returns: {<prefix>: (gene1, gene2, gene3)}
        '''
        if self.args.target is not None:
            ## set None as val; execute_homologue will exit if it detects that val is None
            output = {self.master_prefix: None}
        elif self.args.cluster is None:
            output = {self.master_prefix: tuple(self.args.gene)}
        else:
            output = {**output, **self.parse_cluster()}
        self.gene_sets = output
    
    def parse_genes_for_filter(self, priority = None):
        '''
        To be called after 'mask', 'unmask', 'mask_cluster', 'unmask_cluster', and 'gene' 
            have passed check_filter_args checks
        Returns {<mask/unmask>: (gene1, gene2, gene3)}
        '''
        # genes = (set() if self.args.gene is None else set(self.args.gene))
        genes = set()
        if self.gene_sets: genes |= set(itertools.chain(*self.gene_sets.values()))
        # if self.cluster: genes |= set(itertools.chain(*self.parse_cluster.values()))
        mask_genes = set((tuple() if self.args.mask_gene is None else tuple(self.args.mask_gene)) + \
                         (tuple() if self.args.mask_cluster is None else
                          tuple(itertools.chain(*[self.get_cluster_genes(cluster)
                                                  for cluster in self.args.mask_cluster]))))
        unmask_genes = set((tuple() if self.args.unmask_gene is None else tuple(self.args.unmask_gene)) + \
                           (tuple() if self.args.unmask_cluster is None else
                            tuple(itertools.chain(*[self.get_cluster_genes(cluster)
                                                    for cluster in self.args.unmask_cluster]))))
        ## parse '.' for 'all genes passed to --gene/--cluster'
        def expand_genes(genes_to_expand):
            if REFERENCED_ALL in genes_to_expand and REFERENCED_NONE in genes_to_expand:
                raise MessageError( (f"ERROR: '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                     " are mutually exclusive for --mask-gene and --unmask-gene") )
            if REFERENCED_ALL in genes_to_expand:
                genes_to_expand -= {REFERENCED_ALL}
                genes_to_expand |= genes
            elif REFERENCED_NONE in genes_to_expand:
                genes_to_expand -= {REFERENCED_NONE}
                genes_to_expand -= genes
            return genes_to_expand
        mask_genes = expand_genes(mask_genes)
        unmask_genes = expand_genes(unmask_genes)
        overlap_genes = mask_genes.intersection(unmask_genes)
        ## check for overlaps between mask and unmask if priority is not set
        if priority is None and overlap_genes:
            raise MessageError( ("ERROR: The following genes were provided to both"
                                 " --mask-gene/--mask-cluster and --unmask-gene/--unmask-gene:"
                                 f" {','.join(sorted(overlap_genes))}") )
        ## generate output dictionary
        output = {}
        if priority == "mask":
            output["mask"] = tuple(mask_genes)
            output["unmask"] = tuple(unmask_genes - mask_genes)
        elif priority == "unmask":
            output["mask"] = tuple(mask_genes - unmask_genes)
            output["unmask"] = tuple(unmask_genes)
        else:
            output["mask"] = tuple(mask_genes)
            output["unmask"] = tuple(unmask_genes)
        self.mask_gene_sets = output
        self.genes = list(genes.union(mask_genes))
        return
        
    def parse_PAM(self):
        '''
        To be called after 'pam' and 'length' have passed check_grna_args checks OR
                     after 'pam' has passed check_filter_args (standalone) checks
        Returns: {<prefix>: (gene1, gene2, gene3)}
        '''
        self.PAM = PAM(pam = self.args.pam, gRNA_length = self.args.length)
        self.gRNA_length = self.args.length
    
    #######################
    ##  MASTER ARGPARSE  ##
    #######################
    
    ## implements appropriate parsers for each argument
    def parse_args(self, args, subcmd: str = None):
        self.args = args
        ## set output options
        if self.args.prefix:
            self.master_prefix = self.args.prefix
            self.prefix = self.args.prefix ## this may be updated by target OR full for each cluster
        if self.args.directory:
            self.mv(self.args.directory)
        ## check the appropriate args are used for each subcmd
        ## (no checks for subcmd minimumset)
        if subcmd in ["homologue", "homolog", "target"]:
            self.parse_target_args()
        elif subcmd in ["grna"]:
            self.parse_grna_args()
        elif subcmd in ["filter", "check"]:
            self.parse_filter_args()
        elif subcmd in ["minimumset"]:
            self.parse_minimumset_args()
        elif subcmd in ["full"]:
            self.parse_full_args()
        ## write args to logfile
        self.logfile.args_expanded([self.args, self.params])
        ## iterate through all args and apply appropriate parser based on arg name
        ##  to convert args to required format for MINORg object
        ## call further subcmd-specific argparsers
        pass
    
    def parse_target_args(self):
        self.check_reference_args()
        self.check_target_args(standalone = True)
        self.parse_genes()
        self.subset_annotation(quiet = True)
        ## PASS ARGS TO MINORg
        ## args handled by callbacks: domain
        ## args handled separately by parser: cluster, gene --> self.gene_sets
        ## args handled by check_reference_args: ext_gene, ext_cds, attr_mod
        ## args handled by check_target_args: query, indv, target --> self.query/self.target
        self.copy_args("blastn", "mafft", "rpsblast", "db", "remote_rps", "thread",
                       "target", "feature", "check_recip", "relax_recip",
                       "minid", "mincdslen", "check_id_before_merge", "merge_within")
        return
    
    def parse_grna_args(self):
        self.check_reference_args()
        self.check_grna_args(standalone = True)
        self.parse_genes()
        self.subset_annotation(quiet = True)
        ## PASS ARGS TO MINORg
        ## args handled by callbacks: domain
        ## args handled by check_reference_args: ext_gene, ext_cds, attr_mod
        ## args handled separately by parser: cluster, gene --> self.gene_sets
        self.copy_args("rpsblast", "db", "thread", "pam", "target", "length")
        self.grna_map = self.args.out_map
        self.grna_fasta = self.args.out_fasta
        return
    
    def parse_filter_args(self):
        self.check_reference_args()
        self.check_filter_args(standalone = True)
        self.parse_genes_for_filter()
        self.subset_annotation(quiet = True)
        ## FORMAT ARGS FOR AND PASS ARGS TO MINORg
        ## args handled by check_reference_args: ext_gene, ext_cds, attr_mod
        self.copy_args("pam", "blastn", "mafft", "thread",
                       "mask",
                       "gc_min", "gc_max",
                       "background", "ot_mismatch", "ot_gap", "ot_pamless", "screen_reference",
                       "feature", "max_insertion", "min_within_n", "min_within_fraction")
        ## args that require a little more parsing/have different names
        self.parse_grna_map_from_file(self.args.map)
        if self.args.reset_checks:
            self.grna_hits.clear_checks()
        self.grna_fasta = self.args.grna
        self.grna_map = self.args.out_map
        self.pass_fasta = self.args.out_fasta
        return
    
    def parse_minimumset_args(self):
        if ( self.args.out_map is None or self.args.out_fasta is None ):
            typer.echo(f"Output files will be generated in '{self.out_dir}' with the prefix '{self.prefix}'.")
        self.copy_args("sets", "auto")
        ## args that require a little more parsing/have different names
        self.parse_grna_map_from_file(self.args.map)
        self.grna_fasta = self.args.grna
        self.final_map = self.args.out_map
        self.final_fasta = self.args.out_fasta
        return
    
    def parse_full_args(self):
        self.check_reference_args()
        self.check_target_args(standalone = False)
        self.check_grna_args(standalone = False)
        self.check_filter_args(standalone = False)
        self.parse_genes()
        self.parse_genes_for_filter(priority = "mask")
        self.subset_annotation(quiet = True)
        ## PASS ARGS TO MINORg
        ## args handled by callbacks: domain
        ## args handled separately by parser: cluster, gene --> self.gene_sets
        ## args handled by check_reference_args: ext_gene, ext_cds, attr_mod
        ## args handled by check_target_args: query, indv, target --> self.query/self.target
        self.copy_args("blastn", "mafft", "rpsblast", ## executables
                       "db", "remote_rps", ## database
                       "thread",
                       "target", "feature", "check_recip", "relax_recip", ## target definition
                       "minid", "mincdslen", "check_id_before_merge", "merge_within",
                       "pam", "target", "length", ## gRNA generation
                       "gc_min", "gc_max", ## filtering options
                       "background", "ot_mismatch", "ot_gap", "ot_pamless", "screen_reference",
                       "feature", "max_insertion", "min_within_n", "min_within_fraction",
                       "sets", "auto") ## minimum set options
        return
    
    #######################
    ##  GENE SET MINORg  ##
    #######################

    def set_genes(self, gene_set_prefix):
        if gene_set_prefix not in self.gene_sets:
            print(f"{gene_set_prefix} not in gene_sets ({', '.join(self.gene_sets.keys())}).")
        self.genes = self.gene_sets[gene_set_prefix]
        self.prefix = gene_set_prefix
        
    ##############
    ##  SUBCMD  ##
    ##############
    
    def seq(self):
        for prefix in self.gene_sets:
            self.set_genes(prefix)
            super().seq(quiet = False)
        return
    
    def grna(self):
        for prefix in self.gene_sets:
            self.set_genes(prefix)
            super().grna()
            self.write_all_grna_map()
        return
    
    def filter(self):
        if self.background_check:
            print("Executing background check")
            self.filter_background()
            if self.keep_tmp:
                self.write_mask_report(self.mkfname("mask_report.tsv"))
        if self.feature_check:
            print("Executing feature check")
            self.filter_feature()
        if self.gc_check:
            print("Executing GC check")
            self.filter_gc()
        self.write_all_grna_map()
        self.write_pass_grna_fasta()
        return
    
    ## somewhat tested
    def minimumset(self):
        super().minimumset(report_full_path = False)
        return
    
    def full(self):
        for prefix in self.gene_sets:
            self.set_genes(prefix)
            super().seq(quiet = False)
            ## abort if no targets
            if len(fasta_to_dict(self.target)) == 0:
                raise MessageError("No targets found. Aborting programme.")
            super().grna()
            ## abort if no gRNA
            if len(self.grna_hits) == 0:
                raise MessageError("No gRNA found. Aborting programme.")
            self.filter()
            super().minimumset()
        return
            
