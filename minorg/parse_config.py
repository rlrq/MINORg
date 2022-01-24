import os
import typer
import shutil
import tempfile
import configparser

from enum import Enum
from argparse import Namespace

from minorg.functions import IndexedFasta
from minorg.constants import (
    INDV_GENOMES_ALL,
    INDV_GENOMES_NONE,
    INDV_GENOMES_REF,
    REFERENCED_ALL,
    REFERENCED_NONE
)
from minorg.reference import AnnotatedFasta
# from minorg.annotation import reduce_ann

## ensure that 'config = Config(params, keep_on_crash = True)' is updated to 'config = Config(params, keep_on_crash = False)' when released (keep_on_crash set to False)
## TODO: ensure genetic code specified in reference_set files are used

def get_val_none(val, d: dict, none = None):
    '''
    Attempt to retrieve value of key-value pair in 'd' if val in d.keys().
    Failing that, if val is not '' (empty string), return val.
    Else, return None
    '''
    val = d.get(val, val)
    return val if val else none

def get_val_default(val, default = None, coerce = None):
    '''
    If val is empty (None or empty iterable), return default.
    '''
    if coerce is not None and val is not None:
        try:
            val = coerce(val)
        except ValueError:
            pass
    if val == 0: return val
    elif val: return val
    else: return default

def inverse_dict(d: dict):
    '''
    Converts {<k1>:<v1>, <k2>:<v1>, <k3>:<v2>} into {<v1>:[<k1>,<k2>], <v2>:[<k3>]}
    '''
    return {v1: [k for k, v2 in d.items() if v2 == v1] for v1 in d.values()}

def parse_sep_sdict(s: str, item_sep: str = ';', kv_sep = ':'):
    '''
    Converts raw '<key1><kv_sep><v1>,<v2>,...,<v5><item_sep><key2><kv_sep><v6>,...,<v9>' string into 
    {'<key1>':'<v1>,<v2>,...,<v5>', '<key2>':'<v6>,...,<v9>'} dict.
    '''
    return dict( ( item.split(kv_sep)[0], kv_sep.join(item.split(kv_sep)[1:]) )
                 for item in s.split(item_sep) if item )

def parse_multiline_sdict(s: str, kv_sep: str = ':'):
    '''
    Converts raw '<key1><kv_sep><v1>,<v2>,...,<v5>\n<key2><kv_sep><v6>,...,<v9>' string into 
    {'<key1>':'<v1>,<v2>,...,<v5>', '<key2>':'<v6>,...,<v9>'} dict.
    '''
    return parse_sep_sdict(s, kv_sep = kv_sep, item_sep = '\n')

def parse_multiline_inv_sdict(s: str, v_sep: str = ',', **kwargs):
    '''
    Converts raw '<key>:<v1>,<v2>,...,<v5>' string into {'<v1>':'<key>','<v2>':'<key>', ...,'<v5>':'<key>'} dict.
    '''
    inv_dict = parse_multiline_sdict(s, **kwargs)
    output_dict = {}
    for k, v in inv_dict.items():
        for value in v.split(v_sep):
            output_dict[value] = k
    return output_dict

def parse_multiline_multikey_sdict(s: str, k_sep: str = ';', **kwargs):
    multikey_dict = parse_multiline_sdict(s, **kwargs)
    output = {}
    for k, v in multikey_dict.items():
        output = {**output, **{k2: v for k2 in k.split(k_sep)}}
    return output

def parse_attr_mod_sdict(s: str, attr_sep: str = ',', feature_sep: str = ';', fa_sep: str = ':',
                         alias_sep: str = '='):
    """
    Parse attribute modifications from string to dictionary.

    Arguments:
        s (str): required, string of attribute modifications in format
            '<feature type>:<standard attribute field name>=<nonstandard attribute field name>,<standard attribute field name>=<nonstandard attribute field name>;<feature type>:<standard attribute field name>=<nonstandard attribute field name>'
            (e.g. 'mRNA:Parent=Locus_id')
        attr_sep (str): delimiter for attribute modifications of same feature type (default=',')
        feature_sep (str): delimiter for feature types (default=';')
        fa_sep (str): delimiter between feature type and attribute modifications (default=':')
        alias_sep (str): delimiter between standard attribute field name and non-standard attribute field name
            (default='=')

    Returns
    -------
    dict
        parsed attribute modifications in format
            {<feature>: {<standard attribute field name>: <nonstandard attribute field name>}}
    """
    if not s: return {}
    parse_attr_mod = lambda mappings: dict(mapping.split(alias_sep) for mapping in mappings.split(attr_sep))
    feature_dict = parse_sep_sdict(s, kv_sep = fa_sep, item_sep = feature_sep)
    return {feature: parse_attr_mod(mappings) for feature, mappings in feature_dict.items()}

def generate_autocompletion(param_name, valid_completion_items):    
    def complete(ctx: typer.Context, incomplete: str):
        values = ctx.params.get(param_name) or []
        for value, help_text in valid_completion_items:
            if value.startswith(incompelte) and value not in values:
                yield(value, help_text)
        return
    return complete

def mv_dir_overwrite(src_dir, dst_dir):
    for root, dirs, files in list(os.walk(src_dir))[::-1]:
        out_dir = os.path.join(dst_dir, os.path.relpath(root, src_dir))
        os.makedirs(out_dir, exist_ok = True)
        for fname in files:
            shutil.move(os.path.join(root, fname), os.path.join(out_dir, fname))
        for dirname in dirs:
            os.rmdir(os.path.join(root, dirname))
    return

## parameters defaults
class Param():
    ## names: some (non-generator) iterable
    def __init__(self, default, *names,
                 false_true = ("False", "True"), show_any_default = True,
                 description = None, alias_value_description = None,
                 help_subcmd = {}, **kwargs):
        self.default = default
        self.names = names
        self.false_true = false_true
        self.description = description
        self.alias_value_description = alias_value_description
        self.options = kwargs
        self.help_subcmd = help_subcmd
        if not show_any_default:
            self.options["show_default"] = False
        elif not kwargs.get("show_default", True):
            self.options["help"] = '  '.join([kwargs.get("help", ''), self.custom_default()]).lstrip()
            self.help_subcmd = {k: '  '.join([v, self.custom_default()]).lstrip()
                                for k, v in self.help_subcmd}
    
    ## generate [default, *names] list for use in typer.Options
    ## - if called w/ a pos arg, 1st pos arg will replace self.default as default value
    def __call__(self, *default):
        default = self.default if (isinstance(default, str) or not default) else default[0]
        return [default] + list(self.names)
    
    ## returns help parameter if it exists
    def help(self, subcmd = None):
        return self.help_subcmd.get(subcmd, self.options.get("help", self.custom_default()))
    
    ## some function that converts a "True" to "accept" (or "reject") when displaying default using help param
    ## - use false_true (see __init__) to do this
    ## TODO: integrate this into the help display
    def custom_default(self):
        if type(self.default) is bool:
            return f"[default: {self.false_true[int(self.default)]}]"
        else:
            return f"[default: {self.default}]"
    
    def format_log(self, val):
        if isinstance(val, list) or isinstance(val, tuple):
            val = ','.join(val)
        output = f"{'|'.join(self.names)}:\t{val}"
        if val == self.default:
            output += " (default)"
        return output
    
    ## get shortest name or longest name
    @property
    def short(self):
        shortest_len = min(len(name) for name in self.names)
        shortest_names = [name for name in self.names if len(name) == shortest_len]
        ## tie-breaker: return first result
        return shortest_names[0]
    @property
    def long(self):
        longest_len = max(len(name) for name in self.names)
        longest_names = [name for name in self.names if len(name) == longest_len]
        ## tie-breaker: return first result
        return longest_names[0]

## parameter names namespace
class Params():
    def __init__(self, config_file):
        
        ## read config file
        if config_file and os.path.exists(config_file):
            self.config_file = config_file
            conf = configparser.ConfigParser()
            conf.read(self.config_file)
            def conf_get(*args, type = str, **kwargs):
                if type is str: get = conf.get
                elif type is bool: get = conf.getboolean
                elif type is int: get = conf.getint
                elif type is float: get = conf.getfloat
                return get(*args, **kwargs)
        elif config_file and not os.path.exists(config_file):
            raise InvalidPath(config_file)
        else:
            def conf_get(*args, **kwargs):
                return None
        
        ## several lookups (parsed here because here is where the config file is parsed haha)
        section_lookup = "lookup"
        # get_lookup = lambda x: conf.get(section_lookup, x) ## for strings
        get_lookup = lambda x: get_val_default(conf_get(section_lookup, x), default = '')
        self.domain_mapping = '' ## TODO, analogous to genome_mapping and cluster_mapping
        self.domain_aliases = parse_multiline_inv_sdict(get_lookup("domain alias"))
        self.domain_aliases_inv = {k: sorted(v) for k, v in
                                   inverse_dict(self.domain_aliases).items()} ## {pssid: [<aliases>]}, for help
        self.genome_sets = parse_multiline_multikey_sdict(get_lookup("genome sets"))
        self.indv_genomes = parse_multiline_multikey_sdict(get_lookup("genome alias"))
        self.indv_genomes_raw = parse_multiline_multikey_sdict(get_lookup("genome alias")) ## for --genomes
        self.indv_genomes_inv = {k: sorted(v) for k, v in
                                 inverse_dict(self.indv_genomes).items()} ## {fname: [<aliases>]}, for '-i .')
        self.cluster_sets = parse_multiline_multikey_sdict(get_lookup("cluster sets"))
        self.attr_mod_presets = parse_multiline_multikey_sdict(get_lookup("gff attribute modification presets"))
        self.reference_sets = parse_multiline_multikey_sdict(get_lookup("reference sets"))
        self.assembly_aliases = parse_multiline_multikey_sdict(get_lookup("assembly alias"))
        self.annotation_aliases = parse_multiline_multikey_sdict(get_lookup("annotation alias"))
        self.rps_db_aliases = parse_multiline_multikey_sdict(get_lookup("rps database alias"))
        
        ## general
        section_general = "general"
        self.help = Param(None, "-h", "--help")
        self.version = Param(None, "-v", "--version")
        self.quiet = Param(get_val_default(conf_get(section_general, "quiet", type = bool), default = False),
                           "--quiet/--verbose")
        # self.quiet = Param(get_val_default(conf.getboolean(section_general, "quiet"), False),
        #                    "--quiet/--verbose")
        
        ## user lookups (print to screen)
        self.references = Param(None, "--references")
        self.genomes = Param(None, "--genomes")
        self.clusters = Param(None, "--clusters")
        self.members = Param(None, "--members")
        # self.reference_set = Param(None, "--reference-set")
        # self.genome_set = Param(None, "--genome-set")
        # self.cluster_set = Param(None, "--cluster-set")
        
        ## executables
        section_binary = "binary"
        # get_binary = lambda x: conf.get(section_binary, x)
        get_binary = lambda x, default: get_val_default(conf_get(section_binary, x), default = default)
        help_msg = lambda name: (f"absolute path to {name} executable/binary,"
                                 " OR command name if the executable/binary is accessible"
                                 " via the command line")
        self.rpsblast = Param(get_binary("rpsblast", default = "rpsblast"),
                              "--rpsblast", help = help_msg("path to rpsblast or rpsblast+ executable") )
        self.blastn = Param(get_binary("blastn", default = "blastn"),
                            "--blastn", help = help_msg("path to blastn executable") )
        self.mafft = Param(get_binary("mafft", default = "mafft"),
                           "--mafft", help = help_msg("path to mafft executable") )
        # self.rpsblast = Param(get_val_default(get_binary("rpsblast"), "rpsblast"),
        #                       "--rpsblast", help = help_msg("path to rpsblast or rpsblast+ executable") )
        # self.blastn = Param(get_val_default(get_binary("blastn"), "blastn"),
        #                     "--blastn", help = help_msg("path to blastn executable") )
        # self.mafft = Param(get_val_default(get_binary("mafft"), "mafft"),
        #                    "--mafft", help = help_msg("path to mafft executable") )
        
        ## data
        section_data = "data"
        # get_data = lambda x: conf.get(section_data, x) ## for strings
        get_data = lambda x, **kwargs: conf_get(section_data, x, **kwargs)
        parse_reference = lambda x: None if not x else x.split(',')
        self.reference = Param(parse_reference(get_val_default(get_data("reference"), '')),
                               "-r", "--reference",
                               help = ("comma-separated alias(es) for"
                                       " reference genome AND reference genome annotation"),
                               description = ("reference alias"))
                               # description = ("specifies both reference genome sequences and"
                               #                " annotation in a single parameter"))
        # self.reference = Param(get_val_none(get_data("reference"), self.reference_aliases),
        #                        "-r", "--reference",
        #                        help = "reference genome alias",
        #                        description = "comma-separated reference genome alias(es)")
        self.assembly = Param(get_val_none(get_data("assembly"), self.assembly_aliases),
                              "--assembly",
                              help = "reference assembly alias or path to FASTA file of reference assembly",
                              description = ( "reference assembly alias"
                                              " or path to FASTA file of reference assembly" ))
        self.annotation = Param(get_val_none(get_data("annotation"), self.annotation_aliases),
                                "--ann", "--annotation",
                                help = ( "GFF3 file of reference genome annotation or BED file"
                                         " converted from GFF3 format using bedtools' gff2bed" ))
        self.db = Param(get_val_none(get_data("rps database"), self.rps_db_aliases),
                        "--rps-db",
                        help = "local or remote RPS-BLAST database")
        self.remote_rps = Param(get_val_default(get_data("remote rps", type = bool), False), "--remote-rps",
                            help = "Raise --remote flag when executing RPS-BLAST")
        self.reference_set = Param(get_val_none(get_data("reference set"), self.reference_sets),
                                   "--reference-set",
                                   help = ( "File containing reference alias-filename mapping or"
                                            " alias of file containing genome alias-filename mapping as"
                                            " described by 'reference sets' in the config file."
                                            "\nFormat: <semicolon-separated reference aliases><tab>"
                                            "<path to FASTA file of assembly><tab><path to GFF annotation file>"
                                            "\nOne assembly-annotation pair per line." ))
        self.cluster_set = Param(get_val_none(get_data("cluster set"), self.cluster_sets),
                                    "--cluster-set",
                                    help = ( "File containing cluster alias-members mapping or"
                                             " alias of file containing cluster alias-members mapping as"
                                             " described by 'cluster sets' in the config file."
                                             "\nFormat: <semicolon-separated cluster aliases><tab>"
                                             "<comma-separated gene IDs>"
                                             "\nOne cluster per line." ))
        self.genome_set = Param(get_val_none(get_data("genome set"), self.genome_sets),
                                "--genome-set",
                                help = ( "File containing genome alias-filename mapping or"
                                         " alias of file containing genome alias-filename mapping as"
                                         " described by 'genome sets' in the config file."
                                         "\nFormat: <semicolon-separated genome aliases><tab>"
                                         "<absolute path to FASTA file>"
                                         "\nOne FASTA file per line." ))
        self.attr_mod = Param(get_val_none(get_data("gff attribute modification"),
                                           self.attr_mod_presets, none = {}),
                              "--attr-mod",
                              help = ( "mapping for non-standard attribute fields to standard GFF3 fields."
                                       " Only the standard fields 'ID' and 'Parent' will be used."
                                       "\nFormat: '<feature>:<standard field name>=<nonstandard field name>,"
                                       "<standard field name>=<nonstandard field name>;"
                                       "<feature>:<standard field name>=<nonstandard field name>'"
                                       "\nExample: 'mRNA:Parent=Locus_id,ID=Transcript_id;all:ID=Id'") )
        self.ext_cds = Param([], "--extend-cds",
                             help = ( "FASTA file of CDS sequence;"
                                      " used with --extend-gene to generated extended GFF3 annotation" ))
        self.ext_gene = Param([], "--extend-genome", "--extend-gene",
                              help = ( "FASTA file of genomic gene sequence;"
                                       " used with --extend-cds to generate extended GFF3 annotation" ))
        # self.ext_reference = Param(None, "--extend-reference",
        #                            help = ( "FASTA file of genomic data;"
        #                                     " used with --extend-annotation; molecules must be unique" ))
        # self.ext_annotation = Param(None, "--extend-annotation",
        #                             help = ( "GFF3 file of extended genome(s) annotation or BED file"
        #                                      " converted from GFF3 format using bedtools' gff2bed;"
        #                                      " used with --extend-reference; molecules must be unique" ))
        self.ext_assembly = Param(None, "--extend-assembly",
                                  help = ( "FASTA file of genomic data;"
                                           " used with --extend-annotation; molecules must be unique" ))
        self.ext_annotation = Param(None, "--extend-annotation",
                                    help = ( "GFF3 file of extended genome(s) annotation or BED file"
                                             " converted from GFF3 format using bedtools' gff2bed;"
                                             " used with --extend-assembly; molecules must be unique" ))
        self.sep = Param(get_val_default(get_data("gene-isoform separator"), '.'), "--sep",
                         help = ( "character separating gene ID from isoform number/ID;"
                                  " used with --extend-cds and --extend-gene;"
                                  " for example: with gene ID AT1G01010 and isoform ID AT1G01010.1,"
                                  " the separator is '.'" ))
        self.genetic_code = Param(get_val_default(get_data("genetic code"), 1, coerce = int), "--genetic-code",
                                  description = "NCBI genetic code number or name",
                                  help = ( "NCBI genetic code number or name;"
                                           " used with --domain and --db for translation of DNA to amino acids"
                                           " to search for domain using RPS-BLAST;"
                                           " See: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi") )
        
        ## shared
        section_shared = "shared"
        self.directory = Param(os.getcwd(), "-d", "--dir", "--directory",
                               help = "output directory  [default: <working directory>]",
                               show_any_default = False)
        self.prefix = Param(get_val_default(conf_get(section_shared, "prefix"), "minorg"),
                            "--prefix", help = "output files/directory prefix")
        self.fasta = Param(None, "-f", "--fasta")
        self.thread = Param(1, "--thread")
        ## shared by grna, filter, and mininmumset
        self.out_map = Param(None, "--out-map")
        self.out_fasta = Param(None, "--out-fasta")
        
        ## homolog/homologue
        section_homologue = "homologue"
        get_homologue = lambda x, **kwargs: conf_get(section_homologue, x, **kwargs)
        self.gene = Param(None, "-g", "--gene", "--genes",
                          description = "comma-separated gene ID(s)")
        self.cluster = Param(None, "-c", "--cluster",
                             description = "comma-separated cluster alias(es)",
                             alias_value_description = "cluster members")
        self.indv = Param(["none"], "-i", "--indv",
                          # "-a", "--acc",
                          autocompletion = generate_autocompletion("indv",
                                                                   [(INDV_GENOMES_REF, "<reference genome>")] +
                                                                   list(sorted(self.indv_genomes.items())) +
                                                                   [(INDV_GENOMES_ALL,
                                                                     "<all genomes except reference>")]),
                          description = "comma-separated genome alias(es)",
                          alias_value_description = "the location of their FASTA files")
        self.target = Param(None, "-t", "--target")
        self.feature = Param(["CDS"], "--feature",
                             help = "comma-separated GFF feature(s) to find gRNA in")
        # self.query = Param([conf.get(section_homologue, "query")], "-q", "--query")
        self.query = Param([], "-q", "--query")
        self.domain = Param(get_val_none(get_homologue("domain"),
                                         self.domain_aliases, none = "gene"), "--domain",
                                         # self.domain_aliases, none = None), "--domain",
                            # conf.get(section_homologue, "domain"), "--domain",
                            autocompletion = generate_autocompletion("domain",
                                                                     sorted([( pssmid, ','.join(aliases) )
                                                                             for pssmid, aliases in
                                                                             self.domain_aliases_inv.items()])))
        self.merge_within = Param(get_val_default(get_homologue("merge hits within", type = int), 100),
                                  "--merge-within")
        self.minid = Param(get_val_default(get_homologue("minimum hit id", type = float), 95),
                           "--minid", "--min-id")
        self.minlen = Param(get_val_default(get_homologue("minimum target length", type = int), 1),
                            "--minlen", "--min-len", "--min-target-len")
        self.mincdslen = Param(get_val_default(get_homologue("minimum cds length", type = int), 1),
                               "--mincdslen", "--min-cds-len",
                               help = "minimum bases in target aligned by BLAST with reference CDS sequences")
        self.check_recip = Param(get_val_default(get_homologue("check reciprocal", type = bool), True),
                                 "--check-recip", "--check-reciprocal")
        self.relax_recip = Param(get_val_default(get_homologue("relax reciprocal", type = bool), False),
                                 "--relax-recip", "--relax-reciprocal")
        self.check_id_before_merge = Param(get_val_default(get_homologue("check hit id before merging",
                                                                         type = bool), False),
                                           "--check-id-before-merge")
        
        ## generate_grna
        section_grna = "gRNA"
        get_grna = lambda x, **kwargs: conf_get(section_grna, x, **kwargs)
        self.pam = Param(get_val_default(get_grna("pam"), ".NGG"), "-p", "--pam")
        self.length = Param(get_val_default(get_grna("length", type = int), 20), "-l", "--len", "--length")
        self.span_junction = Param(get_val_default(get_grna("span junction", type = bool), False),
                                   "--span-junction", help = "allow gRNA to span intron-exon boundary")
        
        ## filter
        section_filter = "filter"
        get_filter = lambda x, **kwargs: conf_get(section_filter, x, **kwargs)
        self.check_all = Param(False, "--check-all",
                               false_true = ("execute some checks", "execute all checks"),
                               description = "Filter gRNA by all checks")
        self.gc_check = Param(True, "--gc-check/--skip-gc-check",
                              false_true = ("skip GC content check", "execute GC content check"),
                              description = "Filter gRNA by GC content")
        self.feature_check = Param(False, "--feature-check/--skip-feature-check",
                                   false_true = ("skip within feature check", "execute within feature check"),
                                   description = "Filter gRNA by gene feature")
        self.background_check = Param(False, "--bg-check/--skip-bg-check",
                                      "--background-check/--skip-background-check",
                                      "--ot-check/--skip-ot-check",
                                      false_true = ("skip off-target check", "execute off-target check"),
                                       description = "Filter gRNA by off-target")
        self.reset_checks = Param(False, "--reset-checks/--update-checks",
                                  false_true = ("update only requested checks while retaining the rest",
                                                "reset all checks"),
                                  description = "Reset all checks to NA")
        self.background = Param([], "-b", "--bg", "--background",
                                description = "path to file(s) to be screened for off-targets",
                                help = ("if using '--by-indv', each file should only contain sequences"
                                        "  from a single individual, and all sequences from a single"
                                        "  individual should be contained in a single file"))
        self.alignment = Param(None, "--alignment") ## fname
        self.ot_mismatch = Param(get_val_default(get_filter("minimum off-target mismatch", type = int), 1),
                                 "--ot-mismatch", description = "minimum acceptable off-target mismatch")
        self.ot_gap = Param(get_val_default(get_filter("minimum off-target gap", type = int), 0),
                            "--ot-gap", description = "minimum acceptable off-target gap")
        self.ot_pamless = Param(get_val_default(get_filter("pamless off-target search", type = bool), True),
                                "--ot-pamless", help = ("ignore PAM when searching for off-target"))
        self.ot_indv = Param([get_val_default(get_filter("screen individuals"), REFERENCED_ALL)],
                             "--ot-indv",
                             # autocompletion = generate_autocompletion("indv",
                             #                                          [('.', ("<all genomes passed to '-i',"
                             #                                                  " valid only with full programme"
                             #                                                  " and when '-i' is used>")),
                             #                                           ('-', ("<no genomes passed to '-i',"
                             #                                                  " valid only with full programme"
                             #                                                  " and when '-i' is used>"))] +
                             #                                          list(sorted(self.indv_genomes.items()))),
                             description = "comma-separated genome alias(es)",
                             alias_value_description = "the location of their FASTA files",
                             help_subcmd = {"full" :( f"Use '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                                      " to indicate all and no (respectively)"
                                                      " individuals passed to '-i'."
                                                      f" '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                                      " can be used in combination with genome aliases."
                                                      f" (E.g. '--bg-indv {REFERENCED_NONE},HeLa')")})
        self.exclude = Param(None, "-e", "--exclude") ## fname
        self.gc_min = Param(get_val_default(get_filter("GC minimum", type = float), 0.3), "--gc-min",
                            help = ( "minimum GC content (inclusive),"
                                     " where value should be between 0 (0% GC) and 1 (100% GC)"))
        self.gc_max = Param(get_val_default(get_filter("GC maximum", type = float), 0.7), "--gc-max",
                            help = ( "maximum GC content (inclusive),"
                                     " where value should be between 0 (0% GC) and 1 (100% GC)"))
        self.accept_invalid = Param(get_val_default(get_filter("accept invalid", type = bool), False),
                                    "--accept-invalid/--reject-invalid",
                                    false_true = ("reject", "accept"), show_default = False)
        self.accept_feature_unknown = Param(get_val_default(get_filter("accept feature unknown",
                                                                       type = bool), False),
                                            "--accept-feature-unknown/--reject-feature-unknown",
                                        false_true = ("reject", "accept"), show_default = False)
        self.max_insertion = Param(get_val_default(get_filter("maximum insertion size", type = int), 15),
                                   "--max-insertion",
                                   help = ("maximum allowable insertion size (bp) in feature."
                                           " gRNA overlapping with insertions larger than the specific size"
                                           " will be excluded."))
        self.min_within_n = Param(get_val_default(get_filter(("minimum number of genes which features a"
                                                              " gRNA must fall within"), type = int), 1),
                                  "--min-within-n",
                                  help = ("minimum number of genes (-g) which desired features (-f) a"
                                          " gRNA must fall within to pass the 'within feature' check."))
        self.min_within_fraction = Param(get_val_default(get_filter(("minimum fraction of genes which"
                                                                     " features a gRNA must fall within"),
                                                                    type = float), 0),
                                         "--min-within-fraction",
                                         help = ("minimum fraction (0-1) of genes (-g) which desired"
                                                 " features (-f) a gRNA must fall within to pass"
                                                 " the 'within feature' check."))
        self.skip_bg_check = Param(get_val_default(get_filter("skip background check", type = bool),
                                                   False),
                                   "--skip-bg-check",
                                   false_true = ("execute off-target check", "skip off-target check"),
                                   show_default = False,
                                   help = "skip off-target screen")
        self.screen_reference = Param(get_val_default(get_filter("screen reference", type = bool), False),
                                      "--screen-ref", "--screen-reference", show_any_default = False,
                                      help = "screen for off-targets in reference genome")
        self.unmask_ref = Param(get_val_default(get_filter("unmask gene(s) in reference", type = bool),
                                                False),
                                "--unmask-ref", show_any_default = False,
                                help = "skip masking of genes in reference genome for off-target check")
        self.by_indv = Param(get_val_default(get_filter("screen by individual", type = bool), False),
                             "--by-indv",
                             help = ("if raised, gRNA may have off-target effects in some individuals."
                                     " For example, a gRNA targeting a gene in individual A will pass the"
                                     " the background check in individual A but may fail the background check"
                                     " in individual B. If not raised, all gRNA will pass background checks in"
                                     " all individuals."))
        self.mask = Param([], "--mask")
        self.mask_gene = Param(get_val_default(get_filter("mask"), '.').split(','), "--mask-gene",
                               description = "comma-separated gene ID(s)",
                               help = ("masked genes are hidden from background check so that"
                                       " gRNA hits to them will not be considered off-target."),
                               help_subcmd = {"full": ( "masked genes are hidden from background check so that"
                                                        " gRNA hits to them will not be considered off-target."
                                                        f" Use '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                                        " to indicate all and no (respectively)"
                                                        " genes passed to '-g'."
                                                        f" '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                                        " can be used in combination with gene IDs."
                                                        f" (E.g. '--mask {REFERENCED_NONE},BRC1')")})
        ## the following parameters are only used with the "full" subcmd
        self.unmask_gene = Param(get_val_default(get_filter("unmask"), '-').split(','), "--unmask-gene",
                                 help = ("masked genes are hidden from background check so that"
                                         " gRNA hits to them will not be considered off-target."
                                         f" Use '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                         " to indicate all and no (respectively) genes passed to '-g'."
                                         f" '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                         " can be used in combination with gene IDs."
                                         f" (E.g. '--mask {REFERENCED_NONE},BRC1')"
                                         " --unmask will be processed first, followed by --mask."))
        self.mask_cluster = Param(None, "--mask-cluster",
                                  description = "comma-separated cluster alias(es)",
                                  help = ("genes in masked clusters are hidden from background check so that"
                                          " gRNA hits to them will not be considered off-target."),
                                  alias_value_description = "cluster members")
        self.unmask_cluster = Param(None, "--unmask-cluster",
                                    description = "comma-separated cluster alias(es)",
                                    help = ("genes in masked clusters are hidden from background check so that"
                                            " gRNA hits to them will not be considered off-target."),
                                    alias_value_description = "cluster members")
        self.mask_homologue = Param(True, "--mask-homologue/--unmask-homologue",
                                    false_true = ("leave user-provided non-reference sequences unmasked",
                                                  ("infer homologues in user-provided non-reference"
                                                   " sequences and mask from off-target search")))
        self.bg_indv = Param([get_val_default(get_filter("screen individuals"), '.')], "--bg-indv",
                             autocompletion = generate_autocompletion("indv",
                                                                      [(REFERENCED_ALL,
                                                                        ("<all genomes passed to '-i',"
                                                                         " valid only with full programme"
                                                                         " and when '-i' is used>")),
                                                                       (REFERENCED_NONE,
                                                                        ("<no genomes passed to '-i',"
                                                                         " valid only with full programme"
                                                                         " and when '-i' is used>"))] +
                                                                      list(sorted(self.indv_genomes.items()))),
                             description = "comma-separated genome alias(es)",
                             alias_value_description = "the location of their FASTA files",
                             help_subcmd = {"full" :( f"Use '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                                      " to indicate all and no (respectively)"
                                                      " individuals passed to '-i'."
                                                      f" '{REFERENCED_ALL}' and '{REFERENCED_NONE}'"
                                                      " can be used in combination with genome aliases."
                                                      f" (E.g. '--bg-indv {REFERENCED_NONE},HeLa')" )})
        
        ## minimumset
        section_minimumset = "minimumset"
        get_minimumset = lambda x, **kwargs: conf_get(section_minimumset, x, **kwargs)
        self.grna = Param(None, "--grna")
        self.map = Param(..., "-m", "--map")
        self.in_place = Param(False, "--in-place",
                              false_true = ("writes .map to user-specified new file or default output file",
                                            "overwrites .map file passed to --map"))
        self.sets = Param(get_val_default(get_minimumset("sets", type = int), 1), "-s", "--set", "--sets")
        self.sc_algorithm = Param(get_val_default(get_minimumset("set cover algorithm"), "LAR"),
                                  "--sc-algo", "--sc-algorithm",
                                  help = "algorithm for generating minimum set(s)",
                                  case_sensitive = False)
        self.auto = Param(get_val_default(get_minimumset("auto", type = bool), True), "--auto/--manual",
                          # help = "[default: auto]",
                          false_true = ("manual", "auto"), show_default = False)
        self.input_ver = Param(None, "--input-ver")
        self.output_ver = Param(2, "--output-ver")
        
        return
    
    def parse_attr_mod(self, s: str):
        if s in self.attr_mod_presets:
            s = self.attr_mod_presets[s]
        return parse_attr_mod_sdict(s)
    
    def parse_domain(self, s: str):
        return self.domain_aliases.get(s, s)
    
    def get_param(self, param):
        return getattr(self, param)

class OptionParams():
    def __init__(self):
        self.file_valid = {"resolve_path": True, "exists": True,
                           "file_okay": True, "readable": True}
        self.dir_new = {"resolve_path": True, "exists": False}
        return

# ## TODO: make domains compatible w/ genomes and clusters input methods
# class LookupSets():
#     def __init__(self, name = None):

## parses file/string that contains info about aliases
## ...this is too confusing. Either scrap it or streamline it.
class Lookup:
    def __init__(self, name = None, alias_description = "alias", value_description = "value",
                 param_set = None, param_sets = None, param_alias = None, param_inpt = None,
                 config_sets: str = None, config_aliases = None, fname_set = None, none = None,
                 multikey = True, multiline = True, inverse = False, **parse_kwargs):
        ## descriptions
        self._name = name ## used to print lookup descr. (e.g. <self._name> set lookup --> cluster set lookup)
        self._alias_description = alias_description
        self._value_description = value_description
        ## Param objects, used to print parameter names
        self._param_set = param_set ## (e.g. --cluster-set <set alias or path to file of set>)
        self._param_sets = param_sets ## (e.g. --cluster-sets)
        self._param_alias = param_alias ## (e.g. --clusters)
        self._param_inpt = param_input ## (e.g. --cluster <cluster alias>)
        ## values defined in config file
        self._config_fname = config_fname ## global config_fname
        self._config_sets = parse_multiline_multikey_sdict(config_sets) ## dictionary of <set alias>:<set>
        self._config_aliases = config_aliases ## ?? do I wanna keep this or dynamically retrieve <alias>:<value>?
        self._config_lookup = (config_string if not isinstance(config_string, str) ## probably remove all these?
                               else self.parse_from_string(config_string))
        self._config_lookup_inv = {k: sorted(v) for k, v in inverse_dict(self._config_lookup)}
        self._fname = fname
        self._aliases = {}
        self._aliases_inv = {}
        self._multikey = multikey
        self._multiline = multiline
        self._inverse = inverse
        self._parse_kwargs = parse_kwargs
        self._none = none
    def parse_from_string(self, s):
        if self._inverse: parse_f = parse_multiline_inv_sdict
        elif not self._multiline: parse_f = parse_sep_sdict
        elif self._multikey: parse_f = parse_multiline_multikey_sdict
        else: parse_f = parse_multiline_sdict
        self._lookup = parse_f(s, **self._parse_kwargs)
        self._lookup_inv = {k: sorted(v) for k, v in inverse_dict(self._lookup)}
    def parse_from_file(self, fname = None):
        if fname is None: fname = self.fname
        else: self.fname = fname
        with open(fname, 'r') as f:
            self.parse_from_strings(f.read())
    def fname(self):
        return self._fname if self._fname is not None else self._config_fname
    def lookup(self, k, none = None):
        return self._lookup.get(k, self._config_lookup.get(self._none if none is None else none))
    def lookup_inv(self, k, none = None):
        return self._lookup_inv.get(k, self._config_lookup_inv.get(self._none if none is None else none))
    def print(self, header = None, inv = False, prepend = None):
        if not self.fname():
            print( (f"\nNo {self._name} set lookup file has been specified."
                    f"\nYOu may update the list of files and set a default lookup file by"
                    f" modifying the relevant fields in the config file located at {self._config_fname}"
                    f" or provide a lookup file for a single execution using {self._param.names[0]}.") )
        elif not self._config_lookup and not self._lookup and not (inv and self._lookup_inv):
            print(f"\nNo lookup information found in {self.fname()}")
        else:
            print(f"\nThe following information is retrieved from {self.fname()}:\n")
            if header is not None: print(header)
            if prepend is not None: print(prepend)
            if self.fname() == self._config_fname:
                print(self._config_string)
            else:
                with open(self.fname(), 'r') as f:
                    print(f.read())
            if postpend is not None: print(postpend)
        return        

# class Config:
#     def __init__(self, params, keep_on_crash = False, tmp_on = True, keep_tmp = False):
#         self.verbose = False
#         self.out_dir = None ## final output directory
#         self.new_dirs = []
#         self.keep_on_crash = keep_on_crash
#         self.keep_tmp = keep_tmp
#         # self.move_on_crash = False
#         self.tmp_on = tmp_on
#         self.tmp_files = set()
#         self.prefix = None
#         if self.tmp_on:
#             ## tmp directory; contents to be moved to self.out_dir upon fulfilment of command using self.resolve
#             self.tmpdir = tempfile.mkdtemp()
#             ## self.directory is used outside Config objects.
#             ## Config object will handle cleanup if self.directory points to self.tmpdir
#             self.directory = self.tmpdir
#             print(f"tmpdir: {self.tmpdir}")
#         else:
#             print("Problem!! No temporary directory!!")
#             exit
        
#         ## handle log file
#         self.logfile = None
        
#         ## parameter things
#         ## programme params
#         self.raw_args = None
#         self.subcmd = None
#         self.params = params
#         ## shared params
#         self.bed_red = None ## all bed entries to be collapsed into single file regardless of source
#         self.annotation_red = None ## all bed entries to be collapsed into single file regardless of source
#         self.attr_mod = None
#         ## homologue params
#         self.raw_domain = None
#         self.query_map = []
#         # self.cluster_set = None
#         # self.genome_set = None
#         self.reference_aliases = {}
#         self.cluster_aliases = None
#         self.genome_aliases = None
#         # self.reference_ext = ({} if self.params.reference.default is None
#         #                       else {"Reference": self.params.reference.default})
#         # self.bed_ext = {} if self.params.bed is None else {"Reference": self.params.bed}
#         # self.annotation_ext = {} if self.params.annotation is None else {"Reference": self.params.annotation}
#         self.reference = {}
#         self.reference_ext = {}
#         self.reference_indexed = {}
#         self.bed_ext = {}
#         self.annotation_ext = {}
#         self.genetic_code_ext = {}
#         self.attr_mod_ext = {}
#         # ## transition to using Lookup object instead of a bunch of dictionaries
#         # self.mapping_domain = Lookup(name = "domain mapping", config_string = params.domain_mapping)
#         # self.mapping_cluster = Lookup(name = "cluster mapping", config_string = params.cluster_mapping)
#         # self.mapping_genome = Lookup(name = "genome mapping", config_string = params.cluster_mapping)
#         # self.lookup_domain = Lookup(name = "domain", config_string = params.domain_mapping)
#         # self.lookup_cluster = Lookup(name = "cluster", config_string = params.cluster_mapping)
#         # self.lookup_genome = Lookup(name = "genome", config_string = params.cluster_mapping)
#         # # self.cluster_set = params.cluster_mapping.get(params.cluster_set.default,
#         # #                                                  params.cluster_set.default)
    
#     ## make directory (or directories)
#     def mkdir(self, *directories):
#         output = []
#         for directory in directories:
#             if not os.path.isabs(directory):
#                 directory = os.path.join(self.directory, directory)
#             output.append(directory)
#             if not os.path.exists(directory) and not directory in self.new_dirs:
#                 os.mkdir(directory)
#                 self.new_dirs.append(directory)
#         return output[0] if len(directories) == 1 else output
#     ## generate fname (usage: self.mkfname('ref', 'tmp.txt'))
#     def mkfname(self, *path, tmp = False):
#         if os.path.isabs(os.path.join(*path)): path = os.path.join(*path)
#         else: path = os.path.join(self.directory, *path)
#         if tmp: self.tmp_files.add(path)
#         return path
#     ## mkfname, except it also creates path to directory + empty file
#     def reserve_fname(self, *path, tmp = False):
#         fname = self.mkfname(*path)
#         fdir = os.path.dirname(fname)
#         os.makedirs(fdir, exist_ok = True)
#         open(fname, 'a').close()
#         if tmp: self.tmp_files.add(fname)
#         return fname
#     def rm_tmpfiles(self, *fnames):
#         if not self.keep_tmp:
#             to_remove = set(fnames if fnames else self.tmp_files)
#             files_to_remove = [fname for fname in to_remove if
#                                (os.path.exists(fname) and os.path.isfile(fname))]
#             dir_to_remove = [fname for fname in to_remove if
#                              (os.path.exists(fname) and os.path.isdir(fname))]
#             for fname in files_to_remove:
#                 os.remove(fname)
#             for fname in dir_to_remove:
#                 if not os.listdir(fname):
#                     os.rmdir(fname)
#             ## update tmp_files
#             self.tmp_files -= to_remove
#         return
#     ## remove files upon mid-execution termination or crash etc.
#     def cleanup(self):
#         if self.tmp_on and self.keep_on_crash and os.path.exists(self.tmpdir):
#             self.resolve()
#         if self.tmp_on and not self.keep_on_crash and os.path.exists(self.tmpdir):
#             print(f"cleaning up {self.tmpdir}")
#             shutil.rmtree(self.tmpdir)
#     ## move output files to final output directory
#     def resolve(self):
#         # ## rename log file
#         # shutil.move(self.log_file.path, self.reserve_fname(f"{self.prefix}_minorg.log"))
#         # self.log_file.write(f"{self.prefix}_minorg.log")
#         ## move files to final destination (lol)
#         if self.tmp_on:
#             if self.out_dir:
#                 ## create final output directory if doesn't exist
#                 if not os.path.exists(self.out_dir):
#                     os.makedirs(self.out_dir, exist_ok = True)
#                 ## remove tmp files
#                 self.rm_tmpfiles()
#                 ## copy items
#                 mv_dir_overwrite(self.tmpdir, self.out_dir)
#                 print(f"Output files have been generated in {self.out_dir}")
#             ## remove tmpdir
#             shutil.rmtree(self.tmpdir)
#     def set_reference(self, fasta, annotation, genetic_code = 1, attr_mod = {}):
#         self.reference["Reference"] = AnnotatedFasta(fasta, annotation, attr_mod = attr_mod,
#                                                      genetic_code = genetic_code, name = "Reference")
#         if fasta is not None:
#             self.clear_reference()
#             self.add_reference("Reference", fasta)
#         if annotation is not None:
#             self.clear_annotation()
#             self.add_annotation("Reference", annotation, genetic_code = genetic_code)
#         # self.reference_ext = self.reference_ext if fasta is None else {"Reference": fasta}
#         # self.annotation_ext = self.annotation_ext if annotation is None else {"Reference": annotation}
#         return
#     def clear_reference(self):
#         self.reference_ext = {}
#         self.reference_indexed = {}
#     def clear_annotation(self):
#         self.annotation_ext = {}
#     def add_reference(self, ref_id, fasta):
#         self.reference_ext[ref_id] = fasta
#         self.reference_indexed[ref_id] = IndexedFasta(fasta)
#         return
#     def add_annotation(self, ann_id, annotation, genetic_code = 1, attribute_modification = {}):
#         self.annotation_ext[ann_id] = annotation
#         self.genetic_code_ext[ann_id] = genetic_code
#         self.attr_mod_ext[ann_id] = attribute_modification
#         return
#     def extend_reference(self, ref_id, fasta, annotation, genetic_code = 1, attribute_modification = {}):
#         if ref_id in self.reference_ext:
#             raise Exception(f"ID '{ref_id}' already in use as source id for config.reference_ext.")
#         self.reference[ref_id] = AnnotatedFasta(fasta, annotation, genetic_code = genetic_code,
#                                                 attr_mod = attribute_modification, name = ref_id)
#         self.add_reference(ref_id, fasta)
#         self.add_annotation(ref_id, annotation, genetic_code = genetic_code,
#                             attribute_modification = attribute_modification)
#         return
#     def reduce_annotation(self, ids, fout_fmt = "GFF"):
#         mk_tmpf_name = lambda x: self.reserve_fname("reduced_ann", f"reduced.{x}.gff", tmp = True)
#         self.annotation_red = reduce_ann(gff_beds = self.annotation_ext, ids = tuple(ids), fout_fmt = "GFF",
#                                          mk_tmpf_name = mk_tmpf_name, attr_mods = self.attr_mod_ext)
#         for alias, annotated_ref in self.reference.items():
#             annotated_ref.reduce_annotation(tuple(ids), fout = mk_tmpf_name(alias))
#         return

# ## make log file (TODO)
# class LogFile:
#     def __init__(self, config, args: Namespace, params: Params,
#                  ignore = ["version", "quiet", "members", "clusters", "genomes"]):
#         self.config = config
#         self.params = params
#         self.args = vars(args)
#         self.ignore = set(ignore) ## variables not to be written to log file (e.g. eager ones like --version)
    
#     def format_entry(self, argname):
#         argval = self.args[argname]
#         return params.get_param(argname).format_log(argval)
    
#     def write(self, fout = "minorg.log"):
#         if not os.path.isabs(fout):
#             fout = config.mkfname(fout)
#         output = [self.format_entry(argname) for argname in self.args if argname not in self.ignore]
#         with open(fout, "w+") as f:
#             f.write(' '.join(config.raw_args) + '\n\n' + '\n'.join(output) + '\n')
#             ## write query file mapping (especially useful if -q is used)
#             if config.query_map:
#                 f.write("\nquery_name\tquery_file\n")
#                 for name, fname in config.query_map:
#                     f.write(f"{name}\t{fname}\n")
#         return

## Enum choices
class SetCoverAlgo(str, Enum):
    lar = "LAR"
    greedy = "greedy"

# IndvGenomesAll = Enum(
#     "IndvGenomesAll",
#     {**{k: k for k in params.indv_genomes.keys()}, **{"ref":"ref", "none":"none", 'all':'.'}},
#     type = str)

# IndvGenomesAllClear = Enum(
#     "IndvGenomesAll",
#     {**{k: k for k in params.indv_genomes.keys()}, **{"ref":"ref", "none":"none", 'all':'.', 'clear':'-'}},
#     type = str)

# IndvGenomesAll = Enum(
#     "IndvGenomesAll",
#     {"ref":INDV_GENOMES_REF, "none":INDV_GENOMES_NONE, 'all':INDV_GENOMES_ALL},
#     type = str)

# IndvGenomesAllClear = Enum(
#     "IndvGenomesAll",
#     {"ref":"ref", "none":"none", 'all':'.', 'clear':'-'},
#     type = str)
