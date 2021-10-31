import os
import typer
import shutil
import tempfile
import itertools
import subprocess

## TODO: figure out a way to have _complete.bed report the complete range (i.e. min-max for each gene, not range of all features within the "complete" range)

from functions import (
    # make_make_fname,
    reduce_bed, GFF,
    blast6, extract_features_and_subfeatures,
    fasta_to_dict, dict_to_fasta,
    parse_get_data, write_tsv, splitlines
)

from extract_homologue import (
    find_homologue_multiindv,
    find_homologue_indv,
    get_ref_by_genes_resolve
)

from exceptions import (
    MessageError,
    InputFormatError,
    InvalidPath,
    InvalidFile,
    UnreadableFile
)

from display import (
    make_print_preindent
)

from mafftcommandline_add import MafftCommandline

## These functions will not do any type checking.
## All args are assumed to have been properly parsed.

# def map_query_and_log(prefix, args, config, params):
#     if args.indv:
#         fasta_queries = [[indv, params.indv_genomes[indv]] for indv in args.indv]
#     elif args.query:
#         fasta_queries = [[i, query] for i, query in enumerate(args.query)]    
#     with open(config.log_file as "a+") as f:
#         f.write("\nquery_name\tquery_file")
#         for name, fname in fasta_queries:
#             f.write(f"\n{name}\t{fname}")
#     return fasta_queries

## remove '_seed_' prefix after using --seed for mafft
def remove_pref(fasta, pref):
    import re
    pref_len = len(pref)
    seqs = {(k if not re.search(f"^{pref}", k) else k[pref_len:]): v for k, v in fasta_to_dict(fasta).items()}
    dict_to_fasta(seqs, fasta)
    return

remove_seed = lambda fasta: remove_pref(fasta, "_seed_")

## remove duplicate sequences added so mafft's --seed doesn't throw error
def remove_duplicate(fasta, duplicated = False):
    if duplicated:
        import re
        seqs = {k: v for k, v in fasta_to_dict(fasta).items() if not re.search("^Duplicate\|", k)}
        dict_to_fasta(seqs, fasta)
    return

def execute_homologue(args, config, params, prefix, genes, for_masking = False, keep_bed = False):
    
    ## if no genes provided (i.e. args.target is defined + no need to mask genes)
    if genes is None:
        directory = config.mkdir(config.mkfname(prefix))
        return (directory, prefix, args.target, None, None, None)
    
    ## reduce bed file if none provided
    if config.bed_red is None:
        bed_red = config.mkfname(f"{prefix}_reduced.bed", tmp = True)
        reduce_bed(beds = config.bed_ext.values(), ids = tuple(genes), fout = bed_red,
                   mk_tmpf_name = lambda x: config.mkfname(f"tmp_reduced.{x}.bed", tmp = True))
    #     bed_red = config.mkfname(f"{prefix}_reduced.bed", tmp = True)
    #     extract_features_and_subfeatures(args.bed, tuple(genes), bed_red,
    #                                      quiet = True, fin_fmt = "BED", fout_fmt = "BED")
    else:
        bed_red = config.bed_red
    
    ## generate common functions
    make_fname = lambda *path, **kwargs: config.reserve_fname(prefix, *path, **kwargs)
    def get_seq_ref_genes(genes, feature, fasta_out, out_dir, **kwargs):
        return get_ref_by_genes_resolve(genes = genes, feature = feature, adj_dir = True,
                                        out_dir = out_dir, fout = fasta_out, ## output options
                                        ref_fasta_files = config.reference_ext, ## reference
                                        gff_bed = bed_red, attribute_mod = args.attr_mod, ## reference annotation
                                        minlen = args.minlen, quiet = True,
                                        # seqid_template = "Reference|$source|$domain|$feature|$complete|$gene",
                                        seqid_template = "Reference|$source|$domain|$feature|$complete|$gene|$isoform",
                                        no_bed = (not keep_bed), **kwargs)
    
    ## generate fnames
    fout_dir = os.path.join(config.directory, prefix)
    ref_pref = f"{prefix}_ref_{config.raw_domain}"
    fasta_gene = make_fname("ref", f"{ref_pref}_gene.fasta")
    fasta_cds = make_fname("ref", f"{ref_pref}_CDS.fasta")
    fout_domain_gff_bed = None if not args.domain else make_fname("ref", f"{ref_pref}_domain.gff")
    # ## bed_complete and bed_cds are identical if complete is just well CDS but complete.
    # ##  TODO: figure out a way to have _complete report the complete range (i.e. min-max for each gene)
    # bed_complete = make_fname("ref", f"{ref_pref}_complete.bed") if keep_bed else None
    # bed_cds = make_fname("ref", f"{ref_pref}_CDS.bed") if keep_bed else None
    ## get reference sequences
    get_reference_fa(params = params, config = config, get_ref = get_seq_ref_genes,
                     make_fname = make_fname, genes = genes, out_dir = fout_dir, out_pref = prefix,
                     fout = fasta_gene, fout_cds = fasta_cds, domain = args.domain,
                     fout_domain_gff_bed = fout_domain_gff_bed,
                     # bedout = bed_complete, bedout_cds = bed_cds,
                     db = args.db, rpsblast = str(args.rpsblast))
    
    ## get homologues
    fout_pref = f"{prefix}_{config.raw_domain}"
    fasta_target = make_fname(f"{fout_pref}_targets.fasta" if not for_masking else \
                              f"{fout_pref}_toMask.fasta")
    if set(args.indv) == {"ref"}:
        shutil.copyfile(fasta_gene, fasta_target)
    else:
        print("Finding homologues in non-reference genomes")
        fasta_queries = config.query_map
        ## blast reference sequences to non-ref sequences
        find_homologue_multiindv(fasta_queries = fasta_queries, fout = fasta_target, directory = fout_dir,
                                 fasta_gene = fasta_gene, fasta_cds = fasta_cds, genes = genes,
                                 min_len = args.minlen, min_id = args.minid, min_cds_len = args.mincdslen,
                                 check_reciprocal = args.check_recip, relax = args.relax_recip,
                                 check_id_before_merge = args.check_id_premerge, blastn = args.blastn,
                                 gff_beds = tuple(config.bed_ext.values()),
                                 fasta_ref = tuple(config.reference_ext.values()),
                                 merge_within_range = args.merge_within,
                                 keep_tmp = config.keep_tmp)
    
    if for_masking:
        ## skip alignment step, only return potential homologues (don't return reference seqs)
        return (config.mkfname(prefix), fout_pref, fasta_target, None, None, None)
    else:
        ## align
        ## TODO: make sure not using -g/-c (i.e. no fasta_cds/complete) won't crash this part (2021/05/07)
        ## - tentatively okayed (2021/05/12)
        fasta_aln = make_fname(f"{prefix}_{config.raw_domain}_mafft.fasta")
        tmp_f = make_fname(f"tmp_mafft.fasta", tmp = True)
        duplicated = False
        with open(fasta_aln, "w+") as f:
            stdout, stderr = MafftCommandline(args.mafft, input = fasta_cds, # quiet = True,
                                              thread = args.thread)()
            f.write(stdout)
            # ## if fewer than 2 seqs, duplicate so that --seed doesn't throw error
            # if stdout.count('>') == 1:
            #     if stdout[0] != '\n' and stdout[-1] != '\n':
            #         f.write('\n')
            #     f.write('\n'.join([(x if x[0] != '>' else ">Duplicate|tmp") for x in stdout.split('\n')]))
            #     # # f.write(">Duplicate|" + '|'.join(stdout.split('\n')[0].split('|')[1:]) + '\n')
            #     # f.write(">Duplicate|" + '|'.join(stdout.split('\n')[0].split('|')[1:]) + '\n')
            #     # f.write('\n'.join(stdout.split('\n')[1:]))
            #     duplicated = True
        with open(tmp_f, "w+") as f:
            stdout, stderr = MafftCommandline(args.mafft, add = fasta_gene, # quiet = True,
                                              thread = args.thread, input = fasta_aln)()
            f.write(stdout)
        ## remove '_seed_' prefix
        remove_seed(tmp_f)
        ## remove duplicated seq
        remove_duplicate(tmp_f, duplicated = duplicated)
        if set(args.indv) == {"ref"}:
            shutil.copyfile(tmp_f, fasta_aln)
        else:
            with open(fasta_aln, "w+") as f:
                stdout, stderr = MafftCommandline(args.mafft, add = fasta_target, # quiet = True,
                                                  thread = args.thread, input = tmp_f,
                                                  adjustdirectionaccurately = True)()
                f.write(stdout)
                # subprocess.run(["mafft", "--quiet", "--add", fasta_target, tmp_f], stdout = f)
            ## remove '_seed_' prefix
            remove_seed(fasta_aln)
        config.rm_tmpfiles(tmp_f)
        return (config.mkfname(prefix), fout_pref, fasta_target, fasta_aln, fasta_gene, fasta_cds,
                bed_red, fout_domain_gff_bed)

## collapses identical sequences (arbitrarily selects a seqid to represent each set)
## writes collapsed fasta file and a tsv file mapping identical sequences
def collapse_identical_seqs(fasta, fout_fa, fout_seqid):
    dat = fasta_to_dict(fasta)
    identicals = {k: set(seqid for seqid, seq in dat.items() if str(seq) == str(v))
                  for k, v in dat.items()}
    identical_sets = set(map(lambda x: tuple(sorted(x)), identicals.values()))
    dict_to_fasta({seqids[0]: dat[seqids[0]] for seqids in identical_sets}, fout_fa)
    open(fout_seqid, "w+").write('\n'.join(['\t'.join(seqids) for seqids in identical_sets]))
    return

## bed should be a reduced BED file filtered for relevant entries
def get_reference_fa(params, config, get_ref, make_fname, genes, out_dir, out_pref,
                     fout, fout_cds, db, gff_bed = None, fout_domain_gff_bed = None,
                     bedout = None, bedout_cds = None, rpsblast = "rpsblast",
                     domain = None, domain_aliases = {}, lvl = 0):
    
    ## make common functions
    make_refname = lambda suff, **kwargs: make_fname("ref", f"{out_pref}_ref_{suff}", **kwargs)
    printi = make_print_preindent(initial_lvl = lvl)
    
    ## parse domain
    if domain is not None and domain != "gene":
        if domain not in domain_aliases:
            # printi( ( f"'{domain}' is not a supported domain alias."
            #           f"\nAttempting to parse '{domain}' as CDD PSSM-Id." ) )
            try:
                pssmid = str(int(domain))
            except ValueError:
                raise InputFormatError(message = ( "Unexpected domain input."
                                                   " Please provide a CDD PSSM-Id (numeric)"
                                                   " or one of the following supported domain names:"
                                                   f" {', '.join(domain_aliases.keys())}" ),
                                       hint = params.domain.help() )
        else:
            domain = "gene"
            pssmid = domain_aliases[domain]
        
        ## generate local file names
        fasta_aa = make_refname("pep.fasta")
        fasta_aa_collapsed = make_refname("pep.col.fasta", tmp = True)
        seqid_aa_collapsed = make_refname("pep.col.seqid.txt", tmp = True)
        tsv_domains = make_refname(f"{domain}.tsv", tmp = True)
        
        ## get protein sequences
        printi("Extracting reference peptide sequence(s)", overwrite = True)
        get_ref(genes, "CDS", fasta_aa, out_dir, translate = True, by_gene = False)
        collapse_identical_seqs(fasta_aa, fasta_aa_collapsed, seqid_aa_collapsed)
        
        ## generate file for input as DOMAIN_F to getSeq
        printi("Extracting reference domain range(s)", overwrite = True)
        ## get domain positions in protein
        from Bio.Blast.Applications import NcbirpsblastCommandline
        blast6(blastf = NcbirpsblastCommandline, cmd = rpsblast,
               header = "qseqid,sseqid,pident,length,qstart,qend",
               fout = tsv_domains, query = fasta_aa_collapsed, db = db)
        ## copy tsv_domains somewhere else to keep for posterity
        with open(tsv_domains + ".txt", "w+") as f1:
            with open(tsv_domains, 'r') as f2:
                for line in f2:
                    f1.write(line)
        ## filter blast6 output for relevant entries
        filter_rpsblast_for_domain(fname = tsv_domains, pssmid = pssmid)
        ## expand filtered output of rpsblast+ to sequences w/ identical protein
        expand_rpsblast_to_redundant_peptides(blast_fname = tsv_domains, collapsed_fname = seqid_aa_collapsed)
        
        ## get domain nucleotide sequences
        printi("Extracting reference domain nucleotide sequence(s)", overwrite = True)
        def get_ref_domain(fasta_out, complete, feature = "CDS", **kwargs):
            return get_ref(genes = genes, feature = feature, domain_f = tsv_domains, out_dir = out_dir,
                           qname_dname = ("qseqid", "domain"), qstart_qend = ("qstart", "qend"),
                           domain = domain, verbose = False, fasta_out = fasta_out, complete = complete,
                           by_gene = True, fout_domain_gff_bed = fout_domain_gff_bed, **kwargs)
        ## complete domain sequence
        get_ref_domain(fasta_out = fout, complete = True, bed_out = bedout, feature = "CDS")
        ## CDS-only domain sequence
        get_ref_domain(fasta_out = fout_cds, complete = False, bed_out = bedout_cds, feature = "CDS")
        # ## remove tmp files (tmp files to be removed w/ other tmp files on cleanup)
        # os.remove(tsv_domains)
        # os.remove(fasta_aa_collapsed)
        # os.remove(seqid_aa_collapsed)
        
    ## else if no domain provided
    else:
        def get_ref_gene(fasta_out, complete, feature = "CDS", **kwargs):
            return get_ref(genes = genes, feature = feature, out_dir = out_dir, verbose = False,
                           fasta_out = fout, complete = complete, by_gene = True, **kwargs)
        ## get complete gene sequence
        get_ref(genes, "gene", fout, out_dir, complete = True, by_gene = True, bed_out = bedout)
        ## get CDS-only of gene sequence
        get_ref(genes, "CDS", fout_cds, out_dir, complete = False, by_gene = True, bed_out = bedout_cds)
    return


## overwrites original file
def filter_rpsblast_for_domain(fname, pssmid):
    pssmid = str(pssmid)
    get, dat = parse_get_data(fname, delim = '\t')
    output = [list(x) + [pssmid] for x in dat if pssmid in get(x, "sseqid").split('|')]
    write_tsv(fname, [get(get_cols = True) + ["domain"]] + output)
    return

## overwrites original file
def expand_rpsblast_to_redundant_peptides(blast_fname, collapsed_fname):
    get, dat = parse_get_data(blast_fname, delim = '\t')
    repr_map = [x.split('\t') for x in splitlines(collapsed_fname)]
    ## expand entry for identical peptides previously collapsed
    output = [[seqid] + line[1:] for seqids in repr_map for seqid in seqids for line in dat \
              if line[0] == seqids[0]]
    write_tsv(blast_fname, [get(get_cols = True)] + output)
    return

