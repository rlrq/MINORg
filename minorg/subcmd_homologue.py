import os
import typer
import shutil
import tempfile
import itertools

## TODO: figure out a way to have _complete.bed report the complete range (i.e. min-max for each gene, not range of all features within the "complete" range)

from minorg.functions import (
    # make_make_fname,
    blast6, BlastNR, filter_rpsblast_for_domain,
    parse_get_data, write_tsv, splitlines
)

from minorg.annotation import (
    reduce_ann, extract_features_and_subfeatures
)

from minorg.fasta import (
    fasta_to_dict, dict_to_fasta,
    collapse_identical_seqs
)

from minorg.extract_homologue import (
    find_homologue_multiindv,
    find_homologue_indv# ,
    # get_ref_by_genes_resolve
)

from minorg.exceptions import (
    MessageError,
    InputFormatError,
    InvalidPath,
    InvalidFile,
    UnreadableFile
)

from minorg.display import (
    make_print_preindent
)

from minorg.mafftcommandline_add import MafftCommandline

## These functions will not do any type checking.
## All args are assumed to have been properly parsed.

# ## remove '_seed_' prefix after using --seed for mafft
# def remove_pref(fasta, pref):
#     import re
#     pref_len = len(pref)
#     seqs = {(k if not re.search(f"^{pref}", k) else k[pref_len:]): v for k, v in fasta_to_dict(fasta).items()}
#     dict_to_fasta(seqs, fasta)
#     return

# remove_seed = lambda fasta: remove_pref(fasta, "_seed_")

# ## remove duplicate sequences added so mafft's --seed doesn't throw error
# def remove_duplicate(fasta, duplicated = False):
#     if duplicated:
#         import re
#         seqs = {k: v for k, v in fasta_to_dict(fasta).items() if not re.search("^Duplicate\|", k)}
#         dict_to_fasta(seqs, fasta)
#     return

def execute_homologue(args, config, params, prefix, genes,
                      for_masking = False, keep_bed = False, indv = None):
    
    if indv is None: indv = args.indv
    
    ## if no genes provided (i.e. args.target is defined + no need to mask genes)
    if genes is None:
        directory = config.mkdir(config.mkfname(prefix))
        return (directory, prefix, args.target, None, None, None)
    
    ## reduce BED/GFF file if none provided
    if config.annotation_red is None:
        ann_red = reduce_ann(gff_beds = config.annotation_ext, ids = tuple(genes), fout_fmt = "GFF",
                             mk_tmpf_name = lambda x: config.reserve_fname("reduced_ann",
                                                                           f"reduced.{x}.gff", tmp = True))
        for alias, annotated_ref in config.reference.values():
            annotated_ref.reduce_annotation(tuple(genes),
                                            fout = config.reserve_fname("reduced_ann", f"reduced.{alias}.gff",
                                                                        tmp = True))
    else:
        ann_red = config.annotation_red
    
    ## generate common functions
    seqid_template = "Reference|${source}|${domain}|${feature}|${complete}|${gene}|${isoform}"
    # seqid_source_pattern = f"(?<=^Reference\\|)[^|]+(?=\\|)"
    make_fname = lambda *path, **kwargs: config.reserve_fname(prefix, *path, **kwargs)
    
    ## generate fnames
    # fout_dir = os.path.join(config.directory, prefix)
    ref_pref = f"{prefix}_ref_{config.raw_domain}"
    fasta_gene = make_fname("ref", f"{ref_pref}_gene.fasta")
    fasta_cds = make_fname("ref", f"{ref_pref}_CDS.fasta")
    fout_domain_gff_bed = (None if (args.domain == "gene" or not args.domain)
                           else make_fname("ref", f"{ref_pref}_domain.gff"))
    ## get reference sequences
    print("Extracting reference sequences")
    get_reference_fa(params = params, config = config, make_fname = make_fname, genes = genes,
                     out_pref = prefix, fout = fasta_gene, fout_cds = fasta_cds,
                     domain = args.domain, rpsblast = args.rpsblast, db = args.db,
                     seqid_template = seqid_template, fout_gff = fout_domain_gff_bed)
    
    ## get homologues
    fout_pref = f"{prefix}_{config.raw_domain}"
    fasta_target = make_fname(f"{fout_pref}_targets.fasta" if not for_masking else \
                              f"{fout_pref}_toMask.fasta")
    if set(indv) == {"ref"}:
        shutil.copyfile(fasta_gene, fasta_target)
    else:
        print("Finding homologues in non-reference genomes")
        fasta_queries = config.query_map
        ## blast reference sequences to non-ref sequences
        find_homologue_multiindv(fasta_queries = fasta_queries, fout = fasta_target, directory = fout_dir,
                                 fasta_complete = fasta_gene, fasta_cds = fasta_cds, genes = genes,
                                 min_len = args.minlen, min_id = args.minid, min_cds_len = args.mincdslen,
                                 check_reciprocal = args.check_recip, relax = args.relax_recip,
                                 check_id_before_merge = args.check_id_before_merge, blastn = args.blastn,
                                 gff_beds = config.annotation_ext,
                                 fasta_ref = config.reference_ext,
                                 attribute_mod = config.attr_mod_ext,
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
        ## align reference CDS
        with open(fasta_aln, "w+") as f:
            stdout, stderr = MafftCommandline(args.mafft, input = fasta_cds, # quiet = True,
                                              thread = args.thread)()
            f.write(stdout)
        ## add reference gene to alignment
        with open(tmp_f, "w+") as f:
            stdout, stderr = MafftCommandline(args.mafft, add = fasta_gene, # quiet = True,
                                              thread = args.thread, input = fasta_aln)()
            f.write(stdout)
        ## add non-reference seqs to alignment
        if set(indv) == {"ref"}:
            shutil.copyfile(tmp_f, fasta_aln)
        else:
            with open(fasta_aln, "w+") as f:
                stdout, stderr = MafftCommandline(args.mafft, add = fasta_target, # quiet = True,
                                                  thread = args.thread, input = tmp_f,
                                                  adjustdirectionaccurately = True)()
                f.write(stdout)
        config.rm_tmpfiles(tmp_f)
        return (config.mkfname(prefix), fout_pref, fasta_target, fasta_aln, fasta_gene, fasta_cds,
                ann_red, fout_domain_gff_bed)

## domain ranges will be written to fout_gff
def get_reference_fa(params, config, make_fname, genes, out_pref,
                     fout, fout_cds, db = None, rpsblast = "rpsblast", domain = None, fout_gff = None,
                     seqid_template = "Reference|${source}|${domain}|${feature}|${complete}|${gene}|${isoform}"):
    make_refname = lambda suf: make_fname("ref", "f{out_pref}_ref_{suf}", tmp = True)
    ## get sequences
    seqs_gene = {}
    seqs_cds = {}
    for ref in config.reference.values():
        def get_seqs(fout, feature, complete):
            return ref.get_feature_seq(*genes, fout = fout, feature = feature, fout_gff = fout_gff,
                                       adj_dir = True, by_gene = True, pssmid = domain, complete = complete,
                                       db = db, rpsblast = rpsblast, mktmp = make_fname,
                                       seqid_template = seqid_template, apply_template_to_dict = True)
        ref_seqs_gene = get_seqs(None, None, True)
        ref_seqs_cds = get_seqs(None, "CDS", False)
        seqs_gene = {**seqs_gene, **{seqid: seq for gene_dat in
                                     ref_seqs_gene.values() for seqid, seq in gene_dat.items()}}
        seqs_cds = {**seqs_cds, **{seqid: seq for gene_dat in
                                   ref_seqs_cds.values() for seqid, seq in gene_dat.items()}}
    ## write
    dict_to_fasta(seqs_gene, fout)
    dict_to_fasta(seqs_cds, fout_cds)
    return
