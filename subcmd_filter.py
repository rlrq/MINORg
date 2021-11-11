import itertools

from functions import fasta_to_dict, dict_to_fasta, assign_alias

from filter_grna import (
    # filter_grna
    filter_gc,
    filter_background_gen,
    filter_excluded_seqs,
    # filter_unique_flank,
    filter_in_feature_gen
)

## These functions will not do any type checking.
## All args are assumed to have been properly parsed.

def execute_filter(args, config, directory, prefix, gRNA_hits, grna_map, fasta_grna,
                   fasta_alignment, fasta_target, gff_bed, fasta_exclude = None,
                   domain_gff_bed = None, checks = None, to_mask = None,
                   # checks = ["feature", "background", "GC", "exclude", "flank"],
                   **kwargs):
    mk_fname = lambda *args, **kwargs: config.mkfname(*args, **kwargs, tmp = True)
    ## parse checks
    if checks is None:
        val_or_empty = lambda b, v: [v] if b else []
        checks = list(itertools.chain(*[val_or_empty(b, v)
                                        for b, v in [(args.gc_check, "GC"),
                                                     (args.feature_check, "feature"),
                                                     (args.background_check, "background")]]))
    # ## parse target_names
    # target_names = None if args.target is None else list(fasta_to_dict(args.target).keys())
    ## output file names
    fout_map_pass = config.reserve_fname(directory, f"{prefix}_gRNA_pass.map")
    fout_fasta_pass = config.reserve_fname(directory, f"{prefix}_gRNA_pass.fasta")
    fout_map_all = config.reserve_fname(directory, f"{prefix}_gRNA_all.map")
    ## filter
    checks = set(check.upper() for check in checks)
    if "GC" in checks:
        filter_gc(gRNA_hits, gc_min = args.gc_min, gc_max = args.gc_max)
    if ( "FEATURE" in checks and args.feature
         and gff_bed is not None and fasta_alignment is not None):
        filter_in_feature_gen(gRNA_hits, fasta_alignment, gff_bed,
                              features = args.feature, max_insertion = args.max_insertion,
                              ref = (set(args.indv) == {"ref"}), domain_gff_bed = domain_gff_bed,
                              min_within_n = args.min_within_n, min_within_fraction = args.min_within_fraction)
    if "BACKGROUND" in checks:
        ## parse fasta_background and fasta_reference into {<alias>: <path>} dicts
        fasta_background = assign_alias(args.background, mk_name = lambda i: "bg_{str(i).zfill(3)}")
        fasta_reference = assign_alias(config.reference_ext, mk_name = lambda i: "ref_{str(i).zfill(3)}")
        fasta_query = dict(config.query_map)
        fout_mask = config.mkfname("masked_report.tsv")
        filter_background_gen(gRNA_hits, fasta_grna, fasta_target, {**fasta_background, **fasta_query},
                              gff_bed = gff_bed, fasta_mask = to_mask, fasta_reference = fasta_reference,
                              max_mismatch = max(1, args.ot_mismatch)-1, max_gap = max(1, args.ot_gap)-1,
                              pam = args.pam, pamless_check = args.ot_pamless, blastn = args.blastn,
                              screen_reference = args.screen_ref,
                              mask_reference = (args.screen_ref and not args.unmask_ref), 
                              report_bg = True,
                              fout_mask = fout_mask, mk_fname = mk_fname)
    if "EXCLUDE" in checks and fasta_exclude is not None:
        filter_excluded_seqs(gRNA_hits, args.exclude)
    # if "FLANK" in checks:
    #     filter_unique_flank(gRNA_hits, flank_length, background_fname, out_dir)
    gRNA_screened = gRNA_hits.filter_seqs_all_checks_passed(ignore_invalid = True).filter_hits_all_checks_passed(ignore_invalid = True)
    gRNA_screened.write_mapping(fout_map_pass, version = 2, write_all = True, write_checks = False)
    gRNA_screened.write_fasta(fout_fasta_pass, write_all = True)
    gRNA_hits.write_mapping(fout_map_all, version = 2, write_all = True, write_checks = True)
    return gRNA_screened, fout_map_all
