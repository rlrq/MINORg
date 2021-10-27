from filter_grna import (
    filter_gc,
    # filter_in_background_gen,
    # filter_excluded_seqs,
    # filter_unique_flank,
    filter_in_feature_gen
)

## These functions will not do any type checking.
## All args are assumed to have been properly parsed.

def execute_filter(args, config, directory, prefix, gRNA_hits, grna_map,
                   fasta_alignment, fasta_target, fasta_ref, gff_bed,
                   checks = ["feature", "GC"],
                   # checks = ["feature", "background", "GC", "exclude", "flank"],
                   **kwargs):
    checks = set(check.upper() for check in checks)
    if "GC" in checks:
        filter_gc(gRNA_hits, gc_min = args.gc_min, gc_max = args.gc_max)
    if "FEATURE" in checks:
        filter_in_feature_gen(gRNA_hits, fasta_alignment, gff_bed, features = args.feature, **kwargs)
    # if "BACKGROUND" in checks:
    #     filter_in_background_gen()
    # if "EXCLUDE" in checks and "exclude" in dir(args) and args.exclude is not None:
    #     filter_excluded_seqs(gRNA_hits, args.exclude)
    # if "FLANK" in checks:
    #     filter_unique_flank(gRNA_hits, flank_length, background_fname, out_dir)
    screened_mapping_fname = config.reserve_fname(directory, f"{prefix}_gRNA_final.map")
    all_gRNA_mapping_fname = config.reserve_fname(directory, f"{prefix}_gRNA_all.map")
    screened_gRNA = gRNA_hits.filter_seqs_all_checks_passed(ignore_invalid = True).filter_hits_all_checks_passed(ignore_invalid = True)
    screened_gRNA.write_mapping(screened_mapping_fname, version = 2, write_all = True, write_checks = False)
    gRNA_hits.write_mapping(all_gRNA_mapping_fname, version = 2, write_all = True, write_checks = True)
    return screened_gRNA, screened_mapping_fname
