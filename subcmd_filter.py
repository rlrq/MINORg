from functions import fasta_to_dict, dict_to_fasta

from filter_grna import (
    filter_grna
    # filter_gc,
    # # filter_in_background_gen,
    # # filter_excluded_seqs,
    # # filter_unique_flank,
    # filter_in_feature_gen
)

## These functions will not do any type checking.
## All args are assumed to have been properly parsed.

def execute_filter(args, config, directory, prefix, gRNA_hits, grna_map,
                   fasta_alignment, fasta_target, gff_bed, fasta_exclude = None,
                   domain_gff_bed = None, checks = ["feature", "GC"],
                   # checks = ["feature", "background", "GC", "exclude", "flank"],
                   **kwargs):
    ## parse target_names
    target_names = None if args.target is None else list(fasta_to_dict(args.target).keys())
    ## output file names
    fout_map_pass = config.reserve_fname(directory, f"{prefix}_gRNA_pass.map")
    fout_map_all = config.reserve_fname(directory, f"{prefix}_gRNA_all.map")
    ## filter
    gRNA_screened = filter_grna(gRNA_hits, target_names = target_names,
                                gc_min = args.gc_min, gc_max = args.gc_max, features = args.feature,
                                fasta_alignment = fasta_alignment, fasta_exclude = fasta_exclude,
                                gff_bed = gff_bed, fasta_background = args.background,
                                checks = checks, ref = (set(args.indv) == {"ref"}),
                                domain_gff_bed = domain_gff_bed)
    gRNA_screened.write_mapping(fout_map_pass, version = 2, write_all = True, write_checks = False)
    gRNA_hits.write_mapping(fout_map_all, version = 2, write_all = True, write_checks = True)
    return gRNA_screened, fout_map_pass
