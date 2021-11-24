from minimum_set import get_minimum_sets_from_files_and_write

## These functions will not do any type checking.
## All args are assumed to have been properly parsed.

def execute_minimumset(args, config, prefix = None, directory = None,
                       mapping = None, grna = None, fout_mapping = None, fout_fasta = None):
    if prefix is None:
        prefix = config.prefix
    if directory is None:
        directory = config.directory
    if mapping is None:
        mapping = args.mapping
    if grna is None:
        grna = args.grna
    if fout_mapping is None:
        fout_mapping = args.out_mapping
    if fout_fasta is None:
        fout_fasta = args.out_fasta
    get_minimum_sets_from_files_and_write(mapping = mapping, fasta = grna, exclude_fname = args.exclude,
                                          directory = directory, prefix = prefix,
                                          fout_mapping = fout_mapping, fout_fasta = fout_fasta,
                                          num_sets = args.sets, sc_algorithm = args.sc_algorithm,
                                          output_map_ver = 2, manual_check = (not args.auto),
                                          ignore_invalid = args.accept_invalid,
                                          accept_unknown_within_feature_status = args.accept_feature_unknown)
    return
