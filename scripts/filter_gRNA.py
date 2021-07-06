import os
import re
import sys
import glob
import itertools
# from Bio import SeqIO

sys.path.append("/mnt/chaelab/rachelle/src")
from fasta_manip import fasta_to_dict, dict_to_fasta
sys.path.append(os.path.realpath(__file__))
from find_gRNA_functions import *

## TODO: output
def filter_grna_from_files_and_write(fasta, mapping, input_map_ver = None,
                                     targets = None, background_fname = None, out_dir = None, ## bg filter
                                     gc_min = 0.3, gc_max = 0.7, ## gc filter
                                     alignment_fname = None, complete_fasta = None, cds_fasta = None, ## cds
                                     exclude_fname = None, ## exclude filter
                                     flank_length = 100, ## flank filter
                                     **kwargs):
    ## parse gRNA
    all_gRNA = gRNAHits()
    all_gRNA.parse_from_mapping(mapping, targets = targets, version = input_map_ver)
    
    ## filter against background
    if check_bg and out_dir and targets and len(fasta_to_dict(background_fname)) > 0:
        # ## function for masking that checks if a blast alignment is within a masked region
        # print("Identifying locations of target sequences in background")
        # outside_targets = mask_and_generate_outside(target_fname, background_fname, out_dir = out_dir,
        #                                             ref_genes_fasta = complete_fasta,
        #                                             **{k: v for k, v in kwargs.items()
        #                                                if k in ["mask_reference", "reference_fasta"]})
        # ## filter
        # # print("Filtering background sequences")
        print("Filtering out gRNA with off-target hits in background")
        filter_background(fasta, target, background_fname, all_gRNA,
                          fout_pref = "tmp_background", out_dir = out_dir,
                          # outside_targets = outside_targets,
                          ref_genes_fasta = complete_fasta,
                          **{k: v for k, v in kwargs.items()
                             if k in ["max_mismatch", "max_gap",
                                      "reference_fasta", "screen_reference", "mask_reference"]})
        print("Background filter:", len(all_gRNA.filter_seqs("background", ignore_invalid = False)))
    
    ## filter GC content
    filter_gc(all_gRNA, gc_min, gc_max)
    print("GC filter:", len(all_gRNA.filter_seqs("GC", ignore_invalid = False)))
    
    ## filter within CDS
    if alignment_fname:
        filter_in_cds(all_gRNA, alignment_fname, cds_fasta, complete_fasta,
                      **{k: v for k, v in kwargs.items()
                         if k in ["relax", "min_within_n", "min_within_percentage", "alignment_rvs_pattern"]})
    print("Within CDS filter:", len(all_gRNA.filter_hits("CDS", ignore_invalid = False)))
    
    ## filter against user-specified sequences to exclude
    if exclude_fname:
        filter_excluded_seqs(all_gRNA, exclude_fname)
        print("Exclude filter:", len(all_gRNA.filter_seqs("exclude", ignore_invalid = False)))
    
    ## filter for unique flanking regions
    filter_unique_flank(all_gRNA, flank_length, background_fname, out_dir)
    print("Unique flank filter:", len(all_gRNA.filter_hits("flank", ignore_invalid = False)))
    
    ## write gRNA that pass all of the above filters to file, replace all_gRNA_fname
    screened_gRNA = all_gRNA.filter_seqs_all_checks_passed(ignore_invalid = True).filter_hits_all_checks_passed(ignore_invalid = True)
    screened_gRNA_fname = os.path.join(out_dir, f"{fout_pref}_gRNA_passed.fasta")
    screened_gRNA.write_fasta(screened_gRNA_fname, write_all = True)
    print("All (valid) filters:", len(screened_gRNA))
    
    return
