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
from get_minimum_set import *

arath_getseq_cds_pattern = "^(?=.*AT\dG\d{5})(?!.*\|complete(\|revcomp)?).*"
accs_background_fname_default = "/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta"

# out_dir = "/mnt/chaelab/rachelle/scripts/test_dir/findgrna_v2/findgRNA_B4"; accIDs = ('7273', '6909'); target_fname = out_dir + "/findgRNA_B4_NB-ARC_targets.fasta"; fout_pref = "findgRNA_B4_NB-ARC"; tig_pattern = lambda accID: f"^{accID}\.tig\d"; accs_background_fname = "/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta"; background_usr_fname = ''; pam = "GG"; gRNA_len = 20

def find_cluster_gRNA_in_acc(target_fname, accIDs, out_dir, fout_pref = "findgRNA",
                             tig_pattern = lambda accID: f"^{accID}\.tig\d",
                             accs_background_fname = accs_background_fname_default,
                             background_usr_fname = '', check_bg = True, **kwargs):
    ## read background seqs, filter for accession, and write to file
    background_fname = os.path.join(out_dir, f"{fout_pref}_tmp_background.fasta")
    if check_bg:
        accIDs = set(accIDs)
        if accIDs and accs_background_fname == accs_background_fname_default:
            print("Reading background sequences")
            acc_seqs = {seq_id: seq for seq_id, seq in fasta_to_dict(accs_background_fname).items()
                        for accID in accIDs if re.search(tig_pattern(accID), seq_id)}
            background_exists = write_background_seqs(background_fname, acc_seqs,
                                                      [] if ((not background_usr_fname) or
                                                             background_usr_fname == accs_background_fname)
                                                      else [background_usr_fname])
        else:
            background_fname = background_usr_fname
    ## pass background to find_cluster_gRNA_generic
    find_cluster_gRNA_generic(target_fname, out_dir,
                              fout_pref = fout_pref, background_fname = background_fname,
                              check_bg = check_bg, **kwargs)
    return

def find_cluster_gRNA_accBG(*args, **kwargs):
    find_cluster_gRNA_in_acc(*args, **kwargs)

def find_cluster_gRNA_custom(*args, background_usr_fname = '', **kwargs):
    find_cluster_gRNA_generic(*args, background_fname = background_usr_fname, **kwargs)

def find_cluster_gRNA_generic(target_fname, out_dir, num_sets = 1, manual_check = True, fout_pref = "findgRNA",
                              cds_fasta = '', complete_fasta = '', exclude_fname = '', version = 3,
                              sc_algorithm = "LAR", background_fname = '', alignment_fname = '',
                              gc_min = 0.3, gc_max = 0.8, flank_length = 100, check_bg = True,
                              **kwargs):
    """
    Identifies all possible gRNA sequences from targets (writes this to file)
    Filters gRNA sequences for specificity
    Finds (ideally) the fewest gRNA sequences required to cover all targets (writes this to file)
    Requests user input for exclusion of any sequences and updates list of gRNA candidates (updates file)
    Returns: None
    """
    ## check for target sequences
    os.chdir(out_dir)
    targets = fasta_to_dict(target_fname)
    if not targets:
        print("No target sequences found. Aborting.")
        return
    
    ## get all possible gRNAs from targets
    all_gRNA_fname = os.path.join(out_dir, f"{fout_pref}_gRNA_all.fasta")
    all_gRNA = find_cluster_gRNA(target_fname, **{k: v for k, v in kwargs.items() if k in ("pam", "gRNA_len")})
    all_gRNA.assign_seqid(prefix = "gRNA_")
    all_gRNA.write_fasta(all_gRNA_fname, write_all = True)
    print("Total candidates:", len(all_gRNA))
    
    ## filter against background
    if check_bg:
        # ## function for masking that checks if a blast alignment is within a masked region
        # print("Identifying locations of target sequences in background")
        # outside_targets = mask_and_generate_outside(target_fname, background_fname, out_dir = out_dir,
        #                                             ref_genes_fasta = complete_fasta,
        #                                             **{k: v for k, v in kwargs.items()
        #                                                if k in ["mask_reference", "reference_fasta"]})
        # ## filter
        # # print("Filtering background sequences")
        print("Filtering out gRNA with off-target hits in background")
        filter_background(all_gRNA_fname, target_fname, background_fname, all_gRNA,
                          fout_pref = f"{fout_pref}_background", out_dir = out_dir,
                          # outside_targets = outside_targets,
                          ref_genes_fasta = complete_fasta,
                          **{k: v for k, v in kwargs.items()
                             if k in ["max_mismatch", "max_gap", "report_bg", "pam",
                                      "reference_fasta", "screen_reference", "mask_reference",
                                      "nonref_mask_fname"]})
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
    
    ## write
    screened_mapping_fname = os.path.join(out_dir, f"{fout_pref}_gRNA_passed_targets.txt")
    all_gRNA_mapping_fname = os.path.join(out_dir, f"{fout_pref}_gRNA_all_targets.txt")
    screened_gRNA.write_mapping(screened_mapping_fname, version = version,
                                write_all = True, write_checks = False)
    all_gRNA.write_mapping(all_gRNA_mapping_fname, version = version,
                           write_all = True, write_checks = True)
    
    ## run while manual_check == True AND there are valid combinations
    fout_fasta = os.path.join(out_dir, f"{fout_pref}_gRNA_final.fasta")
    fout_mapping = os.path.join(out_dir, f"{fout_pref}_gRNA_final_targets.txt")
    get_minimum_sets_from_files_and_write(num_sets = num_sets,
                                          mapping = all_gRNA_mapping_fname, targets = target_fname,
                                          input_map_ver = version, output_map_ver = version,
                                          fout_fasta = fout_fasta, fout_mapping = fout_mapping,
                                          ignore_invalid = True, accept_unknown_within_cds_status = False,
                                          sc_algorithm = sc_algorithm, manual_check = manual_check)
    
    ## delete temporary files
    for f in [os.path.join(out_dir, f) for f in os.listdir(out_dir) if re.match(".*tmp.*\..*", f)]:
        os.remove(f)
    return

