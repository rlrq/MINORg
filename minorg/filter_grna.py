import os
import itertools
import regex as re

from Bio import SeqIO
from pathlib import Path

from minorg.functions import (
    assign_alias,
    within_any, gc_content,
    tsv_entries_as_dict, prepend,
    ranges_union, ranges_to_pos, ranges_subtract, ranges_intersect, within, convert_range,
    adjusted_feature_ranges, adjusted_pos,
    imap_progress
)

from minorg.fasta import fasta_to_dict, dict_to_fasta, find_identical_in_fasta
from minorg.blast import blast
from minorg.searchio import searchio_parse

from minorg.index import IndexedFasta
from minorg.annotation import GFF
from minorg.pam import PAM

def filter_background_imap(minorg_ifasta_fout_thread_keep):
    minor_g, ifasta, fout, thread, keep_output = minorg_ifasta_fout_thread_keep
    output = minor_g._offtarget_hits(ifasta, fout = fout, keep_output = keep_output, thread = thread)
    return output

def filter_background_multiindv(args, threads = 1, descr = None):
    if descr:
        msg = lambda curr, last: f"{descr}: {curr}/{last} done."
        return imap_progress(filter_background_imap, args, threads = threads, msg = msg)
    else:
        return imap_progress(filter_background_imap, args, threads = threads)

## TODO: report location of feature(s) for each reference gene

#########################
##  FILTER_BACKGROUND  ##
##   HELPER FUNCTIONS  ##
#########################

class Masked:
    ## takes in a HSP object (produced by Bio.SearchIO)
    def __init__(self, hsp):
        self.hsp = hsp
    @property
    def masked(self): return self.hsp.query_id
    @property
    def molecule(self): return self.hsp.hit_id
    ## Biopython converts ranges to 0-index, end-exclusive! and also sorts them so start is always smaller!
    @property
    def start(self): return self.hsp.hit_start
    @property
    def end(self): return self.hsp.hit_end
    def __eq__(self, other):
        return (isinstance(other, self.__class__) and other.masked == self.masked and
                other.molecule == self.molecule and other.start == self.start and other.end == self.end)
    def __hash__(self):
        return hash((self.masked, self.molecule, self.start, self.end))
    def within(self, r, index = 0, end_incl = False):
        return within(r, convert_range((self.start, self.end), index_in = 0, index_out = index,
                                       start_incl_in = True, start_incl_out = True,
                                       end_incl_in = False, end_incl_out = end_incl))

def mask_identical(to_mask_fname, fasta_fname, fout_fname, **kwargs):
    seqs_to_mask = fasta_to_dict(to_mask_fname)
    ## separate sequences with ambiguous bases (not compatible with BLAST) and those without
    ##  - ambiguous bases typically present in scaffold-level assemblies as runs of 'N's
    standard_bases = {'A','T','G','C','U', 'a', 't', 'g', 'c', 'u'}
    unambig_to_mask = {seqid: seq for seqid, seq in seqs_to_mask.items()
                       if set(seq).issubset(standard_bases)}
    ambig_to_mask = {seqid: seq for seqid, seq in seqs_to_mask.items()
                     if seqid not in unambig_to_mask}
    ## if some (but not all) sequences don't have unambiguous bases
    if (unambig_to_mask and ambig_to_mask):
        import tempfile
        ## replace to_mask_fname with temporary file
        tmp_to_mask_fname = tempfile.mkstemp(suffix = ".fasta")[1]
        ## write unambig_to_mask to file to be used for BLAST
        dict_to_fasta(unambig_to_mask, tmp_to_mask_fname)
    else:
        tmp_to_mask_fname = None
    ## start masking
    masked = []
    ## if at least one sequence doesn't have ambiguous bases
    if unambig_to_mask:
        masked.extend(blast_mask((tmp_to_mask_fname if tmp_to_mask_fname is not None else to_mask_fname),
                                 fasta_fname, fout_fname, **kwargs))
    ## if at least one sequence has ambiguous bases
    if ambig_to_mask:
        masked.extend([Masked(hsp)
                       for query_result in find_identical_in_fasta(ambig_to_mask, fasta_fname)
                       for hit in query_result for hsp in hit.hsps])
    ## remove temporary file if created
    if tmp_to_mask_fname is not None:
        os.remove(tmp_to_mask_fname)
    return masked
    
def blast_mask(to_mask_fname, fasta_fname, fout_fname, blastn = "blastn",
               header = ['qseqid', 'sseqid', 'mismatch', 'gaps',
                         'sstart', 'send', 'qstart', 'qend', 'qlen']):
    ## masking criteria: perfect match from end to end
    col_i = {col_name: i for i, col_name in enumerate(header)}
    is_to_mask = lambda query_result, hsp: (hsp.query_start == 0
                                            and hsp.query_end == query_result.seq_len
                                            and hsp.mismatch_num == 0
                                            and hsp.gap_num == 0)
    ## run blast
    from Bio.Blast.Applications import NcbiblastnCommandline
    blast(NcbiblastnCommandline, cmd = blastn, header = header, fout = fout_fname, task = "megablast",
          query = to_mask_fname, subject = fasta_fname)
    ## get perfect matches
    masked = [Masked(hsp)
              for query_result in searchio_parse(fout_fname, "blast-tab", fields = ' '.join(header))
              for hit in query_result for hsp in hit.hsps
              if is_to_mask(query_result, hsp)]
    return masked

## mask_fnames, background_fname and reference_fasta can be dict of {<alias>: <path>} or [<path>].
##  - if list provided, files will be indexed by well index of position in list.
def mask_and_generate_outside(mask_fnames, background_fnames = None, mask_reference = True,
                              ref_genes_fasta = '', reference_fnames = None, out_dir = None,
                              header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                                        'gaps', 'sstart', 'send', 'qlen', 'slen'],
                              mk_fname = None, blastn = "blastn", fout_mask = None, new_file = True,
                              **kwargs): ## kwargs: sink for irrelevant keyword args
    ## tmp file; reusable, cuz it's not gonna be read
    if mk_fname is not None:
        to_mask_hits = mk_fname("tmp.tsv")
    elif out_dir is not None:
        to_mask_hits = f"{out_dir}/tmp.tsv"
    else:
        import tempfile
        to_mask_hits = tempfile.mkstemp(suffix = ".tsv")[1]
    ## all file paths are assumed to be valid
    ## mask things
    masked = {}
    mask_fnames = assign_alias(mask_fnames, mk_name = lambda i: f"mask_{str(i).zfill(3)}")
    background_fnames = assign_alias(background_fnames, mk_name = lambda i: f"bg_{str(i).zfill(3)}")
    reference_fnames = assign_alias(reference_fnames, mk_name = lambda i: f"ref_{str(i).zfill(3)}")
    ## mask in reference
    if mask_reference and reference_fnames: ## mask in reference genome
        print("Masking gene(s) in reference")
        ## ALT: finish function to define masked regions in reference using GFF(_bed)
        masked = {fname: tuple(itertools.chain(*[mask_identical(mask_fname, fname, to_mask_hits,
                                                                blastn = blastn)
                                                 for mask_fname in mask_fnames.values()]))
                  for fname in reference_fnames.values()}
        ## remove reference FASTA from background_fnames
        ## - ensures that each file is only masked once if reference was also queried for gRNA
        background_fnames = {alias: fname for alias, fname in background_fnames.items()
                             if fname not in reference_fnames.values()}
    ## mask in query or user-provided background
    if background_fnames:
        print("Masking targets (in query or user-provided background)")
        masked = {**masked,
                  **{bg_fname: tuple(itertools.chain(*[mask_identical(mask_fname, bg_fname, to_mask_hits,
                                                                      blastn = blastn)
                                                       for mask_fname in mask_fnames.values()]))
                     for bg_fname in background_fnames.values()}}
    ## write masked ranges if requested
    if fout_mask is not None:
        with open(fout_mask, "w+" if new_file else "a+") as f:
            ## write alias-filename mapping
            inv_fnames = {v: k for k, v in {**background_fnames, **reference_fnames}.items()}
            f.write("#alias\tfname\n")
            for fname, alias in inv_fnames.items():
                f.write(f"#{alias}\t{fname}\n")
            ## header
            f.write('\t'.join(["source", "molecule", "start", "end", "masked"]) + '\n')
            ## write entries
            for fname, masked_in_f in masked.items():
                alias = inv_fnames[fname]
                for entry in masked_in_f:
                    f.write('\t'.join([str(alias),
                                       entry.molecule, str(entry.start+1), str(entry.end), entry.masked]) + '\n')
    ## functions to determine if gRNA sequences is outside masked regions
    ## (to_screen is a dictionary of a blast output entry, where keys are field names)
    def outside_targets(hsp, bg_fname):
        if bg_fname not in masked:
            return True
        for masked_target in masked[bg_fname]:
            same_molecule = hsp.hit_id == masked_target.molecule
            if same_molecule:
                within_target = masked_target.within((hsp.hit_start, hsp.hit_end))
                return not within_target
        return True
    return outside_targets


#########################
##  FILTER_IN_FEATURE  ##
##   HELPER FUNCTIONS  ##
#########################

## min_within_n: minimum number of reference genes for which query_range is within feature_range
## min_within_fraction: minimum proprtion (0-1) of reference genes for which query_range is within feature_range
def within_feature(feature_ranges, seq, query_range, min_within_n = 1, min_within_fraction = 0):
    ## feature_ranges, start, and end should be 0-indexed
    query_range = [adjusted_pos(seq, x + 1) - 1 for x in query_range]
    within_check = [within_any(query_range, ranges) for ranges in feature_ranges.values()]
    within_count = within_check.count(True)
    return ((within_count >= min_within_n) and (within_count/len(within_check) >= min_within_fraction))

def make_target_feature_ranges_function(feature_only_ranges, feature_gaps_ranges,
                                        max_insertion = 15):
    print("Max acceptable insertion length:", max_insertion)
    def get_target_feature_ranges(seq, q_seqid):
        output = {}
        seq_pos = {i for i, c in enumerate(seq) if c != '-'}
        for seqid, gaps_ranges in feature_gaps_ranges.items():
            ## skip if q_seqid is a reference sequence and seqid is NOT <reference + same gene as q_seqid>
            if (q_seqid.split('|')[0] == "Reference" and
                not ((re.search("\|complete\|", seqid)
                      and re.search("\|" + re.escape(re.search("(?<=\|complete\|).+(?=\|\d)", seqid).group(0)) +
                                    "\|\d+-\d+(,|$)", q_seqid))
                     or (re.search("\|gene\|", seqid)
                         and re.search("\|" + '|'.join(seqid.split('|')[5:-1]) + "\|\d+-\d+(,|$)",
                                       q_seqid)))):
                continue
            ## keep gaps if the bases of target in the gap is fewer than max_insertion
            ## (i.e. if the target & ref gene were the only 2 in the alignment,
            ##  the gap would be <= max_insertion)
            acceptable_gaps = []
            for gap_range in gaps_ranges:
                if len(seq_pos.intersection(ranges_to_pos(gaps_ranges))) <= max_insertion:
                    acceptable_gaps.append(gap_range)
            ## generate final acceptable range for target-gene pair
            output[seqid] = ranges_union([feature_only_ranges[seqid] + acceptable_gaps])
        return output
    return get_target_feature_ranges

