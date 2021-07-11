import os
import sys
import itertools
import subprocess
import regex as re
from Bio import SeqIO
sys.path.append("/mnt/chaelab/rachelle/src")
from fasta_manip import fasta_to_dict, dict_to_fasta
from data_manip import splitlines

###############
##  IDK MAN  ##
###############

def infer_full_pam(pam):
    ## use default 3' PAM + 1 base spacer if '.' and 'N' not provided
    pam = pam.upper()
    if '.' not in pam:
        ## assume 3' PAM if not indicated
        if 'N' not in pam: pam = 'N' + pam
        ## 3' PAM
        if pam[0] == 'N': pam = '.' + pam
        ## if 5' PAM
        else: pam = pam + '.'
    return pam

def expand_ambiguous(pam):
    """
    Map ambiguous bases
    """
    from Bio import Seq
    amb_dna = Seq.IUPAC.IUPACData.ambiguous_dna_values
    pam_mapped = ''
    for c in pam:
        mapped = amb_dna.get(c.upper(), c)
        if len(mapped) == 1: pam_mapped += mapped
        else: pam_mapped += f"[{mapped}]"
    return pam_mapped

def make_pam_pattern(pam, gRNA_len):
    """
    Square brackets not allowed in PAM pattern.
    """
    ## infer pam location + spacer if not explicitly described
    pam = infer_full_pam(pam)
    ## map ambiguous bases
    pam_mapped = expand_ambiguous(pam)
    ## generate pattern for compilation
    grna_pre, grna_post = pam_mapped.split('.')
    pam_pattern = ''
    if grna_pre: pam_pattern += f"(?<={grna_pre})"
    pam_pattern += ".{" + str(gRNA_len) + '}'
    if grna_post: pam_pattern += f"(?={grna_post})"
    ## generate pattern for gRNA extraction
    print("PAM pattern:", pam_pattern)
    return re.compile(pam_pattern, re.IGNORECASE)

def find_cluster_gRNA_gen(target_fname, pam = "NGG", gRNA_len = 20): #Let NGG be GG, or NG be G
    """
    Finds all possible target sequences given pam.
    Returns: {<seq>: set(ids of targeted targets)}
    """
    ## specifying 'N' allows users to specify how many wobble/N bases
    ## in the future, maybe use regex?
    pam_pattern = make_pam_pattern(pam, gRNA_len = gRNA_len)
    # TODO: some function to parse NXX or XXN PAM into regex
    dic_target = {}
    parsed = SeqIO.parse(target_fname, "fasta")
    hit_id = 0
    for sequence in parsed:
        fwd_seq, rvs_seq = sequence.seq, sequence.seq.reverse_complement()
        target = Target(fwd_seq, id = sequence.id, strand = '+')
        for strand, seq in {'+': fwd_seq, '-': rvs_seq}.items():
            pam_pos = [x.start() for x in re.finditer(pam_pattern, str(seq), overlapped = True)]
            for start in pam_pos:
                end = start + gRNA_len
                grna_seq = str(seq[start: start + gRNA_len]).upper()
                if len(grna_seq) == gRNA_len:
                    ## gRNAHit's target is set as fwd_seq (variable 'target')
                    ## gRNAHit's range is set relative to fwd_seq (direction-less)
                    gRNAhit = gRNAHit(target,
                                      start if strand == '+' else len(target) - end,
                                      end if strand == '+' else len(target) - start,
                                      strand, hit_id)
                    dic_target[grna_seq] = dic_target.get(grna_seq, set()).union({gRNAhit})
                    hit_id += 1
    # return dic_target
    output_gRNA_hits = gRNAHits()
    output_gRNA_hits.parse_from_dict(dic_target)
    return output_gRNA_hits

###################
##  WRITE FILES  ##
###################

## combines fasta files of background sequences with default background from target search
def write_background_seqs(fout, background_w_targets = {}, background_fnames = []):
    from re import search
    output = background_w_targets
    dup_titles = {}
    for background_fname in background_fnames:
        to_extend = fasta_to_dict(background_fname)
        # ## pre-pending fname is causing some problems. So for now pls ensure that bg molecules don't have the same name(s) as reference molecules.
        # basename = re.search("[^/]+?(?=\.[^./]+?$)", background_fname).group(0)
        # to_extend = {f'{basename}_{seqid}': seq for seqid, seq in fasta_to_dict(background_fname).items()}
        # for title in to_extend:
        #     if title in output:
        #         dup_titles[title] = dup_titles.get(title, 0) + 1
        #         seq = to_extend[title]
        #         del(to_extend[title])
        #         to_extend[f'{title}_dupl{dup_titles[title] + 1}'] = seq
        output = {**output, **to_extend}
    dict_to_fasta(output, fout)
    return (len(output) > 0)


#################
##  FILTERING  ##
##  FUNCTIONS  ##
#################

## filter for gRNA off-target hits in background
def filter_background_gen(seqs_fname, target_fname, background_fname, gRNA_hits, screen_reference = False,
                          fout_pref = "", out_dir = "", max_mismatch = 0, max_gap = 0, mask_reference = True,
                          reference_fasta = '', ref_genes_fasta = False, outside_targets = None, pam = "GG",
                          report_bg = False, nonref_mask_fname = '', pamless_bg_check = False,
                          **kwargs): ## kwargs is for arguments for the mask_and_generate_outside function
    """
    Executes bl2seq with candidate gRNA sequences as query and background as subject.
    Filters for maximum allowable mismatches and gaps.
    Removes candidate gRNA that pass the above criteria but are also found outside of target sequences.
    Inputs:
        seqs_fname: path to FASTA file of gRNA candidate sequences
        target_fname: path to FASTA file of target sequences
        background_fname: path to FASTA file of all background sequences (can include targets)
        gRNA_hits: gRNAHits object
    // DEFUNCT // Returns: filtered dictionary of {<gRNA seq>: set(ids of targeted target)}
    Returns: None. gRNAHits is modified in-place.
    """
    ## parse file names ( only handles a limited number of combinations :/ )
    if not fout_pref:
        fout_pref = os.path.join(out_dir, "proj")
    elif fout_pref and out_dir and not os.path.dirname(fout_pref):
        fout_pref = os.path.join(out_dir, fout_pref)
    ## create function that masks targets and checks if gRNA hit is in or out of targets if not provided
    if not outside_targets:
        ## if nonref_mask_fname (FASTA file of non-targets that are to be masked) is provided
        if nonref_mask_fname and os.path.exists(nonref_mask_fname):
            ## merge nonref_mask_fname & target_fname into nonref_mask_fname
            to_mask = {**{"user_provided_" + k: v for k, v in fasta_to_dict(nonref_mask_fname).items()},
                       **{"pred_target_" + k: v for k, v in fasta_to_dict(target_fname).items()}}
            to_mask_fname = os.path.join(out_dir, "tmp_to_mask.fasta")
            dict_to_fasta(to_mask, to_mask_fname)
        ## else just use target_fname as source of sequences to mask
        else:
            to_mask_fname = target_fname
        ## mask
        outside_targets = mask_and_generate_outside(to_mask_fname, background_fname,
                                                    out_dir = out_dir,
                                                    mask_reference = mask_reference,
                                                    ref_genes_fasta = ref_genes_fasta,
                                                    reference_fasta = reference_fasta)
    ## index fasta files
    indexed_fa = {fname: SeqIO.index(fname, "fasta") for fname in \
                  ( (background_fname, reference_fasta) if reference_fasta else (background_fname,) )}
    ## make exclude function (function that returns True if gRNA fails off-target check)
    exclude_with_pam = make_exclude_function(pam, outside_targets, indexed_fa, max_mismatch, max_gap)
    ## prep for blast
    blast_output_header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                           'gaps', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'sstrand']
    ## function to prepend header if user wants a bg report, else delete blast hits file
    def resolve_bg_blast(fname):
        from file_manip import prepend
        prepend(fname, '\t'.join(blast_output_header) + '\n')
        if report_bg: prepend(fname, '\t'.join(blast_output_header) + '\n')
        else: os.remove(fname)
    ## function to get seqs to exclude
    def get_excl_seqid(seqs_fname, fasta_fname, fout_fname, fully_aligned_check = True):
        def precheck_exclude(x, mismatch_check = True, gap_check = True, unaligned_check = True):
            mismatch_excess = max_mismatch - int(x["mismatch"])
            gap_excess = max_gap - int(x["gaps"])
            unaligned_excess = (int(x["qlen"]) - int(x["qend"])) + (int(x["qstart"]) - 1)
            fully_aligned = (int(x["length"]) + int(x["gaps"]) == int(x["qlen"])) \
                            if fully_aligned_check else True
            mismatch_pass = (mismatch_excess >= 0) if mismatch_check else True
            gap_pass = (gap_excess >= 0) if gap_check else True
            unaligned_pass = unaligned_excess <= (mismatch_excess + gap_excess) if unaligned_check else True
            return fully_aligned and mismatch_pass and gap_pass and unaligned_pass and \
                outside_targets(x, fasta_fname)
        entries = blastn(seqs_fname, fasta_fname, to_dict = True, fout = fout_fname,
                         header = blast_output_header)
        prelim_exclude = [entry for entry in entries if precheck_exclude(entry)]
        final_exclude = prelim_exclude if pamless_bg_check else exclude_with_pam(prelim_exclude, fasta_fname)
        return set(x["qseqid"] for x in final_exclude)
    ## run blast to identify gRNA hits in background
    if background_fname and os.path.exists(background_fname):
        print("Assessing off-target in user-provided sequences")
        nonref_hits = f"{fout_pref}_hit_nonref.tsv"
        excl_seqid = get_excl_seqid(seqs_fname, background_fname, nonref_hits, fully_aligned_check = False)
        resolve_bg_blast(nonref_hits)
    else:
        excl_seqid = set()
    if screen_reference and reference_fasta: ## find in reference
        print("Assessing off-target in reference genome")
        ref_hits = f"{fout_pref}_hit_ref.tsv"
        excl_seqid |= get_excl_seqid(seqs_fname, reference_fasta, ref_hits)
        resolve_bg_blast(ref_hits)
    ## parse bg check status
    print("Filtering gRNA sequences")
    seqs_fasta = fasta_to_dict(seqs_fname)
    seqs_screened = set(map(lambda s: str(s).upper(), seqs_fasta.values()))
    seqs_failed = set(str(seqs_fasta[seqid]).upper() for seqid in excl_seqid)
    seqs_passed = set(seq for seq in seqs_screened if seq not in seqs_failed)
    ## record bg check status
    gRNA_hits.set_seqs_check("background", True, seqs_passed)
    gRNA_hits.set_seqs_check("background", False, seqs_failed)
    return # gRNA_hits

## filter GC content
def filter_gc(gRNA_hits, gc_min = 0, gc_max = 1):
    from seq_manip import gc_content
    for gRNA_seq in gRNA_hits.gRNAseqs().values():
        gRNA_seq.set_gc_check(gc_min = gc_min, gc_max = gc_max)
    return

## filter within CDS
def filter_in_cds_gen(gRNA_hits, alignment_fname, cds_fasta, complete_fasta, max_cds_insertion = 15,
                      alignment_rvs_pattern = '^_R_', min_within_n = 1, min_within_percentage = 0,
                      relax = True, relax_cds_within = 15):
    print("Filtering for sequences within reference CDS")
    ## get CDS & CDS_complete titles
    if not (os.path.exists(alignment_fname) and os.path.exists(cds_fasta) and os.path.exists(complete_fasta)):
        print("No alignment or reference CDS found. Aborting in-CDS filtering step.")
        return # gRNA_hits
    else:
        alignment = fasta_to_dict(alignment_fname)
        cds_titles = sorted(list(fasta_to_dict(cds_fasta).keys()))
        complete_titles = sorted(list(fasta_to_dict(complete_fasta).keys()))
    cds_only_ranges = get_CDS_ranges(alignment_fname, cds_titles = cds_titles) ## dictionary
    cds_gaps_ranges = get_CDS_gaps_ranges(alignment_fname, cds_titles = cds_titles, ## set
                                          complete_titles = complete_titles)
    ref_seq_ids = [x for x in cds_titles + complete_titles if x in alignment]
    ref_plus = len([seq_id for seq_id in ref_seq_ids
                    if re.match(alignment_rvs_pattern, seq_id)]) < len(ref_seq_ids)/2
    get_target_CDS_ranges = make_target_CDS_ranges_function(cds_only_ranges, cds_gaps_ranges,
                                                            max_cds_insertion = max_cds_insertion)
    targets_cds_ranges = {seqid: get_target_CDS_ranges(alignment[seqid]) for seqid in alignment
                          if seqid not in ref_seq_ids}
    ## iterate through all gRNAs
    for gRNA_seq, coverage in gRNA_hits.hits().items():
        new_coverage = []
        ## assess each hit
        for gRNA_hit in coverage:
            seq_id, seq = [(seq_id, seq) for seq_id, seq in alignment.items()
                           if re.search(f"({alignment_rvs_pattern}|^)" + re.escape(gRNA_hit.target_id()) + "$",
                                        seq_id)][0]
            seq_plus = not re.match(alignment_rvs_pattern, seq_id)
            gRNA_hit.set_parent_sense('+' if seq_plus == ref_plus else '-') ## used to tie break set cover
            gRNA_range = gRNA_hit.range() if (not re.match(alignment_rvs_pattern, seq_id)) \
                         else gRNA_hit.reverse_range()
            ## append gRNAHit object if within at least 1 CDS
            if (within_CDS(targets_cds_ranges[seq_id], seq, gRNA_range,
                           min_within_n = min_within_n, min_within_percentage = min_within_percentage)):
                new_coverage.append(gRNA_hit)
                # screened_gRNA[gRNA] = set(new_coverage)
                gRNA_hit.set_cds_check(True)
            else:
                gRNA_hit.set_cds_check(False)
    return # {gRNA: coverage for gRNA, coverage in screened_gRNA.items() if len(coverage) > 0}

## remove gRNA which sequences the user specified to be unwanted
def filter_excluded_seqs(gRNA_hits, exclude_fname):
    print("Filtering user-specified sequences to exclude")
    to_exclude = set(str(x).upper() for x in fasta_to_dict(exclude_fname).values())
    gRNA_hits.set_all_seqs_check_by_function("exclude", lambda gRNA_seq: str(gRNA_seq) not in to_exclude)
    return # {seq: coverage for seq, coverage in gRNAs.items() if str(seq).upper() not in to_exclude}

def filter_unique_flank(gRNA_hits, length, background_fname, out_dir):
    # print("Filtering for unique flanking regions within {length} bp of each gRNA")
    # flanks_fname = f"{out_dir}/tmp_flanks.fasta"
    # flanks_seqs = dict(itertools.chain(*[{f"{gRNA.hit_id()}_front": gRNA.flank(length)[0],
    #                                       f"{gRNA.hit_id()}_back": gRNA.flank(length)[1]}.items()
    #                                      for hits in gRNAs.values() for gRNA in hits]))
    # dict_to_fasta(flanks_seqs, flanks_fname)
    # blast_hits = blastn(flanks_fname, background_fname)
    return # gRNAs ## to update with actual return value after function has been fully written


#########################
##  FILTER_BACKGROUND  ##
##   HELPER FUNCTIONS  ##
#########################


## function to get maximum pattern length
def pattern_max_len(pattern):
    ## count number of sets [], which get reduced to length of 1 each
    nsets = len(re.findall("(?<=\[).+?(?=\])", pattern))
    ## remove sets from pattern
    pattern_no_sets = ''.join(re.findall("(?<=^|\]).+?(?=$|\[)", pattern))
    ## count number of reps {n}, which get reduced to length of n - 1 each
    ## ## n-1 because the char to be repeated itself will count as the 1st rep
    reps = sum(map(lambda x: int(x.split(',')[-1]) - 1, re.findall("(?<=\{).+?(?=\})", pattern_no_sets)))
    ## remove reps from pattern
    nchars = sum(map(lambda x: len(x), re.findall("(?<=^|\}).+?(?=$|\{)", pattern_no_sets)))
    return nchars + nsets + reps ## sum

## helper function to judge if off-target crosses threshold
def make_exclude_function(pam, outside_targets, indexed_fa, max_mismatch, max_gap):
    pam_pattern = expand_ambiguous(infer_full_pam(pam))
    pam_pre, pam_post = pam_pattern.split('.')
    pam_pre_max = pattern_max_len(pam_pre)
    pam_post_max = pattern_max_len(pam_post)
    ## output function
    def exclude(entries, fasta_fname, fully_aligned_check = False, mismatch_check = True, gap_check = True):
        d_entries = {}
        for entry in entries:
            d_entries[entry["sseqid"]] = d_entries.get(entry["sseqid"], []) + [entry]
        output = []
        for molecule, entries in d_entries.items():
            ## check for presence of PAM of potentially excluded gRNA hits; if no PAM, it's not off-target
            seq = indexed_fa[fasta_fname][molecule]
            for x in entries:
                ## note that "pass" in this function means that a gRNA hit entry will be considered off-target
                gaps, qlen, qstart, qend = int(x["gaps"]), int(x["qlen"]), int(x["qstart"]), int(x["qend"])
                sstart, send = int(x["sstart"]), int(x["send"])
                ## get maximum PAM distance (unused gap allowance + uanligned gRNA length)
                gap_excess = max_gap - gaps
                len_excess = qlen - (qend - qstart + 1)
                pre_excess = int(x["qstart"]) - 1
                post_excess = int(x["qlen"]) - int(x["qend"])
                pre_len_max = pam_pre_max + gap_excess + pre_excess
                post_len_max = pam_post_max + gap_excess + post_excess
                ## get the sequences flanking the gRNA hit
                if x["sstrand"] == "plus" or x["sstrand"] == '+':
                    pre_seq = seq[sstart - pre_len_max : sstart]
                    post_seq = seq[send : send + post_len_max]
                else:
                    post_seq = seq[sstart - post_len_max: sstart].reverse_complement()
                    pre_seq = seq[send : send + pre_len_max].reverse_complement()
                ## assess whether pre-gRNA and post-gRNA regions match PAM pattern
                pre_pass = True if not pam_pre else bool(re.search(pam_pre, str(pre_seq)))
                post_pass = True if not pam_post else bool(re.search(pam_post, str(post_seq)))
                if pre_pass and post_pass: output.append(x)
        return output
    return exclude

def blast_mask(to_mask_fname, fasta_fname, fout_fname,
               header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                         'gaps', 'sstart', 'send', 'qlen', 'slen']):
    ## TODO: set word size to length of shortest sequence to mask
    ## masking criteria: perfect match from end to end
    is_to_mask = lambda x: x["pident"][:3] == "100" and x["qlen"] == x["length"] and x["gaps"] == '0'
    ## reformat blast output
    get_pos_dict = lambda x: {k: v for k, v in x.items() if k in ["sseqid", "qseqid", "sstart", "send"]}
    return [get_pos_dict(x)
            for x in blastn(to_mask_fname, fasta_fname, to_dict = True,
                            task = "megablast", fout = fout_fname) if is_to_mask(x)]

def mask_and_generate_outside(to_mask_fname, *background_fnames, mask_reference = True,
                              ref_genes_fasta = '', reference_fasta = '', out_dir,
                              header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                                        'gaps', 'sstart', 'send', 'qlen', 'slen'],
                              **kwargs): ## kwargs: sink for irrelevant keyword args
    ## tmp file; will be cleaned up by master function
    to_mask_hits = f"{out_dir}/tmp.tsv" ## reusable, cuz it's not gonna be read
    ## filter for valid bg fnames
    background_fnames = [fname for fname in background_fnames if
                         (fname and os.path.exists(fname))]
    ## mask things
    masked = {}
    if background_fnames:
        print("Masking targets (in query or user-provided background)")
        masked = {background_fname: blast_mask(to_mask_fname, background_fname, to_mask_hits)
                  for background_fname in background_fnames}
    if mask_reference and ref_genes_fasta and reference_fasta: ## mask in reference genome
        print("Masking gene(s) in reference")
        ## ALT: finish function to define masked regions in reference using GFF(_bed)
        masked = {**masked,
                  **{reference_fasta: blast_mask(ref_genes_fasta, reference_fasta, to_mask_hits)}}

    ## functions to determine if gRNA sequences is outside masked regions
    def outside_targets(to_screen, fasta_fname):
        if fasta_fname not in masked: return True
        for masked_target in masked[fasta_fname]:
            same_molecule = to_screen["sseqid"] == masked_target["sseqid"]
            within_target = within([to_screen[x] for x in ["sstart", "send"]],
                                   [masked_target[x] for x in ["sstart", "send"]])
            if same_molecule and within_target:
                return False
        return True
    return outside_targets


#############
##  BLAST  ##
#############

def blastn(seqs_fname, subject_fname, fout, to_dict = True, task = "blastn-short",
           header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                     'gaps', 'sstart', 'send', 'qlen', 'slen']):
    ## blast
    subprocess.check_output(("blastn", "-task", task, "-query", seqs_fname, "-subject", subject_fname,
                             "-out", f"{fout}", "-outfmt", r' '.join(["6", ' '.join(header)])))
    ## return output
    if not to_dict:
        return fout
    else:
        data = [x.rstrip("\n\r").split('\t') for x in open(fout).readlines()]
        output = [{header[i]: entry[i] for i in range(len(header))} for entry in data]
        return output


#################
##  ALIGNMENT  ##
##   THINGS    ##
#################

## returns a in b
def within(a, b):
    a = [int(x) for x in a]
    b = [int(x) for x in b]
    return min(a) >= min(b) and max(a) <= max(b)

## returns a in any range in ranges
def within_any(a, ranges):
    for r in ranges:
        if within(a, r):
            return True
    return False


def make_target_CDS_ranges_function(cds_only_ranges, cds_gaps_ranges, max_cds_insertion = 15):
    print("Max acceptable CDS insertion length:", max_cds_insertion)
    def get_target_CDS_ranges(seq):
        output = {}
        seq_pos = {i for i, c in enumerate(seq) if c != '-'}
        for seqid, gaps_ranges in cds_gaps_ranges.items():
            ## keep gaps if the bases of target in the gap is fewer than max_cds_insertion
            ## (i.e. if the target & ref gene were the only 2 in the alignment,
            ##  the gap would be <= max_cds_insertion)
            acceptable_gaps = []
            for gap_range in gaps_ranges:
                if len(seq_pos.intersection(ranges_to_pos(gaps_ranges))) <= max_cds_insertion:
                    acceptable_gaps.append(gap_range)
            ## generate final acceptable range for target-gene pair
            output[seqid] = ranges_union([cds_only_ranges[seqid] + acceptable_gaps])
        return output
    return get_target_CDS_ranges

## takes in a query_range, a sequence (aligned, presumably with gaps) to adjust query_range against, and a dictionary of {isoform: [tuple of CDS ranges]} and checks if a given query_range falls within a minimum number or percentage of isoforms' CDS
def within_CDS(cds_ranges, seq, query_range, min_within_n = 1, min_within_percentage = 0):
    ## cds_ranges, start, and end should be 0-indexed
    from fasta_manip import adjusted_pos
    query_range = [adjusted_pos(seq, x + 1) - 1 for x in query_range]
    within_check = [within_any(query_range, ranges) for ranges in cds_ranges.values()]
    within_count = within_check.count(True)
    return ((within_count >= min_within_n) and (within_count/len(within_check) >= min_within_percentage))

## convert ranges of [start, end) into set of positions
## e.g. [(1, 3), (6, 10)] --> {1, 2, 6, 7, 8, 9}
def ranges_to_pos(ranges):
    return set(itertools.chain(*[set(range(x[0], x[1])) for x in ranges]))

## convert set of positions into list of ranges
## e.g. {1, 2, 6, 7, 8, 9} --> [(1, 3), (6, 10)]
def pos_to_ranges(pos):
    if not pos: return []
    pos = sorted(pos)
    output = []
    start_v = pos[0]
    last_v = pos[0] - 1
    for i in range(len(pos)):
        v = pos[i]
        if v - last_v != 1:
            output.append((start_v, last_v + 1))
            start_v = v
        if i == len(pos) - 1:
            output.append((start_v, v + 1))
        last_v = v
    return output

## takes iterable of sets of pos and outputs a single set of pos
def pos_union(pos):
    return set(itertools.chain(*list(pos)))

def ranges_subtract(r1, r2):
    pos_subtract = ranges_to_pos(r1) - ranges_to_pos(r2)
    return pos_to_ranges(pos_subtract)

## e.g. [[(1, 3), (6, 9)], [(2, 3), (6, 10)]] --> {1, 2, 6, 7, 8, 9} --> [(1, 3), (6, 10)]
def ranges_union(ranges):
    pos_set = ranges_to_pos(itertools.chain(*ranges))
    return pos_to_ranges(pos_set)

## defined by CDS==genomic - CDS-only, per gene
def get_CDS_gaps_pos(alignment_fname, cds_titles = [], complete_titles = []):
    ## get reference pos
    cds_pos = get_CDS_pos(alignment_fname, cds_titles = cds_titles)
    shared_pos = get_CDS_pos_shared_with_complete(alignment_fname, cds_titles = cds_titles,
                                                  complete_titles = complete_titles)
    ## get potential cds insertion positions (union of potential insertions of all ref genes)
    insertion_pos = {seqid: shared_pos[seqid] - cds_pos[seqid] for seqid in cds_pos}
    return insertion_pos

## ranges defined CDS==genomic - CDS-only, per gene
def get_CDS_gaps_ranges(alignment_fname, cds_titles = [], complete_titles = []):
    pos = get_CDS_gaps_pos(alignment_fname, cds_titles = cds_titles, complete_titles = complete_titles)
    return {k: pos_to_ranges(v) for k, v in pos.items()}

def get_CDS_pos(alignment_fname, cds_titles = []):
    ranges = get_CDS_ranges(alignment_fname, cds_titles = cds_titles)
    return {k: ranges_to_pos(v) for k, v in ranges.items()}

def get_CDS_pos_shared_with_complete(alignment_fname, cds_titles = [], complete_titles = []):
    ranges = get_CDS_ranges_shared_with_complete(alignment_fname, cds_titles = cds_titles,
                                                 complete_titles = complete_titles)
    return {k: ranges_to_pos(v) for k, v in ranges.items()}

## ranges defined by CDS only
def get_CDS_ranges(alignment_fname, cds_titles = []):
    ## get reference seqs only
    alignment = {seq_id: seq.upper() for seq_id, seq in fasta_to_dict(alignment_fname).items()
                 if (seq_id in cds_titles)}
    cds_ranges = {cds_title: [] for cds_title in cds_titles}
    for cds_title in cds_titles:
        cds_seq = alignment.get(cds_title, '')
        if not cds_seq: continue
        ## get cds ranges
        in_cds = False
        range_start = 0
        for i, c in enumerate(cds_seq):
            if not in_cds and c != '-': ## initiate on CDS non-gap
                in_cds = True
                range_start = i
            elif in_cds and c == '-': ## terminate on CDS gap
                in_cds = False
                cds_ranges[cds_title].append((range_start, i))
            elif in_cds and i == len(cds_seq) - 1: ## terminate on last position
                cds_ranges[cds_title].append((range_start, i + 1))
    return cds_ranges

## previously, this range was defined when the old function has the 'relax' flag raised
def get_CDS_ranges_shared_with_complete(alignment_fname, cds_titles = [], complete_titles = []):
    ## get reference seqs only
    alignment = {seq_id: seq.upper() for seq_id, seq in fasta_to_dict(alignment_fname).items()
                 if (seq_id in cds_titles or seq_id in complete_titles)}
    cds_titles.sort()
    complete_titles.sort()
    cds_ranges = {cds_title: [] for cds_title in cds_titles}
    ## get CDS range for each isoform
    for i in range(len(cds_titles)):
        cds_title, complete_title = cds_titles[i], complete_titles[i]
        cds_seq, complete_seq = alignment.get(cds_title, ''), alignment.get(complete_title, '')
        if not (cds_seq and complete_seq):
            continue
        ## get starting position for isoform
        in_cds = False
        indel = 0 ## keep track of size of insertions in targets relative to reference genome
        range_start = 0
        for i in range(len(cds_seq)):
            cds_c, complete_c = cds_seq[i], complete_seq[i]
            if not in_cds and cds_c == complete_c != '-': ## initiate on non-gap identical
                in_cds = True
                indel = 0
                range_start = i
            elif in_cds and cds_c != complete_c: ## terminate on difference
                in_cds = False
                cds_ranges[cds_title].append((range_start, i - indel))
                indel = 0
            elif in_cds and indel and cds_c == complete_c != '-': ## reset indel if gap is flanked by CDS
                indel = 0
            elif in_cds and cds_c == complete_c == '-': ## if gap shared by CDS & genomic
                indel += 1
            elif in_cds and i == len(cds_seq) -1: ## terminate on last position
                cds_ranges[cds_title].append((range_start, i + 1))
    return cds_ranges

## the indices of cds_titles and their corresponding complete_titles must be identical
def get_CDS_ranges_old(alignment_fname, cds_titles = [], complete_titles = [],
                       relax = False, relax_cds_within = None):
    ## get reference seqs only
    alignment = {seq_id: seq.upper() for seq_id, seq in fasta_to_dict(alignment_fname).items()
                 if (seq_id in cds_titles or seq_id in complete_titles)}
    
    cds_titles.sort()
    complete_titles.sort()
    cds_ranges = {cds_title: [] for cds_title in cds_titles}
    
    ## get CDS range for each isoform
    for i in range(len(cds_titles)):
        cds_title, complete_title = cds_titles[i], complete_titles[i]
        cds_seq, complete_seq = alignment.get(cds_title, ''), alignment.get(complete_title, '')
        if not (cds_seq and complete_seq):
            continue
        ## get starting position for isoform
        in_cds = False
        relax_indel = 0 ## keep track of size of insertions in targets relative to reference genome
        range_start = 0
        for i in range(len(cds_seq)):
            if (not in_cds) and (((not relax) and cds_seq[i] != '-') or ## initiate on CDS non-gap
                                 (relax and \
                                  cds_seq[i] == complete_seq[i] != '-')): ## initiate on non-gap identical
                in_cds = True
                range_start = i
            elif in_cds and (((not relax) and cds_seq[i] == '-') or ## terminate on CDS gap
                             (relax and cds_seq[i] != complete_seq[i])): ## terminate on difference
                in_cds = False
                relax_indel = 0
                cds_ranges[cds_title].append((range_start, i))
            elif in_cds and (relax and relax_cds_within is not None and \
                             ( cds_seq[i] == complete_seq[i] == '-' )): ## if gap shared by CDS & genomic
                relax_indel += 1
                if relax_indel >= relax_cds_within: ## terminate on too-large indel
                    cds_ranges[cds_title].append((range_start, i - relax_indel + 1))
                    relax_indel = 0
                    in_cds = False
            elif in_cds and i == len(cds_seq) -1: ## terminate on last position
                cds_ranges[cds_title].append((range_start, i + 1))
    return cds_ranges


##################
##  SET COVER   ##
##  ALGORITHMS  ##
##################

## note that tie_breaker function should work on dictionaries of {gRNA_seq: {gRNAHit objects}} and return a tuple or list of two values: (gRNA_seq, {gRNAHit objects})
def set_cover(gRNA_hits, target_ids, algorithm = "LAR", id_key = lambda x: x, **kwargs):
    exclude_seqs = set(map(lambda s: str(s).upper(), kwargs["exclude_seqs"]))
    gRNA_coverage = {seq: hits for seq, hits in gRNA_hits.hits().items()
                     if str(seq).upper() not in exclude_seqs}
    ## check if set cover is possible before attempting to solve set cover
    if ( set(id_key(y) for x in gRNA_coverage.values() for y in x) & set(target_ids) ) < set(target_ids):
        print("\nWARNING: The provided gRNA sequences cannot cover all target sequences.\n")
    elif algorithm == "LAR":
        return set_cover_LAR(gRNA_coverage, target_ids, id_key = id_key,
                             **{k: v for k, v in kwargs.items() if k in ["tie_breaker"]})
    elif algorithm == "greedy":
        return set_cover_greedy(gRNA_coverage, target_ids, id_key = id_key,
                                **{k: v for k, v in kwargs.items() if k in ["tie_breaker"]})
    return []

## LAR algorithm
def set_cover_LAR(gRNA_coverage, target_ids, id_key = lambda x: x, tie_breaker = lambda x: tuple(x.items())[0]):
    
    main_sorting_key = lambda k, uncovered_count: len(set(id_key(x) for x in uncovered_count[k]))
    result_cover, covered = {}, set()
    from copy import deepcopy
    uncovered_count = deepcopy(gRNA_coverage)
    ## get set cover
    while {len(v) for v in uncovered_count.values()} != {0}:
        sorting_key = lambda k: main_sorting_key(k, uncovered_count)
        max_val = sorting_key(max(uncovered_count, key = sorting_key))
        max_items = {seq: hits for seq, hits in uncovered_count.items()
                     if sorting_key(seq) == max_val}
        seq, coverage = tie_breaker(max_items)
        coverage = set(id_key(target) for target in coverage)
        covered |= coverage
        result_cover[seq] = coverage
        ## update coverage of unchosen seqs
        uncovered_count = {seq: set(hit for hit in hits if id_key(hit) not in covered)
                           for seq, hits in uncovered_count.items()}
    ## remove redundant sequences
    for seq_id, coverage in result_cover.items():
        coverage_remaining = set().union(*[v for k, v in result_cover.items() if k != seq_id])
        if coverage.issubset(coverage_remaining):
            del result_cover[seq_id]
    
    return set(result_cover.keys())


## greedy algorithm
def set_cover_greedy(gRNA_coverage, target_ids, id_key = lambda x: x, tie_breaker = lambda x: x.items()[0]):
    
    main_sorting_key = lambda k, covered: len(set(id_key(x) for x in gRNA_coverage[k]) - set(covered))
    coverage_set = lambda k: set(id_key(x) for x in gRNA_coverage[k])
    covered, desired = set(), []
    while set(covered) != set(target_ids):
        sorting_key = lambda k: main_sorting_key(k, covered)
        max_val = sorting_key(max(gRNA_coverage, key = sorting_key))
        max_items = {k: v for k, v in gRNA_coverage.items() if sorting_key(k) == max_val}
        subset = tie_breaker(max_items)[0]
        covered |= coverage_set(subset)
        desired.append(subset)
    return set(desired)


###############
##  CLASSES  ##
###############

class CheckObj:
    def __init__(self, *check_names):
        self._checks = {check_name: None for check_name in check_names}
    def set_check(self, check_name, status):
        self._checks[check_name] = (True if status in ("pass", True) else
                                    False if status in ("fail", False) else
                                    None)
    def check(self, check_name, mode = "raw"):
        check_value = self._checks[check_name]
        if mode == "str":
            if check_value == True: return "pass"
            elif check_value == False: return "fail"
            else: return "NA"
        else:
            return check_value
    def some_checks_passed(self, *check_names, **kwargs):
        return all(self.check(check_name, **kwargs) for check_name in check_names)
    def some_valid_checks_passed(self, *check_names):
        return tuple(self.check(check_name) for check_name in check_names).count(False) == 0
    def all_valid_checks_passed(self):
        return self.some_valid_checks_passed(*self._checks.keys())
    def all_checks_passed(self):
        return self.some_checks_passed(*self._checks.keys())

class Target(CheckObj):
    ## seqs are stored in uppercase
    def __init__(self, seq, id = None, strand = None):
        self._seq = str(seq).upper()
        self._id = id
        self._strand = strand ## possible values: None, '+', '-'
        self._sense = None ## possible values: None, '+', '-'
    def __str__(self): return self.seq()
    def __repr__(self): return str(self)
    def __len__(self): return len(str(self))
    def id(self): return self._id
    def seq(self): return self._seq
    def strand(self): return self._strand
    def sense(self): return self._sense
    def set_strand(self, strand): self._strand = strand
    def set_sense(self, sense): self._sense = sense
    def set_sense_by_parent(self, parent_sense):
        parent_sense = ('+' if parent_sense in ('+', "sense") else '-')
        self._sense = (None if self.strand() == None else
                       '+' if parent_sense == self.strand() else '-')
    def parent_sense(self, mode = "raw"):
        parent_sense_value = ( None if (None in (self.sense(), self.strand())) else
                               '+' if (self.strand() == self.sense()) else
                               '-' )
        if mode == "str":
            if parent_sense_value == '+': return "sense"
            elif parent_sense_value == '-': return "antisense"
            else: return "NA"
        else:
            return parent_sense_value
    def valid_len(self): return ( len(self) > 0 )

class gRNASeq(CheckObj):
    ## seqs are stored in uppercase
    def __init__(self, seq):
        super().__init__("background", "exclude", "GC")
        self._seq = str(seq).upper()
        self._id = None
    def __str__(self): return self.seq()
    def __repr__(self): return str(self)
    def id(self): return self._id
    def seq(self): return self._seq
    def set_id(self, id): self._id = id
    def set_bg_check(self, status): self.set_check("background", status)
    def set_exclude_check(self, status): self.set_check("exclude", status)
    def set_gc_check(self, status = None, gc_min = 0, gc_max = 1):
        if status:
            self.set_check("GC", status)
        else:
            from seq_manip import gc_content
            self.set_check("GC", gc_min <= gc_content(self.seq()) <= gc_max)
        return

class gRNAHits:
    ## seqs are stored in uppercase
    def __init__(self, d = {}, gRNA_seqs = {}, gRNA_hits = {}):
        self._gRNAseqs = gRNA_seqs ## dictionary of {seq: <gRNASeq obj>}
        self._hits = gRNA_hits ## dictionary of {seq: [list of <gRNAHit obj>]}
        if d:
            self.parse_from_dict(d)
    def __repr__(self): return self.hits()
    def __len__(self): return len(self.seqs())
    def update_records(self): ## update dictionaries to remove any discrepancies
        self._gRNAseqs = {k: v for k, v in self.gRNAseqs() if k in self.hits()}
        self._hits = {k: v for k, v in self.hits() if k in self.gRNAseqs()}
    def copy(self):
        from copy import deepcopy
        new_obj = gRNAHits()
        new_obj._gRNAseqs = deepcopy(self.gRNAseqs())
        new_obj._hits = deepcopy(self.hits())
        return new_obj
    ################
    ##  BOOLEANS  ##
    ################
    def all_target_len_valid(self): return all(map(lambda hit: hit.target().valid_len(), self.flatten_hits()))
    ###############
    ##  PARSERS  ##
    ###############
    def parse_from_dict(self, d): ## where d = {str(seq): [list of <gRNAHit obj>s]}
        from copy import deepcopy
        self._gRNAseqs = {str(seq).upper(): gRNASeq(seq) for seq in d.keys()}
        self._hits = deepcopy(d)
    def parse_from_mapping(self, fname, targets = None, version = 1):
        ## read data
        seq_targets = {} if not targets else fasta_to_dict(targets)
        raw_mapping = [line.split('\t') for line in splitlines(fname)]
        header_mapping = raw_mapping[0]
        ## determine version where checks start
        if not version: version = (1 if "exclusive" in header_mapping[5] else \
                                   3 if "target length" in header_mapping else 2)
        check_col = header_mapping.index("group") + 1
        # if version == 1: check_col = 7
        # elif version == 2: check_col = 8
        header_checks = header_mapping[check_col:]
        dat_mapping = raw_mapping[1:]
        del raw_mapping
        for i, entry in enumerate(dat_mapping):
            ## note: gRNA_range is relative to + strand (not necessarily sense) when read from ..targets.txt file
            ## note: gRNA_id is probably meaningless as it will be overwritten using ids in the fasta file supplied to the variable 'fasta' in get_minimum_set_from_file
            if version == 1:
                gRNA_id, gRNA_seq, target_id, sense, strand, gRNA_range, group = entry[:check_col]
                gRNA_start, gRNA_end = map(int, re.search("\[(\d+), (\d+)\)", gRNA_range).group(1,2))
                target_len = 0
            elif version == 2:
                ## start, end are 1-indexed, start-inclusive and end-inclusive in v2
                gRNA_id, gRNA_seq, target_id, sense, strand, gRNA_start, gRNA_end, group = entry[:check_col]
                gRNA_start = int(gRNA_start) - 1 ## convert to 0-indexed, start-inclusive
                gRNA_end = int(gRNA_end) ## "convert" to 0-indexed, end-exclusive
                target_len = 0
            elif version == 3:
                ## start, end are 1-indexed, start-inclusive and end-inclusive in v2
                gRNA_id, gRNA_seq, target_id, target_len, sense, strand, gRNA_start, gRNA_end, group = entry[:check_col]
                gRNA_start = int(gRNA_start) - 1 ## convert to 0-indexed, start-inclusive
                gRNA_end = int(gRNA_end) ## "convert" to 0-indexed, end-exclusive
            try:
                ## note: dummy target sequence is used if targets file not provided; target strand assumed '+'
                target = Target(seq_targets.get(target_id, 'N' * int(target_len)), id = target_id, strand = '+')
            except Exception as e:
                print(version, i, header_mapping, entry)
                raise e
            ## create gRNAHit object
            gRNA_hit = gRNAHit(target, gRNA_start, gRNA_end, strand, gRNA_id)
            gRNA_hit.set_parent_sense(sense)
            ## add gRNAHit object to gRNAHits object
            self.add_hit(gRNA_seq, gRNA_hit)
            self.get_gRNAseq_by_seq(gRNA_seq).set_id(gRNA_id)
            ## log checks
            gRNA_checks = entry[check_col:]
            for i, check in enumerate(header_checks):
                if check == "CDS": ## this is the only one that's affected by hit location
                    gRNA_hit.set_check(check, gRNA_checks[i])
                else: ## all other checks apply to all hits with the same gRNA seq so set check to gRNASeq obj
                    self.get_gRNAseq_by_seq(gRNA_seq).set_check(check, gRNA_checks[i])
        return
    #################
    ##  MODIFIERS  ##
    #################
    def assign_gRNAseq_id(self, fasta):
        fasta_inv = {str(v): k for k, v in fasta_to_dict(fasta).items()}
        valid_seqs = set(fasta_inv.keys()) & set(self.seqs())
        if len(valid_seqs) < len(set(self.seqs())):
            print("\nWARNING: The provided FASTA file does not cover all gRNAs.\n")
        for seq in valid_seqs:
            self.get_gRNAseq_by_seq(seq).set_id(fasta_inv[seq])
        return
    def add_hit(self, seq, gRNA_hit):
        self.add_seq(seq)
        self._hits[str(seq)] = self.get_hits(seq) + [gRNA_hit]
    def add_seq(self, seq):
        if not str(seq) in self.gRNAseqs(): self._gRNAseqs[str(seq)] = gRNASeq(seq)
        if not str(seq) in self.hits(): self._hits[str(seq)] = []
    def remove_seqs(self, *seqs):
        if seqs and type(seqs[0]) in (list, tuple, set):
            seqs = list(itertools.chain(*seqs))
        for seq in seqs:
            if str(seq) in self.gRNAseqs(): del self._gRNAseqs[str(seq)]
            if str(seq) in self.hits(): del self._hits[str(seq)]
        return
    ###############
    ##  SETTERS  ##
    ###############
    def set_seqs_check(self, check_name, status, seqs):
        if type(seqs) not in (tuple, list, set):
            seqs = (seqs,)
        for seq in seqs:
            self.get_gRNAseq_by_seq(seq).set_check(check_name, status)
        return
    def set_seqs_check_by_function(self, check_name, func, seqs):
        for seq in seqs:
            self.set_seqs_check(check_name, func(self.get_gRNAseq_by_seq(seq)), [seq])
        return
    def set_all_seqs_check_by_function(self, check_name, func):
        self.set_seqs_check_by_function(check_name, func, self.seqs())
    def rename_seqs(self, fasta):
        seqs_names = {str(v): k for k, v in fasta_to_dict(fasta).items()}
        for seq, name in seqs_names.items():
            gRNA_seq = self.get_gRNAseq_by_seq(seq)
            if gRNA_seq:
                gRNA_seq.set_id(name)
        return
    def assign_seqid(self, prefix = "gRNA_", zfill = 3):
        for i, gRNA_seq in enumerate(self.flatten_gRNAseqs()):
            gRNA_seq.set_id(f"{prefix}{str(i+1).zfill(zfill)}")
        return
    ###############
    ##  GETTERS  ##
    ###############
    def gRNAseqs(self): return self._gRNAseqs
    def seqs(self, output_type = list): return output_type(self.gRNAseqs().keys())
    def hits(self): return self._hits
    def flatten_hits(self, output_type=list): return output_type(itertools.chain(*self.hits().values()))
    def flatten_gRNAseqs(self, output_type=list): return output_type(self.gRNAseqs().values())
    def get_hits(self, seq):
        return self.hits().get(str(seq).upper(), []) ## 'seq' must be able to be coerced using str()
    def get_gRNAseq_by_seq(self, seq):
        return self.gRNAseqs().get(str(seq).upper(), None) ## retrieve gRNASeq obj by seq
    def get_gRNAseq_by_id(self, id):
        output = [gRNA_seq for gRNA_seq in self.flatten_gRNAseqs() if gRNA_seq.id() == id]
        if not output: return None
        elif len(output) == 1: return output[0]
        else:
            print("\nWARNING: There are multiple gRNA sequences with the requested ID. Returning first in list.")
            return output[0]
    def get_gRNAseqs_by_seq(self, *seqs, output_type = list, ignore_invalid = True):
        if not seqs:
            return []
        if type(seqs[0]) in (list, tuple, set):
            seqs = tuple(itertools.chain(*seqs))
        return output_type(self.get_gRNAseq_by_seq(seq) for seq in seqs if self.get_gRNAseq_by_seq(seq) != None)
    def get_gRNAseqs_by_id(self, *ids, output_type = list, ignore_invalid = True):
        if type(ids[0]) in (list, tuple, set):
            ids = tuple(itertools.chain(*ids))
        return output_type(self.get_gRNAseq_by_id(id) for id in ids if self.get_gRNAseq_by_id(seq) != None)
    ######################
    ##      FILTER      ##
    ##  (& return new)  ##
    ######################
    ## filter gRNAHit objects for certain criteria and return new gRNAHits object
    def filter_hits(self, *check_names, exclude_empty_seqs = True, ignore_invalid = True):
        filtered = {seq: [hit for hit in hits
                          if ( (check_names and ( (ignore_invalid and hit.some_valid_checks_passed(*check_names))
                                                  or ( (not ignore_invalid) and
                                                       hit.some_checks_passed(*check_names) ) ) ) or
                               ( (not check_names) and ( (ignore_invalid and hit.all_valid_checks_passed())
                                                         or ( (not ignore_invalid) and
                                                              hit.all_checks_passed() ) ) ) ) ]
                    for seq, hits in self.hits().items()}
        if exclude_empty_seqs:
            filtered = {seq: hits for seq, hits in filtered.items() if hits}
        filtered_seqs = {seq: self.get_gRNAseq_by_seq(seq) for seq in filtered.keys()}
        output = gRNAHits(gRNA_seqs = filtered_seqs, gRNA_hits = filtered)
        return output
    def filter_hits_some_checks_passed(self, *check_names, **kwargs):
        return self.filter_hits(*check_names, **kwargs)
    def filter_hits_all_checks_passed(self, **kwargs):
        return self.filter_hits(**kwargs)
    ## filter gRNASeq objects for certain criteria and return new gRNAHits object
    def filter_seqs(self, *check_names, ignore_invalid = True):
        filtered = {seq: hits for seq, hits in self.hits().items()
                    if ((check_names and ((ignore_invalid and
                                           self.get_gRNAseq_by_seq(seq).some_valid_checks_passed(*check_names))
                                          or (not ignore_invalid and
                                              self.get_gRNAseq_by_seq(seq).some_checks_passed(*check_names)))) or
                         ((not check_names) and ((ignore_invalid and
                                                  self.get_gRNAseq_by_seq(seq).all_valid_checks_passed())
                                                 or (not ignore_invalid and
                                                     self.get_gRNAseq_by_seq(seq).all_checks_passed()))))}
        filtered_seqs = {seq: self.get_gRNAseq_by_seq(seq) for seq in filtered.keys()}
        output = gRNAHits(gRNA_seqs = filtered_seqs, gRNA_hits = filtered)
        return output
    def filter_seqs_some_checks_passed(self, *check_names, **kwargs):
        return self.filter_seqs(*check_names, **kwargs)
    def filter_seqs_all_checks_passed(self, **kwargs):
        return self.filter_seqs(**kwargs)
    #############
    ##  WRITE  ##
    #############
    def write_mapping(self, fout, sets = [], write_all = False,
                      write_checks = False, checks = ["background", "GC", "CDS"],
                      index = 1, start_incl = True, end_incl = True, version = 3,
                      fasta = None ): ## if fasta is provided, it will override default naming behaviour
        if not write_checks: checks = []
        if sets: sets = sets
        elif write_all: sets = [set(self.seqs())]
        else: []
        fasta_inv = {} if not fasta else {str(seq): k for k, seq in fasta_to_dict(fasta).items()}
        if version == 1:
            mapping_header = ["gRNA id", "gRNA sequence", "target id", "target sense",
                              "gRNA strand", "range (0-index, end exclusive)", "group"] + checks
        elif version == 2:
            mapping_header = ["gRNA id", "gRNA sequence", "target id", "target sense",
                              "gRNA strand", "start", "end", "group"] + checks
        else:
            mapping_header = ["gRNA id", "gRNA sequence", "target id", "target length", "target sense",
                              "gRNA strand", "start", "end", "group"] + checks
        mapping_dat = []
        seq_check_names = {"GC", "background"}
        for group, seqs in enumerate(sets):
            for seq in seqs:
                gRNA_seq = self.get_gRNAseq_by_seq(seq)
                gRNA_hits = self.get_hits(seq)
                for gRNA_hit in gRNA_hits:
                    ## put together all required fields
                    if version == 1:
                        mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id()), ## override id if fasta provided
                                            gRNA_seq.seq(), gRNA_hit.target_id(),
                                            gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand(),
                                            "[{}, {})".format(*gRNA_hit.range()), group + 1] +
                                           [(gRNA_seq.check(check_name, mode = "str")
                                             if check_name in seq_check_names else
                                             gRNA_hit.check(check_name, mode = "str"))
                                            for check_name in checks])
                    elif version == 2:
                        ## figure out start and end
                        start, end = map(lambda n: n + index, gRNA_hit.range())
                        if not start_incl: start -= 1
                        if end_incl: end -= 1
                        mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id()), ## override id if fasta provided
                                            gRNA_seq.seq(), gRNA_hit.target_id(),
                                            gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand(),
                                            start, end, group + 1] +
                                           [(gRNA_seq.check(check_name, mode = "str")
                                             if check_name in seq_check_names else
                                             gRNA_hit.check(check_name, mode = "str"))
                                            for check_name in checks])
                    else:
                        ## figure out start and end
                        start, end = map(lambda n: n + index, gRNA_hit.range())
                        if not start_incl: start -= 1
                        if end_incl: end -= 1
                        mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id()), ## override id if fasta provided
                                            gRNA_seq.seq(), gRNA_hit.target_id(), gRNA_hit.target_len(),
                                            gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand(),
                                            start, end, group + 1] +
                                           [(gRNA_seq.check(check_name, mode = "str")
                                             if check_name in seq_check_names else
                                             gRNA_hit.check(check_name, mode = "str"))
                                            for check_name in checks])
        mapping_dat.sort(key = lambda entry: int(entry[mapping_header.index("start")])) ## sort by start
        mapping_dat.sort(key = lambda entry: entry[mapping_header.index("target id")]) ## sort by target id
        mapping_dat.sort(key = lambda entry: entry[mapping_header.index("gRNA id")]) ## sort by gRNA id
        mapping_dat.sort(key = lambda entry: entry[mapping_header.index("group")])
        to_write = '\n'.join(['\t'.join(map(str, entry)) for entry in ([mapping_header] + mapping_dat)]) + '\n'
        open(fout, "w+").write(to_write)
        return
    def write_fasta(self, fout, seqs = [], ids = [], write_all = False, fasta = None):
        ## get relevant gRNA sequences
        if seqs: gRNA_seqs = self.get_gRNAseqs_by_seq(*seqs)
        elif ids: gRNA_seqs = self.get_gRNAseqs_by_id(*ids)
        elif write_all: gRNA_seqs = self.flatten_gRNAseqs()
        else: open(fout, "w+").write('')
        ## rename sequences per fasta file (if fasta file provided)
        fasta_inv = {} if not fasta else {str(seq): k for k, seq in fasta_to_dict(fasta).items()}
        to_write = {fasta_inv.get(gRNA_seq.seq(), gRNA_seq.id()): gRNA_seq.seq() for gRNA_seq in gRNA_seqs}
        if to_write:
            dict_to_fasta(to_write, fout)
        else:
            open(fout, "w+").write('')
        return

class gRNAHit(CheckObj):
    def __init__(self, target, start, end, strand, hit_id):
        ## note: unless something is weird, seq_strand is the same as strand
        super().__init__("background", "GC", "CDS", "exclude", "flank")
        self._target = target ## previously self._seq
        self._range = (start, end) ## relative to original (parent) sequence from which _target was derived
        self._strand = '+' if (strand.lower() == "fwd" or strand == '+') else '-' ## gRNA strand relative to original (parent) sequence from which _target was derived (gRNA strand should be same as _target strand)
        # self._seq_strand = '' # stores _target direction relative to original sequence from which _target was derived
        self._parent_sense = None ## - if original (parent) sequence from which _target was derived is on antisense strand, + if _target is on sense strand
        self._hit_id = hit_id
    def target(self): return self._target
    def start(self): return self._range[0]
    def end(self): return self._range[1]
    def target_id(self): return self.target().id()
    def strand(self): return self._strand
    def hit_id(self): return self._hit_id
    def range(self, output_type = tuple): return output_type(self._range)
    def target_len(self): return len(self.target())
    def target_strand(self): return self.target().strand()
    # def set_target_strand(self, strand): self._target_strand = strand
    def set_parent_sense(self, strand):
        self.target().set_sense_by_parent(strand)
        # self._parent_sense = ('+' if strand in ('+', "sense") else '-')
    def set_bg_check(self, status): self.set_check("background", status)
    def set_gc_check(self, status): self.set_check("GC", status)
    def set_cds_check(self, status): self.set_check("CDS", status)
    def set_exclude_check(self, status): self.set_check("exclude", status)
    def set_flank_check(self, status): self.set_check("flank", status)
    def reverse_range(self): return (self.target_len() - self.end(), self.target_len() - self.start())
    def parent_sense(self, mode = "raw"):
        parent_sense_value = self.target().sense()
        if mode == "str":
            if parent_sense_value == '+':
                return "sense"
            elif parent_sense_value == '-':
                return "antisense"
            else:
                return "NA"
        else:
            return parent_sense_value
    def adj_range(self, mode = "strand"):
        if (not self.strand() or
            (mode == "target" and not self.target_strand()) or
            (mode == "gene" and not self.parent_sense())):
            return (float("nan"), float("nan"))
        elif ((mode == "strand" and (self.strand() == '+')) or \
              (mode == "target" and (self.target_strand() == self.strand())) or \
              (mode == "gene" and (self.parent_sense() == '+'))):
            return self.range()
        return self.reverse_range()
    def flank(self, length = 100):
        target_seq = self.target() if self.strand == '+' else self.target().reverse_complement()
        start, end  = self.adj_range(mode = "target")
        return target_seq[start - length: start], target_seq[end: end + length]


#####################
##  OLD FUNCTIONS  ##
#####################

def find_cluster_gRNA_old(target_fname, pam = "GG", gRNA_len = 20): #Let NGG be GG, or NG be G
    """
    Finds all possible target sequences given pam.
    Returns: {<seq>: set(ids of targeted targets)}
    """
    dic_target = {}
    parsed = SeqIO.parse(target_fname, "fasta")
    hit_id = 0
    pam = pam.upper()
    for sequence in parsed:
        fwd_seq, rvs_seq = sequence.seq, sequence.seq.reverse_complement()
        target = Target(fwd_seq, id = sequence.id, strand = '+')
        for nucleotide in range(gRNA_len,len(sequence)+1):
            for strand, seq in {'+': fwd_seq, '-': rvs_seq}.items():
                start, end = nucleotide - gRNA_len, nucleotide
                if str(seq[nucleotide + 1: nucleotide + len(pam) + 1]).upper() == pam:
                    grna_seq = str(seq[nucleotide - gRNA_len: nucleotide]).upper()
                    if len(grna_seq) == gRNA_len:
                        ## gRNAHit's target is set as fwd_seq (variable 'target')
                        ## gRNAHit's range is set relative to fwd_seq (direction-less)
                        gRNAhit = gRNAHit(target,
                                          start if strand == '+' else len(seq) - end,
                                          end if strand == '+' else len(seq) - start,
                                          strand, hit_id)
                        dic_target[grna_seq] = dic_target.get(grna_seq, set()).union({gRNAhit})
                        hit_id += 1
    # return dic_target
    output_gRNA_hits = gRNAHits()
    output_gRNA_hits.parse_from_dict(dic_target)
    return output_gRNA_hits

## filter for gRNA off-target hits in background
def filter_background_old(seqs_fname, target_fname, background_fname, gRNA_hits, screen_reference = False,
                      fout_pref = "", out_dir = "", max_mismatch = 0, max_gap = 0, mask_reference = True,
                      reference_fasta = '', ref_genes_fasta = False, outside_targets = None, report_bg = False,
                      pam = "GG", nonref_mask_fname = '',
                      **kwargs): ## kwargs is for arguments for the mask_and_generate_outside function
    """
    Executes bl2seq with candidate gRNA sequences as query and background as subject.
    Filters for maximum allowable mismatches and gaps.
    Removes candidate gRNA that pass the above criteria but are also found outside of target sequences.
    Inputs:
        seqs_fname: path to FASTA file of gRNA candidate sequences
        target_fname: path to FASTA file of target sequences
        background_fname: path to FASTA file of all background sequences (can include targets)
        gRNA_hits: gRNAHits object
    // DEFUNCT // Returns: filtered dictionary of {<gRNA seq>: set(ids of targeted target)}
    Returns: None. gRNAHits is modified in-place.
    """
    ## prep for blast
    blast_output_header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                           'gaps', 'sstart', 'send', 'qlen', 'slen']
    ## parse file names ( only handles a limited number of combinations :/ )
    if not fout_pref:
        fout_pref = os.path.join(out_dir, "proj")
    elif fout_pref and out_dir and not os.path.dirname(fout_pref):
        fout_pref = os.path.join(out_dir, fout_pref)
    ## create function that masks targets and checks if gRNA hit is in or out of targets if not provided
    if not outside_targets:
        ## if nonref_mask_fname (FASTA file of non-targets that are to be masked) is provided
        if nonref_mask_fname and os.path.exists(nonref_mask_fname):
            ## merge nonref_mask_fname & target_fname into nonref_mask_fname
            to_mask = {**fasta_to_dict(nonref_mask_fname), **fasta_to_dict(target_fname)}
            to_mask_fname = nonref_mask_fname
            dict_to_fasta(to_mask, nonref_mask_fname)
        ## else just use target_fname as source of sequences to mask
        else:
            to_mask_fname = target_fname
        ## mask
        outside_targets = mask_and_generate_outside(to_mask_fname, background_fname,
                                                    out_dir = out_dir,
                                                    mask_reference = mask_reference,
                                                    ref_genes_fasta = ref_genes_fasta,
                                                    reference_fasta = reference_fasta,
                                                    header = blast_output_header)
    ## helper function to judge if off-target crosses threshold
    def exclude(x, fasta_fname, fully_aligned_check = True, mismatch_check = True, gap_check = True):
        fully_aligned = (int(x["length"]) + int(x["gaps"]) == int(x["qlen"])) \
                        if fully_aligned_check else True
        mismatch_pass = (int(x["mismatch"]) <= max_mismatch) if mismatch_check else True
        gap_pass = (int(x["gaps"]) <= max_gap) if gap_check else True
        return fully_aligned and mismatch_pass and gap_pass and outside_targets(x, fasta_fname)
    ## function to prepend header if user wants a bg report, else delete blast hits file
    def resolve_bg_blast(fname):
        from file_manip import prepend
        if report_bg: prepend(fname, '\t'.join(blast_output_header) + '\n')
        else: os.remove(fname)
    ## function to get seqs to exclude
    def get_excl_seqid(seqs_fname, fasta_fname, fout_fname):
        return set(x["qseqid"]
                   for x in blastn(seqs_fname, fasta_fname, to_dict = True, fout = fout_fname,
                                   header = blast_output_header)
                   if exclude(x, fasta_fname))
    ## run blast to identify gRNA hits in background
    if background_fname and os.path.exists(background_fname):
        nonref_hits = f"{fout_pref}_hit_nonref.tsv"
        excl_seqid = get_excl_seqid(seqs_fname, background_fname, nonref_hits)
        resolve_bg_blast(nonref_hits)
    else:
        excl_seqid = set()
    if screen_reference and reference_fasta: ## find in reference
        ref_hits = f"{fout_pref}_hit_ref.tsv"
        excl_seqid |= get_excl_seqid(seqs_fname, reference_fasta, ref_hits)
        resolve_bg_blast(ref_hits)
    ## parse bg check status
    print("Filtering gRNA sequences")
    seqs_fasta = fasta_to_dict(seqs_fname)
    seqs_screened = set(map(lambda s: str(s).upper(), seqs_fasta.values()))
    seqs_failed = set(str(seqs_fasta[seqid]).upper() for seqid in excl_seqid)
    seqs_passed = set(seq for seq in seqs_screened if seq not in seqs_failed)
    ## record bg check status
    gRNA_hits.set_seqs_check("background", True, seqs_passed)
    gRNA_hits.set_seqs_check("background", False, seqs_failed)
    return # gRNA_hits


## filter within CDS
def filter_in_cds_old(gRNA_hits, alignment_fname, cds_fasta, complete_fasta,
                      relax = False, relax_cds_within = None,
                      alignment_rvs_pattern = '^_R_', min_within_n = 1, min_within_percentage = 0,
                      max_cds_insertion = 0):
    print("Filtering for sequences within reference CDS")
    ## get CDS & CDS_complete titles
    if not (os.path.exists(alignment_fname) and os.path.exists(cds_fasta) and os.path.exists(complete_fasta)):
        print("No alignment or reference CDS found. Aborting in-CDS filtering step.")
        return # gRNA_hits
    else:
        alignment = fasta_to_dict(alignment_fname)
        cds_titles = sorted(list(fasta_to_dict(cds_fasta).keys()))
        complete_titles = sorted(list(fasta_to_dict(complete_fasta).keys()))
    cds_ranges = get_CDS_ranges_old(alignment_fname, cds_titles, complete_titles,
                                    relax = relax, relax_cds_within = relax_cds_within)
    ref_seq_ids = [x for x in cds_titles + complete_titles if x in alignment]
    ref_plus = len([seq_id for seq_id in ref_seq_ids
                    if re.match(alignment_rvs_pattern, seq_id)]) < len(ref_seq_ids)/2
    ## iterate through all gRNAs
    for gRNA_seq, coverage in gRNA_hits.hits().items():
        new_coverage = []
        ## assess each hit
        for gRNA_hit in coverage:
            seq_id, seq = [(seq_id, seq) for seq_id, seq in alignment.items()
                           if re.search(f"({alignment_rvs_pattern}|^)" + re.escape(gRNA_hit.target_id()) + "$",
                                        seq_id)][0]
            seq_plus = not re.match(alignment_rvs_pattern, seq_id)
            gRNA_hit.set_parent_sense('+' if seq_plus == ref_plus else '-') ## used to tie break set cover
            gRNA_range = gRNA_hit.range() if (not re.match(alignment_rvs_pattern, seq_id)) \
                         else gRNA_hit.reverse_range()
            ## append gRNAHit object if within at least 1 CDS
            if (within_CDS(cds_ranges, seq, gRNA_range,
                           min_within_n = min_within_n, min_within_percentage = min_within_percentage)):
                new_coverage.append(gRNA_hit)
                # screened_gRNA[gRNA] = set(new_coverage)
                gRNA_hit.set_cds_check(True)
            else:
                gRNA_hit.set_cds_check(False)
    return # {gRNA: coverage for gRNA, coverage in screened_gRNA.items() if len(coverage) > 0}


################
##  OBSOLETE  ##
################

# def grna_in_background(seqs_fname, target_fname, background_fname, gRNA_hits, screen_reference = False,
#                        fout_pref = "", out_dir = "", max_mismatch = 0, max_gap = 0, mask_reference = True,
#                        reference_fasta = "", ref_genes_fasta = False, outside_targets = None, report_bg = False,
#                        pam = "GG",
#                        **kwargs): ## kwargs is for arguments for the mask_and_generate_outside function
#     ## create function that masks targets and checks if gRNA hit is in or out of targets if not provided
#     if not outside_targets:
#         outside_targets = mask_and_generate_outside(target_fname, background_fname,
#                                                     out_dir = out_dir,
#                                                     mask_reference = mask_reference,
#                                                     ref_genes_fasta = ref_genes_fasta,
#                                                     reference_fasta = reference_fasta)
#     ## prep for exclusion checks
#     from fasta_manip import fasta_to_dict
#     bg_seqs = {fname: fasta_to_dict(fname) for fname in
#                [f for f in (background_fname, reference_fasta) if f]}
#     ## function to check if PAM is present downstream of off-target gRNA match
#     ## ## PAM must be within (+-(max gap + N + PAM length) from end to be considered a match
#     ## ## if hit covers end of gRNA, that's the end.
#     ## ## if hit DOESN'T cover end of gRNA, extend by <length of grna - pos of end of hit> to get end.
#     def has_pam(entry, fasta):
#         pass
#     ## helper function to judge if off-target crosses threshold
#     def exclude(x, fasta_fname, fully_aligned_check = True, mismatch_check = True, gap_check = True):
#         fully_aligned = (int(x["length"]) + int(x["gaps"]) == int(x["qlen"])) \
#                         if fully_aligned_check else True
#         mismatch_pass = (int(x["mismatch"]) <= max_mismatch) if mismatch_check else True
#         gap_pass = (int(x["gaps"]) <= max_gap) if gap_check else True
#         return fully_aligned and mismatch_pass and gap_pass and outside_targets(x, fasta_fname)
#     ## prep for blast
#     blast_output_header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
#                            'gaps', 'sstart', 'send', 'qlen', 'slen']
#     ## function to prepend header if user wants a bg report, else delete blast hits file
#     def resolve_bg_blast(fname):
#         from file_manip import prepend
#         if report_bg: prepend(fname, '\t'.join(blast_output_header) + '\n')
#         else: os.remove(nonref_hits)
#     ## function to get seqs to exclude
#     def get_excl_seqid(seqs_fname, fasta_fname, fout_fname):
#         return set(x["qseqid"]
#                    for x in blastn(seqs_fname, fasta_fname, to_dict = True, fout = fout_fname,
#                                    header = blast_output_header)
#                    if exclude(x, fasta_fname))
#     ## run blast to identify gRNA hits in background
#     if background_fname and os.path.exists(background_fname):
#         nonref_hits = f"{fout_pref}_hit_nonref.tsv"
#         excl_seqid = get_excl_seqid(seqs_fname, background_fname, nonref_hits)
#         resolve_bg_blast(nonref_hits)
#     else:
#         excl_seqid = set()
#     if screen_reference and reference_fasta: ## find in reference
#         ref_hits = f"{fout_pref}_hit_ref.tsv"
#         excl_seqid |= get_excl_seqid(seqs_fname, reference_fasta, ref_hits)
#         resolve_bg_blast(ref_hits)
#     return exlc_seqid


####################################
##  CHOOSE YOUR OWN ADVENTURE :)  ##
####################################

# find_cluster_gRNA = find_cluster_gRNA_old
# filter_background = filter_background_old

find_cluster_gRNA = find_cluster_gRNA_gen
filter_background = filter_background_gen
filter_in_cds = filter_in_cds_gen ## TODO: test filter_in_cds_gen
