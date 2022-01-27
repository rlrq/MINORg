import os
import itertools
import regex as re

from Bio import SeqIO, SearchIO
from pathlib import Path

from minorg.functions import (
    assign_alias, blast,
    within_any, gc_content,
    fasta_to_dict, dict_to_fasta, find_identical_in_fasta,
    tsv_entries_as_dict, prepend,
    ranges_union, ranges_to_pos, ranges_subtract, ranges_intersect, within, convert_range,
    adjusted_feature_ranges, adjusted_pos
)

from minorg.index import IndexedFasta
from minorg.annotation import GFF
from minorg.pam import PAM

## TODO: report location of feature(s) for each reference gene

# ## gRNA_hits is a gRNAHits object
# ## target_names and features are iterables of strings (e.g. ["CDS", "five_prime_UTR"])
# ## - target_names should be list of names of target sequences in fasta_alignment
# ##   - if target_names is not provided, all targets in gRNA_hits will be assumed to be desired targets
# ## fasta_alignment, gff_bed, and fasta_exclude are paths to files.
# ## - in fasta_alignment, all sequences of reference genes should be named according to:
# ##   'Reference|[^|]|[^|]|gene|[^|]|<gene id>'
# ##   - six fields; '|'-delimited
# ##   - 'Reference' in 1st field, 'gene' in 4th field, gene ID in 6th and last field
# ##   - fields 2,3,5 can be anything as long as '|' does not appear anywhere in those fields
# ##   - e.g. Reference|Reference|NB-ARC|gene|complete|AT1G01010
# ## fasta_reference and fasta_background can be either:
# ##  str (path), dictionary ({<alias>: <path>}), or iterable ([<path>, <path>, ...])
# ## - in the case of str, the path will be assigned alias <ref/bg>1
# ## - in the case of non-dict iterable, paths will be assigned alias according to <ref/bg><index in iterable>+1
# ## This function does not perform argument checking so all args are assumed to be in correct format and valid
# def filter_grna(gRNA_hits,
#                 mk_fname = lambda x: x,
#                 checks = ["GC", "feature", "background", "flank", "exclude"],
#                 # ## shared args
#                 # target_names = None,
#                 ## args for filter_gc
#                 gc_min = 0, gc_max = 1,
#                 ## args for filter_in_feature
#                 fasta_alignment: Path = None, features = None, gff_bed: Path = None,
#                 max_insertion = 15, min_within_n = 1, min_within_fraction = 0,
#                 alignment_rvs_pattern = "^_R_$seqid$$", ref: bool = False,
#                 domain_gff_bed: Path = None,
#                 ## args for filter_background
#                 fasta_grna: Path = None, fasta_target: Path = None,
#                 fasta_reference = None, fasta_background = None, fasta_mask: Path = None, 
#                 max_mismatch = 0, max_gap = 0, pam = 'N',
#                 screen_reference = False, mask_reference = True, pamless_check = False, report_bg = False,
#                 ## args for filter_exclude
#                 fasta_exclude: Path = None,
#                 **kwargs): ## kwargs is for arguments for the mask_and_generate_outside function)
#     checks = set(check.upper() for check in checks)
#     if "GC" in checks:
#         filter_gc(gRNA_hits, gc_min = gc_min, gc_max = gc_max)
#     if ( "FEATURE" in checks and features is not None and gff_bed is not None ):
#         filter_in_feature_gen(gRNA_hits, fasta_alignment, gff_bed,
#                               features = features, max_insertion = max_insertion, ref = ref,
#                               min_within_n = min_within_n, min_within_fraction = min_within_fraction)
#     if "BACKGROUND" in checks:
#         ## TODO: parse fasta_reference, fasta_background into dictionaries
#         filter_background_gen(gRNA_hits, fasta_grna, fasta_target, fasta_background, gff_bed = None,
#                               nonref_mask_fname = fasta_mask, reference_fasta = fasta_reference,
#                               max_mismatch = max_mismatch, max_gap = max_gap, pam = pam,
#                               screen_reference = screen_reference, mask_reference = mask_reference,
#                               report_bg = report_bg, pamless_check = pamless_check,
#                               mk_fname = mk_fname)
#     if "EXCLUDE" in checks and fasta_exclude is not None:
#         filter_excluded_seqs(gRNA_hits, args.exclude)
#     # if "FLANK" in checks:
#     #     filter_unique_flank(gRNA_hits, flank_length, background_fname, out_dir)
#     screened_gRNA = gRNA_hits.filter_seqs_all_checks_passed(ignore_invalid = True).filter_hits_all_checks_passed(ignore_invalid = True)
#     ## return only filtered gRNA set. since gRNA_hits is passed as an object, it'll be modified directly.
#     return screened_gRNA

## filter for gRNA off-target hits in background
## fasta_background, fasta_mask, fasta_target, fasta_reference can be dict of {<alias>: <path>} or [<path>].
def filter_background_gen(gRNA_hits, fasta_grna, fasta_target, fasta_background = None, gff_bed = None,
                          fasta_mask = None, fasta_reference = None, fasta_ref_genes = False,
                          max_mismatch = 0, max_gap = 0, pam = "GG",
                          screen_reference = False, mask_reference = True,
                          report_bg = False, pamless_check = False,
                          mk_fname = lambda x: x, blastn = "blastn",
                          ## TODO: all args below this (except kwargs) should be removed/repurposed
                          fout_pref = "", out_dir = "", fout_mask = None,
                          **kwargs): ## kwargs is for arguments for the mask_and_generate_outside function
    """
    Executes bl2seq with candidate gRNA sequences as query and background as subject.
    Filters for maximum allowable mismatches and gaps.
    Removes candidate gRNA that pass the above criteria but are also found outside of target sequences.
    Inputs:
        fasta_grna: path(s) to FASTA file of gRNA candidate sequences
        fasta_target: path(s) to FASTA file of target sequences
        fasta_background: path(s) to FASTA file of all background sequences (can include targets)
        gRNA_hits: gRNAHits object
    // DEFUNCT // Returns: filtered dictionary of {<gRNA seq>: set(ids of targeted target)}
    Returns: None. gRNAHits is modified in-place.
    """
    
    ## parse file names ( only handles a limited number of combinations :/ )
    if not mk_fname:
        if not out_dir:
            import os
            out_dir = os.getcwd()
        if not fout_pref:
            fout_pref = os.path.join(out_dir, "minorg")
        elif fout_pref and out_dir and not os.path.dirname(fout_pref):
            fout_pref = os.path.join(out_dir, fout_pref)
        mk_fname = lambda *x: os.path.join(fout_pref, *x)
    ## parse fasta_background and fasta_reference into {<alias>: <path>} dicts
    fasta_background = assign_alias(fasta_background, mk_name = lambda i: "bg_{str(i).zfill(3)}")
    fasta_reference = assign_alias(fasta_reference, mk_name = lambda i: "ref_{str(i).zfill(3)}")
    reference_fnames = {alias: (fasta if not isinstance(fasta, IndexedFasta) else fasta.filename)
                        for alias, fasta in fasta_reference.items()}
    ## combine fasta_mask, fasta_ref_genes & fasta_target into single dict
    to_mask = {**assign_alias(fasta_mask, mk_name = lambda i: f"mask_{str(i).zfill(3)}"),
               **assign_alias(fasta_target, mk_name = lambda i: f"target_{str(i).zfill(3)}")}
    ## create function that masks targets and checks if gRNA hit is in or out of targets if not provided
    outside_targets = mask_and_generate_outside(to_mask, out_dir = out_dir, blastn = blastn,
                                                mask_reference = mask_reference,
                                                background_fnames = fasta_background,
                                                reference_fnames = reference_fnames,
                                                mk_fname = mk_fname,
                                                fout_mask = fout_mask, new_file = True)
    ## index fasta files
    # indexed_fa = {fname: SeqIO.index(fname, "fasta") for fname in
    #               itertools.chain(*[fnames.values() for fnames in [fasta_background, fasta_reference]])}
    indexed_fa = {fname: IndexedFasta(fname) for fname in
                  itertools.chain(*[fnames.values() for fnames in [fasta_background, fasta_reference]])}
    ## make exclude function (function that returns True if gRNA fails off-target check)
    has_pam = (lambda *args, **kwargs: True) if pamless_check else make_has_pam(pam)
    # exclude_with_pam = make_exclude_function(pam, outside_targets, indexed_fa, max_mismatch, max_gap)
    ## prep for blast
    blast_output_header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                           'gaps', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'sstrand']
    ## function to prepend header if user wants a bg report, else delete blast hits file
    def resolve_bg_blast(fname):
        import os
        if report_bg: prepend(fname, '\t'.join(blast_output_header) + '\n')
        else: os.remove(fname)
    ## function to get seqs to exclude
    def get_excl_seqid(fasta_query, fasta_subject, fout_fname, fully_aligned_check = True):
        def excl_by_loc(hsp):
            return outside_targets(hsp, fasta_subject)
        def excl_by_aln(query_result, hsp, mismatch_check = True, gap_check = True):
            qlen = query_result.seq_len
            mismatch_excess = max_mismatch - (qlen - hsp.ident_num)
            gap_excess = max_gap - hsp.gap_num
            unaligned_excess = qlen - hsp.query_end + hsp.query_start
            mismatch_pass = (mismatch_excess >= 0) if mismatch_check else True
            gap_pass = (gap_excess >= 0) if gap_check else True
            fully_aligned = unaligned_excess == 0 if fully_aligned_check else True
            return fully_aligned and mismatch_pass and gap_pass
        def excl_by_pam(query_result, hsp):
            return has_pam(fasta_subject, query_result, hsp)
        def exclude(query_result, hsp, mismatch_check = True, gap_check = True):
            return (excl_by_loc(hsp) and excl_by_aln(query_result, hsp) and \
                excl_by_pam(query_result, hsp))
        from Bio.Blast.Applications import NcbiblastnCommandline
        blast(NcbiblastnCommandline, header = blast_output_header, fout = fout_fname, outfmt = 5,
              query = fasta_query, subject = fasta_subject, cmd = blastn, task = "blastn-short")
        return set(hsp.query_id
                   for query_result in SearchIO.parse(fout_fname, "blast-xml")
                   for hit in query_result for hsp in hit.hsps
                   if exclude(query_result, hsp))
    ## run blast to identify gRNA hits in background
    excl_seqid = set()
    if screen_reference and fasta_reference: ## find in reference
        print("Assessing off-target in reference genome")
        for i, ref_fname in enumerate(reference_fnames.values()):
            ref_hits = mk_fname(f"hit_ref_{i}.xml")
            excl_seqid |= get_excl_seqid(fasta_grna, ref_fname, ref_hits)
            resolve_bg_blast(ref_hits)
        ## remove reference FASTA from fasta_background
        ## - ensures that each file is screened only once if reference was also queried for gRNA
        fasta_background = {alias: fname for alias, fname in fasta_background.items()
                            if fname not in reference_fnames.values()}
    if fasta_background:
        print("Assessing off-target in user-provided sequences")
        for bg_fname in fasta_background.values():
            nonref_hits = mk_fname("hit_nonref.xml")
            excl_seqid |= get_excl_seqid(fasta_grna, bg_fname, nonref_hits, fully_aligned_check = False)
            resolve_bg_blast(nonref_hits)
    ## parse bg check status
    print("Filtering gRNA sequences")
    grna_seqs = fasta_to_dict(fasta_grna)
    # grna_screened = set(map(lambda s: str(s).upper(), grna_seqs.values()))
    grna_screened = tuple(grna_seqs.values())
    grna_failed = set(str(grna_seqs[seqid]).upper() for seqid in excl_seqid)
    grna_passed = set(str(seq).upper() for seq in grna_screened if str(seq).upper() not in grna_failed)
    ## record bg check status
    gRNA_hits.set_seqs_check("background", True, grna_passed)
    gRNA_hits.set_seqs_check("background", False, grna_failed)
    return # gRNA_hits

## filter GC content
def filter_gc(gRNA_hits, gc_min = 0, gc_max = 1):
    for gRNA_seq in gRNA_hits.gRNAseqs.values():
        gRNA_seq.set_gc_check(gc_min = gc_min, gc_max = gc_max)
    return

## filter within feature
def filter_in_feature_gen(gRNA_hits, fasta_alignment, gff_beds, features = None, target_names = None,
                          max_insertion = 15, min_within_n = 1, min_within_fraction = 0, ref = False,
                          domain_gff_bed = None, alignment_rvs_pattern = "^_R_$seqid$$", memsave = True,
                          seqid_source_pattern = f"(?<=^Reference\\|)[^|]+(?=\\|)"):
    '''
    gff_beds and domain_gff_bed must be dict of {<alias>: <path to file>}
    '''
    print("Filtering for sequences within reference feature")
    if not (os.path.exists(fasta_alignment) and sum(map(os.path.exists, gff_beds.values())) > 0):
        print("No alignment or reference annotations found. Aborting in-feature filtering step.")
        return # gRNA_hits
    else:
        alignment = fasta_to_dict(fasta_alignment)
        # genes = {re.search("^.+(?=\d+-\d+)", '|'.join(seqid.split('|')[6:-1])).group(0):
        genes = {'|'.join(seqid.split('|')[5:-1]):
                 seqid for seqid in alignment
                 if (seqid.split('|')[0] == "Reference" and seqid.split('|')[4] == "complete")}
        ## fmt=None to let GFF obj detect format
        anns = {alias: GFF(fname = gff_bed, quiet = True, memsave = memsave, fmt = None)
                for alias, gff_bed in gff_beds.items()}
        gene_anns = {alias: {gene: ann.get_id(gene, output_list = False) for gene in genes}
                     for alias, ann in anns.items()}
        if domain_gff_bed:
            ## if --domain is used, this function will be supplied with domain_gff_bed, which contains
            ##   feature 'domain'. All fields will be identical to that of the gene with exception
            ##   of feature ('domain'), ranges (whatever the start and end of the domain are),
            ##   and source ('minorg').
            domain_anns = GFF(fname = domain_gff_bed, quiet = True,
                              fmt = ("BED" if domain_gff_bed.split('.')[-1].upper()=="BED" else "GFF3"))
        if features is None:
            feature_ranges = {source: {gene: [(feature.start - 1, feature.end)]
                                       for gene, feature in gene_ann.items()}
                              for source, gene_ann in gene_anns.items()}
        else:
            feature_ranges = {source: {gene: ranges_union([[(x.start - 1, x.end) for x in
                                                            ann.get_subfeatures_full(gene, features)]])
                                       for gene in genes}
                              for source, ann in anns.items() if ann is not None}
    def adjust_feature_ranges(gene, seqid, **kwargs):
        source = re.search(seqid_source_pattern, seqid).group(0)
        gene_ann = gene_anns[source]
        gene_feature_ranges = feature_ranges[source]
        if domain_gff_bed:
            domain_ann = [entry for entry in domain_anns.get_id(gene, output_list = True)
                          if entry.source == source]
            return adjusted_feature_ranges(alignment[seqid],
                                           ranges_intersect([(gene_ann[gene].start-1, gene_ann[gene].end)],
                                                            (ranges_union([[(x.start-1, x.end)] for x in
                                                                           domain_ann])))[0],
                                           gene_feature_ranges[gene],
                                           strand = gene_ann[gene].strand,
                                           **kwargs)
        else:
            return adjusted_feature_ranges(alignment[seqid],
                                           (gene_ann[gene].start-1, gene_ann[gene].end),
                                           gene_feature_ranges[gene],
                                           strand = gene_ann[gene].strand,
                                           **kwargs)
    def is_rvs(aln_seqid, seqid):
        return re.search(Template(alignment_rvs_pattern).substitute(seqid = re.escape(seqid)), aln_seqid)
    feature_only_ranges = {seqid: adjust_feature_ranges(gene, seqid, subtract_gaps = True)
                           for gene, seqid in genes.items()}
    feature_gaps_ranges = {seqid: ranges_subtract(adjust_feature_ranges(gene, seqid, subtract_gaps = False),
                                                  feature_only_ranges[seqid])
                           for gene, seqid in genes.items()}
    ## check orientation of reference genes (should be + but we're double checking to be sure)
    ref_seq_ids = list(genes.values())
    ref_plus = len([seq_id for seq_id in ref_seq_ids
                    if re.match(alignment_rvs_pattern.replace("$seqid", ".+"), seq_id)]) < len(ref_seq_ids)/2
    ## define acceptable ranges in targets
    get_target_feature_ranges = make_target_feature_ranges_function(feature_only_ranges, feature_gaps_ranges,
                                                                    max_insertion = max_insertion)
    targets_feature_ranges = {seqid: get_target_feature_ranges(alignment[seqid], seqid) for seqid in alignment}
    if target_names is not None: target_names = set(target_names)
    from string import Template
    ## iterate through all gRNAs
    for gRNA_seq, coverage in gRNA_hits.hits.items():
        # new_coverage = []
        ## assess each hit
        for gRNA_hit in coverage:
            seq_id_seq = [(seq_id, seq) for seq_id, seq in alignment.items()
                          if (seq_id == gRNA_hit.target_id or is_rvs(seq_id, gRNA_hit.target_id))]
            ## if target not among user-specified targets or not in alignment, skip to next gRNA hit
            if ( (target_names is not None and gRNA_hit.target_id not in target_names) or
                 len(seq_id_seq) == 0 ): continue
            seq_id, seq = [(seq_id, seq) for seq_id, seq in alignment.items()
                           if (seq_id == gRNA_hit.target_id or is_rvs(seq_id, gRNA_hit.target_id))][0]
            seq_plus = seq_id == gRNA_hit.target_id
            gRNA_hit.set_parent_sense('+' if seq_plus == ref_plus else '-') ## used to tie break set cover
            gRNA_range = gRNA_hit.range if (not is_rvs(seq_id, gRNA_hit.target_id)) \
                         else gRNA_hit.reverse_range
            ## append gRNAHit object if within at least 1 feature
            if (within_feature(targets_feature_ranges[seq_id], seq, gRNA_range,
                               min_within_n = min_within_n, min_within_fraction = min_within_fraction)):
                # new_coverage.append(gRNA_hit)
                # # screened_gRNA[gRNA] = set(new_coverage)
                gRNA_hit.set_feature_check(True)
            else:
                gRNA_hit.set_feature_check(False)
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


def make_has_pam(pam, max_gap):
    pam = PAM(pam)
    pam_pre = pam.five_prime()
    pam_post = pam.three_prime()
    pam_pre_max = pam.five_prime_maxlen()
    pam_post_max = pam.three_prime_maxlen()
    ## output function
    def has_pam(fasta, query_result, hsp):
        seq = fasta[hsp.hit_id]
        qlen = query_result.seq_len
        ## get maximum PAM distance (unused gap allowance + uanligned gRNA length)
        gap_excess = max_gap - hsp.gap_num
        len_excess = qlen - (hsp.query_end - hsp.query_start)
        pre_excess = hsp.query_start
        post_excess = qlen - hsp.query_end
        pre_len_max = pam_pre_max + gap_excess + pre_excess
        post_len_max = pam_post_max + gap_excess + post_excess
        ## get the sequences flanking the gRNA hit
        if hsp.hit_strand == 1:
            pre_seq = seq[sstart - pre_len_max : sstart]
            post_seq = seq[send : send + post_len_max]
        else:
            post_seq = seq[sstart - post_len_max: sstart].reverse_complement()
            pre_seq = seq[send : send + pre_len_max].reverse_complement()
        ## assess whether pre-gRNA and post-gRNA regions match PAM pattern
        has_pam_pre = False if not pam_pre else bool(re.search(pam_pre, str(pre_seq)))
        has_pam_post = False if not pam_post else bool(re.search(pam_post, str(post_seq)))
        return pre_pass or post_pass
    return has_pam

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
    def within(self, r, index = 0, end_incl = False):
        return within(r, convert_range((self.start, self.end), index_in = 0, index_out = index,
                                       start_incl_in = True, start_incl_out = True,
                                       end_incl_in = False, end_incl_out = end_incl))

def mask_identical(to_mask_fname, fasta_fname, fout_fname, **kwargs):
    seqs_to_mask = fasta_to_dict(to_mask_fname)
    ## separate sequences with ambiguous bases (not compatible with BLAST) and those without
    ##  - ambiguous bases typically present in scaffold-level assemblies as runs of 'N's
    standard_bases = {'A','T','G','C','U'}
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
    from Bio import SearchIO
    masked = [Masked(hsp)
              for query_result in SearchIO.parse(fout_fname, "blast-tab", fields = ' '.join(header))
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

