import os
import re

import shutil
import tempfile
import itertools
import pybedtools
from datetime import datetime

from minorg.display import make_print_preindent

from minorg.functions import (
    assign_alias,
    parse_get_data,
    get_dat,
    make_custom_get,
    write_table,
    splitlines,
    cat_files,
    fasta_to_dict,
    dict_to_fasta,
    extract_ranges,
    imap_progress,
    blast6,
    blast6multi,
    make_local_print,
    num_lines
)

from minorg.annotation import (
    GFF, Annotation, Attributes
)

from minorg.index import IndexedFasta

from minorg.exceptions import (
    MessageError
)

# from minorg.parse_config import (
#     IndvGenomesAll, IndvGenomesAllClear,
# )



##############################
##      FIND HOMOLOGUE      ##
##  (reused for background  ##
##  if -g is used w/ -t/-q  ##
##############################

# def map_query_and_log(log_file, indvs = None, queries = None, indv_genomes = None):
#     if indvs:
#         fasta_queries = [[indv, indv_genomes[indv]] for indv in indvs]
#     elif queries:
#         fasta_queries = [[i, query] for i, query in enumerate(queries)]
#     with open(log_file as "a+") as f:
#         f.write("\nquery_name\tquery_file")
#         for name, fname in fasta_queries:
#             f.write(f"\n{name}\t{fname}")
#     return fasta_queries

def find_homologue_imap(indv_fout_query_dir_kwargs):
    indv_i, fout, fasta_query, directory, for_find_homologue_indv = indv_fout_query_dir_kwargs
    find_homologue_indv(fout = fout, fasta_query = fasta_query, indv_i = indv_i,
                        directory = directory, **for_find_homologue_indv)
    return fout

def find_homologue_multiindv(fasta_queries, fout, directory, threads = 1,
                             **for_find_homologue_indv):
    '''
    Additional args: 
    fasta_complete (path), fasta_cds (path), genes (tuple), check_reciprocal (bool), quiet (bool), lvl (int), 
    relax (bool), fasta_ref (list of path), gff_beds (list of path)
    Note that fasta_queries is an iterable like such [[<indv_alias/index>, <path to indv's FASTA file>]...]
    '''
    def mkfout(indv_i):
        import os
        return os.path.join(directory, f"{indv_i}_targets.fasta")
    ## generate arguments
    args = [[i, mkfout(i), fasta_query, directory, for_find_homologue_indv]
            for i, fasta_query in fasta_queries]
    print("chkpt1:", [x[0] for x in args])
    ## process
    tmp_fouts = imap_progress(find_homologue_imap, args, threads = threads)
    ## merge files
    cat_files(tmp_fouts, fout, remove = True)
    return

def find_homologue_indv(fout, directory, fasta_complete, fasta_cds, fasta_query, genes = None,
                        check_reciprocal = False, quiet = True, lvl = 0, blastn = "blastn",
                        gff_beds = None, fasta_ref = None, relax = False, keep_tmp = False,
                        attribute_mod = {}, **for_merge_hits_and_filter):
    '''
    Finds homologue in single individual (technically, single fasta file).
    Lots of unnamed args to be passed to other functions.
    '''
    ## execute BLASTN of reference genes against query
    ## note: directory must be valid path
    tsv_blast_ref = os.path.join(directory, f"{fout}_tmp_blastn_ref.tsv")
    tsv_blast_cds = os.path.join(directory, f"{fout}_tmp_blastn_cds.tsv")
    from Bio.Blast.Applications import NcbiblastnCommandline
    blast6(blastf = NcbiblastnCommandline, cmd = blastn,
           header = "qseqid,sseqid,pident,length,sstart,send",
           fout = tsv_blast_ref, query = fasta_complete, subject = fasta_query)
    blast6(blastf = NcbiblastnCommandline, cmd = blastn,
           header = "qseqid,sseqid,pident,length,sstart,send",
           fout = tsv_blast_cds, query = fasta_cds, subject = fasta_query)
    ## check for hits
    if not parse_get_data(tsv_blast_ref)[1]:
        # raise MessageError("No blast hits during homologue search, exiting programme.")
        ## write empty file
        with open(fout, "w+") as f:
            pass
        return
    ## extract homologue
    ## hidden args: accIDs = ('.',), pattern = lambda accID: accID, min_cds_len = 1,
    #                min_len = 1, min_id = 0, merge_within_range = 100,
    #                check_id_before_merge = False
    merge_hits_and_filter(blast6_fname = tsv_blast_ref, fout = fout, fasta = fasta_query,
                          blast6cds_fname = tsv_blast_cds, quiet = quiet, **for_merge_hits_and_filter)
    ## check reciprocal
    if check_reciprocal and genes:
        recip_blast_multiref(fasta_target = fout, directory = directory, genes = genes, blastn = blastn,
                             lvl = lvl, quiet = quiet, gff_beds = gff_beds, fasta_ref = fasta_ref, relax = relax,
                             keep_tmp = keep_tmp, attribute_mod = attribute_mod)
    ## remove tmp files
    if not keep_tmp:
        os.remove(tsv_blast_ref)
        os.remove(tsv_blast_cds)
    return

def merge_hits_and_filter(blast6_fname, fout, fasta, quiet = True, min_cds_len = 0, indv_i = 1,
                          colnames = ("sseqid", "sstart", "send", "pident", "length"),
                          blast6cds_fname = None, lvl = 0, **for_merge_hits):
    printi = make_local_print(quiet, printf = make_print_preindent(initial_lvl = lvl))
    ## get domain pos
    printi("Merging hits")
    header, dat = get_dat(blast6_fname)
    if not dat:
        printi("No domain found, writing empty file.")
        f = open(fout, "w+")
        f.write()
        f.close()
    else:
        merged = merge_hits(dat, header, colnames = colnames, **for_merge_hits)
        if blast6cds_fname is not None:
            merged = filter_min_cds_len(blast6cds_fname = blast6cds_fname, merged = merged,
                                        min_cds_len = min_cds_len, colnames = colnames)
        write_table(merged, fout, header = list(colnames))
        ## get domain seqs
        printi("Writing sequences to file")
        get_merged_seqs(fout, fasta, fout, header = colnames, indv_i = indv_i)
        printi(f"Merged sequences successfully written to {fout}")
    return
    
## merge overlapping ranges if same domain
def merge_hits(data, header = [], merge_within_range = 100, min_id = 90, min_len = 100,
               check_id_before_merge = False, colnames = ("sseqid", "sstart", "send", "pident", "length")):
    get = make_custom_get(header)
    dat = [get(x, *colnames) for x in data]
    get = make_custom_get(colnames)
    dat = [x for x in dat if ((not check_id_before_merge) or get(x, "pident") >= min_id)]
    dat.sort(key = lambda x: (get(x, "sseqid"), min(get(x, "sstart", "send"))))
    output = []
    if len(dat) <= 0:
        return output
    last = {x: get(dat[0], x) for x in colnames}
    s, e = last["sstart"], last["send"]
    last["sstart"], last["send"] = min(s,e), max(s,e)
    plus_minus = [0, 0]
    for i, entry in enumerate(dat):
        curr = {x: get(entry, x) for x in colnames}
        s, e = curr["sstart"], curr["send"]
        curr["sstart"], curr["send"] = min(s,e), max(s,e)
        ## check if current entry overlaps with stored merged entries (same sseqid + overlapping range)
        if curr["sseqid"] != last["sseqid"] or curr["sstart"] > (last["send"] + merge_within_range):
            ## if no overlap, check if last merged entries meet minimum len and id requirement
            if last["pident"] >= min_id and last["length"] >= min_len:
                output.append([last[x] for x in colnames])
            ## update last merged entries
            last = curr
        ## update stored merged entries with data from current entry
        else:
            last["send"] = max(curr["send"], last["send"])
            last["pident"] = max(curr["pident"], last["pident"])
            last["length"] = last["send"] - last["sstart"] + 1
        ## if is last entry and meets minimum id and len requirement, add to output
        if i == len(dat)-1 and last["pident"] >= min_id and last["length"] >= min_len:
            output.append([last[x] for x in colnames])
    return output

def filter_min_cds_len(blast6cds_fname, merged, min_cds_len = 0,
                       colnames = ("sseqid", "sstart", "send", "pident", "length")):
    if min_cds_len <= 0:
        return merged
    else:
        header, dat_cds = get_dat(blast6cds_fname)
        get = make_custom_get(header)
        dat_cds = [get(x, *colnames) for x in dat_cds]
        get = make_custom_get(colnames)
        if len(dat_cds) == 0:
            print("No CDS detected. Returning empty data.")
            return []
        dat_cds.sort(key = lambda x: (get(x, "sseqid"), min(get(x, "sstart", "send"))))
        cds_output = sorted(merge_hits(dat_cds, colnames, merge_within_range = 0, colnames = colnames),
                            key = lambda x: (get(x, "sseqid"), min(get(x, "sstart", "send"))))
        output = []
        def overlap_size(a, b):
            return max(0, min(max(a), max(b)) - max(min(a), min(b)) + 1)
        ## get CDS-complete overlap for each merged complete range
        cds_i_start = 0
        merged.sort(key = lambda x: (get(x, "sseqid"), min(get(x, "sstart", "send"))))
        for entry in merged:
            overlap, sseqid, start, end = 0, *get(entry, "sseqid", "sstart", "send")
            start, end = min(start, end), max(start, end)
            cds_i = cds_i_start
            ## find first overlapping cds entry
            while cds_i < len(cds_output) - 1 and (get(cds_output[cds_i], "sseqid") != sseqid or \
                  (get(cds_output[cds_i], "sseqid") == sseqid and \
                   max(get(cds_output[cds_i], "send", "sstart")) < start)):
                cds_i += 1
            ## update cds_last_start so we don't keep searching earlier cds entries
            cds_i_start = cds_i_start if cds_i >= len(cds_output) - 1 else cds_i
            ## while cds entries, overlap, get total overlap size
            while cds_i < len(cds_output) and \
                  get(cds_output[cds_i], "sseqid") == sseqid and \
                  min(get(cds_output[cds_i], "sstart", "send")) <= end:
                overlap += overlap_size(get(cds_output[cds_i], "sstart", "send"), (start, end))
                cds_i += 1
            if overlap >= min_cds_len:
                output.append(entry + [overlap])
        return output

## read table in and get sequences
def get_merged_seqs(merged_f, fasta, fout, header = [], indv_i = 1):
    ## get domain ranges
    dat = [x.split('\t') for x in splitlines(merged_f)]
    get = make_custom_get(dat[0])
    dat = dat[1:]
    ## parse fasta file
    seqs = fasta_to_dict(fasta)
    ## get sequences
    output = {}
    for i, entry in enumerate(dat):
        ## note: domain_dat is 1-indexed, start and end inclusive, but Python splicing is 0-indexed
        seq = seqs[get(entry, "sseqid")][get(entry, "sstart") - 1:get(entry, "send")]
        key = '|'.join([str(x) for x in \
                        ([get(entry, "sseqid"), i + 1,
                          '-'.join([str(y) for y in get(entry, "sstart", "send")]), indv_i])])
        output[key] = seq
    dict_to_fasta(output, fout)
    return

def recip_blast_multiref(fasta_target, directory, gff_beds, fasta_ref,
                         blastn = "blastn", keep_tmp = False, attribute_mod = {}, **kwargs):
    '''
    Additional args: genes (tuple), quiet (bool), relax (bool), lvl (int)
    gff_beds and fasta_ref must be dictionaries of {<alias>: <path to file>}
    '''
    from pybedtools import BedTool
    from Bio.Blast.Applications import NcbiblastnCommandline
    blast_metrics = ["pident", "bitscore"]
    tsv_blasts = []
    intersect_beds = []
    i = 0
    for alias in gff_beds:
        if not alias in fasta_ref: continue
        tsv_blast = tempfile.mkstemp(prefix = f"tmp_recipblast_{alias}", suffix = "tsv", dir = directory)[1]
        tsv_blasts.append(tsv_blast)
        ## blast
        blast6(blastf = NcbiblastnCommandline, cmd = blastn,
               header = "sseqid,sstart,send,qseqid,qstart,qend," + ','.join(blast_metrics),
               fout = tsv_blast, query = fasta_target, subject = fasta_ref[alias], dir = directory)
        ## convert hits to BED file for bedtools intersect
        get, data = parse_get_data(tsv_blast)
        data = [x for x in data if x]
        ## rearrange coords in cols 2 and 3 so that smaller value is in col 2 and larger in col 3
        data_bed = [x[:1] + (x[1:3] if int(x[1]) < int(x[2]) else x[2:0:-1]) + x[3:] for x in data]
        str_bed = '\n'.join(['\t'.join(x) for x in data_bed])
        ## bedtools intersect
        blast_bed = BedTool(str_bed, from_string = True)
        ## this doesn't work if the -b databases are different formats (e.g. gff and bed :( )
        ##  but yknow what no i'm not doing this i'm either converting gff reference annotations to bed or
        ##  writing the extended annotations as gff. we're gonna standardise.
        ##  of course if we choose gff (which i'm super inclined to actually) we're gonna have to adjust
        ##  the colnames_bed of remove_non_max_bitscore.
        intersect_beds.append(blast_bed.intersect(BedTool(gff_beds[alias]), wao = True))
    remove_non_max_bitscore(fasta_target, intersect_beds, blast_metrics = blast_metrics,
                            **kwargs) ## kwargs: genes, lvl, quiet
    ## remove temporary files
    if keep_tmp:
        intersect_fout = os.path.join(directory, "tmp_recipblast.intersect.tsv")
        with open(intersect_fout, "w+") as f:
            f.write('\n'.join([str(bedtool_obj) for bedtool_obj in intersect_beds]))
    else:
        from os import remove
        for tsv_blast in tsv_blasts:
            remove(tsv_blast)
    return


def recip_blast(fasta_target, directory, gff_beds, fasta_ref, blastn = "blastn", keep_tmp = False, **kwargs):
    '''
    Additional args: genes (tuple), quiet (bool), relax (bool), lvl (int)
    '''
    ## blast
    tsv_blast = os.path.join(directory, "tmp_recipblast.tsv")
    from Bio.Blast.Applications import NcbiblastnCommandline
    blast_metrics = ["pident", "bitscore"]
    blast6multi(blastf = NcbiblastnCommandline, cmd = blastn,
                header = "sseqid,sstart,send,qseqid,qstart,qend," + ','.join(blast_metrics),
                fout = tsv_blast, query = fasta_target, subjects = fasta_ref, dir = directory)
    ## convert hits to BED file for bedtools intersect
    from pybedtools import BedTool
    get, data = parse_get_data(tsv_blast)
    ## rearrange coords in cols 2 and 3 so that smaller value is in col 2 and larger in col 3
    data_bed = [x[:1] + (x[1:3] if int(x[1]) < int(x[2]) else x[2:0:-1]) + x[3:] for x in data]
    str_bed = '\n'.join(['\t'.join(x) for x in data_bed])
    ## bedtools intersect
    blast_bed = BedTool(str_bed, from_string = True)
    ## this doesn't work if the -b databases are different formats (e.g. gff and bed :( )
    ##  but yknow what no i'm not doing this i'm either converting gff reference annotations to bed or
    ##  writing the extended annotations as gff. we're gonna standardise.
    ##  of course if we choose gff (which i'm super inclined to actually) we're gonna have to adjust
    ##  the colnames_bed of remove_non_max_bitscore.
    intersect_bed = blast_bed.intersect([BedTool(gff_bed) for gff_bed in gff_beds], wao = True)
    remove_non_max_bitscore(fasta_target, intersect_bed, blast_metrics = blast_metrics,
                            **kwargs) ## kwargs: genes, lvl, quiet
    ## remove temporary files
    if not keep_tmp:
        from os import remove
        remove(tsv_blast)
    else:
        intersect_fout = os.path.join(directory, "tmp_recipblast.intersect.tsv")
        with open(intersect_fout, "w+") as f:
            f.write('\n'.join([str(bedtool_obj) for bedtool_obj in intersect_bed]))
    return

## if relax == True, candidate targets with equal bitscore to a target gene and a non-target gene are retained
## elif relax == False, candidate targets are only retained if all max bitscore genes are subset of target genes
def remove_non_max_bitscore(fasta, bedtool, genes, relax = False, lvl = 0, quiet = True,
                            colnames_blast=["chrom", "start", "end", "candidate", "cstart", "cend"],
                            blast_metrics = ["bitscore"],
                            colnames_bed=["bed_chrom", "bed_start", "bed_end", "id", "score", "strand", "source",
                                          "feature", "phase", "attributes", "overlap"],
                            colnames_gff=["bed_chrom", "source", "feature", "bed_start", "bed_end", "score",
                                          "strand", "phase", "attributes", "overlap"],
                            attribute_mod = {}):
    import itertools
    from pybedtools import BedTool
    printi = make_local_print(quiet = quiet, printf = make_print_preindent(initial_lvl = lvl))
    genes = set(genes)
    cols_bed = colnames_blast + blast_metrics + colnames_bed
    cols_gff = colnames_blast + blast_metrics + colnames_gff
    ## make get function
    get_bed = make_custom_get(cols_bed)
    get_gff = make_custom_get(cols_gff)
    def get(data, *args, **kwargs):
        if not data:
            helper_get = get_bed
        else:
            entry_len = len(data) if not isinstance(data[0], (list, tuple)) else len(data[0])
            if entry_len == len(cols_bed): helper_get = get_bed
            elif entry_len == len(cols_gff): helper_get = get_gff
            else: return get_bed([], "dummy", suppress_print = True)
        return helper_get(data, *args, **kwargs)
    data = {}
    for entry in (str(bedtool).split('\n') if isinstance(bedtool, BedTool)
                  else tuple(itertools.chain(*[str(bedtool_obj).split('\n') for bedtool_obj in bedtool]))):
        entry = entry.split('\t')
        if get(entry, "feature") in ("gene", "pseudogene", '.'):
            data[get(entry, "candidate")] = (data.get(get(entry, "candidate"), []) +
                                             [{"ann": Annotation(get(entry, *colnames_gff), None,
                                                                  attr_mod = attribute_mod),
                                               "bitscore": get(entry, "bitscore")}])
    ## get largest bitscore for each candidate target
    max_bitscore = {candidate: max([entry["bitscore"] for entry in data[candidate]]) for candidate in data}
    ## identify sequences to discard and print warnings if candidate has max bitscore with target and non-target
    ## note that due to the algorithm, candidates that don't overlap with ANY features will also be kept
    throw = []
    for candidate in data:
        # max_bitscore_genes = set(get(entry, "id") for entry in data[candidate]
        #                          if get(entry, "bitscore") == max_bitscore[candidate])
        max_bitscore_genes = set(entry["ann"].get_attr("ID", fmt = str)
                                 for entry in data[candidate]
                                 if entry["bitscore"] == max_bitscore[candidate])
        if max_bitscore_genes.issubset(genes): ## if max score genes are subset of target genes
            continue
        else:
            if max_bitscore_genes.isdisjoint(genes): ## if no target genes have max score
                throw.append(candidate)
            else: ## if overlapping but not subset
                if relax:
                    printi( (f"Warning: candidate target '{candidate}' has hit(s) with bitscore"
                             f" {max_bitscore[candidate]} that overlap(s) with target gene(s)"
                             f" {genes & max_bitscore_genes} and non-target gene(s)"
                             f" {max_bitscore_genes - genes}."
                             " This sequence will be retained as 'relax' has been set to True.") )
                else:
                    throw.append(candidate)
                    printi( (f"Warning: candidate target '{candidate}' has hit(s) with bitscore"
                             f" {max_bitscore[candidate]} that overlap(s) with target gene(s)"
                             f" {genes & max_bitscore_genes} and non-target gene(s)"
                             f" {max_bitscore_genes - genes}."
                             " This sequence will be removed from the list of candidate targets"
                             " as 'relax' has been set to False.") )
    ## read original candidate targets, filter, and write
    seqs = fasta_to_dict(fasta)
    dict_to_fasta({seq_id: seq for seq_id, seq in seqs.items() if seq_id not in throw}, fasta)
    return
