import os
import re

import shutil
import tempfile
import itertools
import subprocess
import pybedtools
from datetime import datetime

from display import make_print_preindent

from functions import (
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

from exceptions import (
    MessageError
)

from parse_config import (
    IndvGenomesAll, IndvGenomesAllClear,
)



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
                        **for_merge_hits_and_filter):
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
        recip_blast(fasta_target = fout, directory = directory, genes = genes, blastn = blastn,
                    lvl = lvl, quiet = quiet, gff_beds = gff_beds, fasta_ref = fasta_ref, relax = relax,
                    keep_tmp = keep_tmp)
    ## remove tmp files
    if not keep_tmp:
        os.remove(tsv_blast_ref)
        os.remove(tsv_blast_cds)
    return

def merge_hits_and_filter(blast6_fname, fout, fasta, quiet = True, min_cds_len = 0, indv_i = 1,
                          colnames = ("sseqid", "sstart", "send", "pident", "length"),
                          blast6cds_fname = None, lvl = 0, **for_merge_hits):
    printi = make_local_print(quiet, printf = make_print_preindent(initial_lvl = lvl))
    indv_i = indv_i if not (isinstance(indv_i, IndvGenomesAll) or isinstance(indv_i, IndvGenomesAllClear)) \
             else indv_i.value
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
    # with open(merged_f, 'r') as f:
    #     dat = [x.replace('\n', '').split('\t') for x in f.readlines()]
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
    # intersect_bed = []
    # for gff_bed in gff_beds:
    #     bed = BedTool(gff_bed)
    #     ## get chromosomes in the gff_bed file
    #     chrom = set(x.split('\t')[0] for x in open(bed.groupby(o = "first", c = 1, g = 1).fn, 'r').readlines())
    #     ## filter for blast entries with chromosomes from the gff_bed file
    #     dat_filt = [x for x in data_bed if x[0] in chrom]
    #     ## skip gff_bed if no blast entries after filtering
    #     if not dat_filt: continue
    #     ## intersect
    #     dat_intersect = BedTool('\n'.join(['\t'.join(x) for x in dat_filt])).intersect(bed, wao = True)
    #     ## if gff_bed is a gff file, reorder columns
    #     if len(str(dat_intersect).split('\n')[0].split('\t')) == (6 + len(blast_metrics) + 9):
    #         str_intersect = ''
    #         tmp_gff = GFF(fmt = "GFF")
    #         for entry in dat_intersect.split('\n'):
    #             entry = entry.split('\t')
    #             if entry[6+len(blast_metrics)+4] == '.':
    #                 str_intersect += '\t'.join(entry) + '\n'
    #             else:
    #                 pass ## convert
    #         dat_intersect = BedTool(str_intersect, from_string = True)
    #     intersect_bed.append(dat_intersect)
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
                                          "feature", "phase", "attributes", "overlap"]):
    import itertools
    printi = make_local_print(quiet = quiet, printf = make_print_preindent(initial_lvl = lvl))
    genes = set(genes)
    colnames = colnames_blast + blast_metrics + colnames_bed
    get = make_custom_get(colnames)
    data = {}
    for entry in (x.split('\t') for x in
                  # tuple(itertools.chain(*[str(bedtool_obj).split('\n') for bedtool_obj in bedtool]))
                  str(bedtool).split('\n')
                  if ( x.count('\t') == len(colnames) - 1 \
                       and get(x.split('\t'), "feature") in ("gene", "pseudogene", '.') )):
        data[get(entry, "candidate")] = data.get(get(entry, "candidate"), []) + [entry]
    ## get largest bitscore for each candidate target
    max_bitscore = {candidate: max(get(data[candidate], "bitscore")) for candidate in data}
    ## identify sequences to discard and print warnings if candidate has max bitscore with target and non-target
    ## note that due to the algorithm, candidates that don't overlap with ANY features will also be kept
    throw = []
    for candidate in data:
        max_bitscore_genes = set(get(entry, "id") for entry in data[candidate]
                                 if get(entry, "bitscore") == max_bitscore[candidate])
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

#################
##  ORIGINAL   ##
##   GET_REF   ##
##  FUNCTIONS  ##
#################

def get_ref_by_genes_resolve(genes, fout = None, bed_out = None, no_bed = True, directory = None,
                             show_progress_bar = False, **kwargs):
    
    def mv_dir_overwrite(src_dir, dst_dir):
        for root, dirs, files in list(os.walk(src_dir))[::-1]:
            out_dir = os.path.join(dst_dir, os.path.relpath(root, src_dir))
            os.makedirs(out_dir, exist_ok = True)
            for fname in files:
                shutil.move(os.path.join(root, fname), os.path.join(out_dir, fname))
            for dirname in dirs:
                os.rmdir(os.path.join(root, dirname))
        return
    
    with tempfile.TemporaryDirectory() as tmpdir:
        
        store_fa = os.path.join(tmpdir, "fa_fnames.txt")
        store_bed = os.path.join(tmpdir, "bed_fnames.txt")
        
        ## get seqs
        if show_progress_bar:
            from typer import progressbar
            with progressbar(genes) as progress:
                for gene in progress:
                    get_ref_by_gene(gene, **kwargs, acc = "ref", directory = tmpdir,
                                    store_fasta = store_fa, store_bed = store_bed)
        else:
            for gene in genes:
                get_ref_by_gene(gene, **kwargs, acc = "ref", directory = tmpdir,
                                store_fasta = store_fa, store_bed = store_bed)
        
        ## merge fasta files if requested
        if fout is not None and os.path.exists(store_fa):
            cat_files(splitlines(store_fa), fout, remove = True)
        
        ## merge bed if requested
        if bed_out is not None and os.path.exists(store_bed):
            cat_files(splitlines(store_bed), bed_out, remove = True)
        
        ## remove bed files + bed dir if requested
        if no_bed and os.path.exists(store_bed):
            f_beds = splitlines(store_bed)
            if f_beds:
                for f_bed in splitlines(store_bed):
                    if os.path.exists(f_bed):
                        os.remove(f_bed)
                try:
                    os.rmdir(os.path.dirname(f_beds[0]))
                except OSError:
                    pass
        
        ## mv remaining files if directory provided
        if directory is not None:
            os.remove(store_fasta) ## delete tmp files first
            os.remove(store_bed)
            mv_dir_overwrite(tmp_dir, directory)
        
    return

############################
##  SELECTED SECTIONS OF  ##
##   GET_SEQS_FUNCTIONS   ##
############################

chrom_pref="Chr"
bed_path = None
ref_fasta = None

fields={"gff3": {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
                         "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
                         "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
                         "Is_circular": "Is_circular"}},
        "gtf2": {"all": {"ID": "transcript_id", "Parent": "gene_id"}}}

get_gffbed = make_custom_get(["chrom", "start", "end", "ID", "score", "strand", "source", "type", "phase", "attributes"])

def has_overlap(r1, r2):
    r1 = sorted(r1)
    r2 = sorted(r2)
    return (r2[0] <= r1[0] <=r2[1]) or (r1[0] <= r2[0] <= r1[1])

def has_any_overlap(l1, l2):
    any_overlap = [has_overlap(r1, r2) for r1 in l1 for r2 in l2]
    return True in any_overlap

def has_cat_overlap(l1, l2):
    return set(l1).intersection(set(l2)) != set()

def merge_ranges(*l):
    all_ranges = list(sorted(set(itertools.chain(*l))))
    if len(all_ranges) <= 1:
        return(all_ranges)
    final_ranges = [tuple(all_ranges.pop(0))]
    while all_ranges:
        if has_overlap(final_ranges[-1], all_ranges[0]):
            r1, r2 = tuple(final_ranges.pop(-1)), tuple(all_ranges.pop(0))
            final_ranges.append((min(*r1, *r2), max(*r1, *r2)))
        else:
            final_ranges.append(all_ranges.pop(0))
    return(final_ranges)

def store_fname(fout, fname):
    if fout:
        with open(fout, "a+") as f:
            f.write(fname + '\n')
    return

def make_get_attribute(attribute_fields, attribute_mod = {}, **kwargs):
    def parse_gff_attributes(s, field_sep_inter = ';', field_sep_intra = ','):
        ## if multiple separate entries for same field (e.g. 'Parent=abc.1;Parent=abc.2'), parse properly
        attributes_l = [x.split('=') for x in re.findall(f"[^{field_sep_intra}{field_sep_inter}=]+=.+?(?=[^{field_sep_intra}{field_sep_inter}=]+=.+?|$)", s)]
        attributes = {}
        for attribute in attributes_l:
            attributes[attribute[0]] = attributes.get(attribute[0], []) + \
                                       re.search(f'^.+?(?=[{field_sep_intra}{field_sep_inter}]?$)',
                                                 attribute[1]).group(0).split(field_sep_intra)
        return attributes
    def get_attribute(entry, a, feature = ''):
        ## map the field name to what's used in the file
        if not isinstance(entry, str):
            attributes = parse_gff_attributes(get_gffbed(entry, "attributes"), **kwargs)
            feature = get_gffbed(entry, "type")
        else:
            attributes = parse_gff_attributes(entry, **kwargs)
        ## get attributes
        if feature in attribute_mod and a in attribute_mod[feature]:
            mapped_field = attribute_mod[feature][a]
        elif feature in attribute_fields and a in attribute_fields[feature]:
            mapped_field = attribute_fields[feature][a]
        elif a in attribute_fields["all"]:
            mapped_field = attribute_fields["all"][a]
        else:
            mapped_field = a
        return attributes.get(mapped_field, [])
    return get_attribute

def grep_bedmerge(gene, bed, feature, out_dir, merge = False, encoding = "utf-8",
                  attribute_fields = fields["gff3"], field_sep_inter = ';', field_sep_intra = ',',
                  attribute_mod = {}, store_bed = None):
    ## get chrom, start, end of gene's features
    data = [x.split('\t') for x in splitlines(bed) if len(x) > 1]
    ## make local get_attribute function based on gff/gtf format provided
    get_attribute = make_get_attribute(attribute_fields, attribute_mod = attribute_mod)
    ## get relevant feature entries
    if feature == "gene":
        bed_raw = [x for x in data if get_gffbed(x, "type") == "gene" and gene in get_attribute(x, "ID")]
    elif feature == "mRNA":
        bed_raw = [x for x in data if get_gffbed(x, "type") == "mRNA" and (gene in get_attribute(x, "ID") or
                                                                           gene in get_attribute(x, "Parent"))]
    elif feature == "CDS":
        tmp_mrna = set(itertools.chain(*[get_attribute(x, "ID") for x in data
                                         if get_gffbed(x, "type") == "mRNA"
                                         and (gene in get_attribute(x, "ID") or
                                              gene in get_attribute(x, "Parent"))]))
        bed_raw = [x for x in data if (get_gffbed(x, "type") == "CDS" and
                                       has_cat_overlap(get_attribute(x, "Parent"), tmp_mrna))]
    else:
        bed_raw = [x for x in data if gene in get_gffbed(x, "attributes")]
    ## coerce into final columns (and merge overlapping entries if so required)
    if merge:
        from pybedtools import BedTool
        bedtoolobj = BedTool('\n'.join(['\t'.join(x) for x in bed_raw]), from_string = True)
        output = [entry.split('\t') for entry in \
                  str(bedtoolobj.merge(c = "6,8,10", o = "collapse,collapse,collapse")).split('\n') if entry]
        # tmp_f = os.path.join(out_dir, "bed", gene + '_' + feature + ".bed.tmp")
        # with open(tmp_f, "w+") as f:
        #     f.write('\n'.join(['\t'.join(x) for x in bed_raw]))
        #     output = [entry.split('\t') for entry in subprocess.check_output(("bedtools","merge","-c","6,8,10","-o","collapse,collapse,collapse"), stdin=f).decode(encoding).split('\n') if entry]
        # os.remove(tmp_f)
    else:
        output = [[x[i] for i in (0,1,2,5,7,9)] for x in bed_raw]
    ## write bed file of regions used (0-indexed)
    fout = os.path.join(out_dir, "bed", gene + '_' + feature + ".bed")
    os.makedirs(os.path.dirname(fout), exist_ok=True)
    with open(fout, "w+") as f:
        f.write('\n'.join(['\t'.join(x) for x in output]) + '\n')
    store_fname(store_bed, fout)
    return {"fout": fout, "data": output}

###################
## get reference ##
###################

def get_ref_raw(chrom, start, end, fasta_out, encoding="utf-8", ref_fasta_files=ref_fasta,
                store_fasta = None, src = None, **kwargs):
    import fileinput
    from Bio import SeqIO
    ## convert fasta input to dictionary ({source_id: path_to_file}) if not already a dictionary
    if not isinstance(ref_fasta_files, dict):
        import collections
        import six
        ## if some kind of (non-generator) iterable
        if (isinstance(ref_fasta_files, collections.Iterable)
            and not isinstance(ref_fasta_files, six.string_types)):
            ref_fasta_files = {i: fname for fname in ref_fasta_files}
        ## else if string (or even Path)
        else:
            ref_fasta_files = {"Reference": ref_fasta_files}
            src = "Reference"
    ## get relevant files and create inverse mapping ({path_to_file: source_id})
    file_map = {fname: src_id for src_id, fname in ref_fasta_files.items() if src is None or src_id == src}
    ## stream all files
    output = {}
    with fileinput.input(list(file_map.keys())) as f:
        for record in SeqIO.parse(f, "fasta"):
            if record.id == chrom:
                ## extract sequence. Prepend seqid with source_id
                output[f"{file_map[f.filename()]}|{chrom}:{start+1}..{end}"] = record.seq[start:end]
                ## exit loop and stream immediately after finding first match?
                break
    dict_to_fasta(output, fasta_out)
    store_fname(store_fasta, fasta_out)
    return

## note: setting by_gene to True will collapse identical entries from all isoforms
def get_ref_by_gene(gene, feature, out_dir, bed=bed_path, encoding="utf-8",
                    ref_fasta_files=ref_fasta, complete=False, domain="", domain_f="",
                    start_inc=True, end_inc=True, merge=False, translate=False, adj_dir=False,
                    attribute_fields = fields["gff3"], by_gene=False, attribute_mod={},
                    store_fasta = None, store_bed = None, quiet = True, lvl = 0,
                    seqid_template = "$source_id|$gene|$feature|$isoform|$domain|$n|$complete|$revcomp",
                    **kwargs):
    
    ## display function
    printi = (make_print_preindent(initial_lvl = lvl) if not quiet else lambda *x, **kwargs: None)
    
    ## get relevant gffbed entries
    if domain_f:
        feature = "CDS"
    data = grep_bedmerge(gene, bed, feature, out_dir, encoding = encoding, merge = merge, attribute_fields = attribute_fields, attribute_mod = attribute_mod, store_bed = store_bed)["data"]
    if not data:
        return
    ## get chrom, strand, and max boundaries
    chrom = data[0][0]
    start = min([int(x[1]) for x in data])
    end = max([int(x[2]) for x in data])
    strand = data[0][3]
    ## extract sequences from fasta file
    get_attribute = make_get_attribute(attribute_fields)
    if merge or feature == "gene":
        isoforms = {gene: data}
    else:
        isoforms = {isoform: [x for x in data if isoform in get_attribute(x[-1], "Parent")]
                    for isoform in set(itertools.chain(*[get_attribute(x[-1], "Parent") for x in data]))}
    isoforms = {k: sorted(v, key = lambda x: get_gffbed(x, "start"))
                for k, v in isoforms.items()}
    fasta_out_l = []
    seq_ranges = {}
    ## get fasta file of sequence data
    fasta_raw = os.path.join(out_dir, (gene + "_ref_" + feature + ("_complete" if complete else '') + \
                                       (('_' + ("domain" if not domain else domain)) if domain_f else '') + \
                                       ("_protein" if (translate and (feature=="CDS")) else '') + \
                                       ".fasta").replace('|', '_'))
    get_ref_raw(chrom, start, end, fasta_raw, encoding = encoding, ref_fasta_files = ref_fasta_files)
    ref_seq_original = list(fasta_to_dict(fasta_raw).values())[0]
    source_id = list(fasta_to_dict(fasta_raw).keys())[0].split('|')[0]
    ## create function to make seqid
    from string import Template
    def mk_seqid(**kwargs):
        return Template(seqid_template).substitute(source = source_id, gene = gene, feature = feature,
                                                   complete = ("complete" if complete else "stitched"),
                                                   domain = ("domain" if not domain else domain),
                                                   revcomp = ("revcomp" if adj_dir and strand == '-'
                                                              else "sense"), **kwargs)
    ## iterate through isoforms
    for isoform, isoform_dat in isoforms.items():
        ref_seq = ref_seq_original
        ## get fasta file of sequence data
        fasta_out = os.path.join(out_dir, (isoform + "_ref_" + feature + ("_complete" if complete else '') + \
                                           (('_' + ("domain" if not domain else domain)) if domain_f else '') + \
                                           ("_protein" if (translate and (feature=="CDS")) else '') + \
                                           ".fasta").replace('|', '_'))
        ## if extracting only domain-specific range
        if domain_f:
            domain_data = get_domain_in_genome_coords(gene, domain, domain_f, out_dir,
                                                      bed=bed, encoding=encoding,
                                                      isoform = isoform, attribute_fields = attribute_fields,
                                                      attribute_mod = attribute_mod,
                                                      start_inc = start_inc, end_inc = end_inc,
                                                      store_bed = store_bed,
                                                      **{k: v for k, v
                                                         in kwargs.items()
                                                         if k in ["qname_dname", "qstart_qend"]})
            if (not domain_data):
                continue
        else:
            domain_data = [(start, end)]
        seqs_to_write = {}
        for i, domain_range in enumerate(domain_data):
            d_start, d_end = domain_range
            ranges = [(d_start, d_end)]
            ## trim sequence if complete flag not raised or if domain required
            if (not complete) or domain_f:
                if complete and domain_f and d_start and d_end:
                    # ranges = [(max(start, d_start) - start, min(end, d_end) - start)]
                    ranges = [(max(start, d_start), min(end, d_end))]
                elif domain_f and d_start and d_end:
                    # ranges = [(max(int(x[1]), d_start) - start, min(int(x[2]), d_end) - start) \
                    #           for x in isoform_dat if has_overlap((int(x[1]), int(x[2])), (d_start, d_end))]
                    ranges = [(max(int(x[1]), d_start), min(int(x[2]), d_end)) \
                              for x in isoform_dat if has_overlap((int(x[1]), int(x[2])), (d_start, d_end))]
                else:
                    # ranges = [(int(x[1])-start, int(x[2])-start) for x in isoform_dat]
                    ranges = [(int(x[1]), int(x[2])) for x in isoform_dat]
                ref_seq = extract_ranges(ref_seq_original, [(x[0] - start, x[1] - start) for x in ranges])
            if (adj_dir or translate) and strand == '-':
                ref_seq = ref_seq.reverse_complement()
            ## translate sequence if translate flag raised AND feature is CDS
            if translate:
                if feature == "CDS" and not complete:
                    ref_seq = ref_seq.translate(to_stop = True)
                else:
                    printi("Translation is only possible when the selected feature is 'CDS' and the flag 'complete' is not raised.")
            seq_name = mk_seqid(isoform = isoform, n = i+1)
            seqs_to_write[seq_name] = ref_seq
            ## for by_gene
            if by_gene:
                overlap_ranges = []
                overlap_seq_names = []
                ## iterate through processed ranges
                for logged_ranges, logged_seq_names in seq_ranges.items():
                    ## if new range overlaps with already processed ranges
                    if has_any_overlap(ranges, logged_ranges):
                        ## note processed ranges that overlap with new range
                        overlap_ranges.append(logged_ranges)
                        overlap_seq_names.extend(logged_seq_names)
                ## if the new range overlaps w/ any of the processed ranges
                if overlap_ranges:
                    ## replace overlapped processed ranges w/ new merged range
                    for logged_ranges in overlap_ranges:
                        del(seq_ranges[logged_ranges])
                    ranges = merge_ranges(*overlap_ranges, ranges)
                # ranges = [(x[0] - start, x[1] - start) for x in ranges]
                seq_ranges[tuple(sorted(ranges))] = seq_ranges.get(tuple(sorted(ranges)), []) + [seq_name] + \
                                                    overlap_seq_names
            else:
                # ranges = [(x[0] - start, x[1] - start) for x in ranges]
                seq_ranges[tuple(sorted(ranges))] = seq_ranges.get(tuple(sorted(ranges)), []) + [seq_name]
        if seqs_to_write:
            dict_to_fasta(seqs_to_write, fasta_out)
            fasta_out_l.append(fasta_out)
    fasta_out_final = os.path.join(out_dir, (gene + "_ref_" + feature + \
                                             ("_complete" if complete else '') + \
                                             (('_' + ("domain" if not domain else domain))
                                              if domain_f else '') + \
                                             ("_protein" if (translate and (feature=="CDS") and
                                                             not complete) else '')+\
                                             ".fasta").replace('|', '_'))
    if fasta_out_l:
        if (not by_gene) and (fasta_out_l[0] != fasta_out_final):
            cat_files(sorted(fasta_out_l), fasta_out_final)
        if by_gene:
            isoform_seqs = fasta_to_dict(fasta_out_final)
            final_seqs = {}
            i = 0
            ## get fasta file of sequence data
            fasta_out = os.path.join(out_dir, (isoform + "_ref_" + feature +
                                               ("_complete" if complete else '') + \
                                               (('_' + ("domain" if not domain else domain))
                                                if domain_f else '') + \
                                               ("_protein" if (translate and (feature=="CDS")) else '') + \
                                               ".fasta").replace('|', '_'))
            get_ref_raw(chrom, start, end, fasta_out, encoding = encoding, ref_fasta_files = ref_fasta_files)
            fasta_out_l.append(fasta_out)
            ## format seqid
            for ranges, seq_names in sorted(seq_ranges.items()):
                ranges = [(s - start, e - start) for s, e in ranges]
                ## $isoform field is replaced with ranges
                seq_name = mk_seqid(isoform = ','.join(f'{r[0]}-{r[1]}' for r in ranges),
                                    n = (i+1 if domain_f else "NA"))
                # final_seqs[seq_name] = isoform_seqs[seq_names[0]]
                seq_out = ''.join([str(ref_seq_original[r[0]:r[1]]) for r in sorted(ranges)])
                if (adj_dir or translate) and strand == '-':
                    from Bio import Seq
                    seq_out = str(Seq.Seq(seq_out).reverse_complement())
                final_seqs[seq_name] = seq_out
                i += 1
            dict_to_fasta(final_seqs, fasta_out_final)
        ## remove intermediate files
        for fasta_out in fasta_out_l:
            if fasta_out != fasta_out_final:
                if os.path.exists(fasta_out):
                    os.remove(fasta_out)
        printi("{}\tSequences were successfully written to: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), fasta_out_final))
    elif not fasta_out_l:
        f = open(fasta_out_final, "w+")
        f.write('')
        f.close()
        printi("{}\t{} is an empty file".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), fasta_out_final))
    ## remove fasta_raw file
    if fasta_raw != fasta_out_final and os.path.exists(fasta_raw):
        os.remove(fasta_raw)
    store_fname(store_fasta, fasta_out_final)

# def get_ref_by_range(chrom, start, end, out_dir, encoding="utf-8", ref_fasta_files=ref_fasta,
#                      store_fasta = None):
#     fasta_out = os.path.join(out_dir, "chr{}_p{}-{}_ref.fasta".format(chrom, start + 1, end))
#     get_ref_raw(chrom_pref + chrom, start, end, fasta_out, encoding = encoding,
#                 ref_fasta_files = ref_fasta_files, store_fasta = store_fasta)
#     store_fname(store_fasta, fasta_out)
#     printi("Sequences were successfully written to: {}".format(fasta_out))

## output is...0-indexed, I think...it seems like it follows .bed conventions...
def get_domain_in_genome_coords(gene, domain, domain_f, out_dir, pos_type="aa", isoform='',
                                bed=bed_path, encoding="utf-8", start_inc = True, end_inc = True,
                                attribute_fields = fields["gff3"], attribute_mod = {}, store_bed = None,
                                qname_dname=("name", "domain name"), qstart_qend=("start", "end")):
    
    qname, dname = qname_dname
    qstart, qend = qstart_qend
    domain_raw = [x.split('\t') for x in splitlines(domain_f)]
    domain_header = domain_raw[0]
    domain_data = [x for x in domain_raw[1:] if len(domain_raw) > 1 and len(x) == len(domain_header) and \
                   (((not domain) or domain == x[domain_header.index(dname)]) and \
                    (gene if not isoform else isoform) in x[domain_header.index(qname)].split('|'))]
    output = []
    for domain_dat in domain_data:
        ## convert to 0-index, start inclusive, stop exclusive + adjust for unit (e.g. aa or nt)
        def adj_pos(start, end):
            account_unit = lambda x: x * (3 if pos_type == "aa" else 1)
            return (account_unit(int(start) - (1 if start_inc else 0)),
                    account_unit(int(end) - (0 if end_inc else 1)))
        domain_start, domain_end = adj_pos(int(domain_dat[domain_header.index(qstart)]),
                                           int(domain_dat[domain_header.index(qend)]))
        ## get CDS to define boundaries of translated nucleotides
        cds_dat = grep_bedmerge((gene if not isoform else isoform), bed, "CDS", out_dir, encoding = encoding,
                                attribute_fields = attribute_fields, attribute_mod = attribute_mod,
                                store_bed = store_bed)["data"]
        if len(cds_dat) == 0:
            continue
        ## extract boundaries
        chrom = cds_dat[0][0]
        plus_strand = (cds_dat[0][3] == '+')
        bounds = sorted([(int(x[1]), int(x[2])) for x in cds_dat],
                        key = lambda x: x[0], reverse = (not plus_strand))
        last_end = 0
        genome_start, genome_end = None, None
        for i, v in enumerate(bounds):
            curr_cds_start = last_end
            curr_cds_end = last_end + (v[1] - v[0])
            get_g_pos = lambda x: ((v[0] + (x - curr_cds_start)) if plus_strand else\
                                   (v[1] - (x - curr_cds_start)))
            if curr_cds_start <= domain_start < curr_cds_end:
                genome_start = get_g_pos(domain_start)
            if curr_cds_start < domain_end <= curr_cds_end: ## recall that end is exclusive
                genome_end = get_g_pos(domain_end)
            last_end = curr_cds_end
        # if genome_start and genome_end:
        #     output.append((min(genome_start, genome_end), max(genome_start, genome_end)))
        output.append((min(genome_start, genome_end), max(genome_start, genome_end)))
    return output
