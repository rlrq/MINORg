import os
import re

import shutil
import tempfile
import itertools
import pybedtools
from datetime import datetime
from pybedtools import BedTool

from minorg import MINORgError

from minorg.display import make_print_preindent

from minorg.functions import (
    assign_alias,
    parse_get_data,
    get_dat,
    make_custom_get,
    write_table,
    splitlines,
    is_empty_file,
    cat_files,
    imap_progress,
    make_local_print,
    num_lines
)

from minorg.fasta import (
    fasta_to_dict,
    dict_to_fasta,
    extract_ranges
)

from minorg.blast import (
    BlastResult,
    blast,
    blast6,
    blast6multi
)

from minorg.annotation import (
    GFF, Annotation, Attributes
)

from minorg.index import IndexedFasta

# from minorg.parse_config import (
#     IndvGenomesAll, IndvGenomesAllClear,
# )

## TODO: test merge_hits_and_filter


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
    """
    Wrapper of find_homologue_indv function for pickling.
    """
    indv_i, fout, fasta_query, directory, for_find_homologue_indv = indv_fout_query_dir_kwargs
    find_homologue_indv(fout = fout, fasta_query = fasta_query, indv_i = indv_i,
                        directory = directory, **for_find_homologue_indv)
    return fout

# def find_homologue_multiindv(fasta_queries, fout, directory, threads = 1,
#                              **for_find_homologue_indv):
#     """
#     Parallel discovery of homologues in multiple query FASTA.
    
#     Arguments:
#         fasta_queries (list): require, list of paths(s) to FASTA file(s) in which to find homologues
#         fout (str): required, path to output FASTA file in which to write homologues
#         directory (str): required, path in which to write temporary files
#         threads (int): maximum number of parallel processes
#         **for_find_homologue_indv: additional arguments for find_homologue_indv
    
#     Additional args: 
#     fasta_complete (path), fasta_cds (path), genes (tuple), check_reciprocal (bool), quiet (bool), lvl (int), 
#     relax (bool), fasta_ref (list of path), gff (list of path)
#     Note that fasta_queries is an iterable like such [[<indv_alias/index>, <path to indv's FASTA file>]...]
#     """
#     def mkfout(indv_i):
#         import os
#         return os.path.join(directory, f"{indv_i}_targets.fasta")
#     ## generate arguments
#     args = [[i, mkfout(i), fasta_query, directory, for_find_homologue_indv]
#             for i, fasta_query in fasta_queries]
#     ## process
#     tmp_fouts = imap_progress(find_homologue_imap, args, threads = threads)
#     ## merge files
#     cat_files(tmp_fouts, fout, remove = True)
#     return

# def find_homologue_indv(fout, directory, fasta_complete, fasta_cds, fasta_query,
#                         quiet = True, lvl = 0, keep_tmp = False,
#                         check_reciprocal = False, relax = False,
#                         genes = None, blastn = "blastn", bedtools = '',
#                         fasta_ref = None, gff = None,
#                         attribute_mod = {}, **for_merge_hits_and_filter):
#     """
#     Find homologue in single FASTA file.
    
#     1. Execute BLASTN of reference sequences (``fasta_complete`` and ``fasta_cds``) against ``fasta_query``.
#     2. :func:`~minorg.extract_homologue.merge_hits_and_filter`: Merge hits based on proximity to each other, filtering by length and % identity, to generate candidate homologues.
#     3. :func:`~minorg.extract_homologue.recip_blast_multiref`: If check_reciprocal=True, execute BLASTN of candidate homologues to reference genome.
#         - :func:`~minorg.extract_homologue.recip_blast_multiref`: Remove candidate homolougues if the hit with the best bitscore is NOT to a gene in ``gene``.
#         - :func:`~minorg.extract_homologue.recip_blast_multiref`: If relax=False, candidate homologues which best bitscore hit overlaps with gene in ``gene`` AND ALSO a gene NOT IN ``gene`` will be removed.
    
#     Arguments:
#         fout (str): required, path to output FASTA file in which to write homologues
#         directory (path): required, directory in which to write temporary files
#         fasta_complete (str): required, path to FASTA file containing genomic reference sequences
#         fasta_cds (str): required, path to FASTA file containing CDS reference sequences
#         fasta_query (str): required, path to FASTA file in which to search for homologues
#         quiet (bool): print only essential messages
#         lvl (int): optional, indentation level of printed messages
#         keep_tmp (bool): keep temporary files (default=False)
#         check_reciprocal (bool): execute reciprocal check
#         genes (tuple): optional, required only if check_reciprocal=True, on-target gene IDs
#         blastn (str): optional, required only if check_reciprocal=True, 
#             blastn command (e.g. 'blastn') if available at CLI else path to blastn executable
#         bedtools (str): path to directory contaiing BEDTools executables if bedtools is not
#             in command-search path
#         fasta_ref (dict): optional, required only if check_reciprocal=True,
#             dictionary of paths to reference genome FASTA files in format {'<alias>': '/path/to/FASTA'}.
#             Each reference genome should have the same alias in ``fasta_ref`` and ``gff``.
#         gff (dict): optional, required only if check_reciprocal=True,
#             dictionary of paths to reference genome GFF3 files in format {'<alias>': '/path/to/GFF3'}.
#             Each reference genome should have the same alias in ``fasta_ref`` and ``gff``.
#         attribute_mod (dict): optional, 
#             required only if non-standard attriute field names are present in GFF3 files.
#             Dictionary describing attribute modification.
#         **for_merge_hits_and_filter: additional arguments for 
#             :func:`~minorg.extract_homologue.for_merge_hits_and_filter`
#     """
#     ## execute BLASTN of reference genes against query
#     ## note: directory must be valid path
#     tsv_blast_ref = os.path.join(directory, f"{fout}_tmp_blastn_ref.tsv")
#     tsv_blast_cds = os.path.join(directory, f"{fout}_tmp_blastn_cds.tsv")
#     from Bio.Blast.Applications import NcbiblastnCommandline
#     blast(blastf = NcbiblastnCommandline, cmd = blastn,
#           header = None, ## default fields
#           # header = "qseqid,sseqid,pident,length,sstart,send",
#           fout = tsv_blast_ref, query = fasta_complete, subject = fasta_query)
#     blast(blastf = NcbiblastnCommandline, cmd = blastn,
#           header = None, ## default fields
#           # header = "qseqid,sseqid,pident,length,sstart,send",
#           fout = tsv_blast_cds, query = fasta_cds, subject = fasta_query)
#     ## check for hits
#     # if not parse_get_data(tsv_blast_ref)[1]:
#     if is_empty_file(tsv_blast_ref):
#         # raise MINORgError("No blast hits during homologue search, exiting programme.")
#         ## write empty file
#         with open(fout, "w+") as f:
#             pass
#         return
#     ## extract homologue
#     ## hidden args: accIDs = ('.',), pattern = lambda accID: accID, min_cds_len = 1,
#     #                min_len = 1, min_id = 0, merge_within_range = 100,
#     #                check_id_before_merge = False
#     merge_hits_and_filter(blast6_fname = tsv_blast_ref, fout = fout, fasta = fasta_query,
#                           blast6cds_fname = tsv_blast_cds, quiet = quiet, **for_merge_hits_and_filter)
#     ## check reciprocal
#     if check_reciprocal and genes:
#         recip_blast_multiref(fasta_target = fout, directory = directory, genes = genes, blastn = blastn,
#                              lvl = lvl, quiet = quiet, gff = gff, fasta_ref = fasta_ref, relax = relax,
#                              keep_tmp = keep_tmp, attribute_mod = attribute_mod, bedtools = bedtools)
#     ## remove tmp files
#     if not keep_tmp:
#         os.remove(tsv_blast_ref)
#         os.remove(tsv_blast_cds)
#     return

def merge_hits_and_filter(blast6_fname, fout, fasta, quiet = True, min_cds_len = 0, indv_i = 1,
                          # colnames = ("sseqid", "sstart", "send", "pident", "length"),
                          colnames = ("molecule", "start", "end", "max_pident", "length"),
                          blast6cds_fname = None, lvl = 0, **for_merge_hits):
    """
    Merge hits to generate homologues and filter by minimum inferred CDS length.
    
    1. :func:`~minorg.extract_homologue.merge_hits`: Merge hits in proximity to each other to generate candidate homologues.
    2. [optional] Filter candidate homologues by minimum CDS length (``min_cds_len``), as determined by how many bases in the candidate are covered by hits from reference CDS sequences.
    
    Arguments:
        blast6_fname (str): path to BLASTN output file (outfmt 6) prepended with header,
            where query was reference genomic sequences and subject was ``fasta``
        fout (str): path to output FASTA file in which to write homologues
        fasta (str): path to FASTA file in which to find homologues
        quiet (bool): print only essential messages
        min_cds_len (int): minimum CDS length in homologues (default=0)
        indv_i (str or int): alias for ``fasta``, used for generating sequence names
        colnames (tuple or list): fields in ``blast6_fname`` and ``blast6cds_fname``
        blast6cds_fname (str): optional, 
            if provided then candidate homologues will be filtere by minimum CDS length,
            path to BLASTN output file (outfmt 6) prepended with header,
            where query was reference CDS sequences and subject was ``fasta``
        lvl (int): optional, indentation level of printed messages
    """
    printi = make_local_print(quiet, printf = make_print_preindent(initial_lvl = lvl))
    ## get domain pos
    printi("Merging hits")
    dat = BlastResult(blast6_fname, "blast-tab")
    if len(dat) == 0:
        printi("No domain found, writing empty file.")
        f = open(fout, "w+")
        f.write()
        f.close()
    else:
        merged = merge_hits(dat, colnames = colnames, **for_merge_hits)
        if blast6cds_fname is not None:
            merged = filter_min_cds_len(blast6cds_fname = blast6cds_fname, merged = merged,
                                        min_cds_len = min_cds_len, colnames = colnames)
            colnames = list(colnames) + ["overlap"]
        ## use fout as temporary file to store merged ranges
        write_table([tuple(x[colname] for colname in colnames) for x in merged], fout, header = list(colnames))
        ## get domain seqs
        printi("Writing sequences to file")
        get_merged_seqs(fout, fasta, fout, header = colnames, indv_i = indv_i)
        printi(f"Merged sequences successfully written to {fout}")
    return
    
## merge overlapping ranges if same domain
def merge_hits(data, merge_within_range = 100, min_id = 90, min_len = 100,
               check_id_before_merge = False, colnames = ("sseqid", "sstart", "send", "pident", "length")):
    """
    Merge BLAST hits to generate candidate homologues.
    
    1. [optional] If check_id_before_merge=True, discard hits with % identity < ``min_id``
    2. Merge hits within ``merge_within_range`` bp of each other
        - Discard if no hit included in the merge has a % identity >= ``min_id``
        - Discard if resultant merged range is shorter than ``min_len`` bp
    
    Arguments:
        data (:class:`minorg.functions.BlastResult`): BLASTN result
        merge_within_range (int): maximum number of bases between hits to be merged
        min_id (float): minimum hit % identity
        min_len (int): minimum candidate homologue length
        check_id_before_merge (bool): discard hits with % identity < ``min_id`` before merging
        colnames (tuple): BLASTN fields to retain in output tsv of merged hits
    
    Returns
    -------
    list
        Of dict of merged hit data, where fields for each entry are: molecule, start, end, length
    """
    dat = [x for x in data if x.hsp.ident_pct >= min_id]
    dat.sort(key = lambda x: (x.subject.id, x.hsp.hit_start))
    output = []
    if len(dat) <= 0:
        return output
    first = dat[0]
    curr_merged = {"molecule": first.subject.id, "start": first.hsp.hit_start, "end": first.hsp.hit_end,
                   "max_pident": first.hsp.ident_pct, "length": first.hsp.hit_end - first.hsp.hit_start}
    for i, entry in enumerate(dat):
        ## check if current entry overlaps with stored merged entries (same sseqid + overlapping range)
        if (entry.subject.id != curr_merged["molecule"]
            or entry.hsp.hit_start > (curr_merged["end"] + merge_within_range)):
            ## if no overlap, check if last merged entries meet minimum len and id requirement
            if (curr_merged["max_pident"] >= min_id
                and (curr_merged["end"] - curr_merged["start"]) >= min_len):
                output.append(curr_merged)
            ## update last merged entries
            curr_merged = {"molecule": entry.subject.id, "start": entry.hsp.hit_start, "end": entry.hsp.hit_end,
                           "max_pident": entry.hsp.ident_pct, "length": entry.hsp.hit_end - entry.hsp.hit_start}
        ## update stored merged entries with data from current entry
        else:
            curr_merged["end"] = max(curr_merged["end"], entry.hsp.hit_end)
            curr_merged["max_pident"] = max(curr_merged["max_pident"], entry.hsp.ident_pct)
            curr_merged["length"] = curr_merged["end"] - curr_merged["start"]
        ## if is last entry and meets minimum id and len requirement, add to output
        if i == len(dat)-1 and curr_merged["max_pident"] >= min_id and curr_merged["length"] >= min_len:
            output.append(curr_merged)
    return output

def filter_min_cds_len(blast6cds_fname, merged, min_cds_len = 0,
                       colnames = ("sseqid", "sstart", "send", "pident", "length")):
    """
    Filter candidate homologues by minimum CDS length.
    
    Arguments:
        blast6cds_fname (str): required, path to BLASTN output file (outfmt 6) prepended with header,
            where query was reference CDS sequences and subject was ``fasta``
        merged (list): list of candidate homologue data output by :func:`~minorg.extract_homologue.merge_hits`
        min_cds_len (int): minimum CDS length in homologues (default=0)
    
    Returns
    -------
    list
        Of dict of merged hit data, where fields for each entry are: molecule, start, end, length
    """
    if min_cds_len <= 0:
        for entry in merged:
            entry["cds_overlap"] = "NA"
        return merged
    else:
        dat_cds = BlastResult(blast6cds_fname, "blast-tab")
        if len(dat_cds) == 0:
            print("No CDS detected. Returning empty data.")
            return []
        cds_output = sorted(merge_hits(dat_cds, merge_within_range = 0, colnames = colnames),
                            key = lambda x: (x["molecule"], x["start"]))
        output = []
        def overlap_size(a, b):
            return max(0, min(a[1], b[1]) - max(a[0], b[0]))
        ## get CDS-complete overlap for each merged complete range
        cds_i_start = 0
        merged.sort(key = lambda x: (x["molecule"], x["start"]))
        for entry in merged:
            overlap, sseqid, start, end = 0, entry["molecule"], entry["start"], entry["end"]
            cds_i = cds_i_start
            ## find first overlapping cds entry
            while cds_i < len(cds_output) - 1 and \
                  (cds_output[cds_i]["molecule"] != sseqid or \
                   (cds_output[cds_i]["molecule"] == sseqid and \
                    cds_output[cds_i]["end"] <= start)):
                cds_i += 1
            ## update cds_last_start so we don't keep searching earlier cds entries
            cds_i_start = cds_i_start if cds_i >= len(cds_output) - 1 else cds_i
            ## while cds entries overlap, get total overlap size
            while cds_i < len(cds_output) and \
                  cds_output[cds_i]["molecule"] == sseqid and \
                  cds_output[cds_i]["start"] < end:
                overlap += overlap_size((cds_output[cds_i]["start"], cds_output[cds_i]["end"]),
                                        (start, end))
                cds_i += 1
            if overlap >= min_cds_len:
                entry["overlap"] = overlap
                output.append(entry)
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
        seq = seqs[get(entry, "molecule")][get(entry, "start"):get(entry, "end")]
        key = '|'.join([str(x) for x in \
                        ([get(entry, "molecule"), i + 1,
                          f"{get(entry, 'start') + 1}-{get(entry, 'end')}", indv_i])])
        output[key] = seq
    dict_to_fasta(output, fout)
    return

def recip_blast_multiref(fasta_target, directory, gff, fasta_ref,
                         blastn = "blastn", bedtools = '', keep_tmp = False, attribute_mod = {}, **kwargs):
    """
    Additional args: genes (tuple), quiet (bool), relax (bool), lvl (int)
    
    gff and fasta_ref must be dictionaries of {<alias>: <path to file>}
    
    Arguments:
        fasta_target (str): path to FASTA file containing query sequences for reciprocal BLAST
        directory (str): path to directory in which to write temporary files
        gff (dict): dictionary of GFF3 files for subjects of reciprocal BLAST
        fasta_ref (dict): dictionary of FASTA files containing subject sequences for reciprocal BLAST,
        blastn (str): optional, required only if check_reciprocal=True, 
            blastn command (e.g. 'blastn') if available at CLI else path to blastn executable
        bedtools (str): optional, required only if bedtool is not in command-search path;
            path to directory contaiing BEDTools executables
        attribute_mod (dict): optional, 
            required only if non-standard attriute field names are present in GFF3 files.
            Dictionary describing attribute modification.
    """
    pybedtools.helpers.set_bedtools_path(path = bedtools)
    from Bio.Blast.Applications import NcbiblastnCommandline
    blast_metrics = ["pident", "bitscore"]
    tsv_blasts = []
    intersect_beds = []
    i = 0
    for alias in gff:
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
        intersect_beds.append(blast_bed.intersect(BedTool(gff[alias]), wao = True))
    remove_non_max_bitscore(fasta_target, intersect_beds, blast_metrics = blast_metrics, bedtools = bedtools,
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


def recip_blast(fasta_target, directory, gff, fasta_ref,
                blastn = "blastn", bedtools = '', keep_tmp = False, **kwargs):
    '''
    Additional args: genes (tuple), quiet (bool), relax (bool), lvl (int)
    '''
    pybedtools.helpers.set_bedtools_path(path = bedtools)
    ## blast
    tsv_blast = os.path.join(directory, "tmp_recipblast.tsv")
    from Bio.Blast.Applications import NcbiblastnCommandline
    blast_metrics = ["pident", "bitscore"]
    blast6multi(blastf = NcbiblastnCommandline, cmd = blastn,
                header = "sseqid,sstart,send,qseqid,qstart,qend," + ','.join(blast_metrics),
                fout = tsv_blast, query = fasta_target, subjects = fasta_ref, dir = directory)
    ## convert hits to BED file for bedtools intersect
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
    intersect_bed = blast_bed.intersect([BedTool(fname) for fname in gff], wao = True)
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
                            bedtools = '', attribute_mod = {}) -> None:
    """
    Remove query sequences for which the subject feature in the query-subject hit with the max bitscore is
    not a target gene/feature. This occurs in-place.
    
    Arguments:
        fasta (str): path to FASTA file of query sequences to reduce
        bedtool (:class:`BedTool`): BedTool object where BLAST hits have been intersected with
            subject GFF3 files
        genes (list): gene/feature IDs of targets
        relax (bool): retain query sequences even if max bitscore hit overlaps with non-target feature so long
            as it also overlaps with a target feature
        lvl (int): printing indentation
        quiet (bool): silence non-essential messages
        colnames_blast (list): column names of BLAST output
        blast_metrics (list): additional column names of metrics in BLAST output
        colnames_bed (list): column names if annotation intersected with is BED format
        colnames_gff (list): column names if annotation intersected with is GFF3 format
        bedtools (str): path to directory contaiing BEDTools executables if bedtool is not
            in command-search path
        attribute_mod (dict): optional, 
            required only if non-standard attriute field names are present in GFF3 files.
            Dictionary describing attribute modification.
    """
    import itertools
    pybedtools.helpers.set_bedtools_path(path = bedtools)
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
