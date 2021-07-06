import os
import shutil
import tempfile

## for get_ref_by_gene
import re
import sys
import itertools
import subprocess
import pybedtools
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from datetime import datetime
from functions import make_custom_get, splitlines
from display import make_print_preindent

# from get_seqs_functions import get_ref_by_gene

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
        
        ## merge fasta files if required
        if fout is not None:
            seqs_out = {}
            from functions import fasta_to_dict, dict_to_fasta, splitlines
            for fname in splitlines(store_fa):
                seqs_out = {**seqs_out, **fasta_to_dict(fname)}
            dict_to_fasta(seqs_out, fout)
        
        ## merge bed if required
        if bed_out is not None:
            from functions import cat_files, splitlines
            cat_files(splitlines(store_bed), bed_out)
        
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
    # print("l:", l)
    all_ranges = list(sorted(set(itertools.chain(*l))))
    # print("all_ranges:", all_ranges)
    if len(all_ranges) <= 1:
        return(all_ranges)
    final_ranges = [tuple(all_ranges.pop(0))]
    while all_ranges:
        if has_overlap(final_ranges[-1], all_ranges[0]):
            r1, r2 = tuple(final_ranges.pop(-1)), tuple(all_ranges.pop(0))
            # final_ranges.append(tuple(min(*r1, *r2), max(*r1, *r2)))
            final_ranges.append((min(*r1, *r2), max(*r1, *r2)))
        else:
            final_ranges.append(all_ranges.pop(0))
    # print("final_ranges:", final_ranges)
    return(final_ranges)

def store_fname(fout, fname):
    if fout:
        open(fout, "a+").write(fname + '\n')
    return

get_gffbed = make_custom_get(["chrom", "start", "end", "ID", "score", "strand", "source", "type", "phase", "attributes"])

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
    data = [x[:-1].split('\t') for x in open(bed, 'r').readlines() if len(x) > 1]
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
        bed_raw = [x for x in data if get_gffbed(x, "type") == "CDS" and has_cat_overlap(get_attribute(x, "Parent"), tmp_mrna)]
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
                store_fasta = None, **kwargs):
    from functions import fasta_to_dict, dict_to_fasta
    if isinstance(ref_fasta_files, dict):
        ref_seq = list(fasta_to_dict(ref_fasta_files[chrom]).values())[0][start:end] ## 0-indexed
        dict_to_fasta({"Reference|{}:{}..{}".format(chrom, start + 1, end): ref_seq}, fasta_out)
    else:
        dict_to_fasta({"Reference|{}:{}..{}".format(chrom, start + 1, end): fasta_to_dict(ref_fasta_files)[chrom][start:end]}, fasta_out)
    store_fname(store_fasta, fasta_out)

## note: setting by_gene to True will collapse identical entries from all isoforms
def get_ref_by_gene(gene, feature, out_dir, bed=bed_path, encoding="utf-8",
                    ref_fasta_files=ref_fasta, complete=False, domain="", domain_f="",
                    start_inc=True, end_inc=True, merge=False, translate=False, adj_dir=False,
                    attribute_fields = fields["gff3"], by_gene=False, attribute_mod={},
                    store_fasta = None, store_bed = None, verbose = True, lvl = 0, **kwargs):
    
    ## display function
    printi = (make_print_preindent(initial_lvl = lvl) if verbose else lambda *x, **kwargs: None)
    
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
    fasta_out = os.path.join(out_dir, (gene + "_ref_" + feature + ("_complete" if complete else '') + \
                                       (('_' + ("domain" if not domain else domain)) if domain_f else '') + \
                                       ("_protein" if (translate and (feature=="CDS")) else '') + \
                                       ".fasta").replace('|', '_'))
    get_ref_raw(chrom, start, end, fasta_out, encoding=encoding, ref_fasta_files=ref_fasta_files)
    from function import fasta_to_dict
    ref_seq_original = list(fasta_to_dict(fasta_out).values())[0]
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
                    ranges = [(max(start, d_start) - start, min(end, d_end) - start)]
                elif domain_f and d_start and d_end:
                    ranges = [(max(int(x[1]), d_start) - start, min(int(x[2]), d_end) - start) \
                              for x in isoform_dat if has_overlap((int(x[1]), int(x[2])), (d_start, d_end))]
                else:
                    ranges = [(int(x[1])-start, int(x[2])-start) for x in isoform_dat]
                from functions import extract_ranges
                ref_seq = extract_ranges(ref_seq_original, ranges)
            if (adj_dir or translate) and strand == '-':
                ref_seq = ref_seq.reverse_complement()
            ## translate sequence if translate flag raised AND feature is CDS
            if translate:
                if feature == "CDS" and not complete:
                    ref_seq = ref_seq.translate(to_stop = True)
                else:
                    printi("Translation is only possible when the selected feature is 'CDS' and the flag 'complete' is not raised.")
            seq_name = "{}|{}|{}|{}".format("Reference", gene, feature, isoform) + \
                       (('|' + ("domain" if not domain else domain) + f"|{i+1}") if domain_f else '') + \
                       ("|complete" if complete else '') + \
                       ("|revcomp" if adj_dir and strand == '-' else '')
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
                seq_ranges[tuple(sorted(ranges))] = seq_ranges.get(tuple(sorted(ranges)), []) + [seq_name] + overlap_seq_names
            else:
                seq_ranges[tuple(sorted(ranges))] = seq_ranges.get(tuple(sorted(ranges)), []) + [seq_name]
        if seqs_to_write:
            from functions import dict_to_fasta
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
            from functions import cat_files
            cat_files(sorted(fasta_out_l), fasta_out_final)
            # for fasta_out in fasta_out_l:
            #     os.remove(fasta_out)
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
            get_ref_raw(chrom, start, end, fasta_out, encoding=encoding, ref_fasta_files=ref_fasta_files)
            for ranges, seq_names in sorted(seq_ranges.items()):
                seq_name_l = seq_names[0].split('|')
                seq_name_l[3] = ','.join(f'{r[0]}-{r[1]}' for r in ranges)
                if domain_f:
                    seq_name_l[5] = str(i + 1)
                seq_name = '|'.join(seq_name_l)
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
                os.remove(fasta_out)
        printi("{}\tSequences were successfully written to: {}".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), fasta_out_final))
    elif not fasta_out_l:
        f = open(fasta_out_final, "w+")
        f.write('')
        f.close()
        printi("{}\t{} is an empty file".format(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), fasta_out_final))
    store_fname(store_fasta, fasta_out_final)

def get_ref_by_range(chrom, start, end, out_dir, encoding="utf-8", ref_fasta_files=ref_fasta,
                     store_fasta = None):
    fasta_out = os.path.join(out_dir, "chr{}_p{}-{}_ref.fasta".format(chrom, start + 1, end))
    get_ref_raw(chrom_pref + chrom, start, end, fasta_out, encoding=encoding, ref_fasta_files=ref_fasta_files)
    store_fname(store_fasta, fasta_out)
    printi("Sequences were successfully written to: {}".format(fasta_out))

## output is...0-indexed, I think...it seems like it follows .bed conventions...
def get_domain_in_genome_coords(gene, domain, domain_f, out_dir, pos_type="aa", isoform='',
                                bed=bed_path, encoding="utf-8", start_inc = True, end_inc = True,
                                attribute_fields = fields["gff3"], attribute_mod = {}, store_bed = None,
                                qname_dname=("name", "domain name"), qstart_qend=("start", "end")):
    
    qname, dname = qname_dname
    qstart, qend = qstart_qend
    with open(domain_f, 'r') as f:
        domain_raw = [x[:-1].split('\t') for x in f.readlines()]
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
