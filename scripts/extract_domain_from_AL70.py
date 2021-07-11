import os
import re
import sys
sys.path.append("/mnt/chaelab/rachelle/src")
from basics import *

annaLenaContigs = "/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta"

def main(blast6_fname, accIDs, fout, fasta = annaLenaContigs,
         colnames = ("sseqid", "sstart", "send", "pident", "length"), **for_merge_hits):
    ## get domain pos
    print("Merging hits")
    header, dat = get_dat(blast6_fname)
    if not dat:
        print("No domain found, writing empty file.")
        f = open(fout, "w+")
        f.write()
        f.close()
    else:
        merged = merge_hits(dat, accIDs, header, colnames = colnames, **for_merge_hits)
        write_table(merged, fout, header = list(colnames))
        ## get domain seqs
        print("Writing sequences to file")
        get_merged_seqs(fout, fasta, fout, header = colnames)
        print(f"Merged sequences successfully written to {fout}")
    return
    
def get_dat(fname):
    with open(fname, 'r') as f:
        data = [x[:-1].split('\t') for x in f.readlines()]
    header = data[0]
    data = data[1:]
    return (header, data)

def make_custom_get(header, parse_num = True):
    def get_col(colname, data):
        return [get_col_in_row(x, colname) for x in data]
    def get_col_in_row(row, colname):
        output = row[header.index(colname)]
        if isinstance(output, (list, tuple)):
            if [str(x).isdigit() for x in output].count(True) == len(output):
                output = [int(x) for x in output]
            elif [str(x).replace('.','',1).replace('-','',1).isdigit() \
                  for x in output].count(True) == len(output):
                output = [float(x) for x in output]
            return output
        else:
            return output if not str(output).replace('.','',1).replace('-','',1).isdigit() else \
                float(output) if not str(output).isdigit() else int(output)
    def helper(data, *colnames):
        if isinstance(data[0], (list, tuple)):
            output = [get_col(colname, data) for colname in colnames]
            return output[0] if len(output) == 1 else [[output[r][c] for r in range(len(output))] \
                                                       for c in range(len(output[0]))]
        else:
            output = [get_col_in_row(data, colname) for colname in colnames]
            return output[0] if (len(output) == 1) else output
    return helper

    
## merge overlapping ranges if same domain
def merge_hits(data, accIDs, header = [], merge_within_range = 100, min_id = 90, min_len = 100,
               pattern = lambda accID: f"{accID}.tig", check_id_before_merge = False,
               colnames = ("sseqid", "sstart", "send", "pident", "length")):
    get = make_custom_get(header)
    dat = [get(x, *colnames) for x in data]
    get = make_custom_get(colnames)
    dat = [x for x in dat
           if (((not check_id_before_merge) or get(x, "pident") >= min_id)
               and True in [bool(re.search(pattern(accID), get(x, "sseqid"))) for accID in accIDs])]
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
        ## check if current entry overlaps with stored merged entries (same contig + overlapping range)
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

def write_table(data, fout, header = [], sep = '\t'):
    if header:
        data = [header] + data
    to_write = '\n'.join(['\t'.join([str(y) for y in x]) for x in data])
    f = open(fout, "w+")
    f.write(to_write)
    f.close()

## read table in and get sequences
def get_merged_seqs(merged_f, fasta, fout, header = []):
    ## get domain ranges
    with open(merged_f, 'r') as f:
        dat = [x[:-1].split('\t') for x in f.readlines()]
    get = make_custom_get(dat[0])
    dat = dat[1:]
    ## parse fasta file
    from fasta_manip import fasta_to_dict
    seqs = fasta_to_dict(fasta)
    ## get sequences
    output = {}
    for i, entry in enumerate(dat):
        ## note: domain_dat is 1-indexed, start and end inclusive, but Python splicing is 0-indexed
        seq = seqs[get(entry, "sseqid")][get(entry, "sstart") - 1:get(entry, "send")]
        key = '|'.join([str(x) for x in \
                        ([get(entry, "sseqid"), i + 1,
                          '-'.join([str(y) for y in get(entry, "sstart", "send")])])])
        output[key] = seq
    from fasta_manip import dict_to_fasta
    dict_to_fasta(output, fout)
