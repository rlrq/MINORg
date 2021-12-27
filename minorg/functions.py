import itertools
from minorg.index import IndexedFile, IndexedFasta

def set_overlap(a, b):
    for x in a:
        if x in b:
            return True
    return False

def assign_alias(val, mk_name = lambda i: i, start_index = 1):
    if not val:
        return {}
    elif isinstance(val, str):
        return {mk_name(start_index): val}
    elif isinstance(val, dict):
        return val
    else:
        return {mk_name(i+start_index): v for v in val}

###############
##  DISPLAY  ##
###############

from display import (print_indent, make_print_preindent,
                     print_overwrite_multi, make_print_overwrite_multi_preindent)

def make_local_print(quiet, printf = print):
    def local_print(*args, **kwargs):
        if not quiet: printf(*args, **kwargs)
    return local_print

# def make_make_fname(config): ## config is a Config object
#     def make_fname(*path):
#         fname = os.path.join(config.directory, *path)
#         config.reserve_fname(fname)
#         return fname
#     return make_fname

## display for imap_unordered
import multiprocessing as mp
def imap_progress(f, args, threads = 1, overwrite = True, overwrite_last = True, return_output = True,
                  msg = lambda curr, last: f"{curr}/{last} done.", lvl = 0, quiet = False):
    printi = make_local_print(quiet = False, printf = make_print_preindent(lvl + 1))
    pool = mp.Pool(processes = threads)
    total = len(args)
    output = []
    for i, result in enumerate(pool.imap_unordered(f, args), 1):
        output.append(result)
        if overwrite_last or (overwrite and i < len(fnames)):
            printi(msg(i, total), overwrite = True)
        else:
            printi(msg(i, total), overwrite = False)
    pool.close()
    return output if return_output else None

########################
##  BLAST FORMATTING  ##
########################

def blast6multi(blastf, header, fout, subjects,
                threshold = None, n = None, metric = None, dir = None, **kwargs):
    ## execute blast6 for first subject only
    blast6(blastf, header, fout, subject = subjects[0], **kwargs)
    ## execute blast for subsequent subjects and append to output file of first blast(6)
    for subject in subjects[1:]:
        blast(blastf, header, fout, append = True, dir = dir, subject = subject, **kwargs)
    ## execute filters, with in-place modification of fout
    if threshold and metric:
        filter_blast6_threshold(fout, fout, threshold = threshold, metric = metric)
    if n and metric:
        filter_blast6_topn(fout, fout, n = n, metric = metric)
    return

def blast(blastf, header, fout, append = False, dir = None, **kwargs):
    import os
    ## parse header
    if isinstance(header, str):
        header = header.split(',')
    ## execute blast
    if append and os.path.exists(fout):
        import tempfile
        tmpf = tempfile.mkstemp(dir = dir)[1]
        blast_cline = blastf(out = tmpf, outfmt = f"'6 {' '.join(header)}'", **kwargs)
        blast_cline()
        with open(fout, 'a+') as f:
            ## add newline if file to append to doesn't already end with newline
            if not file_ends_with_newline(fout):
                f.write('\n')
            ## append tmpf
            with open(tmpf, 'r') as f2:
                for line in f2:
                    f.write(line)
        os.remove(tmpf)
    else:
        blast_cline = blastf(out = fout, outfmt = f"'6 {' '.join(header)}'", **kwargs)
        blast_cline()
    return

## executes blast w/ outfmt 6 and prepends colnames to file
## header should either be an ordered iterable of fields or str where fields are comma-separated
def blast6(blastf, header, fout, **kwargs):
    blast(blastf, header, fout, **kwargs)
    ## reformat output file
    prepend_line(fout, ('\t'.join(header) if not isinstance(header, str) else header.replace(',', '\t')))
    return

## retains only top n hits by specified metric for each query sequence
## - hits with same value as one of the top n will also be retained
## reads file twice
def blast6topn(blastf, header, fout, n = 1, metric = "bitscore", **kwargs):
    blast6(blastf, header, fout, **kwargs)
    filter_blast6_topn(fout, fout, n = n, metric = metric)
    return

## reads file just once
def blast6threshold(blastf, header, fout, threshold = 0, metric = "bitscore", **kwargs):
    blast6(blastf, header, fout, **kwargs)
    filter_blast6_threshold(fout, fout, threshold = threshold, metric = metric)
    return

## filter blast6 output by top n
def filter_blast6_topn(fname, fout, n = 1, metric = "bitscore"):
    with open(fname, 'r') as f:
        ## parse header
        header = f.readline().replace('\n', '').split('\t')
        i_metric, i_qseqid = header.index(metric), header.index("qseqid")
        ## filter
        values = {}
        for i, entry in enumerate(f):
            dat = entry.replace('\n', '').split('\t')
            qseqid, value = dat[i_qseqid], float(dat[i_metric])
            values[qseqid] = values.get(qseqid, []) + [(value, i)]
        thresholds = {qseqid: sorted(vals)[-n][0] for qseqid, vals in values.items()}
        ## note: i+1 because the first line (header) has already been read
        ##   so line 1 is indexed 0 during enumerate(f)
        lines_to_retain = set(i+1 for qseqid, threshold in thresholds.items()
                              for val, i in values[qseqid] if val >= threshold)
    with open(fname, 'r') as f:
        output = [line for i, line in enumerate(f) if i in lines_to_retain]
    with open(fout, 'w+') as f:
        ## we're not joining the entries w/ '\n' since we didn't strip them to begin with
        f.write(''.join(['\t'.join(header) + '\n'] + output))
    return

## filter blast6 output by threshold
def filter_blast6_threshold(fname, fout, threshold = 0, metric = "bitscore"):
    with open(fname, 'r') as f:
        ## parse header
        header = f.readline().replace('\n', '').split('\t')
        i_metric = header.index(metric)
        ## filter
        output = [entry for entry in f
                  if float(entry.replace('\n', '').split('\t')[i_metric]) >= threshold]
    with open(fout, 'w+') as f:
        ## we're not joining the entries w/ '\n' since we didn't strip them to begin with
        f.write(''.join(['\t'.join(header) + '\n'] + output))
    return

##################
##  DATA_MANIP  ##
##################

def make_custom_get(header, parse_num = True):
    def get_col(colname, data):
        return [get_col_in_row(x, colname) for x in data]
    def get_col_in_row(row, colname):
        if not colname in header: print(f"Column '{colname}' is not found in headers {header}")
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
    def helper(data = None, *colnames, get_cols = False, suppress_print = False, ncol = False):
        if get_cols: return header
        if ncol: return len(header)
        if len(data) == 0:
            if not suppress_print:
                print("No data found; returning empty list")
            return []
        if isinstance(data[0], (list, tuple)):
            output = [get_col(colname, data) for colname in colnames]
            return output[0] if len(output) == 1 else [[output[r][c] for r in range(len(output))] \
                                                       for c in range(len(output[0]))]
        else:
            output = [get_col_in_row(data, colname) for colname in colnames]
            return output[0] if (len(output) == 1) else output
    return helper

def splitlines(fname, ignore_empty_lines = True):
    with open(fname, 'r') as f:
        data = f.read().split('\n')
    if ignore_empty_lines:
        return [line for line in data if line]
    else:
        return data

def get_dat(fname):
    data = [x.split('\t') for x in splitlines(fname)]
    header = data[0]
    data = data[1:]
    return (header, data)

def parse_get_data(fname, delim = None, detect = True):
    if detect:
        ext = fname.split('.')[-1]
        if ext == "csv":
            delim = ','
        elif ext == "tsv":
            delim = '\t'
    elif delim is None:
        delim = '\t'
    data = [line.split(delim) for line in splitlines(fname)]
    get = make_custom_get(data[0])
    return get, data[1:]

def write_tsv(fname, dat):
    with open(fname, "w+") as f:
        f.write('\n'.join(['\t'.join([str(y) for y in x]) for x in dat]))
    return

def write_table(data, fout, header = [], sep = '\t'):
    if header:
        data = [header] + data
    write_tsv(fout, data)
    return

def get_count_dict(iterable):
    output = {}
    for e in iterable:
        output[e] = output.get(e, 0) + 1
    return output


##################
##  FILE_MANIP  ##
##################

import fileinput

def cat_files(fnames, fout, remove = False):
    with open(fout, "w+") as f, fileinput.input(fnames) as fin:
        for line in fin:
            f.write(line)
    if remove:
        import os
        for fname in fnames:
            if os.path.exists(fname) and fname != fout:
                os.remove(fname)
    return fout

def prepend(fname, string, insert_newline = False):
    with open(fname, 'r') as f:
        dat = f.read()
    with open(fname, "w+") as f: ## prepend header
        f.write(string + ('\n' if insert_newline else '') + dat)
    return

def prepend_line(fname, string):
    prepend(fname, string, insert_newline = True)
    return

def last_line(fname):
    import os
    with open(fname, 'rb') as f:
        ## catch OSError in case of one line file
        try:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        return f.readline().decode()

def file_ends_with_newline(fname):
    import re
    return bool(re.search("\n$", last_line(fname)))

def N_trailing_newlines(fname):
    import re
    return len(re.search("\n*$", last_line(fname)))

def append(fname, string, insert_newline_if_absent = True, insert_trailing_newline = True):
    import re
    if insert_newline_if_absent and file_ends_with_newline(fname):
        string = '\n' + string
    if insert_trailing_newline and not re.search("\n$", string):
        string = string + '\n'
    with open(fname, "a+") as f:
        f.write(string)
    return

def num_lines(fname):
    i = -1 ## in case of empty files
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i+1

## header should be some kind of reusable iterable (e.g. list or tuple) of field names
def tsv_entries_as_dict(fname, header, f_filter = lambda x: True, fields = None):
    if not fields: fields = set(header)
    output = []
    with open(fname, 'r') as f:
        for line in f:
            entry = line.replace('\n', '').split('\t')
            if f_filter(entry):
                output.append({col_name: entry[i] for i, col_name in enumerate(header) if col_name in fields})
    return output

###################
##  FASTA_MANIP  ##
###################

def fasta_to_dict(fname):
    """
    returns dictionary of sequences indexed by sequence name
    """
    from Bio import SeqIO
    seqs = {}
    for seq_record in SeqIO.parse(fname, "fasta"):
        seqs[seq_record.id] = seq_record.seq
    return seqs

def dict_to_SeqRecordList(d, description = '', seq_id_func = lambda x:x, iupac_letters = None,
                          seq_type = "DNA", gap_char = '-', gapped = False):
    out_l = []
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC, Gapped
    for seq_id,seq in d.items():
        out_l.append(SeqRecord(seq if isinstance(seq, Seq) else \
                               Seq(str(seq),
                                   Gapped(iupac_letters if iupac_letters else \
                                          IUPAC.ambiguous_dna if seq_type == "DNA" else \
                                          IUPAC.ambiguous_protein if seq_type == "protein" else \
                                          IUPAC.ambiguous_rna, gap_char)) if gapped else \
                               Seq(str(seq)),
                               id = seq_id_func(seq_id), description = description))
    return out_l

def dict_to_fasta(d, fout, seq_type = "detect", gap_char = '-', gapped = False):
    def detect_iupac_letters(iupac_alphabet):
        char_set = set(itertools.chain(*[set(str(seq)) for seq in d.values()]))
        iupac_set = set(str(iupac_alphabet) + str(iupac_alphabet).lower() + gap_char)
        return char_set - iupac_set == set()
    if seq_type == "detect":
        from Bio.Alphabet import IUPAC
        iupac_letters = IUPAC.unambiguous_dna if detect_iupac_letters(IUPAC.unambiguous_dna) else \
                        IUPAC.unambiguous_rna if detect_iupac_letters(IUPAC.unambiguous_rna) else \
                        IUPAC.extended_dna if detect_iupac_letters(IUPAC.extended_dna) else \
                        IUPAC.ambiguous_dna if detect_iupac_letters(IUPAC.ambiguous_dna) else \
                        IUPAC.ambiguous_rna if detect_iupac_letters(IUPAC.ambiguous_rna) else \
                        IUPAC.protein if detect_iupac_letters(IUPAC.protein) else \
                        IUPAC.extended_protein if detect_iupac_letters(IUPAC.extended_protein) else \
                        None
    else:
        iupac_letters = None
    from Bio import SeqIO
    SeqIO.write(dict_to_SeqRecordList(d, gap_char = gap_char, gapped = gapped, iupac_letters = iupac_letters),
                fout, "fasta")

def extract_ranges(seq, ranges, strand = '+'):
    ranges_sorted = sorted(ranges, key = lambda x: int(x[0]), reverse = (strand == '-'))
    output = seq[:0]
    for start, end in ranges_sorted:
        output += seq[int(start):int(end)] if strand == '+' else \
                  seq[int(end)-1:int(start)-1:-1]
    return output

def find_identical_in_fasta(query, subject):
    """
    Searches for exact matches.
    DOES NOT EXPAND AMBIGUOUS BASES. (i.e. 'N' matches ONLY 'N' character, not any base)
    """
    from Bio import SeqIO
    subject_seqs = SeqIO.parse(open(subject, 'r'), "fasta")
    if isinstance(query, dict):
        from Bio import Seq
        query_fwd = {seqid: (seq if isinstance(seq, Seq.Seq) else Seq.Seq(seq))
                     for seqid, seq in query.items()}
    else:
        query_fwd = fasta_to_dict(query)
    query_rvs = {seqid: seq.reverse_complement() for seqid, seq in query_fwd.items()}
    output = []
    for subject_seq in subject_seqs:
        for query_id in query_fwd:
            fwd_pos = subject_seq.seq.find(query_fwd[query_id])
            rvs_pos = subject_seq.seq.find(query_rvs[query_id])
            if fwd_pos >= 0:
                output.append((query_id, subject_seq.id, fwd_pos, fwd_pos + len(query_fwd[query_id])))
            if rvs_pos >= 0:
                output.append((query_id, subject_seq.id, rvs_pos, rvs_pos + len(query_rvs[query_id])))
    return output


###################
##  RANGE ARITH  ##
###################

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

def ranges_subtract(r1, r2):
    pos_subtract = ranges_to_pos(r1) - ranges_to_pos(r2)
    return pos_to_ranges(pos_subtract)

## e.g. [[(1, 3), (6, 9)], [(2, 3), (6, 10)]] --> {1, 2, 6, 7, 8, 9} --> [(1, 3), (6, 10)]
def ranges_union(ranges):
    pos_set = ranges_to_pos(itertools.chain(*ranges))
    return pos_to_ranges(pos_set)

## e.g. [(1, 3), (4, 20)], [(2, 10)] --> {2, 4, 5, 6, 7, 8, 9} --> [(2, 3), (4, 10)]
def ranges_intersect(r1, r2):
    pos_set_1 = set(ranges_to_pos(r1))
    pos_set_2 = set(ranges_to_pos(r2))
    pos_intersection = pos_set_1.intersection(pos_set_2)
    return pos_to_ranges(pos_intersection)

## convert ranges of [start, end) into set of positions
## e.g. [(1, 3), (6, 10)] --> {1, 2, 6, 7, 8, 9}
def ranges_to_pos(ranges):
    return set(itertools.chain(*[set(range(x[0], x[1])) for x in ranges]))

## convert set of positions into list of ranges of [start, end)
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

def adjusted_pos(seq, pos, gap_char = '-'):
    """
    returns position, adjusted for gaps in alignment
    0-indexed
    """
    last_pos = 0
    while True:
        curr_gaps = seq[last_pos:pos+1].count(gap_char)
        if curr_gaps == 0:
            return pos
        last_pos = pos + 1
        pos += curr_gaps

def adjusted_ranges(seq, *ranges, subtract_gaps = True, gap_char = '-'):
    ## note that non-inclusive end of range pos has to be converted to inclusive so there aren't trailing gaps
    adj_ranges = [(adjusted_pos(seq, r[0]), adjusted_pos(seq, r[1]-1) + 1) for r in ranges]
    ## return adj_ranges if gaps do not need to be accounted for
    if not subtract_gaps: return adj_ranges
    ## remove gaps from range
    non_gap_pos = set(i for i, char in enumerate(seq) if char != gap_char)
    feature_pos = ranges_to_pos(adj_ranges)
    non_gap_feature_pos = non_gap_pos.intersection(feature_pos)
    non_gap_feature_ranges = pos_to_ranges(non_gap_feature_pos)
    return non_gap_feature_ranges

# ## test adjusted_pos
# seq = "0123---45678--90"
# x = adjusted_ranges(seq, (2,6), (7, 10))
# join_ranges = lambda seq, ranges: ''.join(seq[r[0]:r[1]] for r in ranges)
# join_ranges(seq, x) == "2345789"

## takes potentially gapped sequence, sequence range in genome, feature range(s) in genome, strand
## outputs adjusted feature range(s) in the gapped sequence as list of tuples
## ranges are 0-indexed, start inclusive and end exclusive
def adjusted_feature_ranges(seq, seq_range, feature_ranges, strand = '+', sort = True, **kwargs):
    """
    returns feature ranges, adjusted for strand and gaps in alignment
    0-indexed, start inclusive, end exclusive
    """
    if strand == '+':
        seq_start = seq_range[0]
        start_adj_ranges = [tuple([x - seq_start for x in r]) for r in feature_ranges]
    else:
        seq_start = seq_range[1]
        start_adj_ranges = [tuple(sorted([abs(x - seq_start) for x in r])) for r in feature_ranges]
    adj_ranges = adjusted_ranges(seq, *start_adj_ranges, **kwargs)
    return list(sorted(adj_ranges)) if sort else adj_ranges

###########
##  PAM  ##
###########

import regex as re
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

#################
##  SEQ_MANIP  ##
#################

def gc_content(seq):
    seq = str(seq).upper().replace('-', '')
    gc_count = seq.count('G') + seq.count('C')
    return gc_count/len(seq)
