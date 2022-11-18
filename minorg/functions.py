import os
import itertools
from minorg.index import IndexedFile, IndexedFasta
from minorg.fasta import (
    fasta_to_dict,
    dict_to_fasta,
    extract_ranges,
    find_identical_in_fasta,
    collapse_identical_seqs
)

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

def non_string_iter(val):
    import collections
    import six
    ## if some kind of (non-generator) iterable
    return (not isinstance(val, dict)
            and isinstance(val, collections.Iterable)
            and not isinstance(val, six.string_types))

def fill_template(template, **kwargs):
    return Template

##############
##  SPLITS  ##
##############

def split_none(iterable):
    return None if ( iterable is None or len(iterable) == 0 ) else \
        (iterable if isinstance(iterable, str) else ','.join(iterable)).split(',')

## previously: split_callback_str
def split_str(val):
    return ','.join(split_none(val))

## previously: split_callback_list
def split_list(val):
    return split_none(val)

###############
##  DISPLAY  ##
###############

from minorg.display import (print_indent, make_print_preindent,
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
from multiprocess import Pool
def imap_progress(f, args, threads = 1,
                  overwrite = True, overwrite_last = True, return_output = True,
                  msg = lambda curr, last: f"{curr}/{last} done.", lvl = 0, quiet = False):
    printi = make_local_print(quiet = False, printf = make_print_preindent(lvl + 1))
    total = len(args)
    output = []
    def print_update(curr):
        if overwrite_last or (overwrite and curr < total):
            printi(msg(curr, total), overwrite = True)
        else:
            printi(msg(curr, total), overwrite = False)
    print_update(0)
    ## don't use multiprocessing Pool if in docker environment. It simply does not work.
    ## environment variable MINORG_IN_DOCKER is defined as 'Yes' in Dockerfile
    if os.environ.get("MINORG_IN_DOCKER", False) == "Yes":
        for i, arg in enumerate(args):
            output.append(f(arg))
            print_update(i)
    else:
        pool = Pool(threads)
        for i, result in enumerate(pool.imap_unordered(f, args), 1):
            output.append(result)
            print_update(i)
        pool.close()
    return output if return_output else None

##################
##  DATA_MANIP  ##
##################

def is_empty_file(fname):
    return os.stat(fname).st_size == 0

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
    with open(fname, 'w') as f:
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

def empty_file(fname):
    with open(fname, 'w') as f:
        f.write()
    return

def cat_files(fnames, fout, remove = False):
    with open(fout, 'w') as f, fileinput.input(fnames) as fin:
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
    with open(fname, 'w') as f: ## prepend header
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
##  RANGE ARITH  ##
###################

def is_range(r):
    return len(r) == 2 and {type(r[0]), type(r[1])}.issubset({int, float})

def convert_range(ranges, index_in = 0, index_out = 1, start_incl_in = True, start_incl_out = True,
                  end_incl_in = False, end_incl_out = True):
    """
    Convert range(s) from one indexing system to another.
    
    Default values converts from 0-index [start, end) ranges to 1-index [start, end] ranges.
    
    Arguments:
        ranges (list or tuple): range(s) to convert.
            Ranges may nested and discontiuous (e.g. [(1, 3), (5, 8)]) or single level (e.g. (1, 3)).
        index_in (int): index of input range(s)
        index_out (int): index of output range(s)
        start_incl_in (bool): whether input range(s) is/are start-inclusive
        start_incl_out (bool): whether output range(s) is/are start-inclusive
        end_incl_in (bool): whether input range(s) is/are end_inclusive
        end_incl_out (bool): whether output range(s) is/are end_inclusive
    
    Returns
    -------
    tuple or list
        In same structure as input.
    """
    ## compensate index
    index_offset = index_out - index_in
    start_offset = int(start_incl_out) - int(start_incl_in)
    end_offset = int(end_incl_in) - int(end_incl_out)
    if is_range(ranges):
        start, end = ranges
        return (start + index_offset + start_offset, end + index_offset + end_offset)
    else:
        return [convert_range(r, index_in = index_in, index_out = index_out,
                              start_incl_in = start_incl_in, start_incl_out = start_incl_out,
                              end_incl_in = end_incl_in, end_incl_out = end_incl_out) for r in ranges]

## returns a in b
def within(a, b) -> bool:
    """
    Return whether range ``a`` is within or equal to range ``b``.
    
    Arguments:
        a (tuple or list): single-level range (e.g. (4, 6))
        b (tuple or list): single-level range (e.g. (2, 7))
    
    Returns
    -------
    bool
        Whether range ``a`` is within or equal to range ``b``.
    """
    a = [int(x) for x in a]
    b = [int(x) for x in b]
    return min(a) >= min(b) and max(a) <= max(b)

## returns a in any range in ranges
def within_any(a, ranges) -> bool:
    """
    Return whether range ``a`` is within at least one range in ``ranges``.
    
    Arguments:
        a (tuple or list): single-level range (e.g (4, 6))
        ranges (tuple or list): list of ranges (e.g. [(-1, 1), (2, 7), (5, 12)])
    
    Returns
    -------
    bool
        Whether range ``a`` is within any range in ``ranges``
    """
    for r in ranges:
        if within(a, r):
            return True
    return False

def ranges_subtract(r1, r2):
    """
    Subtract second range(s) from first range(s). Ranges should only be of integer values.
    
    Ranges are converted to sets of integer values so set operations can be applied.
    Resultant set is then converted back to a list of range(s).
    
    Arguments:
        r1 (list): range(s) to be subtracted from.
            Range(s) may nested and discontiuous (e.g. [(1, 3), (5, 8)]) or single level (e.g. (1, 3)).
        r2 (list): range(s) to subtract from ``r1``.
            Range(s) may nested and discontiuous (e.g. [(1, 3), (5, 8)]) or single level (e.g. (1, 3)).
    
    Returns
    -------
    list
        Of range(s)
    """
    pos_subtract = ranges_to_pos(r1) - ranges_to_pos(r2)
    return pos_to_ranges(pos_subtract)

## e.g. [[(1, 3), (6, 9)], [(2, 3), (6, 10)]] --> {1, 2, 6, 7, 8, 9} --> [(1, 3), (6, 10)]
def ranges_union(ranges):
    pos_set = ranges_to_pos(list(itertools.chain(*ranges)))
    return pos_to_ranges(pos_set)

## e.g. [(1, 3), (4, 20)], [(2, 10)] --> {2, 4, 5, 6, 7, 8, 9} --> [(2, 3), (4, 10)]
def ranges_intersect(r1, r2):
    pos_set_1 = set(ranges_to_pos(r1))
    pos_set_2 = set(ranges_to_pos(r2))
    pos_intersection = pos_set_1.intersection(pos_set_2)
    return pos_to_ranges(pos_intersection)

## convert ranges of [start, end) into set of positions
## e.g. [(1, 3), (6, 10)] --> {1, 2, 6, 7, 8, 9}
def ranges_to_pos(r):
    """
    Convert integer start-inclusive & end-exclusive range(s) to set of integer values in the range.
    
    Arguments:
        r (tuple or list): range(s) to convert.
            Ranges may nested and discontiuous (e.g. [(1, 3), (5, 8)]) or single level (e.g. (1, 3)).
    
    Returns
    -------
    set
        Of integer values within range
    """
    if is_range(r):
        return set(range(r[0], r[1]))
    else:
        return set(itertools.chain(*[ranges_to_pos(x) for x in r]))

## convert set of positions into list of ranges of [start, end)
## e.g. {1, 2, 6, 7, 8, 9} --> [(1, 3), (6, 10)]
def pos_to_ranges(pos) -> list:
    """
    Convert an iterable of integer positions to range(s).
    E.g. pos_to_range({1, 2, 6, 7, 8, 9}) --> [(1, 3), (6, 10)]
    
    Arguments:
        pos (iterable): iterable of integer positions
    
    Returns
    -------
    list
        Of 0-index, start-inclusive, end-exclusive ranges
    """
    if not pos: return []
    pos = sorted(set(pos))
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
    Return 0-indexed position adjusted for gaps in alignment
    
    Arguments:
        seq (str or Biopython Seq.Seq): gapped sequence
        pos (int): ungapped position
        gap_char (str): gap character
    
    Returns
    -------
    int
        ``pos`` adjusted for gaps in sequence
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

#################
##  SEQ_MANIP  ##
#################

def gc_content(seq):
    seq = str(seq).upper().replace('-', '')
    gc_count = seq.count('G') + seq.count('C')
    return gc_count/len(seq)
