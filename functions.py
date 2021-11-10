import itertools

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
##  GFF_MANIP  ##
##  (CLASSES)  ##
#################
    
class GFF:
    
    def __init__(self, fname = None, data = [], attr_mod = {}, fmt = "GFF3", quiet = False, **kwargs):
        self._fmt = fmt
        self._fname = fname
        self._data = data
        self._attr_mod = attr_mod
        self._attr_fields = {"all": {"ID": "ID", "Name": "Name", "Alias": "Alias", "Parent": "Parent",
                                     "Target": "Target", "Gap": "Gap", "Derives_from": "Derives_from",
                                     "Note": "Note", "Dbxref": "Dbxref", "Ontology_term": "Ontology_term",
                                     "Is_circular": "Is_circular"}}
        self._attr_fields_inv = {}
        self._quiet = quiet
        self._kwargs = kwargs
        self._molecules = {}
        self.update_attr_fields()
        self.parse()
    
    def __iter__(self):
        for entry in self._data:
            yield entry
    
    def __len__(self):
        return len(self._data)
    
    def parse(self):
        if self._fname is None: return
        ## create add_entry function
        if self._fmt.upper() == "BED":
            def add_entry(entry):
                gff_fmt = [entry[0], entry[6], entry[7], str(int(entry[1]) + 1), entry[2],
                           entry[4], entry[5], entry[8], entry[9]]
                self.add_entry(Annotation(gff_fmt, self, **self._kwargs))
        else: ## GFF3
            def add_entry(entry):
                self.add_entry(Annotation(entry, self, **self._kwargs))
        ## start parsing
        self._data = []
        with open(self._fname, 'r') as f:
            for entry in f:
                if entry[:1] == '#' or entry == '\n': continue
                entry = entry.replace('\n', '').split('\t')
                add_entry(entry)
        self.index_molecules()
    
    def read_file(self, fname):
        self._fname = fname
        self.parse()
    
    def index_molecules(self, reset = True):
        if reset:
            self._molecules = {}
        index = 0
        curr = None
        for entry in self._data:
            if entry.molecule != curr and entry.molecule not in self._molecules:
                curr = entry.molecule
                index += 1
                self._molecules[entry.molecule] = index
        return
    
    def sort(self):
        self._data.sort(key = lambda entry: (self._molecules[entry.molecule], entry.start, -entry.end,
                                             entry.source, entry.feature, entry.attributes))
    
    def add_entry(self, gff_entry):
        self._data.append(gff_entry)
    
    def update_attr_fields(self, attr_mod = None):
        '''
        Update attribute fields w/ user-provided attribute field modification dictionary.
        (Otherwise, use self._attr_fields)
        '''
        if attr_mod is None: attr_mod = self._attr_mod
        for feature, mods in attr_mod.items():
            if feature in self._attr_fields:
                for canon, atypical in mods.items():
                    self._attr_fields[feature][canon] = atypical
            else:
                self._attr_fields[feature] = mods
        self._attr_fields_inv = self.invert_attr_fields()
        return
    
    def invert_attr_fields(self):
        output = {feature: {v: k for k, v in attr_map.items()}
                 for feature, attr_map in self._attr_fields.items()}
        return output
    
    def get_subfeatures(self, feature_ids, *features, index = False):
        '''
        Get subfeatures of features w/ user-provided feature_ids.
        '''
        features = set(features)
        if type(feature_ids) is str:
            feature_ids = {feature_ids}
        else:
            feature_ids = set(feature_ids)
        indices = [i for i, entry in enumerate(self._data) if
                   ((not features or entry.feature in features) and
                    entry.has_attr("Parent", feature_ids))]
        if index: return indices
        else: return [self._data[i] for i in indices]
    
    def get_subfeatures_full(self, feature_ids, *features, index = False):
        '''
        Get all features that are subfeatures of user-provided feature_ids
        AND subfeatures of those subfeatures, until there are no sub-sub...sub-features left.
        '''
        printi = make_local_print(quiet = self._quiet)
        features = set(features)
        output_indices = []
        curr_indices = []
        iteration_n = 0
        while True:
            iteration_n += 1
            printi(f"Executing iteration {iteration_n}")
            curr_indices = self.get_subfeatures(feature_ids, index = True)
            feature_ids = set(itertools.chain(*[self.get_i(i).get_attr("ID") for i in curr_indices]))
            output_indices.extend(curr_indices)
            if not feature_ids: break
        if features:
            output_indices = [i for i in output_indices if self.get_i(i).feature in features]
        ## prepare output
        output_indices = sorted(output_indices)
        if index: return output_indices
        else: return self.get_i(sorted(output_indices))
    
    def get_features_and_subfeatures(self, feature_ids, index = False, full = True):
        '''
        Gets features w/ feature_ids AND subfeatures of those features.
        If full = True, executes get_subfeatures_full for subfeature discovery, else get_subfeatures
        '''
        features = self.get_id(feature_ids, index = True, output_list = True)
        if full: subfeatures = self.get_subfeatures_full(feature_ids, index = True)
        else: subfeatures = self.get_subfeatures(feature_ids, index = True)
        ## prepare output
        final_features = sorted(features + subfeatures)
        if index: return final_features
        else: return self.get_i(final_features)
    
    def get_i(self, indices):
        '''
        Get GFF entry by index
        '''
        if type(indices) is int: return self._data[indices]
        else: return [self._data[i] for i in indices]
    
    def get_id(self, feature_ids, index = False, output_list = False):
        ## if even if output_list is only used when len(indices) == 1 or == 0 AND type(feature_ids) is str.
        ##  Always returns list otherwise.
        '''
        Get GFF entry by feature ID
        '''
        if type(feature_ids) is str:
            indices = [i for i, entry in enumerate(self._data) if entry.has_attr("ID", [feature_ids])]
        else:
            indices = [i for i, entry in enumerate(self._data) if entry.has_attr("ID", feature_ids)]
        if type(feature_ids) and not output_list and len(indices) <= 1:
            if not indices: return None
            elif index: return indices[0]
            else: return self.get_i(indices[0])
        elif index: return indices
        else: return [self.get_i(i) for i in indices]
    
    def write(self, fout, entries, **kwargs):
        '''
        Writes entries to file
        '''
        with open(fout, "w+") as f:
            f.write('\n'.join([entry.generate_str(**kwargs) for entry in entries]) + '\n')
    
    def write_i(self, fout, indices, **kwargs):
        '''
        Executes get_i, then writes output to file
        '''
        self.write(fout, self.get_i(indices), **kwargs)
        return
    
    def write_id(self, fout, feature_ids, **kwargs):
        '''
        Executes get_id, then writes output to file
        '''
        self.write(fout, self.get_id(feature_ids), **kwargs)
        
class Annotation:
    def __init__(self, entry, gff, **kwargs):
        self._gff = gff
        self.molecule = entry[0]
        self.source = entry[1]
        self.feature = entry[2]
        self.start = int(entry[3])
        self.end = int(entry[4])
        self.score = entry[5] if entry[5] == '.' else float(entry[5])
        self.strand = entry[6]
        self.phase = entry[7] if entry[7] == '.' else int(entry[7])
        self.attributes = Attributes(entry[8], gff, self, **kwargs)
        self._fields = {"molecule": self.molecule,
                        "source": self.source,
                        "feature": self.feature,
                        "start": self.start,
                        "end": self.end,
                        "score": self.score,
                        "strand": self.strand,
                        "phase": self.phase,
                        "attributes": self.attributes}
    def generate_attr(self, original = True, fields = None):
        if original: return self.attributes._raw
        else: return self.attributes.standardise_fields()
    def generate_str(self, fmt = "GFF"):
        if fmt.upper() in {"GFF", "GFF3"}:
            output = self.generate_gff()
        elif fmt.upper() in {"BED"}:
            output = self.generate_bed()
        return '\t'.join(map(str, output))
    def generate_gff(self):
        return list(map(str, [self.molecule, self.source, self.feature, self.start, self.end,
                              self.score, self.strand, self.phase, self.generate_attr()]))
    def generate_bed(self):
        return list(map(str, [self.molecule, self.start - 1, self.end, self.attributes.get("ID", fmt = str),
                              self.score, self.strand, self.source, self.feature, self.phase,
                              self.generate_attr()]))
    def get(self, *fields):
        return [self.f_dict[field] for field in fields]
    
    ## wrappers for Attribute methods
    def get_attr(self, a, **kwargs): return self.attributes.get(a, **kwargs)
    def is_attr(self, a, val, **kwargs): return self.attributes.is_attr(a, val, **kwargs)
    def has_attr(self, a, vals, **kwargs): return self.attributes.has_attr(a, vals, **kwargs)

class Attributes:
    def __init__(self, val, gff, entry, field_sep_inter = ';', field_sep_intra = ','):
        self._entry = entry
        self._gff = gff
        self._raw = val
        self._sep_inter = field_sep_inter
        self._sep_intra = field_sep_intra
        self._data = self._parse()
    def __repr__(self):
        return self._raw
    def __str__(self):
        return self._raw
    def _parse(self):
        import re
        ## if multiple separate entries for same field (e.g. 'Parent=abc.1;Parent=abc.2'), parse properly
        attributes_l = [x.split('=') for x in re.findall(f"[^{self._sep_intra}{self._sep_inter}=]+=.+?(?=[^{self._sep_intra}{self._sep_inter}=]+=.+?|$)", self._raw)]
        attributes = {}
        for attribute in attributes_l:
            attributes[attribute[0]] = attributes.get(attribute[0], []) + \
                                       re.search(f'^.+?(?=[{self._sep_intra}{self._sep_inter}]?$)',
                                                 attribute[1]).group(0).split(self._sep_intra)
        return attributes
    def get(self, a, fmt = list):
        feature = self._entry.feature
        attr_fields = self._gff._attr_fields
        if feature in attr_fields and a in attr_fields[feature]:
            mapped_field = attr_fields[feature][a]
        elif a in attr_fields["all"]:
            mapped_field = attr_fields["all"][a]
        else:
            mapped_field = a
        iterable = self._data.get(mapped_field, [])
        ## format return
        if fmt is str and iterable: return self._sep_intra.join(iterable)
        elif fmt is str and not iterable: return '.'
        else: return fmt(iterable)
    def is_attr(self, a, val):
        '''
        Checks if 'val' is a value of attribute 'a'
        '''
        return val in self.get(a)        
    def has_attr(self, a, vals):
        '''
        Checks if at least 1 value in 'vals' is also a value of attribute 'a'
        '''
        return (set(self.get(a)) & set(vals)) ## checks intersection != empty set
    def standardise_fields(self):
        def get_mod(field):
            feature_mod = get_recursively(self._gff._attr_fields_inv, field, self.feature, field)
            if feature_mod != field: return feature_mod
            else: return get_recursively(self._gff._attr_fields_inv, field, "all", field)
        return ';'.join([f"{get_mod(k)}={'.'.join(v)}" for k, v in self._data.items()])


#############################
##  EXTRACT_GFF_FEAUTURES  ##
#############################

def extract_features_and_subfeatures(gff_fname, feature_ids, fout, quiet = False,
                                     fin_fmt = "GFF", fout_fmt = "GFF"):
    printi = make_local_print(quiet = quiet)
    printi("Reading files")
    gff = GFF(fname = gff_fname, fmt = fin_fmt, quiet = quiet)
    # feature_ids = splitlines(feature_id_fname)
    printi("Retrieving features")
    output_features = gff.get_features_and_subfeatures(feature_ids, index = True, full = True)
    printi("Writing features")
    gff.write_i(fout, output_features, fmt = fout_fmt)
    return

#################
##  SEQ_MANIP  ##
#################

def gc_content(seq):
    seq = str(seq).upper().replace('-', '')
    gc_count = seq.count('G') + seq.count('C')
    return gc_count/len(seq)

#################
##  ANN MANIP  ##
#################

def reduce_ann(gff_beds, ids, fout, mk_tmpf_name = None):
    if mk_tmpf_name is None:
        import tempfile
        mk_tmpf_name = lambda x: tempfile.mkstemp()[1]
    bed_reds = []
    for i, gff_bed in enumerate(gff_beds):
        bed_red = mk_tmpf_name(i)
        bed_reds.append(bed_red)
        fmt = "BED" if gff_bed.split('.')[-1].upper() == "BED" else "GFF"
        extract_features_and_subfeatures(gff_bed, ids, bed_red, quiet = True,
                                         fin_fmt = fmt, fout_fmt = "BED")
    cat_files(bed_reds, fout, remove = True)
    return


####################
##  gRNA CLASSES  ##
####################

class CheckObj:
    def __init__(self, *check_names):
        self._checks = {check_name: None for check_name in check_names}
    def clear_checks(self):
        for check_name in self._checks.keys():
            self._checks[check_name] = None
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
                if check == "feature": ## this is the only one that's affected by hit location
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
    def clear_checks(self):
        for gRNA_seq in self.flatten_gRNAseqs():
            gRNA_seq.clear_checks()
        for hit in self.flatten_hits():
            hit.clear_checks()
        return
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
                      write_checks = False, checks = ["background", "GC", "feature"],
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
        super().__init__("background", "GC", "feature", "exclude", "flank")
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
    def set_feature_check(self, status): self.set_check("feature", status)
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
