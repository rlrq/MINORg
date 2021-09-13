def set_overlap(a, b):
    for x in a:
        if x in b:
            return True
    return False

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

## executes blast w/ outfmt 6 and prepends colnames to file
## header should either be an ordered iterable of fields or str where fields are comma-separated
def blast6(blastf, header, fout, **kwargs):
    ## parse header
    if isinstance(header, str):
        header = header.split(',')
    ## execute blast
    blast_cline = blastf(out = fout, outfmt = f"'6 {' '.join(header)}'", **kwargs)
    blast_cline()
    ## reformat output file
    with open(fout, 'r') as f:
        dat = f.read()
    with open(fout, "w+") as f: ## prepend header
        f.write('\t'.join(header) + '\n' + dat)
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
            if os.path.exists(fname):
                os.remove(fname)
    return fout


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
        import itertools
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
        self.update_attr_fields()
        self.parse()
    
    def parse(self):
        if self._fname is None: return
        self._data = []
        raw_data = [x.split('\t') for x in splitlines(self._fname) if x[:1] != '#']
        if self._fmt.upper() == "BED":
            def add_entry(entry):
                gff_fmt = [entry[0], entry[6], entry[7], str(int(entry[1]) + 1), entry[2],
                           entry[4], entry[5], entry[8], entry[9]]
                self.add_entry(Annotation(gff_fmt, self, **self._kwargs))
        else: ## GFF3
            def add_entry(entry):
                self.add_entry(Annotation(entry, self, **self._kwargs))
        for entry in [x for x in raw_data]:
            add_entry(entry)
    
    def read_file(self, fname):
        self._fname = fname
        self.parse()
    
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
        import itertools
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
        features = self.get_id(feature_ids, index = True)
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
    
    def get_id(self, feature_ids, index = False):
        '''
        Get GFF entry by feature ID
        '''
        indices = [i for i, entry in enumerate(self._data) if entry.has_attr("ID", feature_ids)]
        if not indices and type(feature_ids) is str: return None
        elif type(feature_ids) is str: return self.get_i(indices[0])
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
    def get_attr(self, a): return self.attributes.get(a)
    def is_attr(self, a, val): return self.attributes.is_attr(a, val)
    def has_attr(self, a, vals): return self.attributes.has_attr(a, vals)

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
