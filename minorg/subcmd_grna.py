import os
import itertools
import regex as re

from functions import (
    splitlines,
    fasta_to_dict,
    dict_to_fasta,
    infer_full_pam,
    expand_ambiguous,
    make_pam_pattern,
    Target,
    gRNAHits,
    gRNAHit
)

## These functions will not do any type checking.
## All args are assumed to have been properly parsed.

def execute_grna(args, config, directory, prefix, target):
    all_gRNA_fname = config.mkfname(directory, f"{prefix}_gRNA_all.fasta")
    all_gRNA = find_multi_gRNA(target, pam = args.pam, gRNA_len = args.length)
    all_gRNA.assign_seqid(prefix = "gRNA_")
    all_gRNA.write_fasta(all_gRNA_fname, write_all = True)
    return all_gRNA, all_gRNA_fname, all_gRNA_fname

#################
##  FIND gRNA  ##
#################

def find_multi_gRNA(target_fname, pam = "NGG", gRNA_len = 20): #Let NGG be GG, or NG be G
    """
    Finds all possible target sequences given pam.
    Returns: {<seq>: set(ids of targeted targets)}
    """
    ## specifying 'N' allows users to specify how many wobble/N bases
    pam_pattern = make_pam_pattern(pam, gRNA_len = gRNA_len)
    # TODO: some function to parse NXX or XXN PAM into regex (done)
    dic_target = {}
    from Bio import SeqIO
    parsed = SeqIO.parse(target_fname, "fasta")
    hit_id = 0
    for sequence in parsed:
        fwd_seq, rvs_seq = sequence.seq, sequence.seq.reverse_complement()
        target = Target(fwd_seq, id = sequence.id, strand = '+')
        for strand, seq in {'+': fwd_seq, '-': rvs_seq}.items():
            pam_pos = [x.start() for x in re.finditer(pam_pattern, str(seq), overlapped = True)]
            for start in pam_pos:
                end = start + gRNA_len
                grna_seq = str(seq[start: start + gRNA_len]).upper()
                if len(grna_seq) == gRNA_len:
                    ## gRNAHit's target is set as fwd_seq (variable 'target')
                    ## gRNAHit's range is set relative to fwd_seq (direction-less)
                    gRNAhit = gRNAHit(target,
                                      start if strand == '+' else len(target) - end,
                                      end if strand == '+' else len(target) - start,
                                      strand, hit_id)
                    dic_target[grna_seq] = dic_target.get(grna_seq, set()).union({gRNAhit})
                    hit_id += 1
    # return dic_target
    output_gRNA_hits = gRNAHits()
    output_gRNA_hits.parse_from_dict(dic_target)
    return output_gRNA_hits


# ###############
# ##  CLASSES  ##
# ###############

# class CheckObj:
#     def __init__(self, *check_names):
#         self._checks = {check_name: None for check_name in check_names}
#     def set_check(self, check_name, status):
#         self._checks[check_name] = (True if status in ("pass", True) else
#                                     False if status in ("fail", False) else
#                                     None)
#     def check(self, check_name, mode = "raw"):
#         check_value = self._checks[check_name]
#         if mode == "str":
#             if check_value == True: return "pass"
#             elif check_value == False: return "fail"
#             else: return "NA"
#         else:
#             return check_value
#     def some_checks_passed(self, *check_names, **kwargs):
#         return all(self.check(check_name, **kwargs) for check_name in check_names)
#     def some_valid_checks_passed(self, *check_names):
#         return tuple(self.check(check_name) for check_name in check_names).count(False) == 0
#     def all_valid_checks_passed(self):
#         return self.some_valid_checks_passed(*self._checks.keys())
#     def all_checks_passed(self):
#         return self.some_checks_passed(*self._checks.keys())

# class Target(CheckObj):
#     ## seqs are stored in uppercase
#     def __init__(self, seq, id = None, strand = None):
#         self._seq = str(seq).upper()
#         self._id = id
#         self._strand = strand ## possible values: None, '+', '-'
#         self._sense = None ## possible values: None, '+', '-'
#     def __str__(self): return self.seq()
#     def __repr__(self): return str(self)
#     def __len__(self): return len(str(self))
#     def id(self): return self._id
#     def seq(self): return self._seq
#     def strand(self): return self._strand
#     def sense(self): return self._sense
#     def set_strand(self, strand): self._strand = strand
#     def set_sense(self, sense): self._sense = sense
#     def set_sense_by_parent(self, parent_sense):
#         parent_sense = ('+' if parent_sense in ('+', "sense") else '-')
#         self._sense = (None if self.strand() == None else
#                        '+' if parent_sense == self.strand() else '-')
#     def parent_sense(self, mode = "raw"):
#         parent_sense_value = ( None if (None in (self.sense(), self.strand())) else
#                                '+' if (self.strand() == self.sense()) else
#                                '-' )
#         if mode == "str":
#             if parent_sense_value == '+': return "sense"
#             elif parent_sense_value == '-': return "antisense"
#             else: return "NA"
#         else:
#             return parent_sense_value
#     def valid_len(self): return ( len(self) > 0 )

# class gRNASeq(CheckObj):
#     ## seqs are stored in uppercase
#     def __init__(self, seq):
#         super().__init__("background", "exclude", "GC")
#         self._seq = str(seq).upper()
#         self._id = None
#     def __str__(self): return self.seq()
#     def __repr__(self): return str(self)
#     def id(self): return self._id
#     def seq(self): return self._seq
#     def set_id(self, id): self._id = id
#     def set_bg_check(self, status): self.set_check("background", status)
#     def set_exclude_check(self, status): self.set_check("exclude", status)
#     def set_gc_check(self, status = None, gc_min = 0, gc_max = 1):
#         if status:
#             self.set_check("GC", status)
#         else:
#             from seq_manip import gc_content
#             self.set_check("GC", gc_min <= gc_content(self.seq()) <= gc_max)
#         return

# class gRNAHits:
#     ## seqs are stored in uppercase
#     def __init__(self, d = {}, gRNA_seqs = {}, gRNA_hits = {}):
#         self._gRNAseqs = gRNA_seqs ## dictionary of {seq: <gRNASeq obj>}
#         self._hits = gRNA_hits ## dictionary of {seq: [list of <gRNAHit obj>]}
#         if d:
#             self.parse_from_dict(d)
#     def __repr__(self): return self.hits()
#     def __len__(self): return len(self.seqs())
#     def update_records(self): ## update dictionaries to remove any discrepancies
#         self._gRNAseqs = {k: v for k, v in self.gRNAseqs() if k in self.hits()}
#         self._hits = {k: v for k, v in self.hits() if k in self.gRNAseqs()}
#     def copy(self):
#         from copy import deepcopy
#         new_obj = gRNAHits()
#         new_obj._gRNAseqs = deepcopy(self.gRNAseqs())
#         new_obj._hits = deepcopy(self.hits())
#         return new_obj
#     ################
#     ##  BOOLEANS  ##
#     ################
#     def all_target_len_valid(self): return all(map(lambda hit: hit.target().valid_len(), self.flatten_hits()))
#     ###############
#     ##  PARSERS  ##
#     ###############
#     def parse_from_dict(self, d): ## where d = {str(seq): [list of <gRNAHit obj>s]}
#         from copy import deepcopy
#         self._gRNAseqs = {str(seq).upper(): gRNASeq(seq) for seq in d.keys()}
#         self._hits = deepcopy(d)
#     def parse_from_mapping(self, fname, targets = None, version = 1):
#         ## read data
#         seq_targets = {} if not targets else fasta_to_dict(targets)
#         raw_mapping = [line.split('\t') for line in splitlines(fname)]
#         header_mapping = raw_mapping[0]
#         ## determine version where checks start
#         if not version: version = (1 if "exclusive" in header_mapping[5] else \
#                                    3 if "target length" in header_mapping else 2)
#         check_col = header_mapping.index("group") + 1
#         # if version == 1: check_col = 7
#         # elif version == 2: check_col = 8
#         header_checks = header_mapping[check_col:]
#         dat_mapping = raw_mapping[1:]
#         del raw_mapping
#         for i, entry in enumerate(dat_mapping):
#             ## note: gRNA_range is relative to + strand (not necessarily sense) when read from ..targets.txt file
#             ## note: gRNA_id is probably meaningless as it will be overwritten using ids in the fasta file supplied to the variable 'fasta' in get_minimum_set_from_file
#             if version == 1:
#                 gRNA_id, gRNA_seq, target_id, sense, strand, gRNA_range, group = entry[:check_col]
#                 gRNA_start, gRNA_end = map(int, re.search("\[(\d+), (\d+)\)", gRNA_range).group(1,2))
#                 target_len = 0
#             elif version == 2:
#                 ## start, end are 1-indexed, start-inclusive and end-inclusive in v2
#                 gRNA_id, gRNA_seq, target_id, sense, strand, gRNA_start, gRNA_end, group = entry[:check_col]
#                 gRNA_start = int(gRNA_start) - 1 ## convert to 0-indexed, start-inclusive
#                 gRNA_end = int(gRNA_end) ## "convert" to 0-indexed, end-exclusive
#                 target_len = 0
#             elif version == 3:
#                 ## start, end are 1-indexed, start-inclusive and end-inclusive in v2
#                 gRNA_id, gRNA_seq, target_id, target_len, sense, strand, gRNA_start, gRNA_end, group = entry[:check_col]
#                 gRNA_start = int(gRNA_start) - 1 ## convert to 0-indexed, start-inclusive
#                 gRNA_end = int(gRNA_end) ## "convert" to 0-indexed, end-exclusive
#             try:
#                 ## note: dummy target sequence is used if targets file not provided; target strand assumed '+'
#                 target = Target(seq_targets.get(target_id, 'N' * int(target_len)), id = target_id, strand = '+')
#             except Exception as e:
#                 print(version, i, header_mapping, entry)
#                 raise e
#             ## create gRNAHit object
#             gRNA_hit = gRNAHit(target, gRNA_start, gRNA_end, strand, gRNA_id)
#             gRNA_hit.set_parent_sense(sense)
#             ## add gRNAHit object to gRNAHits object
#             self.add_hit(gRNA_seq, gRNA_hit)
#             self.get_gRNAseq_by_seq(gRNA_seq).set_id(gRNA_id)
#             ## log checks
#             gRNA_checks = entry[check_col:]
#             for i, check in enumerate(header_checks):
#                 if check == "CDS": ## this is the only one that's affected by hit location
#                     gRNA_hit.set_check(check, gRNA_checks[i])
#                 else: ## all other checks apply to all hits with the same gRNA seq so set check to gRNASeq obj
#                     self.get_gRNAseq_by_seq(gRNA_seq).set_check(check, gRNA_checks[i])
#         return
#     #################
#     ##  MODIFIERS  ##
#     #################
#     def assign_gRNAseq_id(self, fasta):
#         fasta_inv = {str(v): k for k, v in fasta_to_dict(fasta).items()}
#         valid_seqs = set(fasta_inv.keys()) & set(self.seqs())
#         if len(valid_seqs) < len(set(self.seqs())):
#             print("\nWARNING: The provided FASTA file does not cover all gRNAs.\n")
#         for seq in valid_seqs:
#             self.get_gRNAseq_by_seq(seq).set_id(fasta_inv[seq])
#         return
#     def add_hit(self, seq, gRNA_hit):
#         self.add_seq(seq)
#         self._hits[str(seq)] = self.get_hits(seq) + [gRNA_hit]
#     def add_seq(self, seq):
#         if not str(seq) in self.gRNAseqs(): self._gRNAseqs[str(seq)] = gRNASeq(seq)
#         if not str(seq) in self.hits(): self._hits[str(seq)] = []
#     def remove_seqs(self, *seqs):
#         if seqs and type(seqs[0]) in (list, tuple, set):
#             seqs = list(itertools.chain(*seqs))
#         for seq in seqs:
#             if str(seq) in self.gRNAseqs(): del self._gRNAseqs[str(seq)]
#             if str(seq) in self.hits(): del self._hits[str(seq)]
#         return
#     ###############
#     ##  SETTERS  ##
#     ###############
#     def set_seqs_check(self, check_name, status, seqs):
#         if type(seqs) not in (tuple, list, set):
#             seqs = (seqs,)
#         for seq in seqs:
#             self.get_gRNAseq_by_seq(seq).set_check(check_name, status)
#         return
#     def set_seqs_check_by_function(self, check_name, func, seqs):
#         for seq in seqs:
#             self.set_seqs_check(check_name, func(self.get_gRNAseq_by_seq(seq)), [seq])
#         return
#     def set_all_seqs_check_by_function(self, check_name, func):
#         self.set_seqs_check_by_function(check_name, func, self.seqs())
#     def rename_seqs(self, fasta):
#         seqs_names = {str(v): k for k, v in fasta_to_dict(fasta).items()}
#         for seq, name in seqs_names.items():
#             gRNA_seq = self.get_gRNAseq_by_seq(seq)
#             if gRNA_seq:
#                 gRNA_seq.set_id(name)
#         return
#     def assign_seqid(self, prefix = "gRNA_", zfill = 3):
#         for i, gRNA_seq in enumerate(self.flatten_gRNAseqs()):
#             gRNA_seq.set_id(f"{prefix}{str(i+1).zfill(zfill)}")
#         return
#     ###############
#     ##  GETTERS  ##
#     ###############
#     def gRNAseqs(self): return self._gRNAseqs
#     def seqs(self, output_type = list): return output_type(self.gRNAseqs().keys())
#     def hits(self): return self._hits
#     def flatten_hits(self, output_type=list): return output_type(itertools.chain(*self.hits().values()))
#     def flatten_gRNAseqs(self, output_type=list): return output_type(self.gRNAseqs().values())
#     def get_hits(self, seq):
#         return self.hits().get(str(seq).upper(), []) ## 'seq' must be able to be coerced using str()
#     def get_gRNAseq_by_seq(self, seq):
#         return self.gRNAseqs().get(str(seq).upper(), None) ## retrieve gRNASeq obj by seq
#     def get_gRNAseq_by_id(self, id):
#         output = [gRNA_seq for gRNA_seq in self.flatten_gRNAseqs() if gRNA_seq.id() == id]
#         if not output: return None
#         elif len(output) == 1: return output[0]
#         else:
#             print("\nWARNING: There are multiple gRNA sequences with the requested ID. Returning first in list.")
#             return output[0]
#     def get_gRNAseqs_by_seq(self, *seqs, output_type = list, ignore_invalid = True):
#         if not seqs:
#             return []
#         if type(seqs[0]) in (list, tuple, set):
#             seqs = tuple(itertools.chain(*seqs))
#         return output_type(self.get_gRNAseq_by_seq(seq) for seq in seqs if self.get_gRNAseq_by_seq(seq) != None)
#     def get_gRNAseqs_by_id(self, *ids, output_type = list, ignore_invalid = True):
#         if type(ids[0]) in (list, tuple, set):
#             ids = tuple(itertools.chain(*ids))
#         return output_type(self.get_gRNAseq_by_id(id) for id in ids if self.get_gRNAseq_by_id(seq) != None)
#     ######################
#     ##      FILTER      ##
#     ##  (& return new)  ##
#     ######################
#     ## filter gRNAHit objects for certain criteria and return new gRNAHits object
#     def filter_hits(self, *check_names, exclude_empty_seqs = True, ignore_invalid = True):
#         filtered = {seq: [hit for hit in hits
#                           if ( (check_names and ( (ignore_invalid and hit.some_valid_checks_passed(*check_names))
#                                                   or ( (not ignore_invalid) and
#                                                        hit.some_checks_passed(*check_names) ) ) ) or
#                                ( (not check_names) and ( (ignore_invalid and hit.all_valid_checks_passed())
#                                                          or ( (not ignore_invalid) and
#                                                               hit.all_checks_passed() ) ) ) ) ]
#                     for seq, hits in self.hits().items()}
#         if exclude_empty_seqs:
#             filtered = {seq: hits for seq, hits in filtered.items() if hits}
#         filtered_seqs = {seq: self.get_gRNAseq_by_seq(seq) for seq in filtered.keys()}
#         output = gRNAHits(gRNA_seqs = filtered_seqs, gRNA_hits = filtered)
#         return output
#     def filter_hits_some_checks_passed(self, *check_names, **kwargs):
#         return self.filter_hits(*check_names, **kwargs)
#     def filter_hits_all_checks_passed(self, **kwargs):
#         return self.filter_hits(**kwargs)
#     ## filter gRNASeq objects for certain criteria and return new gRNAHits object
#     def filter_seqs(self, *check_names, ignore_invalid = True):
#         filtered = {seq: hits for seq, hits in self.hits().items()
#                     if ((check_names and ((ignore_invalid and
#                                            self.get_gRNAseq_by_seq(seq).some_valid_checks_passed(*check_names))
#                                           or (not ignore_invalid and
#                                               self.get_gRNAseq_by_seq(seq).some_checks_passed(*check_names)))) or
#                          ((not check_names) and ((ignore_invalid and
#                                                   self.get_gRNAseq_by_seq(seq).all_valid_checks_passed())
#                                                  or (not ignore_invalid and
#                                                      self.get_gRNAseq_by_seq(seq).all_checks_passed()))))}
#         filtered_seqs = {seq: self.get_gRNAseq_by_seq(seq) for seq in filtered.keys()}
#         output = gRNAHits(gRNA_seqs = filtered_seqs, gRNA_hits = filtered)
#         return output
#     def filter_seqs_some_checks_passed(self, *check_names, **kwargs):
#         return self.filter_seqs(*check_names, **kwargs)
#     def filter_seqs_all_checks_passed(self, **kwargs):
#         return self.filter_seqs(**kwargs)
#     #############
#     ##  WRITE  ##
#     #############
#     def write_mapping(self, fout, sets = [], write_all = False,
#                       write_checks = False, checks = ["background", "GC", "CDS"],
#                       index = 1, start_incl = True, end_incl = True, version = 3,
#                       fasta = None ): ## if fasta is provided, it will override default naming behaviour
#         if not write_checks: checks = []
#         if sets: sets = sets
#         elif write_all: sets = [set(self.seqs())]
#         else: []
#         fasta_inv = {} if not fasta else {str(seq): k for k, seq in fasta_to_dict(fasta).items()}
#         if version == 1:
#             mapping_header = ["gRNA id", "gRNA sequence", "target id", "target sense",
#                               "gRNA strand", "range (0-index, end exclusive)", "group"] + checks
#         elif version == 2:
#             mapping_header = ["gRNA id", "gRNA sequence", "target id", "target sense",
#                               "gRNA strand", "start", "end", "group"] + checks
#         else:
#             mapping_header = ["gRNA id", "gRNA sequence", "target id", "target length", "target sense",
#                               "gRNA strand", "start", "end", "group"] + checks
#         mapping_dat = []
#         seq_check_names = {"GC", "background"}
#         for group, seqs in enumerate(sets):
#             for seq in seqs:
#                 gRNA_seq = self.get_gRNAseq_by_seq(seq)
#                 gRNA_hits = self.get_hits(seq)
#                 for gRNA_hit in gRNA_hits:
#                     ## put together all required fields
#                     if version == 1:
#                         mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id()), ## override id if fasta provided
#                                             gRNA_seq.seq(), gRNA_hit.target_id(),
#                                             gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand(),
#                                             "[{}, {})".format(*gRNA_hit.range()), group + 1] +
#                                            [(gRNA_seq.check(check_name, mode = "str")
#                                              if check_name in seq_check_names else
#                                              gRNA_hit.check(check_name, mode = "str"))
#                                             for check_name in checks])
#                     elif version == 2:
#                         ## figure out start and end
#                         start, end = map(lambda n: n + index, gRNA_hit.range())
#                         if not start_incl: start -= 1
#                         if end_incl: end -= 1
#                         mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id()), ## override id if fasta provided
#                                             gRNA_seq.seq(), gRNA_hit.target_id(),
#                                             gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand(),
#                                             start, end, group + 1] +
#                                            [(gRNA_seq.check(check_name, mode = "str")
#                                              if check_name in seq_check_names else
#                                              gRNA_hit.check(check_name, mode = "str"))
#                                             for check_name in checks])
#                     else:
#                         ## figure out start and end
#                         start, end = map(lambda n: n + index, gRNA_hit.range())
#                         if not start_incl: start -= 1
#                         if end_incl: end -= 1
#                         mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id()), ## override id if fasta provided
#                                             gRNA_seq.seq(), gRNA_hit.target_id(), gRNA_hit.target_len(),
#                                             gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand(),
#                                             start, end, group + 1] +
#                                            [(gRNA_seq.check(check_name, mode = "str")
#                                              if check_name in seq_check_names else
#                                              gRNA_hit.check(check_name, mode = "str"))
#                                             for check_name in checks])
#         mapping_dat.sort(key = lambda entry: int(entry[mapping_header.index("start")])) ## sort by start
#         mapping_dat.sort(key = lambda entry: entry[mapping_header.index("target id")]) ## sort by target id
#         mapping_dat.sort(key = lambda entry: entry[mapping_header.index("gRNA id")]) ## sort by gRNA id
#         mapping_dat.sort(key = lambda entry: entry[mapping_header.index("group")])
#         to_write = '\n'.join(['\t'.join(map(str, entry)) for entry in ([mapping_header] + mapping_dat)]) + '\n'
#         open(fout, "w+").write(to_write)
#         return
#     def write_fasta(self, fout, seqs = [], ids = [], write_all = False, fasta = None):
#         ## get relevant gRNA sequences
#         if seqs: gRNA_seqs = self.get_gRNAseqs_by_seq(*seqs)
#         elif ids: gRNA_seqs = self.get_gRNAseqs_by_id(*ids)
#         elif write_all: gRNA_seqs = self.flatten_gRNAseqs()
#         else: open(fout, "w+").write('')
#         ## rename sequences per fasta file (if fasta file provided)
#         fasta_inv = {} if not fasta else {str(seq): k for k, seq in fasta_to_dict(fasta).items()}
#         to_write = {fasta_inv.get(gRNA_seq.seq(), gRNA_seq.id()): gRNA_seq.seq() for gRNA_seq in gRNA_seqs}
#         if to_write:
#             dict_to_fasta(to_write, fout)
#         else:
#             open(fout, "w+").write('')
#         return

# class gRNAHit(CheckObj):
#     def __init__(self, target, start, end, strand, hit_id):
#         ## note: unless something is weird, seq_strand is the same as strand
#         super().__init__("background", "GC", "CDS", "exclude", "flank")
#         self._target = target ## previously self._seq
#         self._range = (start, end) ## relative to original (parent) sequence from which _target was derived
#         self._strand = '+' if (strand.lower() == "fwd" or strand == '+') else '-' ## gRNA strand relative to original (parent) sequence from which _target was derived (gRNA strand should be same as _target strand)
#         # self._seq_strand = '' # stores _target direction relative to original sequence from which _target was derived
#         self._parent_sense = None ## - if original (parent) sequence from which _target was derived is on antisense strand, + if _target is on sense strand
#         self._hit_id = hit_id
#     def target(self): return self._target
#     def start(self): return self._range[0]
#     def end(self): return self._range[1]
#     def target_id(self): return self.target().id()
#     def strand(self): return self._strand
#     def hit_id(self): return self._hit_id
#     def range(self, output_type = tuple): return output_type(self._range)
#     def target_len(self): return len(self.target())
#     def target_strand(self): return self.target().strand()
#     # def set_target_strand(self, strand): self._target_strand = strand
#     def set_parent_sense(self, strand):
#         self.target().set_sense_by_parent(strand)
#         # self._parent_sense = ('+' if strand in ('+', "sense") else '-')
#     def set_bg_check(self, status): self.set_check("background", status)
#     def set_gc_check(self, status): self.set_check("GC", status)
#     def set_cds_check(self, status): self.set_check("CDS", status)
#     def set_exclude_check(self, status): self.set_check("exclude", status)
#     def set_flank_check(self, status): self.set_check("flank", status)
#     def reverse_range(self): return (self.target_len() - self.end(), self.target_len() - self.start())
#     def parent_sense(self, mode = "raw"):
#         parent_sense_value = self.target().sense()
#         if mode == "str":
#             if parent_sense_value == '+':
#                 return "sense"
#             elif parent_sense_value == '-':
#                 return "antisense"
#             else:
#                 return "NA"
#         else:
#             return parent_sense_value
#     def adj_range(self, mode = "strand"):
#         if (not self.strand() or
#             (mode == "target" and not self.target_strand()) or
#             (mode == "gene" and not self.parent_sense())):
#             return (float("nan"), float("nan"))
#         elif ((mode == "strand" and (self.strand() == '+')) or \
#               (mode == "target" and (self.target_strand() == self.strand())) or \
#               (mode == "gene" and (self.parent_sense() == '+'))):
#             return self.range()
#         return self.reverse_range()
#     def flank(self, length = 100):
#         target_seq = self.target() if self.strand == '+' else self.target().reverse_complement()
#         start, end  = self.adj_range(mode = "target")
#         return target_seq[start - length: start], target_seq[end: end + length]
