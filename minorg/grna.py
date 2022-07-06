import itertools
import warnings

from typing import Union

from minorg import MINORgWarning

from minorg.functions import (
    gc_content,
    splitlines
)

from minorg.fasta import (
    fasta_to_dict,
    dict_to_fasta
)

####################
##  gRNA CLASSES  ##
####################

class CheckObj:
    """
    Object with checks that can be set.
    
    Attributes:
        _checks (dict): stores check values (format: {'<check name>': <check value (True, False, or None)>}
    """
    def __init__(self, *check_names):
        self._checks = {check_name: None for check_name in check_names}
    @property
    def check_names(self) -> list:
        """
        Get all check names
        
        Returns
        -------
        list of str
            Of check names
        """
        return list(self._checks.keys())
    def clear_checks(self) -> None:
        """
        Set values of all checks to None.
        """
        for check_name in self._checks.keys():
            self._checks[check_name] = None
    def set_check(self, check_name, status) -> None:
        """
        Set check value.
        
        Arguments:
            check_name (str): check name
            status (str or bool or None): check status
        
        Check value will be set to:
            - True: If status = "pass" OR status = True
            - False: If status = "fail" OR status = False
            - None: If status is any other value
        """
        self._checks[check_name] = (True if status in ("pass", True) else
                                    False if status in ("fail", False) else
                                    None)
    def check(self, check_name, mode = "bool") -> Union[str, bool, None]:
        """
        Get check status.
        
        Arguments:
            check_name (str): check name
            mode (str): return type (valid values: "bool", "str")
        
        Returns
        -------
        str
            If ``mode = 'str'``.
            If check status is True, returns 'pass'.
            If check status is False, returns 'fail'.
            If check status is None, returns 'NA'.
        bool
            If ``mode = 'bool'`` AND check status is not None
        None
            If ``mode = 'bool'`` AND check status is None
        """
        check_value = self._checks[check_name]
        if mode == "str":
            if check_value == True: return "pass"
            elif check_value == False: return "fail"
            else: return "NA"
        elif mode == "bool":
            return check_value
        else:
            raise Exception(f"Invalid return type: '{mode}'")
    def check_exists(self, check_name) -> bool:
        """
        Whether a given check exists.
        
        Arguments:
            check_name (str): check name
        
        Returns
        -------
        bool
            Whether ``check_name`` is a valid check name
        """
        return check_name in self._checks
    def some_checks_passed(self, *check_names) -> bool:
        """
        Whether self passes all of some combination of checks.
        
        Arguments:
            *check_names (str): check names
        
        Returns
        -------
        bool
            Whether self passed all of the specified checks
        """
        return all(self.check(check_name) for check_name in check_names)
    def some_valid_checks_passed(self, *check_names):
        """
        Whether self passes all of the valid checks in some combination of checks.
        
        Only returns True if none of the values of the specified checks is False.
        
        Arguments:
            *check_names (str): check names
        
        Returns
        -------
        bool
            Whether self's values for all specified checks are True OR None BUT NOT False
        """
        return tuple(self.check(check_name) for check_name in check_names).count(False) == 0
    def all_valid_checks_passed(self) -> bool:
        """
        Whether self passes all valid checks.
        
        Return True if the values of all checks are True OR None BUT NOT False.
        
        Returns
        -------
        bool
        """
        return self.some_valid_checks_passed(*self._checks.keys())
    def all_checks_passed(self) -> bool:
        """
        Whether self passes all checks.
        
        Return True only if the values of all checks are True.
        
        Returns
        -------
        bool
        """
        return self.some_checks_passed(*self._checks.keys())

class Target(CheckObj):
    """
    gRNA target sequence.
    
    Attributes:
        _seq (str): sequence, in uppercase
        _id (str): sequence name
        _strand (str or None): strand (possible values: None, '+', '-')
        _sense (str or None): sense
            (possible values: None, '+', '-'; where '+' means sense and '-' means antisense)
    """
    ## seqs are stored in uppercase
    def __init__(self, seq, id = None, strand = None):
        """
        Create a Target object.
        
        Arguments:
            seq (str or Bio.Seq.Seq): sequence
            id (str): sequence name
            strand (str): sequence strand (values: '+', '-')
        """
        self._seq = str(seq).upper()
        self._id = id
        self._strand = strand ## possible values: None, '+', '-'
        self._sense = None ## possible values: None, '+', '-'
    def __str__(self): return self.seq
    def __repr__(self): return str(self)
    def __len__(self): return len(self._seq)
    @property
    def id(self) -> str: return self._id
    @property
    def seq(self) -> str: return self._seq
    @property
    def strand(self) -> str: return self._strand
    @property
    def sense(self) -> str: return self._sense
    def set_sense(self, sense) -> None:
        """
        Set sense.
        
        Arguments:
            sense (str): sense of sequence. Valid values: '+', 'sense', '-', 'antisense'
        """
        self._sense = sense
    def set_sense_by_parent(self, parent_sense) -> None:
        """
        Set sense using parent sequence's sense AND self's strand.
        
        Sense if
            - ``parent_sense='+'`` AND self's strand is '+'
            - OR ``parent_sense='-'`` AND self's strand is '-'
        Else antisense.
        
        Arguments:
            parent_sense (str): sense of parent sequence. Valid values: '+', 'sense', '-', 'antisense'
        """
        parent_sense = (None if parent_sense == "NA"
                        else '+' if parent_sense in ('+', "sense") else '-')
        self._sense = (None if (self.strand is None or parent_sense is None) else
                       '+' if parent_sense == self.strand else '-')
    def parent_sense(self, mode = "raw") -> Union[str, None]:
        """
        Return parent sequence's sense status.
        
        Arguments:
            mode (str): return string format (valid values: "raw", "str")
        
        Returns
        -------
        str
            If ``mode='str'`` ('sense', 'antisense', 'NA') OR parent sense has been set ('+', '-')
        None
            If ``mode='raw'`` AND parent sense has not been set
        """
        parent_sense_value = ( None if (None in (self.sense, self.strand)) else
                               '+' if (self.strand == self.sense) else
                               '-' )
        if mode == "str":
            if parent_sense_value == '+': return "sense"
            elif parent_sense_value == '-': return "antisense"
            else: return "NA"
        else:
            return parent_sense_value
    def valid_len(self) -> bool:
        """
        Whether length of sequence is equal to or greater than 1.
        """
        return ( len(self) > 0 )

class gRNASeq(CheckObj):
    """
    gRNA sequence.
    
    Tracks checks applicable to gRNA sequence: off-target (check name: background), GC content (check name: GC)
    
    Attributes:
        _seq (str): sequence, in uppercase
        _id (str): sequence name
    """
    ## seqs are stored in uppercase
    def __init__(self, seq):
        """
        Create a gRNASeq object.
        
        Arguments:
            seq (str or Bio.Seq.Seq): sequence
        """
        # super().__init__("background", "exclude", "GC")
        super().__init__("background", "GC")
        self._seq = str(seq).upper()
        self._id = None
    def __str__(self): return self.seq
    def __repr__(self): return str(self)
    @property
    def id(self) -> str: return self._id
    @property
    def seq(self) -> str: return self._seq
    @id.setter
    def id(self, id) -> None: self._id = id
    def set_bg_check(self, status) -> None:
        """
        Set status for check name 'background' (False if gRNA has off-target effects).
        
        Arguments:
            status (str or bool or None): status of check.
                Valid values: 'pass', 'fail', 'NA', True, False, None
        """
        self.set_check("background", status)
    def set_exclude_check(self, status) -> None:
        """
        Set status for check name 'exclude' (False if user wishes to exclude gRNA with this sequence).
        
        Arguments:
            status (str or bool or None): status of check (valid values: 'pass', 'fail', 'NA', True, False, None)
        """
        self.set_check("exclude", status)
    def set_gc_check(self, status = None, gc_min = 0, gc_max = 1) -> None:
        """
        Set status for check name 'GC' (False if gRNA GC content is not within desired range).
        
        Arguments:
            status (str or bool or None): status of check (valid values: 'pass', 'fail', 'NA', True, False, None)
                If ``status`` is provided, ``gc_min`` and ``gc_max`` will be ignored.
            gc_min (float): minimum GC content (betweew 0 and 1, where 0 is no GC content and 1 is all GC)
            gc_max (float): maximum GC content (betweew 0 and 1, where 0 is no GC content and 1 is all GC)
        """
        if status:
            self.set_check("GC", status)
        else:
            self.set_check("GC", gc_min <= gc_content(self.seq) <= gc_max)
        return

class gRNAHits:
    """
    Tracks multiple gRNASeq and gRNAHit objects.
    
    Attributes:
        _gRNAseqs (dict): stores gRNASeq objects by sequence (format: {'<seq>': <gRNASeq object>})
        _hits (dict): stores gRNAHit objects by sequence (format: {'<seq>': [<gRNAHit objects>]})
    """
    ## seqs are stored in uppercase
    def __init__(self, d = None, gRNA_seqs = None, gRNA_hits = None):
        """
        Create a gRNAHits object.
        
        Arguments:
            d (dict): dictionary of gRNAHit objects in same format as :attr:`~minorg.grna.gRNAHits._hits`.
                If provided, overrides ``gRNA_seqs`` and ``gRNA_hits``.
            gRNA_seqs (dict): dictionary of gRNASeq objects in same format as :attr:`~minorg.grna.gRNAHits._gRNAseqs`
            gRNA_hits (dict): dictionary of gRNAHit objects in same format as :attr:`~minorg.grna.gRNAHits._hits`
        """
        self._gRNAseqs = {} if gRNA_seqs is None else gRNA_seqs ## dictionary of {seq: <gRNASeq obj>}
        self._hits = {} if gRNA_hits is None else gRNA_hits ## dictionary of {seq: [list of <gRNAHit obj>]}
        if d:
            self.parse_from_dict(d)
    def __repr__(self): return f"gRNAHits(gRNA = {len(self)})"
    def __len__(self): return len(self.seqs)
    @property
    def gRNAseqs(self) -> dict: return self._gRNAseqs
    @property
    def hits(self) -> dict: return self._hits
    @property
    def seqs(self) -> list:
        """
        List of gRNA sequences.
        
        :type: list of str
        """
        return list(self.gRNAseqs.keys())
    @property
    def check_names(self) -> list:
        """
        List of check names.
        
        :type: list of str
        """
        return list(set(self.check_names_seqs + self.check_names_hits))
    @property
    def check_names_seqs(self) -> list:
        """
        List of check names for gRNA sequences (gRNASeq objects).
        
        :type: list of str
        """
        return list(set(itertools.chain(*[grna_seq.check_names for grna_seq in self.flatten_gRNAseqs()])))
    @property
    def check_names_hits(self) -> list:
        """
        List of check names for gRNA hits (gRNAHit objects).
        
        :type: list of str
        """
        return list(set(itertools.chain(*[grna_hit.check_names for grna_hit in self.flatten_hits()])))
    def update_records(self) -> None: ## update dictionaries to remove any discrepancies
        """
        Remove gRNASeq objects from :attr:`~minorg.grna.gRNAHits._gRNAseqs` 
        if their sequences are not also in :attr:`~minorg.grna.gRNAHits._hits`.
        Also remove gRNAHit objects from :attr:`~minorg.grna.gRNAHits._hits`
        if their sequences are not also in :attr:`~minorg.grna.gRNAHits._gRNAseqs`.
        """
        self._gRNAseqs = {k: v for k, v in self.gRNAseqs if k in self.hits}
        self._hits = {k: v for k, v in self.hits if k in self.gRNAseqs}
    def copy(self) -> 'gRNAHits':
        """
        Deepcopy self to new :class:`~minorg.grna.gRNAHits` object.
        
        Returns
        -------
        :class:`~minorg.grna.gRNAHits`
        """
        from copy import deepcopy
        new_obj = gRNAHits()
        new_obj._gRNAseqs = deepcopy(self.gRNAseqs)
        new_obj._hits = deepcopy(self.hits)
        return new_obj
    ################
    ##  BOOLEANS  ##
    ################
    def all_target_len_valid(self) -> bool:
        """
        Whether all targets of gRNA have valid length (i.e. is greater than 0).
        
        Returns
        -------
        bool
        """
        return all(map(lambda hit: hit.target.valid_len(), self.flatten_hits()))
    def valid_seq_check(self, check_name) -> bool:
        """
        Whether a given check has been set for at least one gRNA sequence.
        
        Arguments:
            check_name (str): check name
        
        Returns
        -------
        bool
        """
        return set(seq.check(check_name) for seq in self.flatten_gRNAseqs()) != {None}
    def valid_hit_check(self, check_name) -> bool:
        """
        Whether a given check has been set for at least one gRNA hit.
        
        Arguments:
            check_name (str): check name
        
        Returns
        -------
        bool
        """
        return set(hit.check(check_name) for hit in self.flatten_hits()) != {None}
    def unused_seq_check(self, check_name) -> bool:
        """
        Whether a given check is NOT set for all gRNA sequences.
        
        Arguments:
            check_name (str): check name
        
        Returns
        -------
        bool
        """
        return set(seq.check(check_name) for seq in self.flatten_gRNAseqs()) == {None}
    def unused_hit_check(self, check_name) -> bool:
        """
        Whether a given check is NOT set for all gRNA hits.
        
        Arguments:
            check_name (str): check name
        
        Returns
        -------
        bool
        """
        return set(hit.check(check_name) for hit in self.flatten_hits()) == {None}
    def seq_check_exists(self, check_name):
        """
        Whether a given check exists for all gRNA sequences.
        
        Arguments:
            check_name (str): check name
        
        Returns
        -------
        bool
        """
        return all(seq.check_exists(check_name) for seq in self.flatten_gRNAseqs())
    def hit_check_exists(self, check_name):
        """
        Whether a given check exists for all gRNA hits.
        
        Arguments:
            check_name (str): check name
        
        Returns
        -------
        bool
        """
        return all(hit.check_exists(check_name) for hit in self.flatten_hits())
    ###############
    ##  PARSERS  ##
    ###############
    def parse_from_dict(self, d) -> None: ## where d = {str(seq): [list of <gRNAHit obj>s]}
        """
        Read data from dictionary of gRNAHit objects and deep copy it to self.
        
        Arguments:
            d (dict): dictionary of gRNAHit objects in same format as :attr:`~minorg.grna.gRNAHits._hits`.
        """
        from copy import deepcopy
        self._gRNAseqs = {str(seq).upper(): gRNASeq(seq) for seq in d.keys()}
        self._hits = deepcopy(d)
        check_names = set()
        for hit in self.flatten_hits():
            check_names |= set(hit.check_names)
        for seq in self.flatten_gRNAseqs():
            check_names |= set(hit.check_names)
    def parse_from_mapping(self, fname, targets = None) -> None:
        """
        Read gRNA data from MINORg .map file.
        
        Arguments:
            fname (str): required, path to file
            targets (str): optional, path to file containing target sequences.
                Used to get target sequence length for tie breaking by favouring gRNA hits closer to 5' end.
        """
        ## read data
        seq_targets = ({} if not targets
                       else {seqid: Target(seq, id = seqid, strand = '+')
                             for seqid, seq in fasta_to_dict(targets).items()})
        # print("parse_from_mapping, seq_targets:", seq_targets)
        raw_mapping = [line.split('\t') for line in splitlines(fname)]
        header_mapping = raw_mapping[0]
        ## determine version where checks start
        version = (1 if "exclusive" in header_mapping[5] else \
                   3 if "target length" in header_mapping else \
                   2 if "group" in header_mapping else 4)
        group_set_colname = "set" if version == 4 else "group"
        check_col = header_mapping.index(group_set_colname) + 1 ## group is last column before any checks
        header_checks = header_mapping[check_col:]
        dat_mapping = raw_mapping[1:]
        del raw_mapping
        for i, entry in enumerate(dat_mapping):
            ## note: gRNA_range is relative to + strand (not necessarily sense) when read from ..targets.txt file
            ## note: gRNA_id is probably meaningless as it can be overwritten using ids in the fasta file supplied to the variable 'fasta' in get_minimum_set_from_file
            if version == 1:
                gRNA_id, gRNA_seq, target_id, sense, strand, gRNA_range, group = entry[:check_col]
                gRNA_start, gRNA_end = map(int, re.search("\[(\d+), (\d+)\)", gRNA_range).group(1,2))
                target_len = 0
            elif version == 2 or version == 4:
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
                target = seq_targets.get(target_id, Target('N' * int(target_len), id = target_id, strand = '+'))
            except Exception as e:
                print(version, i, header_mapping, entry)
                raise e
            ## create gRNAHit object
            gRNA_hit = gRNAHit(target, gRNA_start, gRNA_end, strand, gRNA_id)
            gRNA_hit.set_parent_sense(sense)
            ## add gRNAHit object to gRNAHits object
            self.add_hit(gRNA_seq, gRNA_hit)
            self.get_gRNAseq_by_seq(gRNA_seq).id = gRNA_id
            ## log checks
            gRNA_checks = entry[check_col:]
            for i, check in enumerate(header_checks):
                ## these checks apply to all hits with the same gRNA seq so set check to gRNASeq obj
                if check in {"background", "exclude", "GC"}:
                    self.get_gRNAseq_by_seq(gRNA_seq).set_check(check, gRNA_checks[i])
                ## set all other checks by hit, including standard "feature" check
                ##  checks in columns added by user are also set by hit since we don't know anything
                else:
                    gRNA_hit.set_check(check, gRNA_checks[i])
        return
    #################
    ##  MODIFIERS  ##
    #################
    def assign_gRNAseq_id(self, fasta) -> str:
        """
        (Re)name gRNA according to FASTA file.
        
        Arguments:
            fasta (str): required, path to FASTA file
        """
        fasta_inv = {str(v): k for k, v in fasta_to_dict(fasta).items()}
        valid_seqs = set(fasta_inv.keys()) & set(self.seqs)
        if len(valid_seqs) < len(set(self.seqs)):
            print("\nWARNING: The provided FASTA file does not cover all gRNAs.\n")
        for seq in valid_seqs:
            self.get_gRNAseq_by_seq(seq).id = fasta_inv[seq]
        return
    def add_hit(self, seq, gRNA_hit):
        """
        Add gRNA hit.
        
        If entry for the gRNA's sequence does not exist in self._gRNAseqs,
        call :meth:`~minorg.grna.gRNAHits.add_seq` to add it too.
        
        Arguments:
            seq (str or Bio.Seq.Seq): gRNA sequence
            gRNA_hit (:class:`~minorg.grna.gRNAHit`): gRNAHit object
        """
        self.add_seq(seq)
        self._hits[str(seq)] = self.get_hits(seq) + [gRNA_hit]
    def add_seq(self, seq) -> None:
        """
        Add gRNA sequence if it doesn't already in self._gRNAseqs 
        and create empty entry for it in self._hits.
        
        Arguments:
            seq (str or Bio.Seq.Seq): gRNA sequence
        """
        if not str(seq) in self.gRNAseqs: self._gRNAseqs[str(seq)] = gRNASeq(seq)
        if not str(seq) in self.hits: self._hits[str(seq)] = []
    def remove_seqs(self, *seqs) -> None:
        """
        Remove gRNA sequence and associated hits.
        
        Arguments:
            *seqs (str or Bio.Seq.Seq): gRNA sequence to remove
        """
        if seqs and type(seqs[0]) in (list, tuple, set):
            seqs = list(itertools.chain(*seqs))
        for seq in seqs:
            if str(seq) in self.gRNAseqs: del self._gRNAseqs[str(seq)]
            if str(seq) in self.hits: del self._hits[str(seq)]
        return
    ###############
    ##  SETTERS  ##
    ###############
    def clear_checks(self) -> None:
        """
        Clear checks for both gRNA sequences and hits.
        """
        for gRNA_seq in self.flatten_gRNAseqs():
            gRNA_seq.clear_checks()
        for hit in self.flatten_hits():
            hit.clear_checks()
        return
    def set_seqs_check(self, check_name, status, seqs) -> None:
        """
        Set a given check to a given value for multiple gRNA sequences.
        
        Arguments:
            check_name (str): check name
            status (str or bool or None): check status
            seqs (list of str): list of gRNA sequences for which to set the check
        """
        if type(seqs) not in (tuple, list, set):
            seqs = (seqs,)
        for seq in seqs:
            self.get_gRNAseq_by_seq(seq).set_check(check_name, status)
        return
    def set_seqs_check_by_function(self, check_name, func, seqs) -> None:
        """
        Set a given check for multiple gRNA sequences using a function that accepts str of gRNA sequence.
        
        Arguments:
            check_name (str): check name
            func (func): function that accepts str of gRNA sequence
            seqs (list of str): list of gRNA sequences for which to set the check
        """
        for seq in seqs:
            self.set_seqs_check(check_name, func(self.get_gRNAseq_by_seq(seq)), [seq])
        return
    def set_all_seqs_check_by_function(self, check_name, func) -> None:
        """
        Set a given check for all gRNA sequences using a function that accepts str of gRNA sequence.
        
        Arguments:
            check_name (str): check name
            func (func): function that accepts str of gRNA sequence
        """
        self.set_seqs_check_by_function(check_name, func, self.seqs)
    def rename_seqs(self, fasta) -> None:
        """
        Rename gRNA according to FASTA file.
        Functionally identical to :meth:`~minorg.grna.gRNAHits.assign_gRNAseq_id`
        except it doesn't check whether the FASTA file covers all gRNA sequences.
        
        Arguments:
            fasta (str): required, path to FASTA file
        """
        seqs_names = {str(v): k for k, v in fasta_to_dict(fasta).items()}
        for seq, name in seqs_names.items():
            gRNA_seq = self.get_gRNAseq_by_seq(seq)
            if gRNA_seq:
                gRNA_seq.id = name
        return
    def assign_seqid(self, prefix = "gRNA_", zfill = 3, assign_all = True) -> None:
        """
        Assign sequence ID to gRNA sequences using format <prefix><unique gRNA number>.
        
        Arguments:
            prefix (str): prefix for gRNA ID (default='gRNA\\_')
            zfill (int): number of leading zeroes (default=3)
            assign_all (bool): assign new sequence IDs to all gRNA sequences
                regardless of whether they already have a sequence ID
        """
        for i, gRNA_seq in enumerate(self.flatten_gRNAseqs()):
            if gRNA_seq.id is None or assign_all:
                gRNA_seq.id = f"{prefix}{str(i+1).zfill(zfill)}"
        return
    ###############
    ##  GETTERS  ##
    ###############
    def flatten_hits(self) -> list:
        """
        Get list of gRNAHit objects.
        """
        return list(itertools.chain(*self.hits.values()))
    def flatten_gRNAseqs(self) -> list:
        """
        Get list of gRNASeq objects.
        """
        return list(self.gRNAseqs.values())
    def get_hits(self, seq) -> list:
        """
        Get hits of gRNA with given sequence.
        
        Arguments:
            seq (str or Bio.Seq.Seq): gRNA sequence
        
        Returns
        -------
        list
            Of gRNAHit objects
        """
        return self.hits.get(str(seq).upper(), []) ## 'seq' must be able to be coerced using str()
    def get_gRNAseq_by_seq(self, seq) -> Union['gRNASeq', None]:
        """
        Get gRNASeq object with given sequence.
        
        Arguments:
            seq (str or Bio.Seq.Seq): gRNA sequence
        
        Returns
        -------
        gRNASeq
            If exists
        None
            If doesn't exist
        """
        return self.gRNAseqs.get(str(seq).upper(), None) ## retrieve gRNASeq obj by seq
    def get_gRNAseq_by_id(self, id) -> Union['gRNASeq', None]:
        """
        Get gRNASeq object with given sequence ID.
        
        Arguments:
            seq (str): gRNA sequence ID
        
        Returns
        -------
        gRNASeq
            If exists
        None
            If doesn't exist
        """
        output = [gRNA_seq for gRNA_seq in self.flatten_gRNAseqs() if gRNA_seq.id == id]
        if not output: return None
        elif len(output) == 1: return output[0]
        else:
            print("\nWARNING: There are multiple gRNA sequences with the requested ID. Returning first in list.")
            return output[0]
    def get_gRNAseqs_by_seq(self, *seqs) -> list:
        """
        Get multiple gRNASeq objects by sequence.
        Sequences without associated gRNASeq objects will be skipped.
        
        Arguments:
            *seqs (str or Bio.Seq.Seq): gRNA sequence(s)
        
        Returns
        -------
        list
            Of gRNASeq objects
        """
        if not seqs:
            return []
        # if type(seqs[0]) in (list, tuple, set):
        #     seqs = tuple(itertools.chain(*seqs))
        return [self.get_gRNAseq_by_seq(seq) for seq in seqs if self.get_gRNAseq_by_seq(seq) != None]
    def get_gRNAseqs_by_id(self, *ids) -> list:
        """
        Get multiple gRNASeq objects by sequence ID.
        Sequence IDs without associated gRNASeq objects will be skipped.
        
        Arguments:
            *ids (str): gRNA sequence ID(s)
        
        Returns
        -------
        list
            Of gRNASeq objects
        """
        if type(ids[0]) in (list, tuple, set):
            ids = tuple(itertools.chain(*ids))
        return [self.get_gRNAseq_by_id(id) for id in ids if self.get_gRNAseq_by_id(seq) != None]
    ######################
    ##      FILTER      ##
    ##  (& return new)  ##
    ######################
    ## filter gRNAHit objects for certain criteria and return new gRNAHits object
    def filter_hits(self, *check_names, exclude_empty_seqs = True, accept_invalid = False,
                    accept_invalid_field = True, all_checks = False, quiet = False,
                    report_invalid_field = False):
        """
        Filter gRNA hits by checks and return new gRNAHits object.
        
        Arguments:
            *check_names (str): check name(s) (not required if ``all_checks=True``)
            exclude_empty_seqs (bool): exclude gRNA sequences from new gRNAHits object
                if none of their hits pass the filter(s)
            accept_invalid (bool): score unset checks as pass
            accept_invalid_field (bool): score unset checks as pass if a given check is not set for ALL hits
            all_checks (bool): filter using all checks
            quiet (bool): print only essential messages
        
        Returns
        -------
        :class:`~minorg.grna.gRNAHits`
        """
        if check_names:
            ## remove checks that don't exist
            check_names = [check_name for check_name in check_names if self.hit_check_exists(check_name)]
        elif all_checks:
            ## get all possible check names
            check_names = self.check_names_hits
            # check_names = list(set(itertools.chain(*[hit.check_names for hit in self.flatten_hits()])))
        else:
            raise Exception("Either *check_names OR all_checks is required.")
        ## warn about any invalid fields
        if ( (report_invalid_field or not quiet) and
             (not all(self.valid_hit_check(check_name) for check_name in check_names)) ):
            warnings.warn( ("The following hit check(s) have not been set: " +
                            ','.join([check_name for check_name in check_names
                                      if not self.valid_hit_check(check_name)])),
                           MINORgWarning)
        ## remove any invalid fields
        if accept_invalid_field:
            check_names = [check_name for check_name in check_names if self.valid_hit_check(check_name)]
        if not check_names:
            if not quiet:
                print("No valid hit check names remaining. Returning new gRNAHits object with all hits.")
            filtered_hits = {seq: hits for seq, hits in self.hits.items()}
        else:
            filtered_hits = {seq: [hit for hit in hits
                                   if ( (accept_invalid and
                                         hit.some_valid_checks_passed(*check_names))
                                        or ( (not accept_invalid) and
                                             hit.some_checks_passed(*check_names)) )]
                             for seq, hits in self.hits.items()}
            if exclude_empty_seqs:
                filtered_hits = {seq: hits for seq, hits in filtered_hits.items() if hits}
        filtered_seqs = {seq: self.get_gRNAseq_by_seq(seq) for seq in filtered_hits.keys()}
        output = gRNAHits(gRNA_seqs = filtered_seqs, gRNA_hits = filtered_hits)
        return output
    def filter_hits_some_checks_passed(self, *check_names, **kwargs):
        """
        Wrapper for :meth:`~minorg.grna.gRNAHits.filter_hits`.
        Filter hits by some checks.
        
        Arguments:
            *check_names (str): required, check_names
            **kwargs: other arguments passed to :meth:`~minorg.grna.gRNAHits.filter_hits`
        """
        if not check_names:
            raise Exception("At least one check name is required.")
        return self.filter_hits(*check_names, **kwargs)
    def filter_hits_all_checks_passed(self, **kwargs):
        """
        Wrapper for :meth:`~minorg.grna.gRNAHits.filter_hits`.
        Filter hits by ALL checks.
        
        Arguments:
            **kwargs: other arguments passed to :meth:`~minorg.grna.gRNAHits.filter_hits`
        """
        return self.filter_hits(all_checks = True, **kwargs)
    ## filter gRNASeq objects for certain criteria and return new gRNAHits object
    def filter_seqs(self, *check_names, accept_invalid = True, accept_invalid_field = True,
                    all_checks = False, quiet = False, report_invalid_field = False):
        """
        Filter gRNA sequences by checks and return new gRNAHits object.
        
        Arguments:
            *check_names (str): check name(s) (not required if ``all_checks=True``)
            accept_invalid (bool): score unset checks as pass
            accept_invalid_field (bool): score unset checks as pass if a given check is not set for ALL hits
            all_checks (bool): filter using all checks
            quiet (bool): print only essential messages
        
        Returns
        -------
        :class:`~minorg.grna.gRNAHits`
        """
        if check_names:
            ## remove checks that don't exist
            check_names = [check_name for check_name in check_names if self.seq_check_exists(check_name)]
        elif all_checks:
            ## get all possible check names
            check_names = list(set(itertools.chain(*[seq.check_names for seq in self.flatten_gRNAseqs()])))
        else:
            raise Exception("Either *check_names OR all_checks is required.")
        ## warn about any invalid fields
        if ( (report_invalid_field or not quiet) and
             (not all(self.valid_seq_check(check_name) for check_name in check_names)) ):
            warnings.warn( ("The following seq check(s) have not been set: " +
                            ','.join([check_name for check_name in check_names
                                      if not self.valid_seq_check(check_name)])),
                           MINORgWarning)
        ## remove any invalid fields
        if accept_invalid_field:
            check_names = [check_name for check_name in check_names if self.valid_seq_check(check_name)]
        if not check_names:
            if not quiet:
                print("No valid seq check names remaining. Returning new gRNAHits object with all sequences.")
            filtered_hits = {seq: hits for seq, hits in self.hits.items()}
        else:
            filtered_hits = {seq: hits for seq, hits in self.hits.items()
                             if ((accept_invalid
                                  and self.get_gRNAseq_by_seq(seq).some_valid_checks_passed(*check_names))
                                 or (not accept_invalid and
                                     self.get_gRNAseq_by_seq(seq).some_checks_passed(*check_names)))}
        filtered_seqs = {seq: self.get_gRNAseq_by_seq(seq) for seq in filtered_hits.keys()}
        output = gRNAHits(gRNA_seqs = filtered_seqs, gRNA_hits = filtered_hits)
        return output
    def filter_seqs_some_checks_passed(self, *check_names, **kwargs):
        """
        Wrapper for :meth:`~minorg.grna.gRNAHits.filter_seqs`.
        Filter gRNA sequences by some checks.
        
        Arguments:
            *check_names (str): required, check_names
            **kwargs: other arguments passed to :meth:`~minorg.grna.gRNAHits.filter_seqs`
        """
        if not check_names:
            raise Exception("At least one check name is required")
        return self.filter_seqs(*check_names, **kwargs)
    def filter_seqs_all_checks_passed(self, **kwargs):
        """
        Wrapper for :meth:`~minorg.grna.gRNAHits.filter_seqs`.
        Filter hits by ALL checks.
        
        Arguments:
            **kwargs: other arguments passed to :meth:`~minorg.grna.gRNAHits.filter_seqs`
        """
        return self.filter_seqs(all_checks = True, **kwargs)
    #############
    ##  WRITE  ##
    #############
    def write_mapping(self, fout, sets = [], write_all = False,
                      write_checks = False, checks = ["background", "GC", "feature"],
                      index = 1, start_incl = True, end_incl = True, version = 4,
                      fasta = None ) -> None: ## if fasta is provided, it will override default naming behaviour
        """
        Write MINORg .map file for gRNA mapping.
        
        A MINORg .map file is tab-delimited and includes a header line.
        In version 4, columns are:
        
            - gRNA id: gRNA sequence ID
            - gRNA sequence: gRNA sequence
            - target id: sequence ID of gRNA target
            - target sense: sense of gRNA target ('sense' or 'antisense')
            - gRNA strand: strand of gRNA (relative to target) ('+' or '-')
            - start: position of start of gRNA in target
            - end: position of end of gRNA in target
            - set: gRNA set number (set to 1 for all gRNA if ``write_all=True``)
        
        Any columns after 'set' are check statuses.
        
        Arguments:
            fout (str): required, path to output file
            sets (list): list of grouped gRNA sequences (e.g. [('<seq1>', '<seq2>'), ('<seq3>',)]).
                If provided, only gRNA sequences included in ``sets`` will be written.
                All sets must be mutually exclusive.
                Used to assign set numbers in 'set' column.
                Overrides ``write_all``.
            write_all (bool): write all gRNA sequences
            write_checks (bool): write check statuses
            checks (list): check names of checks to write (default=['background', 'GC', 'feature'])
            index (int): index of written gRNA range (default=1)
            start_incl (bool): whether written gRNA range is start-inclusive (default=True)
            end_incl (bool): whether written gRNA range is end-inclusive (default=True)
            fasta (str): optional, path to FASTA file, used for renaming gRNA
            version (int): .map file version
        """
        if not write_checks: checks = []
        if sets: sets = sets
        elif write_all: sets = [set(self.seqs)]
        else:
            print("Either 'sets' OR 'write_all' is required. Writing empty file.")
            open(fout, "w+").write('')
            return
        fasta_inv = {} if not fasta else {str(seq): k for k, seq in fasta_to_dict(fasta).items()}
        if version == 1:
            mapping_header = ["gRNA id", "gRNA sequence", "target id", "target sense",
                              "gRNA strand", "range (0-index, end exclusive)", "group"] + checks
        elif version == 2:
            mapping_header = ["gRNA id", "gRNA sequence", "target id", "target sense",
                              "gRNA strand", "start", "end", "group"] + checks
        elif version == 3:
            mapping_header = ["gRNA id", "gRNA sequence", "target id", "target length", "target sense",
                              "gRNA strand", "start", "end", "group"] + checks
        else:
            mapping_header = ["gRNA id", "gRNA sequence", "target id", "target sense",
                              "gRNA strand", "start", "end", "set"] + checks
        mapping_dat = []
        seq_check_names = {"GC", "background"}
        for group, seqs in enumerate(sets):
            for seq in seqs:
                gRNA_seq = self.get_gRNAseq_by_seq(seq)
                gRNA_hits = self.get_hits(seq)
                for gRNA_hit in gRNA_hits:
                    ## put together all required fields
                    if version == 1:
                        mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id), ## override id if fasta provided
                                            gRNA_seq.seq, gRNA_hit.target_id,
                                            gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand,
                                            "[{}, {})".format(*gRNA_hit.range), group + 1] +
                                           [(gRNA_seq.check(check_name, mode = "str")
                                             if check_name in seq_check_names else
                                             gRNA_hit.check(check_name, mode = "str"))
                                            for check_name in checks])
                    elif version == 2 or version == 4:
                        ## figure out start and end
                        start, end = map(lambda n: n + index, gRNA_hit.range)
                        if not start_incl: start -= 1
                        if end_incl: end -= 1
                        mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id), ## override id if fasta provided
                                            gRNA_seq.seq, gRNA_hit.target_id,
                                            gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand,
                                            start, end, group + 1] +
                                           [(gRNA_seq.check(check_name, mode = "str")
                                             if check_name in seq_check_names else
                                             gRNA_hit.check(check_name, mode = "str"))
                                            for check_name in checks])
                    else:
                        ## figure out start and end
                        start, end = map(lambda n: n + index, gRNA_hit.range)
                        if not start_incl: start -= 1
                        if end_incl: end -= 1
                        mapping_dat.append([fasta_inv.get(seq, gRNA_seq.id), ## override id if fasta provided
                                            gRNA_seq.seq, gRNA_hit.target_id, gRNA_hit.target_len,
                                            gRNA_hit.parent_sense(mode = "str"), gRNA_hit.strand,
                                            start, end, group + 1] +
                                           [(gRNA_seq.check(check_name, mode = "str")
                                             if check_name in seq_check_names else
                                             gRNA_hit.check(check_name, mode = "str"))
                                            for check_name in checks])
        mapping_dat.sort(key = lambda entry: int(entry[mapping_header.index("start")])) ## sort by start
        mapping_dat.sort(key = lambda entry: entry[mapping_header.index("target id")]) ## sort by target id
        mapping_dat.sort(key = lambda entry: entry[mapping_header.index("gRNA id")]) ## sort by gRNA id
        mapping_dat.sort(key = lambda entry: entry[mapping_header.index("set" if version == 4
                                                                        else "group")])
        to_write = '\n'.join(['\t'.join(map(str, entry)) for entry in ([mapping_header] + mapping_dat)]) + '\n'
        open(fout, "w+").write(to_write)
        return
    def write_fasta(self, fout, ids = [], seqs = [], write_all = False, fasta = None) -> None:
        """
        Write gRNA sequences to FASTA file.
        
        Arguments:
            fout (str): required, path to output file
            ids (list): list of IDs (str) of gRNA to write
            seqs (list): list of sequences (str) of gRNA to write, overrides ``ids``
            write_all (bool): write all gRNA seqences, overrides ``seqs`` and ``ids``
            fasta (str): optional, path to FASTA file, used for renaming gRNA
        """
        ## get relevant gRNA sequences
        if write_all: gRNA_seqs = self.flatten_gRNAseqs()
        elif seqs: gRNA_seqs = self.get_gRNAseqs_by_seq(*seqs)
        elif ids: gRNA_seqs = self.get_gRNAseqs_by_id(*ids)
        else:
            print("Either 'ids' OR 'seqs' OR 'write_all' is required. Writing empty file.")
            open(fout, "w+").write('')
            return
        ## rename sequences per fasta file (if fasta file provided)
        fasta_inv = {} if not fasta else {str(seq): k for k, seq in fasta_to_dict(fasta).items()}
        self.assign_seqid(assign_all = False)
        to_write = {fasta_inv.get(gRNA_seq.seq, gRNA_seq.id): gRNA_seq.seq for gRNA_seq in gRNA_seqs}
        if to_write:
            dict_to_fasta(to_write, fout)
        else:
            open(fout, "w+").write('')
        return

class gRNAHit(CheckObj):
    """
    gRNA hit
    
    Attributes:
        _target (:class:`~minorg.grna.Target`): Target object of sequence that gRNA targets
        _range (tuple): position of gRNA in target
        _strand (str): strand of gRNA in target
        _hit_id: unique hit identifier
    """
    def __init__(self, target, start, end, strand, hit_id):
        """
        Create a gRNAHit object.
        
        Arguments:
            target (:class:`~minorg.grna.Target`): Target object of sequence that gRNA targets
            start (int): start position of gRNA in target (relative to parent sequence of target)
            end (int): end position of gRNA in target (relative to parent sequence of target)
            hit_id: unique hit identifier
        """
        ## note: unless something is weird, seq_strand is the same as strand
        # super().__init__("background", "GC", "feature", "exclude", "flank")
        super().__init__("background", "GC", "feature", "exclude")
        self._target = target ## previously self._seq
        self._range = (start, end) ## relative to original (parent) sequence from which _target was derived
        self._strand = '+' if (strand.lower() == "fwd" or strand == '+') else '-' ## gRNA strand relative to original (parent) sequence from which _target was derived (gRNA strand should be same as _target strand)
        # self._seq_strand = '' # stores _target direction relative to original sequence from which _target was derived
        # self._parent_sense = None ## - if original (parent) sequence from which _target was derived is on antisense strand, + if _target is on sense strand
        self._hit_id = hit_id
    @property
    def target(self) -> 'Target': return self._target
    @property
    def start(self) -> int: return self._range[0]
    @property
    def end(self) -> int: return self._range[1]
    @property
    def target_id(self) -> str: return self.target.id
    @property
    def strand(self) -> str: return self._strand
    @property
    def hit_id(self): return self._hit_id
    @property
    def range(self) -> tuple: return tuple(self._range)
    @property
    def reverse_range(self) -> tuple: return (self.target_len - self.end, self.target_len - self.start)
    @property
    def target_len(self) -> int: return len(self.target)
    @property
    def target_strand(self) -> str: return self.target.strand
    # def set_target_strand(self, strand): self._target_strand = strand
    def set_parent_sense(self, strand) -> None:
        self.target.set_sense_by_parent(strand)
        # self._parent_sense = ('+' if strand in ('+', "sense") else '-')
    def set_bg_check(self, status): self.set_check("background", status)
    def set_gc_check(self, status): self.set_check("GC", status)
    def set_feature_check(self, status): self.set_check("feature", status)
    def set_exclude_check(self, status): self.set_check("exclude", status)
    # def set_flank_check(self, status): self.set_check("flank", status)
    def parent_sense(self, mode = "raw"):
        parent_sense_value = self.target.sense
        if mode == "str":
            if parent_sense_value == '+':
                return "sense"
            elif parent_sense_value == '-':
                return "antisense"
            else:
                return "NA"
        else:
            return parent_sense_value
    def adj_range(self, mode = "strand") -> tuple:
        """
        Return range based on strand, target, or gene.
        
        ``mode='strand'``:
        - Return range relative to plus strand. Uses self.strand.
        ``mode='target'``:
        - Return range relative to target sequence. Uses self.strand == self.target_strand.
        ``mode='gene'``:
        - Return range relative to sense strand. Uses self.target.sense == self.strand.
        
        Arguments:
            mode (str): reference for range. Valid values: 'strand', 'target', 'gene'
        
        Returns
        -------
        tuple
            Of (start, end)
        tuple
            Of (float('nan'), float('nan')) if the reference cannot be determined.
        """
        if (not self.strand or
            (mode == "target" and not self.target_strand) or
            (mode == "sense" and not self.parent_sense())):
            return (float("nan"), float("nan"))
        elif ((mode == "strand" and (self.strand == '+')) or \
              (mode == "target" and (self.target_strand == self.strand)) or \
              (mode == "sense" and (self.target.sense == self.strand))):
            return self.range
        return self.reverse_range
    def flank(self, length = 100):
        target_seq = self.target if self.strand == '+' else self.target.reverse_complement()
        start, end  = self.adj_range(mode = "target")
        return target_seq[start - length: start], target_seq[end: end + length]
