import regex as re

from minorg.fasta import fasta_to_dict, dict_to_fasta
from minorg.functions import prepend_line, parse_get_data, splitlines, write_tsv
from minorg.searchio import searchio_parse
from minorg.ot_regex import (
    CHAR_UNALIGNED,
    CHAR_DEL,
    CHAR_INS,
    CHAR_MISMATCH,
    CHAR_MATCH,
    CHAR_GAP,
    BTOP_GAP
)

## TODO: test BlastHSP

########################
##  BLAST FORMATTING  ##
########################

class BlastHSP:
    """
    Class that binds HSP object with QueryResult and Hit as query and subject respectively.
    """
    import copy
    def __init__(self, hsp, query, subject):
        self.hsp = hsp
        self.query = query
        self.subject = subject
        self.qlen = None if "seq_len" not in dir(self.query) else self.query.seq_len
        self.sbtop = None
        self.tbtop = None
        self.tbtop_rvs = None
        self.parse_btop()
    
    @property
    def btop(self):
        return self.hsp.btop
    
    def parse_btop(self):
        if "btop" in dir(self.hsp):
            self.sbtop = self.expand_btop_str(include_unaligned = True)
            self.tbtop = self.expand_btop_tuple(include_unaligned = True)
            self.tbtop_rvs = self.expand_btop_tuple(include_unaligned = True, rvs_index = True)
    
    def expand_btop_str(self, include_unaligned = False):
        """
        Expands btop to a string of characters of length equal to alignment, where:
            '.' is a match,
            'm' is a mismatch,
            'i' is an insertion in the query,
            'd' is a deletion in the query.
        
        Arguments:
            include_unaligned (bool): add space character for each unaligned position at 3' and 5' ends
        
        Returns
        -------
        str
            Expanded btop pattern
        """
        output = ''
        btop = self.btop
        ## processed characters are removed from btop until none are left
        while btop:
            num = re.search('^\d+', btop)
            if num: ## if match
                num = num.group(0)
                output += CHAR_MATCH*int(num)
                btop = btop[len(num):]
            else: ## if gap or mismatch
                unit = btop[:2]
                if unit[0] == BTOP_GAP: output += CHAR_DEL ## deletion in query
                elif unit[1] == BTOP_GAP: output += CHAR_INS ## insertion in query
                else: output += CHAR_MISMATCH ## mismatch
                btop = btop[2:]
        if include_unaligned:
            output = (CHAR_UNALIGNED*self.hsp.query_start) + \
                     output + \
                     (CHAR_UNALIGNED*(self.qlen-self.hsp.query_end))
        return output
    
    def expand_btop_tuple(self, rvs_index = False, include_unaligned = False):
        """
        Calls str_expand_btop and groups expanded btop to a tuple of characters
        of length equal to aligned query, where:
            '.' is a match,
            'm' is a mismatch,
            'i' is an insertion in the query,
            'd' is a deletion in the query (deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N).

        Arguments:
            rvs_index (bool): default=False. A deletion between positions N and N+1 will be grouped with:
                position N+1 if rvs_index=False;
                position N, if rvs_index=True.

        Returns
        -------
        tuple of str
            Expanded btop pattern
        """
        str_expanded = self.expand_btop_str()
        output = []
        non_del = {CHAR_MATCH, CHAR_MISMATCH, CHAR_INS}
        ## processed characters are removed from str_expanded until none are left
        while str_expanded:
            c = str_expanded[0]
            if c in non_del: ## if not deletion
                output.append(c)
                str_expanded = str_expanded[1:]
            elif rvs_index: ## if deletion in reversed index
                output[-1] = output[-1] + c
                str_expanded = str_expanded[1:]
            else: ## if deletion in index counting up
                output.append(str_expanded[:2])
                str_expanded = str_expanded[2:]
        if include_unaligned:
            output = ([CHAR_UNALIGNED]*self.hsp.query_start) + output + \
                     ([CHAR_UNALIGNED]*(self.qlen-self.hsp.query_end))
        return tuple(output)
    
    def splice_expanded_btop(self, start = None, end = None, fmt = tuple,
                             rvs_index = False, include_unaligned = True):
        """
        Returns expanded btop spliced to specified range and in specified format
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
            fmt (type): str or tuple; return expanded btop type
            rvs_index (bool): used when fmt=tuple;
                deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N
            include unaligned (bool): add space character for each unaligned position at 3' and 5' ends
        
        Returns
        -------
        str or tuple
            Expanded btop spliced to specified range
        """
        if start is None: start = 0
        if end is None: end = self.qlen
        ## originally: incr (int): default=1; splice increment (use -1 when start < end; use 1 when start > end)
        incr = 1 if (end >= start) else -1
        btop = self.tbtop_rvs if rvs_index else self.tbtop
        btop = btop[start:end:incr]
        if not include_unaligned:
            btop = type(btop)(x for x in btop if x != CHAR_UNALIGNED)
        if fmt == str: return ''.join(btop)
        else: return btop
        
    def _num_char_in_range(self, *char, start = None, end = None, rvs_index = False):
        """
        Returns number of times a given set of characters appears within a given range
        
        Arguments:
            *char (str): characters to search for
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
            rvs_index (bool): used when fmt=tuple;
                deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N
        
        Returns
        -------
        int
            Number of times characters appear within range
        """
        sbtop = self.splice_expanded_btop(start, end, fmt = str, rvs_index = rvs_index,
                                          include_unaligned = True)
        output = 0
        for c in char:
            output += sbtop.count(c)
        return output
    
    def num_insertion(self, start = None, end = None):
        """
        Returns number of insertions relative to subject
        (that is, one or more bases that is/are present in the subject but not in the query).
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
        
        Returns
        -------
        int
            Number of insertions within range
        """
        return self._num_char_in_range(CHAR_INS, start = start, end = end)
    
    def num_deletion(self, start = None, end = None, rvs_index = False):
        """
        Returns number of deletions relative to subject
        (that is, one or more bases that is/are present in the query but not in the subject).
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
            rvs_index (bool): used when fmt=tuple;
                deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N
        
        Returns
        -------
        int
            Number of deletions within range
        """
        return self._num_char_in_range(CHAR_DEL, start = start, end = end, rvs_index = rvs_index)
    
    def num_mismatch(self, start = None, end = None):
        """
        Returns number of mismatches.
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
        
        Returns
        -------
        int
            Number of mismatches within range
        """
        return self._num_char_in_range(CHAR_MISMATCH, start = start, end = end)
    
    def num_gap(self, start = None, end = None, rvs_index = False):
        """
        Returns number of gaps.
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
            rvs_index (bool): used when fmt=tuple;
                deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N
        
        Returns
        -------
        int
            Number of gaps within range
        """
        return self._num_char_in_range(CHAR_DEL, CHAR_INS, start = start, end = end, rvs_index = rvs_index)
    
class BlastResult:
    """
    Generator that reads blast-tab format and yields BlastHSP,
    which stores a HSP object with its associated QueryResult and Hit objects 
    as attributes query and subject respectively.
    """
    def __init__(self, filename, fmt, **kwargs):
        """
        Create a BlastResult object.
        
        Arguments:
            filename (str): path to file
            fmt (str): file format (e.g. 'blast-tab', 'blast-xml')
            **kwargs: additional arguments for SearchIO (notably, fields)
        """
        self.filename = filename
        self.fmt = fmt
        self.kwargs = kwargs
    def __iter__(self):
        for query_result in searchio_parse(self.filename, self.fmt, **self.kwargs):
            for hit in query_result:
                for hsp in hit:
                    yield BlastHSP(hsp, query_result, hit)

# class Test:
#     def __init__(self):
#         return
#     def __iter__(self):
#         for i in range(3):
#             for j in range(i):
#                 yield (i, j)

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

def blast(blastf, header, fout, append = False, dir = None, outfmt = 6, **kwargs):
    import os
    ## parse header
    if isinstance(header, str):
        header = header.split(',')
    blast_outfmt = f"{outfmt}{'' if (not header or outfmt not in [6,7,10]) else (' ' + ' '.join(header))}"
    ## execute blast
    if append and os.path.exists(fout):
        import tempfile
        tmpf = tempfile.mkstemp(dir = dir)[1]
        blast_cline = blastf(out = tmpf, outfmt = blast_outfmt,
                             **kwargs)
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
        blast_cline = blastf(out = fout, outfmt = blast_outfmt, **kwargs)
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
    with open(fout, 'w') as f:
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
    with open(fout, 'w') as f:
        ## we're not joining the entries w/ '\n' since we didn't strip them to begin with
        f.write(''.join(['\t'.join(header) + '\n'] + output))
    return

## overwrites original file
def filter_rpsblast_for_domain(fname, *pssmids, fout = None):
    pssmids = set(str(pssmid) for pssmid in pssmids)
    pssmid_str = '|'.join(sorted(pssmids))
    get, dat = parse_get_data(fname, delim = '\t')
    output = [list(x) + [pssmid_str] for x in dat if (get(x, "sseqid").split('|')[-1] in pssmids or
                                                      get(x, "sseqid").split(':')[-1] in pssmids)]
    if fout:
        write_tsv(fout, [get(get_cols = True) + ["domain"]] + output)
    else:
        ## return dictionary of hits {<qseqid>: [(hit1_start, hit1_end), (hit2_start, hit2_end)...]}
        output_dict = {}
        for entry in output:
            qseqid = get(entry, "qseqid")
            output_dict[qseqid] = output_dict.get(qseqid, []) + [tuple(get(entry, "qstart", "qend"))]
        return output_dict

class BlastNR:
    def __init__(self, query, blastf, header, fout, directory = None, keep_tmp = False, **kwargs):
        self.directory = directory
        self.query = query
        self.blastf = blastf ## e.g. NcbirpsblastCommandline
        self.header = header
        self.fout = fout
        self.kwargs = kwargs
        ## tmp files
        self.query_nr = None
        self.query_nr_map = None
        self.collapse_query()
        self.blast()
        self.expand()
        if not keep_tmp:
            import os
            os.remove(self.query_nr)
            os.remove(self.query_nr_map)
        return
    def collapse_query(self, fout_fasta = None, fout_map = None):
        import tempfile
        self.query_nr = fout_fasta if fout_fasta is not None else tempfile.mkstemp(dir = self.directory)[1]
        self.query_nr_map = fout_map if fout_map is not None else tempfile.mkstemp(dir = self.directory)[1]
        dat = fasta_to_dict(self.query)
        identicals = {k: set(seqid for seqid, seq in dat.items() if str(seq) == str(v))
                      for k, v in dat.items()}
        identical_sets = set(map(lambda x: tuple(sorted(x)), identicals.values()))
        ## write nr sequences
        dict_to_fasta({seqids[0]: dat[seqids[0]] for seqids in identical_sets}, self.query_nr)
        ## write nr mapping
        with open(self.query_nr_map, 'w') as f:
            f.write('\n'.join(['\t'.join(seqids) for seqids in identical_sets]))
        return
    def blast(self):
        header = self.header if "qseqid" in self.header.split(',') else f"qseqid,{header}"
        if self.kwargs.get("remote", False):
            incompatible = {"num_threads"}
            kwargs = {k: v for k, v in self.kwargs.items() if k not in incompatible}
        else:
            kwargs = self.kwargs
        # print("blast kwargs:", kwargs)
        blast6(blastf = self.blastf, header = header, fout = self.fout,
               query = self.query_nr, **kwargs)
    def expand(self, write = True):
        get, dat = parse_get_data(self.fout, delim = '\t')
        repr_map = [x.split('\t') for x in splitlines(self.query_nr_map)]
        ## expand entry for identical peptides previously collapsed
        qseqid_i = get(get_cols = True).index("qseqid")
        output = [line[:qseqid_i] + [seqid] + line[qseqid_i+1:]
                  for seqids in repr_map for seqid in seqids for line in dat \
                  if line[0] == seqids[0]]
        if write:
            write_tsv(self.fout, [get(get_cols = True)] + output)
        else:
            return output
