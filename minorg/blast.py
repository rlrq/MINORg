from minorg.fasta import fasta_to_dict, dict_to_fasta
from minorg.functions import prepend_line, parse_get_data, splitlines, write_tsv
from minorg.searchio import searchio_parse

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

class BlastResult(list):
    """
    Read blast-tab format and flattens it into a list of BlastHSP,
    which stores a HSP object with it associated QueryResult and Hit objects 
    as attributes query and subject respectively.
    """
    def __init__(self, filename, fmt, **kwargs):
        """
        Create a BlastTab object.
        
        Arguments:
            filename (str): path to file
            fmt (str): file format (e.g. 'blast-tab', 'blast-xml')
            **kwargs: additional arguments for SearchIO (notably, fields)
        """
        self.filename = filename
        super().__init__()
        for query_result in searchio_parse(filename, fmt, **kwargs):
            for hit in query_result:
                for hsp in hit:
                    self.append(BlastHSP(hsp, query_result, hit))

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
    blast_outfmt = f"'{outfmt}{'' if (not header or outfmt not in [6,7,10]) else (' ' + ' '.join(header))}'"
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
    output = [list(x) + [pssmid_str] for x in dat if get(x, "sseqid").split('|')[-1] in pssmids]
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
        blast6(blastf = self.blastf, header = header, fout = self.fout, query = self.query_nr, **self.kwargs)
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
