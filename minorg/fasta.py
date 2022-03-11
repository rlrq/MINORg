import tempfile
import itertools

from minorg.index import IndexedFasta
from minorg.searchio import searchio_parse

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

def dict_to_SeqRecordList(d, description = '', seq_id_func = lambda x:x,
                          seq_type = "DNA", gap_char = '-', gapped = False):
    out_l = []
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    for seq_id,seq in d.items():
        out_l.append(SeqRecord(seq if isinstance(seq, Seq) else \
                               Seq(str(seq)),
                               id = seq_id_func(seq_id), description = description))
    return out_l

def dict_to_fasta(d, fout, seq_type = "detect", gap_char = '-', gapped = False):
    from Bio import SeqIO
    SeqIO.write(dict_to_SeqRecordList(d, gap_char = gap_char, gapped = gapped),
                fout, "fasta")

def extract_ranges(seq, ranges, strand = '+'):
    ranges_sorted = sorted(ranges, key = lambda x: int(x[0]), reverse = (strand == '-'))
    output = seq[:0]
    for start, end in ranges_sorted:
        output += seq[int(start):int(end)] if strand == '+' else \
                  seq[int(end)-1:int(start)-1:-1]
    return output

def find_identical_in_fasta(query, subject, chunk = 100000):
    """
    Searches for exact matches in memory saving way using 
    sliding within with size ``2*(max(chunk, len(query_seq)))``
    and overlap ``max(chunk, len(query_seq))``.
    DOES NOT EXPAND AMBIGUOUS BASES. (i.e. 'N' matches ONLY 'N' character, not any base)
    
    Arguments:
        query (dict): {'<seqid>': '<str of seq or Bio.Seq.Seq obj>'}
        subject (str): path to FASTA file
        chunk (int): window overlap size (bp) for search
    
    Returns:
        list: list of SearchIO QueryResult objects as if read from blast-tab file
              with fields "qseqid sseqid sstart send"
    """
    isubject = IndexedFasta(subject)
    output = []
    def find_overlapping(query, subject, query_id):
        overlap = max(len(query), chunk)
        query_rvs = query.reverse_complement()
        for i in range(0, len(subject), overlap):
            subject_window = subject[i:(i+(2*overlap))]
            ## fwd query
            pos = subject_window.find(query)
            while pos >= 0:
                ## BLAST tab fmt: 1-index, end inclusive
                output.append((query_id, subject.name, i + pos + 1, i + pos + len(query)))
                pos = subject_window.find(query, start = pos + 1)
            ## rvs query
            pos = subject_window.find(query_rvs)
            while pos >= 0:
                ## BLAST tab fmt: 1-index, end inclusive
                output.append((query_id, subject.name, i + pos + 1, i + pos + len(query)))
                pos = subject_window.find(query_rvs, start = pos + 1)
        return
    from Bio import Seq
    for subject_seq in isubject.values():
        for query_id, seq in query.items():
            query_seq = (seq if isinstance(seq, Seq.Seq) else Seq.Seq(seq))
            find_overlapping(query_seq, subject_seq, query_id)
    ## write to file and parse with SearchIO into generator
    import tempfile
    tmp_f = tempfile.mkstemp(suffix = ".tsv")[1]
    with open(tmp_f, 'w') as f:
        for entry in output:
            f.write('\t'.join(map(str, entry)) + '\n')
    searchio = list(searchio_parse(tmp_f, "blast-tab", fields = "qseqid sseqid sstart send"))
    import os
    os.remove(tmp_f)
    return searchio

## collapses identical sequences (arbitrarily selects a seqid to represent each set)
## writes collapsed fasta file and a tsv file mapping identical sequences
def collapse_identical_seqs(fasta, fout_fa, fout_seqid):
    dat = fasta_to_dict(fasta)
    identicals = {k: set(seqid for seqid, seq in dat.items() if str(seq) == str(v))
                  for k, v in dat.items()}
    identical_sets = set(map(lambda x: tuple(sorted(x)), identicals.values()))
    dict_to_fasta({seqids[0]: dat[seqids[0]] for seqids in identical_sets}, fout_fa)
    open(fout_seqid, "w+").write('\n'.join(['\t'.join(seqids) for seqids in identical_sets]))
    return
