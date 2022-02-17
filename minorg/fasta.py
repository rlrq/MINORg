import tempfile
import itertools

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
    ## find all identical seqs (even overlapping ones)
    def find_overlapping(query, subject, query_id):
        ## Bio.Seq.Seq.find only returns first match. loop until no matches left.
        pos = subject.seq.find(query)
        while pos >= 0:
            output.append((query_id, subject.id, pos + 1, pos + len(query)))
            pos = subject.seq.find(query, start = pos + 1)
        return
    for subject_seq in subject_seqs:
        for query_id in query_fwd:
            find_overlapping(query_fwd[query_id], subject_seq, query_id)
            find_overlapping(query_rvs[query_id], subject_seq, query_id)
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
