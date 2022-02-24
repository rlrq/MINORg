import os
import itertools
import regex as re

from minorg.functions import (
    splitlines,
    fasta_to_dict,
    dict_to_fasta
)

from minorg.grna import (
    Target,
    gRNAHits,
    gRNAHit
)

from minorg.pam import PAM

#################
##  FIND gRNA  ##
#################

def find_multi_gRNA(target_fname, pam = "NGG", gRNA_len = 20): #Let NGG be GG, or NG be G
    """
    Finds all possible target sequences given pam.
    Returns: {<seq>: set(ids of targeted targets)}
    """
    ## specifying 'N' allows users to specify how many wobble/N bases
    pam = PAM(pam = pam, gRNA_length = gRNA_len)
    pam_pattern = pam.regex()
    # pam_pattern = make_pam_pattern(pam, gRNA_len = gRNA_len)
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

