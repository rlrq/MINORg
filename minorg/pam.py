###########
##  PAM  ##
###########

import regex as re
from minorg.constants import PAM_PATTERNS

DEFAULT_gRNA_LENGTH = 20
DEFAULT_PAM = "NGG"

class PAM():
    def __init__(self, pam = "NGG", gRNA_length = DEFAULT_gRNA_LENGTH):
        """
        pam: PAM object or str of PAM pattern (if PAM object provided, attributes are copied)
        gRNA_length: int length of gRNA (bp)
        """
        if isinstance(pam, PAM):
            self.raw_pam = pam.raw_pam
            self.pam = pam.pam
            self.gRNA_length = (1 if gRNA_length is None
                                else gRNA_length if gRNA_length != DEFAULT_gRNA_LENGTH
                                else pam.gRNA_length)
        else:
            if pam in PAM_PATTERNS:
                pam = PAM_PATTERNS[pam]
            self.raw_pam = pam
            self.pam = pam
            self.gRNA_length = gRNA_length
        self.parse()
    
    def __repr__(self):
        return self.pam
    
    ## previously: infer_full_pam
    def infer_full(self):
        ## use default 3' PAM + 1 base spacer if '.' and 'N' not provided
        pam = self.pam.upper()
        if '.' not in pam:
            ## assume 3' PAM if not indicated
            if 'N' not in pam: pam = 'N' + pam
            ## 3' PAM
            if pam[0] == 'N': pam = '.' + pam
            ## if 5' PAM
            else: pam = pam + '.'
        self.pam = pam
        return
    
    def expand_ambiguous(self):
        """
        Map ambiguous bases
        """
        from Bio.Data import IUPACData
        amb_dna = IUPACData.ambiguous_dna_values
        pam_mapped = ''
        for c in self.pam:
            mapped = amb_dna.get(c.upper(), c)
            if len(mapped) == 1: pam_mapped += mapped
            else: pam_mapped += f"[{mapped}]"
        self.pam = pam_mapped
        return
    
    def parse(self):
        ## infer pam location + spacer if not explicitly described
        self.infer_full()
        ## map ambiguous bases
        self.expand_ambiguous()
        return
    
    ## previously: make_pam_pattern
    def regex(self, gRNA_length: int = None):
        """
        Square brackets not allowed in PAM pattern unless they contain only A, T, G, C, or U.
        """
        if gRNA_length is None:
            gRNA_length = self.gRNA_length
        ## generate pattern for compilation
        grna_pre, grna_post = self.pam.split('.')
        pam_pattern = ''
        if grna_pre: pam_pattern += f"(?<={grna_pre})"
        pam_pattern += ".{" + str(gRNA_length) + '}'
        if grna_post: pam_pattern += f"(?={grna_post})"
        ## generate pattern for gRNA extraction
        print("PAM pattern:", pam_pattern)
        return re.compile(pam_pattern, re.IGNORECASE)
    
    ## function to get maximum pattern length
    def _pattern_max_len(self, pattern):
        ## count number of sets [], which get reduced to length of 1 each
        nsets = len(re.findall("(?<=\[).+?(?=\])", pattern))
        ## remove sets from pattern
        pattern_no_sets = ''.join(re.findall("(?<=^|\])[^\]\[]+?(?=$|\[)", pattern))
        ## count number of reps {n}, which get reduced to length of n - 1 each
        ## ## n-1 because the char to be repeated itself will count as the 1st rep
        reps = sum(map(lambda x: int(x.split(',')[-1]) - 1, re.findall("(?<=\{).+?(?=\})", pattern_no_sets)))
        ## remove reps from pattern
        nchars = sum(map(lambda x: len(x), re.findall("(?<=^|\}).+?(?=$|\{)", pattern_no_sets)))
        return nchars + nsets + reps ## sum
    
    ## functions to get pam patterns on 5' and 3' sides of gRNA
    def five_prime(self):
        self.parse()
        return self.pam.split('.')[0]
    def three_prime(self):
        self.parse()
        return self.pam.split('.')[-1]
    def five_prime_maxlen(self):
        return self._pattern_max_len(self.five_prime())
    def three_prime_maxlen(self):
        return self._pattern_max_len(self.three_prime())
