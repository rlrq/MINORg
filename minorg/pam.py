###########
##  PAM  ##
###########

import regex as re
# from minorg.constants import PAM_PATTERNS

DEFAULT_gRNA_LENGTH = 20
DEFAULT_PAM = "NGG"

class PAM():
    """
    PAM object.
    
    Tracks PAM pattern, gRNA length, as well as PAM regex.
    
    Attributes
    ----------
    raw_pam: str
        original PAM pattern used to create this PAM object
    pam: str
        expanded PAM pattern
    gRNA_length: int
        gRNA length
    """
    def __init__(self, pam = "NGG", gRNA_length = DEFAULT_gRNA_LENGTH):
        """
        Create a PAM object.
        
        Arguments:
            pam (str/:class:`~minorg.pam.PAM`): PAM object or str of PAM pattern (if PAM object provided, attributes are copied)
            gRNA_length (int): int length of gRNA (bp)
        
        Returns:
            :class:`~minorg.pam.PAM`
        """
        if isinstance(pam, PAM):
            self.raw_pam = pam.raw_pam
            self.pam = pam.pam
            self.gRNA_length = (1 if gRNA_length is None
                                else gRNA_length if gRNA_length != DEFAULT_gRNA_LENGTH
                                else pam.gRNA_length)
        else:
            # if pam in PAM_PATTERNS:
            #     pam = PAM_PATTERNS[pam]
            self.raw_pam = pam
            self.pam = pam
            self.gRNA_length = gRNA_length
        self.parse()
    
    def __repr__(self):
        return self.pam
    
    ## previously: infer_full_pam
    def infer_full(self) -> None:
        """
        Expand PAM pattern to full pattern by inserting 'N' and '.' where necessary.
        
        In the absence of 'N' in the PAM pattern, MINORg will assume 3' PAM with
        1 spacer base (such as in the 3' 'NGG' of SpCas9).
        If a pattern includes an 'N' at either end, MINORg will assume that the
        gRNA is directly adjacent to the 'N' base of the pattern. 
        To specify a 5' PAM in the absence of 'N' in the PAM pattern, 
        '.' should be inserted where the gRNA is.
        """
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
    
    def expand_ambiguous(self) -> None:
        """
        Map ambiguous bases.
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
    
    def parse(self) -> None:
        """
        Wrapper for :meth:`~minorg.pam.PAM.infer_full` and :meth:`~minorg.pam.PAM.expand_ambiguous`.
        """
        ## infer pam location + spacer if not explicitly described
        self.infer_full()
        ## map ambiguous bases
        self.expand_ambiguous()
        return
    
    ## previously: make_pam_pattern
    def regex(self, gRNA_length: int = None) -> None:
        """
        Generate regex based on PAM and gRNA length.
        
        Square brackets not allowed in PAM pattern unless they contain only A, T, G, C, or U.
        
        Arguments:
            gRNA_length (int): optional, gRNA lenth (bp); if not provided, self.gRNA_length is used.
        
        Returns: 
            _regex.Pattern: regex pattern that can be used to match gRNA
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
        """
        Get maximum pattern length.
        
        Arguments:
            pattern (str): regex pattern for searching for gRNA using PAM
        
        Returns: 
            int: Maximum length of match (including look aheads and look behinds)
        """
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
        """
        Get PAM pattern on 5' side of gRNA
        
        Returns:
            str
        """
        self.parse()
        return self.pam.split('.')[0]
    def three_prime(self):
        """
        Get PAM pattern on 3' side of gRNA
        
        Returns:
            str
        """
        self.parse()
        return self.pam.split('.')[-1]
    def five_prime_maxlen(self):
        """
        Get maximum length of PAM pattern on 5' side of gRNA.
        
        Calls :meth:`~minorg.pam.PAM._pattern_max_length` on output of :meth:`~minorg.pam.PAM.five_prime`.
        
        Returns:
            int
        """
        return self._pattern_max_len(self.five_prime())
    def three_prime_maxlen(self):
        """
        Get maximum length of PAM pattern on 3' side of gRNA.
        
        Calls :meth:`~minorg.pam.PAM._pattern_max_length` on output of :meth:`~minorg.pam.PAM.three_prime`.
        
        Returns:
            int
        """
        return self._pattern_max_len(self.three_prime())



## PAM

"""
Reference for PAM patterns: https://www.synthego.com/guide/how-to-use-crispr/pam-sequence
"""

## if these variables are updated, ensure that the dictionary in constants.py is updated as well
SpCas9 = ".NGG"
"""str: SpCas9 (Streptococcus pyogenes) PAM
   
   3' NGG"""
spcas9 = ".NGG"
"""str: SpCas9 (Streptococcus pyogenes) PAM
   
   3' NGG"""

SaCas9N = ".NGRRN"
"""str: SaCas9N (Staphylococcus aureus) PAM
   
   3' NGRRN"""
sacas9n = ".NGRRN"
"""str: SaCas9N (Staphylococcus aureus) PAM
   
   3' NGRRN"""
SaCas9T = ".NGRRT"
"""str: SaCas9T (Staphylococcus aureus) PAM

   3' NGRRT"""
sacas9t = ".NGRRT"
"""str: SaCas9T (Staphylococcus aureus) PAM

   3' NGRRT"""

NmeCas9 = ".NNNNGATT"
"""str: NmeCas9 (Neisseria meningitidis) PAM

   3' NNNNGATT"""
nmecas9 = ".NNNNGATT"
"""str: NmeCas9 (Neisseria meningitidis) PAM

   3' NNNNGATT"""

CjCas9 = ".NNNNRYAC"
"""str: CjCas9 (Campylobacter jejuni) PAM

   3' NNNNRYAC"""
cjcas9 = ".NNNNRYAC"
"""str: CjCas9 (Campylobacter jejuni) PAM

   3' NNNNRYAC"""

StCas9 = ".NNAGAAW"
"""str: StCas9 (Streptococcus thermophilus) PAM

   3' NNAGAAW"""
stcas9 = ".NNAGAAW"
"""str: StCas9 (Streptococcus thermophilus) PAM

   3' NNAGAAW"""

Cas12a = "TTTV."
"""str: Cas12a (Cpf-1) PAM

   5' TTTV"""
cas12a = "TTTV."
"""str: Cas12a (Cpf-1) PAM

   5' TTTV"""

AacCas12b = "TTN."
"""str: AacCas12b (Alicyclobacillus acidiphilus) PAM

   5' TTN"""
aaccas12b = "TTN."
"""str: AacCas12b (Alicyclobacillus acidiphilus) PAM

   5' TTN"""

BhCas12b = "DTTN."
"""str: BhCas12b (Bacillus hisashii) PAM

   5' DTTN"""
bhcas12b = "DTTN."
"""str: BhCas12b (Bacillus hisashii) PAM

   5' DTTN"""

