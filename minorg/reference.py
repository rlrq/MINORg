"""Classes that integrate FASTA and GFF3 files for retrieval of feature sequences"""

import sys
sys.path.append("/mnt/chaelab/rachelle/scripts/minorgpy")

import Bio
import copy
import tempfile

from Bio import Seq

from minorg.exceptions import MessageError

from minorg.functions import (
    is_range,
    convert_range,
    ranges_to_pos,
    ranges_union,
    ranges_intersect,
    adjusted_ranges,
    adjusted_pos,
    empty_file,
    parse_get_data,
    write_tsv
)

from minorg.blast import (
    BlastNR,
    filter_rpsblast_for_domain
)

from minorg.fasta import (
    dict_to_fasta,
    collapse_identical_seqs
)

from minorg.annotation import Annotation, GFF# , subset_ann
from minorg.index import IndexedFasta

from minorg.display import (
    make_print_preindent
)

# ## test
# gid = "AT3G44400"
# gid2 = "AT3G44480"
# pssm = 366375
# rblast = 'rpsblast+'
# rdb = '/mnt/chaelab/shared/blastdb/Cdd.v3.18/Cdd'
# mktmp = lambda fname: "/mnt/chaelab/rachelle/scripts/minorgpy/test/homologue/tmp/" + fname
# ref = AnnotatedFasta(fasta = "/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta", gff = "/mnt/chaelab/rachelle/scripts/minorgpy/test/homologue/test_DM2_2/reduced_ann/reduced.TAIR10.gff", name = "TAIR10")
# g = ref.annotated_feature(gid)
# g._get_domain_ranges(pssm, mktmp = mktmp, rpsblast = rblast, db = rdb, by_gene = True)
# g._get_domain_seq(pssm, by_gene = True, feature = "CDS", adj_dir = True, mktmp = mktmp)
# ## {'AT3G44400': {(813, 1524): Seq('GCTCATATGGAAAGGACGGAACAATTGTTACGCCTAGATTTGGATGAAGTAAGG...GCC')}}
# g._get_domain_seq(pssm, by_gene = False, feature = "CDS", adj_dir = True, mktmp = mktmp)
# # {'AT3G44400.1': {(780, 1491): Seq('GCTCATATGGAAAGGACGGAACAATTGTTACGCCTAGATTTGGATGAAGTAAGG...GCC')}, 'AT3G44400.2': {(813, 1524): Seq('GCTCATATGGAAAGGACGGAACAATTGTTACGCCTAGATTTGGATGAAGTAAGG...GCC')}}
# ref.get_feature_seq(gid, feature = "CDS", pssm_id = pssm, complete = False, adj_dir = True,
#                     translate = False, mktmp = mktmp)
# ref.get_feature_seq(gid2, feature = "CDS", pssm_id = pssm, complete = False, adj_dir = True,
#                     translate = False, mktmp = mktmp, db = rdb, rpsblast = rblast)

class InvalidFeatureID(MessageError):
    """
    Raised when feature_id does not exist.
    """
    def __init__(self, feature_id):
        super().__init__( f"Error: Invalid feature ID. {feature_id} does not exist." )

#######################
##  REFERENCE_MANIP  ##
#######################

class AnnotatedFeature(Annotation):
    """
    Representation of GFF3 feature annotation, including its subfeatures.
    
    Attributes:
        annotated_fasta (:class:`AnnotatedFasta`): parent :class:`AnnotatedFasta` object
        annotation (:class:`minorg.annotation.GFF`): 
    """
    def __init__(self, id, annotated_fasta, gff = None, rps_hits = None):
        """
        Create an AnnotatedFeature object.
        
        Annotation data for feature and its subfeatures will be extracted from
        either ``gff`` or ``annotated_fasta`` 's :attr:`minorg.reference.AnnotatedFasta.annotation` attribute.
        
        Arguments:
            id (str): required, ID of feature (typically gene name or isoform name)
            annotated_fasta (:class:`AnnotatedFasta`): required, :class:`AnnotatedFasta` object
            gff (str or :class:`minorg.annotation.GFF`): optional, 
                path to GFF file or :class:`minorg.annotation.GFF` object
            rps_hits (str): optional, path to RPS-BLAST output (outfmt 6, with header)
        """
        self.annotated_fasta = annotated_fasta
        self.assembly = self.annotated_fasta.assembly
        self.rps_hits = rps_hits
        self.id = id
        ## subset original GFF data & copy attributes from original GFF obj
        annotation = (self.annotated_fasta.annotation if gff is None else
                      gff if isinstance(gff, GFF) else GFF(fname = gff))
        data = annotation.get_features_and_subfeatures(self.id, full = True, index = False)
        if not data:
            raise InvalidFeatureID(id)
        self.annotation = annotation.empty_copy()
        self.annotation._data = data
        ## generate Annotation obj attributes
        try:
            super().__init__(entry = self.annotation.get_id(self.id, output_list=False).generate_gff(),
                             gff = self.annotation)
        except Exception as e:
            print(self.id, self.annotation.get_id(self.id, output_list = True))
            raise e
    
    @property
    def length(self):
        """
        Length of feature
        """
        return self.end - self.start

    def copy_annotation(self) -> Annotation:
        """
        Generate new :class:`minorg.annotation.Annotation` object with self's annotation data.
        
        Used to make new Annotation for domains that will be modified with domain range.
        
        Returns
        -------
        :class:`minorg.annotation.Annotation`
        """
        return Annotation(entry = self.generate_gff(), gff = self.annotation)
    
    def convert_range(self, ranges, index = 0, start_incl = True, end_incl = False) -> list:
        """
        Convert range(s) from 0-index [start, end).
        
        Arguments:
            index (int): output range index
            start_incl (bool): whether output range is start-inclusive
            end_incl (bool): whether output range is end-inclusive
        
        Returns
        -------
        tuple or list
            In same structure as input.
        """
        return convert_range(ranges, index_in = 0, index_out = index,
                             start_incl_in = True, start_incl_out = start_incl,
                             end_incl_in = False, end_incl_out = end_incl)
    
    # def range(self, **conversion_kwargs) -> list:
    #     """
    #     Get genomic range of feature (not adusted for sense)
    #     """
    #     return convert_range([(self.start, self.end)], **conversion_kwargs)[0]
    
    def subset(self, id) -> 'AnnotatedFeature':
        return AnnotatedFeature(id = id, annotated_fasta = self.annotated_fasta, gff = self.annotation)
    
    def range_from_genomic(self, ranges, **conversion_kwargs) -> list:
        """
        Output ranges relative to sense strand of self, where first base on sense strand of feature is 0
        
        Arguments:
            ranges (list): genomic range(s), 0-indexed, start-inclusive, end-exclusive
                (e.g. [(start1, end1), (start2, end2), ...])
            index (int): output range index
            start_incl (bool): whether output range is start-inclusive
            end_incl (bool): whether output range is end-inclusive
        
        Returns
        -------
        list of tuples
            Of ranges
        """
        if self.plus:
            ranges = [(start - self.start, end - self.start) for start, end in ranges]
        else:
            ranges = [(abs(end - self.end), abs(start - self.end)) for start, end in ranges]
        return self.convert_range(sorted(ranges), **conversion_kwargs)
    
    def range_to_genomic(self, ranges, **conversion_kwargs) -> list:
        """
        Convert ranges relative to sense strand of feature to their genomic positions.
        
        Arguments:
            ranges (list): range(s) relative to sense strand of self, 0-indexed, start-inclusive, end-exclusive
                (e.g. [(start1, end1), (start2, end2), ...])
            index (int): output range index
            start_incl (bool): whether output range is start-inclusive
            end_incl (bool): whether output range is end-inclusive
        
        Returns
        -------
        list of tuples
            Of ranges
        """
        if is_range(ranges):
            start, end = ranges
            if self.plus: r = (start + self.start, end + self.start)
            else: r = (abs(end - self.end), abs(start - self.end))
            return self.convert_range([r], **conversion_kwargs)[0]
        else:
            return [self.range_to_genomic(r, **conversion_kwargs) for r in ranges]
    
    def feature_range(self, feature_type, genomic = False, union = True) -> list:
        """
        Return ranges of entries of a given feature type (e.g. "CDS")
        
        Arguments:
            feature_type (str): GFF3 feature type
            genomic (bool): output range in genome instead of relative to sense strand of current feature
            union (bool): merge overlapping ranges
        
        Returns
        -------
        list of tuples
            Of ranges
        """
        feature_ranges = [(ann.start, ann.end) for ann in
                          self.annotation.get_features(feature_type, index = False)]
        if not genomic:
            feature_ranges = self.range_from_genomic(feature_ranges)
        if union:
            return ranges_union([feature_ranges])
        return feature_ranges
    
    def subfeature_range(self, subfeature_id, feature_type = None, genomic = False, union = True) -> list:
        """
        Return range of given feature type (e.g. "CDS") of a sub-feature.
        (e.g. If self.type == "gene", self.subfeature_range("<mRNA isoform ID>", "CDS") 
        will return CDS range(s) of the specified mRNA isoform)
        
        Arguments:
            subfeature_id (str): subfeature_id
            feature_type (str): GFF3 feature type
            genomic (bool): output range in genome instead of relative to sense strand of current feature
            union (bool): merge overlapping ranges
        
        Returns
        -------
        list of tuples
            Of ranges
        """
        subfeature_ann = self.subset(subfeature_id)
        if feature_type is None: feature_type = subfeature_ann.type
        output = subfeature_ann.feature_range(feature_type, union = union, genomic = genomic)
        if not genomic:
            offset = (subfeature_ann.start - self.start) if self.plus else (self.end - subfeature_ann.end)
            output = [(start + offset, end + offset) for start, end in output]
        return output
    
    def adj_feature_range(self, feature_type, ranges, complete = True, **feature_range_kwargs) -> list:
        """
        Convert range in subfeatures of feature type to range in current feature.
        (For example, convert the range of a domain in a concatenated CDS sequence to position in mRNA.)
        
        Arguments:
            feature_type (str): feature type
            ranges (list of tuple): ranges relative to subfeatures of feature type
            complete (bool): merge output range(s) together into single range.
                Output will stll be a list of tuple. (e.g. [(<smallest start>, <largest end>)])
        
        Returns
        -------
        list of tuples
            Of ranges
        """
        feature_bounds = self.feature_range(feature_type, **feature_range_kwargs)
        feature_pos = ranges_to_pos(feature_bounds)
        ## coopt adjusted_ranges (meant to account for gaps in alignment)
        feature_dummy_seq = ''.join([('A' if p in feature_pos else '-')
                                     for p in range(max(feature_pos))])
        return adjusted_ranges(feature_dummy_seq, *ranges, subtract_gaps = (not complete))
    
    def adj_subfeature_range(self, subfeature_id, ranges, complete = True, **subfeature_range_kwargs) -> list:
        """
        Convert range in a subfeature to range in current feature.
        (For example, convert the range in mRNA to position in gene.)
        
        Arguments:
            subfeature_id (str): subfeature ID
            ranges (list of tuple): ranges relative to sub-subfeatures of feature type in subfeature
            complete (bool): merge output range(s) together into single range.
                Output will stll be a list of tuple. (e.g. [(<smallest start>, <largest end>)])
        
        Returns
        -------
        list of tuples
            Of ranges
        """
        subfeature_bounds = self.subfeature_range(subfeature_id, genomic = False, **subfeature_range_kwargs)
        subfeature_pos = ranges_to_pos(subfeature_bounds)
        ## coopt adjusted_ranges (meant to account for gaps in alignment)
        subfeature_dummy_seq = ''.join([('A' if p in subfeature_pos else '-')
                                        for p in range(max(subfeature_pos))])
        return adjusted_ranges(subfeature_dummy_seq, *ranges, subtract_gaps = (not complete))
    
    def range_in_genomic_coords(self, ranges, index=0, start_incl=True, end_incl=False) -> list:
        """
        Convert range in feature to genomic coordinates.
        
        Arguments:
            ranges (list of tuple): [(start1, end1), (start2, end2), ...]
                - relative to sense strand of feature, where first base on sense strand of feature is 0
                - 0-indexed; start-inclusive, end-exclusive
            index (int): output index
            start_incl (bool): start inclusive output ranges
            end_incl (bool): end inclusive output ranges
        
        Returns
        -------
        list of tuples
            Of ranges
        """
        if self.plus: ranges = [(start + self.start, end + self.start) for start, end in ranges]
        else: ranges = [(self.length - end, self.length - start) for start, end in ranges]
        if not start_incl: ranges = [(start - 1, end) for start, end in ranges]
        if end_incl: ranges = [(start, end - 1) for start, end in ranges]
        return [(start + index, end + index) for start, end in ranges]
    
    def _get_seq(self, feature_type=None, ranges=[],
                 adj_dir=True, complete=False, translate=False) -> Bio.Seq.Seq:
        """
        Get sequence of ``feature_type`` by concatenating subfeatures of feature_type.
        Final range from which sequence will be extracted is the intersection of ``ranges`` and 
        ranges of the relevant subfeatures.
        
        Arguments:
            feature_type (str): GFF feature type (e.g. "CDS")
            ranges (list of tuple): 0-indexed ranges, relative to sense strand of self
            adj_dir (bool): output sense strand
            complete (bool): merge output range(s) together into single range.
                Output will stll be a list of tuple. (e.g. [(<smallest start>, <largest end>)])
            translate (bool): translate sequence.
                For biologically sound results, ``translate=True`` only be used with
                ``feature_type="CDS"``, ``adj_dir=True``, and ``complete=False``.
        
        Returns
        -------
        Bio.Seq.Seq
        """
        if feature_type is None: feature_type = self.type
        feature_ranges = self.feature_range(feature_type, union = True, genomic = True)
        if ranges:
            output_ranges = sorted(ranges_intersect(feature_ranges, self.range_to_genomic(ranges)))
        else:
            output_ranges = sorted(feature_ranges)
        if complete:
            output_ranges = [(min(r[0] for r in output_ranges), max(r[1] for r in output_ranges))]
        ## start extracting seq
        seq = Bio.Seq.Seq('')
        for start, end in output_ranges:
            seq += self.assembly[self.seqid][start:end]
        if adj_dir and not self.plus:
            seq = seq.reverse_complement()
        if translate:
            seq = seq.translate(table = self.annotated_fasta.genetic_code)
        return seq
    
    def get_seq(self,
                ## shared options
                fout=None, feature_type=None, adj_dir=True, complete=True, translate=False,
                ## _get_seq options
                ranges=[],
                ## _get_domain_seq options
                pssm_id=None, by_gene=False, db=None, remote_rps=False, rpsblast="rpsblast",
                rps_hits=None, thread=None,
                mktmp=None, quiet=False, fout_gff=None,
                ## seqid options
                seqid_template = "$source|$gene|$feature|$isoform|$domain|$n|$complete|$strand|$sense",
                apply_template_to_dict = False) -> dict:
        """
        Get sequence of feature or subfeature.
        All arguments are strictly speaking optional.
        If no arguments are specified, self's sense strand sequence will be returned.
        For more complicated queries, 
        in addition to modifying flags ``adj_dir``, ``complete``, and/or ``translate``,
        users may optionally specify a feature type, Pssm-Id, and/or
        restrict range relative to self.
        
        Arguments:
            fout (str): path to output file; if not provided, dictionary of seqs will be returned
            feature_type (str): GFF feature type (e.g. "CDS")
            adj_dir (bool): output sense strand
            complete (bool): merge output range(s) together into single range.
                Output will stll be a list of tuple. (e.g. [(<smallest start>, <largest end>)])
            translate (bool): translate sequence.
                For biologically sound results, ``translate=True`` only be used with
                ``feature_type="CDS"``, ``adj_dir=True``, and ``complete=False``.
            ranges (list of tuple): 0-indexed ranges, relative to sense strand of self
            pssm_id (str or list): Pssm-Id of domain(s) to restrict output sequence range to.
                If provided and ``rps_hits=None``, 
                RPS-BLAST will be conducted on concatenated CDS to search for domain(s).
            by_gene (bool): merge overlapping domain ranges from different isoforms
            db (str): path to RPS-BLAST database
            remote_rps (bool): use remote RPS-BLAST; if True, ``db`` is not required
            rpsblast (str): path to rpsblast or rpsblast+ executable
                (or name of command, if available at command line)
            rps_hits (str): path to RPS-BLAST output, outfmt 6, with header. 
                If provided, ``db`` and ``rpsblast`` are not required.
            thread (int): number of threads for RPS-BLAST
            mktmp (func): function that takes a file name and generates an absolute path, 
                used for temporary files.
                ``tempfile.mkstemp`` will be used if not provided.
            quiet (bool): print only essential messages
            fout_gff (str): path to output GFF3 file in which to write domain features, entirely optional
            seqid_template (str): optional, template for output sequence name.
                Template will be parsed by strings.Template.
                The default template is "$source|$gene|$feature|$isoform|$domain|$n|$complete|$strand|$sense".
                    
                    - $source: reference genome alias
                    - $gene: gene/isoform ID
                    - $feature: GFF3 feature type
                    - $isoform: isoform ID if by_gene = False, else same as $gene
                    - $domain: Pssm-Id
                    - $n: if multiple domains are present, they will be numbered according to 
                      proximity to 5' of sense strand
                    - $complete: 'complete' if ``complete=True`` else 'stitched'
                    - $strand: 'minus' if ``adj_dir=False`` and not self.plus else 'plus'
                    - $sense: 'NA' if self.strand is not set
                      else 'antisense' if ``adj_dir=False`` and not self.plus else 'sense'
            
            apply_template_to_dict (bool): flatten output dict and assign sequence names
        
        Returns
        -------
        dict of dict of Bio.Seq.Seq
            If ``apply_template_to_dict=False``. 
            Format if ``pssm_id=None``: {self.id: {(0, self.length): <sequence>}}
            Format if pssm_id is provided:
            {<isoform_id>: {(start pos of domain in self, end pos of domain in self): <sequence>}}
        dict of Bio.Seq.Seq
            If ``apply_template_to_dict=True``
        """
        feature_type = feature_type if feature_type is not None else self.type
        ## select appropriate get_seq function
        ## seqs fmt: {id: {(<start>, <end>): <seq>}}
        if pssm_id and pssm_id != "gene":
            if isinstance(pssm_id, str): pssm_ids = pssm_id.split(',')
            elif isinstance(pssm_id, int): pssm_ids = [str(pssm_id)]
            else: pssm_ids = pssm_id
            seqs = self._get_domain_seq(*pssm_ids, by_gene = by_gene, feature_type = feature_type,
                                        adj_dir = adj_dir, complete = complete, translate = translate,
                                        db = db, remote_rps = remote_rps, rpsblast = rpsblast,
                                        rps_hits = rps_hits, thread = thread,
                                        fout_gff = fout_gff, mktmp = mktmp, quiet = quiet)
        else:
            seq = self._get_seq(adj_dir = adj_dir, complete = complete, translate = translate,
                                feature_type = feature_type, ranges = ranges)
            seqs = {self.id: {(0, self.length): seq}}
        ## format seqid
        if seqid_template and (fout or apply_template_to_dict):
            ## create function to make seqid
            from string import Template
            def mk_seqid(**kwargs):
                return Template(seqid_template).substitute(source = self.annotated_fasta.name,
                                                           gene = self.id, feature = feature_type,
                                                           complete = ("complete" if complete else "stitched"),
                                                           domain = (','.join(pssm_ids) if pssm_id else '.'),
                                                           strand = ("minus" if (adj_dir and not self.plus)
                                                                     else "plus"),
                                                           sense = ("NA"
                                                                    if (not self.strand or self.strand== '.')
                                                                    else "antisense"
                                                                    if (not adj_dir and not self.plus)
                                                                    else "sense"), **kwargs)
            ## flatten seqs
            seqs = {mk_seqid(isoform = f"{r[0]+1}-{r[1]}", n = n+1): seq
                    for isoform, isoform_dat in seqs.items()
                    for n, seq_dat in enumerate(sorted(isoform_dat.items()))
                    for r, seq in [seq_dat]}
        ## output
        if fout:
            seqs = {seqid if isinstance(seqid, str) else '|'.join(seqid): seq
                    for seqid, seq in seqs.items()}
            dict_to_fasta(seqs, fout)
        else:
            return seqs
    
    def _get_domain_seq(self, *pssm_ids, by_gene=False, feature_type=None,
                        ## get_seq options
                        adj_dir=True, complete=False, translate=False,
                        ## get_domain_ranges options
                        db=None, remote_rps=False, rpsblast="rpsblast", rps_hits=None, thread=None,
                        mktmp=None, quiet=False, fout_gff=None) -> dict:
        """
        Get sequence of domain in feature or subfeature.
        Except for ``pssm_ids``, all arguments are strictly speaking optional.
        If no arguments are specified, 
        the spliced sense strand sequence of domain in each mRNA isoform will be returned.
        For more complicated queries, 
        users may modify the following flags ``adj_dir``, ``complete``, and/or ``translate`` and/or
        and optionally specify a feature type.
        
        Arguments:
            *pssm_ids (str): required, Pssm-Ids of desired domain(s)
            by_gene (bool): merge overlapping domain ranges from different isoforms
            feature_type (str): GFF feature type (e.g. "CDS")
            adj_dir (bool): output sense strand
            complete (bool): merge output range(s) together into single range.
                Output will stll be a list of tuple. (e.g. [(<smallest start>, <largest end>)])
            translate (bool): translate sequence.
                For biologically sound results, ``translate=True`` only be used with
                ``feature_type="CDS"``, ``adj_dir=True``, and ``complete=False``.
            db (str): path to RPS-BLAST database
            remote_rps (bool): use remote RPS-BLAST database; if True, ``db`` is not required
            rpsblast (str): path to rpsblast or rpsblast+ executable
            rps_hits (str): path to RPS-BLAST ouput, fmt 6, with header
            thread (int): number of threads for RPS-BLAST
            mktmp (func): function that takes a file name and generates an absolute path, 
                used for temporary files.
                ``tempfile.mkstemp`` will be used if not provided.
            quiet (bool): print only essential messages
            fout_gff (str): path to output GFF3 file in which to write domain features, entirely optional
            seqid_template (str): optional, template for output sequence name.
                Template will be parsed by strings.Template.
                The default template is "$source|$gene|$feature|$isoform|$domain|$n|$complete|$strand|$sense".
                    
                    - $source: reference genome alias
                    - $gene: gene/isoform ID
                    - $feature: GFF3 feature type
                    - $isoform: isoform ID if by_gene = False, else same as $gene
                    - $domain: Pssm-Id
                    - $n: if multiple domains are present, they will be numbered according to 
                      proximity to 5' of sense strand
                    - $complete: 'complete' if ``complete=True`` else 'stitched'
                    - $strand: 'minus' if ``adj_dir=False`` and not self.plus else 'plus'
                    - $sense: 'NA' if self.strand is not set
                      else 'antisense' if ``adj_dir=False`` and not self.plus else 'sense'
            
            apply_template_to_dict (bool): flatten output dict and assign sequence names
        
        Returns
        -------
        dict of dict of Bio.Seq.Seq
            {<isoform_id>: {(start pos of domain in self, end pos of domain in self): <Bio.Seq.Seq sequence>}}
        """
        ## get domain ranges (complete) (relative to each isoform)
        domain_ranges = self._get_domain_ranges(*pssm_ids, db = db, remote_rps = remote_rps,
                                                rpsblast = rpsblast, rps_hits = rps_hits, thread = thread,
                                                by_gene = by_gene, genomic = False,
                                                mktmp = mktmp, quiet = quiet, fout_gff = fout_gff)
        ## get seqs; discontiguous domains are NOT concatenated
        domain_seqs = {isoform: {tuple(hit_range): self.subset(isoform)._get_seq(feature_type,
                                                                                 ranges = [hit_range],
                                                                                 complete = complete,
                                                                                 translate = translate)
                                 for hit_range in hit_ranges}
                       for isoform, hit_ranges in domain_ranges.items()}
        return domain_seqs ## {isoform: {<range of first domain>: <seq>, <range of second domain>: <seq>}...}
    
    def _get_domain_ranges(self, *pssm_ids, db=None, remote_rps=False, rpsblast="rpsblast",
                           rps_hits=None, thread=None,
                           genomic=False, by_gene=False, mktmp=None, quiet=False, fout_gff=None) -> dict:
        """
        Retrieve ranges of domain(s). Overlapping domain hits are merged.
        
        Except for ``pssm_ids``, all arguments are strictly speaking optional.
        If no arguments are specified, 
        the complete range of domain(s) (i.e. including any intervening introns)
        in each mRNA isoform relative to self's sense strand will be returned.
        
        If self contains mRNA w/ CDS, CDS will be combined and translated for each mRNA separately,
            searched for the relevant domain (``pssm_ids``) against a domain database (``db``),
            and coordinates will be mapped back to the sense strand of self
        If ``rps_hits`` not provided AND self.rps_hits is None, rpsblast will be executed to generate rps_hits file
            - rps_hits file is basically blast6 with header
            - required fields: qseqid, sseqid. qseqid must be mRNA ID.
        
        Arguments:
            *pssm_ids (str): required, Pssm-Ids of desired domain(s)
            db (str): path to RPS-BLAST database, used with ``rpsblast``
            remote_rps (bool): use remote RPS-BLAST, used with ``rpsblast``; if True, ``db`` is not required
            rpsblast (str): path to rpsblast or rpsblast+ executable, used with ``db``
                (or name of command, if available at command line)
            rps_hits (str): path to RPS-BLAST output, outfmt 6, with header. 
                If provided, ``db`` and ``rpsblast`` are not required.
            thread (int): number of threads for RPS-BLAST, used with ``db`` and ``rpsblast`` (default=1)
            genomic (bool): output range in genomic coordinates instead of relative to current feature
            by_gene (bool): merge overlapping ranges of different isoforms
            mktmp (func): function that takes a file name and generates an absolute path, 
                used for temporary files.
                ``tempfile.mkstemp`` will be used if not provided.
            quiet (bool): print only essential messages
            fout_gff (str): path to output GFF3 file in which to write domain features, entirely optional
        
        Returns
        -------
        dict
            Of ranges, grouped by isform if by_gene=False.
            (e.g. {"geneA.1": [(50, 100)], "geneA.2": [(20, 70), (220, 270)]})
        dict
            Of ranges, single entry under self's ID if by_gene=True.
            (e.g. {"geneA": [(50, 100), (250, 300)]})
        """
        ## some functions
        mktmpfname = ((lambda suf: tempfile.mkstemp()[1]) if mktmp is None
                      else (lambda suf: mktmp(f"{self.get_attr('ID',fmt=str)}_{suf}")))
        def printi(msg):
            if not quiet: print(msg)
        ## if rps_hits doesn't alreay exist, execute rps-blast
        if rps_hits is None and self.rps_hits is None:
            ## get isoform sequences
            printi("Extracting peptide sequence(s)")
            isoforms = {entry.get_attr("ID", fmt=list)[0]: self.subset(id=entry.get_attr("ID", fmt=list)[0])
                        for entry in self.annotation if entry.type == "mRNA"}
            isoform_seqs = {isoform: isoform_ann._get_seq(feature_type = "CDS", adj_dir = True, translate = True)
                            for isoform, isoform_ann in isoforms.items()}
            ## collapse identical sequences into non-redundant set, write to file
            tmp_pep = mktmpfname("pep.fasta")
            dict_to_fasta(isoform_seqs, tmp_pep)
            ## execute rpsblast
            printi("Extracting domain range(s)")
            tmp_rpsblast = mktmpfname(f"{'-'.join(map(str, pssm_ids))}.tsv")
            from Bio.Blast.Applications import NcbirpsblastCommandline
            BlastNR(query = tmp_pep, blastf = NcbirpsblastCommandline,
                    header = "qseqid,sseqid,pident,length,qstart,qend",
                    fout = tmp_rpsblast, keep_tmp = False,
                    db = db, cmd = rpsblast, remote_rps = remote_rps,
                    # mt_mode = (1 if thread is not None else 0),
                    num_threads = (1 if thread is None else thread))
            ## store path to output
            self.rps_hits = tmp_rpsblast
        ## filter hits for desired domain
        pssm_ids = list(map(str, pssm_ids))
        rps_hits = rps_hits if rps_hits is not None else self.rps_hits
        raw_hits = filter_rpsblast_for_domain(rps_hits, *pssm_ids)
        ## get coords (converted from 1-index [start, end] peptide to 0-index [start, end) genome)
        aa_to_nt = lambda r: (r[0]*3, (r[1])*3)
        nt_hits = {isoform: ranges_union([[aa_to_nt((start-1, end)) for start, end in isoform_hits]])
                   for isoform, isoform_hits in raw_hits.items()}
        adj_hits = {isoform: self.subset(isoform).adj_feature_range("CDS", hit_ranges, complete = True,
                                                                    genomic = False, union = True)
                    for isoform, hit_ranges in nt_hits.items()}
        ## combine ranges if by_gene or fout_gff is requested
        if by_gene or fout_gff:
            by_gene_hits = {self.id:
                            ranges_union([self.adj_subfeature_range(isoform, hit_ranges, complete = True)
                                          for isoform, hit_ranges in adj_hits.items()])}
        ## write domain ranges to file
        if fout_gff:
            domain_gff = GFF(fname = fout_gff, attr_mod = self.annotated_fasta.attr_mod, quiet = True)
            domain_ranges = by_gene_hits[self.id]
            for r in domain_ranges:
                start, end = self.range_to_genomic(r, index = 0, end_incl = False)
                domain_ann = self.copy_annotation()
                domain_ann.start = start
                domain_ann.end = end
                domain_ann.type = "domain"
                domain_ann.source = "minorg" if self.annotated_fasta.name is None else self.annotated_fasta.name
                domain_gff.add_entry(domain_ann, duplicate_check = True)
            if len(domain_gff) != 0:
                domain_gff.write(domain_gff._fname, domain_gff._data, fmt = domain_gff._fmt)
        ## output
        if by_gene:
            return by_gene_hits
        else:
            return adj_hits

## inherits from IndexedFasta so AnnotatedFasta_obj["chr"][:50] works
class AnnotatedFasta(IndexedFasta):

    """
    Representation of FASTA file with GFF3 annotations
    
    Attributes:
        name (str): name/alias of reference genome
        fasta (str): path to genome FASTA file
        _gff (str): path to GFF3 file
        gff (str): path to GFF3 file.
            If :meth:`~minorg.reference.AnnotatedFasta.subset_annotation` has been called,
            this points to the file of subsetted annotations.
        genetic_code (str or int): NCBI genetic code
        attr_mod (dict): attribute modification 
        memsave (bool): whether to index GFF3 file instead of reading to memory
        assembly (:class:`~minorg.reference.AnnotatedFasta`): self, stores parsed sequence data.
            As this class inherits from IndexedFasta,
            ``self.assembly`` provides an analogous way of retrieving sequence data 
            as ``self.annotation`` retrieves annotation data.
        annotation (:class:`~minorg.annotation.GFF`): GFF object, stores parsed annotation data.
            If :meth:`~minorg.reference.AnnotatedFasta.subset_annotation` has been called,
            this stores data of subsetted annotations.
        _annotated_features (dict): stores recently retrieved AnnotatedFeature objects
    """
    
    def __init__(self, fasta, gff, name = "reference", attr_mod = {}, genetic_code = 1,
                 memsave = True):
        """
        Create an AnnotatedFasta object.
        
        Arguments:
            fasta (str): path to FASTA file of genome assembly
            gff (str): path to GFF3 file of genome annotation
            name (str): name of genome (default='refrence')
            attr_mod (dict): attribute modifications (default={})
                (fmt:{'<feature>':{'<standard attribute field name>':'<nonstandard attribute field name used>'}})
                (e.g.: {'mRNA': {'Parent': 'Locus_id'}})
            genetic_code (int or str): NCBI genetic code (default=1)
                (see: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
                (1 --> Standard Code; 2 --> Vertebrate Mitochondrial Code; etc.)
            memsave (bool): index GFF3 file instead of reading data to memory
        """
        super().__init__(fasta)
        self.name = name
        self.fasta = fasta
        self._gff = gff
        self.gff = gff
        self.genetic_code = genetic_code if genetic_code else 1
        self.attr_mod = attr_mod if attr_mod else {}
        self.memsave = memsave
        self.assembly = self
        self.annotation = None
        self.parse_gff()
        self._annotated_features = {}
    
    def __repr__(self):
        return f"AnnotatedFasta(name = {self.name}, fasta = {self.fasta}, gff = {self.gff})"
        
    def parse_gff(self) -> None:
        """Parse GFF3 file stored :attr:`~minorg.reference.AnnotatedFasta.gff`"""
        self.annotation = GFF(self.gff, attr_mod = self.attr_mod, memsave = self.memsave)
    
    def subset_annotation(self, *ids, fout = None, memsave = None, preserve_order = True) -> None:
        """
        Subset GFF3 annotations by feature ID (including subfeatures) 
        and generate new file of subsetted features.
        
        Arguments:
            *ids (str): required, list or tuple of feature IDs
            fout (str): optional, path to output file.
                If not provided, ``tempfile.mkstemp`` will be used to generate a temporary file.
        """
        if len(ids) == 0:
            print("At least 1 feature ID is required for AnnotatedFasta.subset_annotation")
        if fout is None:
            fout = tempfile.mkstemp()[1]
        subset_annotation = self.annotation.subset(feature_ids = ids, subfeatures = True,
                                                   preserve_order = preserve_order)
        subset_annotation.write(fout = fout)
        # subset = subset_ann(gff_beds = {self.name: self._gff}, ids = ids, fout_fmt = "GFF",
        #                     mk_tmpf_name = lambda x: fout, attr_mods = {self.name: self.attr_mod},
        #                     memsave = (memsave if memsave is not None else self.memsave),
        #                     sort = sort)
        # self.gff = fout
        # self.parse_gff()
        self.gff = fout
        self.annotation = subset_annotation
        return
    
    def reduce_annotation(self, ids, fout = None, preserve_order = True) -> None:
        """
        Wrapper for ``self.subset_annotation`` so functions written before the renaming won't break.
        """
        return self.subset_annotation(ids, fout = fout, preserve_order = preserve_order)
    
    def feature_range(self, feature_id) -> tuple:
        """
        Get range of feature.
        
        Arguments:
            feature_id (str): feature ID
        
        Returns
        -------
        tuple
            Of range. E.g (<start>, <end>)
        """
        ann = self.annotation.get_id(feature_id, output_list = False)
        return (ann.start, ann.end)
    
    def feature_pos(self, feature_id) -> set:
        """
        Get range of feature and convert to set of positions.
        
        Arguments:
            feature_id (str): feature ID
        
        Returns
        -------
        set
            Of integer value coordinates within feature range
        """
        return ranges_to_pos([self.feature_range(feature_id)])
    
    def annotated_feature(self, feature_id) -> Annotation:
        """
        Retrieve annotated feature.
        
        Also store feature's annotation in ``self._annotated_features`` for future quick retrieval.
        If this method has been called before with the same feature ID, 
        output will be retrieved from ``self._annotated_features``.
        
        Returns
        -------
        :class:`minorg.annotation.Annotation`
            Annotation object of feature
        """
        annotated_feature = self._annotated_features.get(feature_id, None)
        if annotated_feature is None:
            annotated_feature = AnnotatedFeature(id = feature_id, annotated_fasta = self)
            self._annotated_features[feature_id] = annotated_feature
        return annotated_feature
    
    def _get_feature_seq(self, feature_id, fout=None, feature_type=None,
                         ## get_seq options
                         complete=False, adj_dir=False, translate=False, by_gene=False,
                         ## get_domain_seq options
                         pssm_id=None, db=None, remote_rps=False, rpsblast="rpsblast",
                         rps_hits=None, thread=None,
                         mktmp=None, fout_gff=None, quiet=True,
                         seqid_template = "$source|$gene|$feature|$isoform|$domain|$n|$complete|$strand|$sense",
                         apply_template_to_dict = False) -> dict:
        """
        Get sequences of a single feature.
        Except for ``feature_ids``, all arguments are strictly speaking optional.
        If no arguments are specified, the sense strand of the specified features will be returned.
        For more complicated queries, 
        in addition to modifying flags ``adj_dir``, ``complete``, and/or ``translate``,
        users may optionally specify a feature type and/or Pssm-Id.
        
        Arguments:
            fout (str): path to output file; if not provided, dictionary of seqs will be returned
            feature_type (str): GFF feature type (e.g. "CDS")
            adj_dir (bool): output sense strand
            complete (bool): merge output range(s) together into single range.
                Output will stll be a list of tuple. (e.g. [(<smallest start>, <largest end>)])
            translate (bool): translate sequence.
                For biologically sound results, ``translate=True`` only be used with
                ``feature_type="CDS"``, ``adj_dir=True``, and ``complete=False``.
            pssm_id (str or list): Pssm-Id of domain(s) to restrict output sequence range to.
                If provided and ``rps_hits=None``, 
                RPS-BLAST will be conducted on concatenated CDS to search for domain(s).
            by_gene (bool): merge overlapping domain ranges from different isoforms
            db (str): path to RPS-BLAST database
            remote_rps (bool): use remote RPS-BLAST; if True, ``db`` is not required
            rpsblast (str): path to rpsblast or rpsblast+ executable
                (or name of command, if available at command line)
            rps_hits (str): path to RPS-BLAST output, outfmt 6, with header. 
                If provided, ``db`` and ``rpsblast`` are not required.
            thread (int): number of threads for RPS-BLAST
            mktmp (func): function that takes a file name and generates an absolute path, 
                used for temporary files.
                ``tempfile.mkstemp`` will be used if not provided.
            quiet (bool): print only essential messages
            fout_gff (str): path to output GFF3 file in which to write domain features, entirely optional
            seqid_template (str): optional, template for output sequence name.
                Template will be parsed by strings.Template.
                The default template is "$source|$gene|$feature|$isoform|$domain|$n|$complete|$strand|$sense".
                
                    - $source: reference genome alias
                    - $gene: gene/isoform ID
                    - $feature: GFF3 feature type
                    - $isoform: isoform ID if by_gene = False, else same as $gene
                    - $domain: Pssm-Id
                    - $n: if multiple domains are present, they will be numbered according to 
                      proximity to 5' of sense strand
                    - $complete: 'complete' if ``complete=True`` else 'stitched'
                    - $strand: 'minus' if ``adj_dir=False`` and not self.plus else 'plus'
                    - $sense: 'NA' if self.strand is not set
                      else 'antisense' if ``adj_dir=False`` and not self.plus else 'sense'
            
            apply_template_to_dict (bool): flatten output dict and assign sequence names
        
        Returns
        -------
        dict of dict of Bio.Seq.Seq
            If ``apply_template_to_dict=False``. 
            Format if ``pssm_id=None``: {self.id: {(0, self.length): <sequence>}}
            Format if pssm_id is provided:
            {<isoform_id>: {(start pos of domain in self, end pos of domain in self): <sequence>}}
        dict of Bio.Seq.Seq
            If ``apply_template_to_dict=True``
        """
        ## get relevant GFF entries
        def invalid_feature_output():
            if fout: empty_file(fout)
            else: return {}
        ## terminate if feature doesn't exist or feature has no data
        try:
            ann = self.annotated_feature(feature_id)
            if not ann.length:
                return invalid_feature_output()
        except InvalidFeatureID:
            return invalid_feature_output()
        ## get sequences
        seqs = ann.get_seq(adj_dir = adj_dir, complete = complete, translate = translate,
                           feature_type = feature_type,
                           pssm_id = pssm_id, by_gene = by_gene, db = db, remote_rps = remote_rps,
                           mktmp = mktmp, rpsblast = rpsblast, rps_hits = rps_hits, quiet = quiet,
                           fout_gff = fout_gff,
                           seqid_template = seqid_template,
                           apply_template_to_dict = (bool(fout) or apply_template_to_dict),
                           thread = thread)
        ## output
        if fout:
            ## write
            dict_to_fasta(seqs, fout)
        else:
            return seqs
    
    def get_feature_seq(self, *feature_ids, fout=None,                        
                        ## shared options
                        feature_type=None, complete=False, adj_dir=False, translate=False, by_gene=False,
                        ## get_domain options
                        pssm_id=None, db=None, remote_rps=False, rpsblast="rpsblast", rps_hits=None, thread=None,
                        mktmp=None, fout_gff=None, quiet=True,
                        seqid_template = "$source|$gene|$feature|$isoform|$domain|$n|$complete|$strand|$sense",
                        apply_template_to_dict = False) -> dict:
        """
        Get sequences of multiple features.
        Except for ``feature_ids``, all arguments are strictly speaking optional.
        If no arguments are specified, the sense strand of the specified features will be returned.
        For more complicated queries, 
        in addition to modifying flags ``adj_dir``, ``complete``, and/or ``translate``,
        users may optionally specify a feature type and/or Pssm-Id.
        
        Arguments:
            fout (str): path to output file; if not provided, dictionary of seqs will be returned
            feature_type (str): GFF feature type (e.g. "CDS")
            adj_dir (bool): output sense strand
            complete (bool): merge output range(s) together into single range.
                Output will stll be a list of tuple. (e.g. [(<smallest start>, <largest end>)])
            translate (bool): translate sequence.
                For biologically sound results, ``translate=True`` only be used with
                ``feature_type="CDS"``, ``adj_dir=True``, and ``complete=False``.
            pssm_id (str or list): Pssm-Id of domain(s) to restrict output sequence range to.
                If provided and ``rps_hits=None``, 
                RPS-BLAST will be conducted on concatenated CDS to search for domain(s).
            by_gene (bool): merge overlapping domain ranges from different isoforms
            db (str): path to RPS-BLAST database
            remote_rps (bool): use remote RPS-BLAST; if True, ``db`` is not required
            rpsblast (str): path to rpsblast or rpsblast+ executable
                (or name of command, if available at command line)
            rps_hits (str): path to RPS-BLAST output, outfmt 6, with header. 
                If provided, ``db`` and ``rpsblast`` are not required.
            thread (int): number of threads for RPS-BLAST
            mktmp (func): function that takes a file name and generates an absolute path, 
                used for temporary files.
                ``tempfile.mkstemp`` will be used if not provided.
            quiet (bool): print only essential messages
            fout_gff (str): path to output GFF3 file in which to write domain features, entirely optional
            seqid_template (str): optional, template for output sequence name.
                Template will be parsed by strings.Template.
                The default template is "$source|$gene|$feature|$isoform|$domain|$n|$complete|$strand|$sense".
                
                    - $source: reference genome alias
                    - $gene: gene/isoform ID
                    - $feature: GFF3 feature type
                    - $isoform: isoform ID if by_gene = False, else same as $gene
                    - $domain: Pssm-Id
                    - $n: if multiple domains are present, they will be numbered according to 
                      proximity to 5' of sense strand
                    - $complete: 'complete' if ``complete=True`` else 'stitched'
                    - $strand: 'minus' if ``adj_dir=False`` and not self.plus else 'plus'
                    - $sense: 'NA' if self.strand is not set
                      else 'antisense' if ``adj_dir=False`` and not self.plus else 'sense'
            
            apply_template_to_dict (bool): flatten output dict and assign sequence names
        
        Returns
        -------
        dict of dict of Bio.Seq.Seq
            If ``apply_template_to_dict=False``. 
            Format if ``pssm_id=None``: {self.id: {(0, self.length): <sequence>}}
            Format if pssm_id is provided:
            {<isoform_id>: {(start pos of domain in self, end pos of domain in self): <sequence>}}
        dict of Bio.Seq.Seq
            If ``apply_template_to_dict=True``
        """
        ## get all seqs
        seqs = {feature_id: self._get_feature_seq(feature_id, feature_type = feature_type, pssm_id = pssm_id,
                                                  complete = complete, adj_dir = adj_dir, translate = translate,
                                                  by_gene = by_gene, db = db, remote_rps = remote_rps,
                                                  rpsblast = rpsblast, rps_hits = rps_hits, thread = thread,
                                                  mktmp = mktmp, fout_gff = fout_gff,
                                                  quiet = quiet,
                                                  seqid_template = seqid_template,
                                                  apply_template_to_dict = (bool(fout) or
                                                                            apply_template_to_dict))
                for feature_id in feature_ids}
        ## output
        if fout:
            seqs_flat = {seqid: seq
                         for feature_id, feature_seqs in seqs.items()
                         for seqid, seq in feature_seqs.items()}
            ## write blank file if no seqs
            if not seqs_flat:
                empty_file(fout)
            else:
                ## write
                dict_to_fasta(seqs_flat, fout)
        else:
            return seqs
