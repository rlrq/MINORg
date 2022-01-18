import sys
sys.path.append("/mnt/chaelab/rachelle/scripts/minorgpy")

import Bio
import copy
import tempfile

from Bio import Seq

from minorg.exceptions import MessageError

from minorg.functions import (
    BlastNR,
    filter_rpsblast_for_domain,
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

from minorg.fasta import (
    dict_to_fasta,
    collapse_identical_seqs
)

from minorg.annotation import GFF, subset_ann
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
# ref.get_feature_seq(gid, feature = "CDS", pssmid = pssm, complete = False, adj_dir = True,
#                     translate = False, mktmp = mktmp)
# ref.get_feature_seq(gid2, feature = "CDS", pssmid = pssm, complete = False, adj_dir = True,
#                     translate = False, mktmp = mktmp, db = rdb, rpsblast = rblast)

class InvalidFeatureID(MessageError):
    def __init__(self, feature_id):
        super().__init__( f"Error: Invalid feature ID. {feature_id} does not exist." )

#######################
##  REFERENCE_MANIP  ##
#######################

class AnnotatedFeature(GFF):
    """
    Representation of GFF3 feature annotation
    """
    def __init__(self, id, annotated_fasta, gff = None, rps_hits = None):
        """
        Creat an AnnotatedFeature object.
        
        Arguments:
            id (str): required, ID of feature (typically gene name or isoform name)
            annotated_fasta (AnnotatedFasta): required, :class:`AnnotatedFasta` object
            gff (str or GFF): optional, path to GFF file or :class:`minorg.annotation.GFF` object
            rps_hits (str): optional, path to RPS-BLAST output (outfmt 6, with header)
        """
        self.annotated_fasta = annotated_fasta
        self.annotation = (self.annotated_fasta.annotation if gff is None else
                           gff if isinstance(gff, GFF) else GFF(gff))
        self.assembly = self.annotated_fasta.assembly
        self.rps_hits = rps_hits
        self.id = id
        ## subset original GFF data & copy attributes from original GFF obj
        data = self.annotation.get_features_and_subfeatures(self.id, full = True, index = False)
        if not data:
            raise InvalidFeatureID(id)
        super().__init__(data = data)
        self.annotation.empty_copy(other = self)
        ## get feature-specific attributes
        self.ann = self.get_id(self.id, output_list = False)
        self.chrom = self.ann.molecule
        self.start = self.ann.start ## gff is 1-index
        self.start0 = self.start - 1 ## python's 0-index is start-inclusive
        self.end = self.ann.end ## gff's 1-index is end-inclusive
        self.end0 = self.end ## python's 0-index is end-exclusive
        self.feature = self.ann.feature
        self.strand = self.ann.strand
        self.plus = self.strand == '+'
        self.length = self.end0 - self.start0
    
    def convert_range(self, ranges, index = 0, start_incl = True, end_incl = False):
        '''
        Converts from 0-index [start, end)
        '''
        return convert_range(ranges, index_in = 0, index_out = index,
                             start_incl_in = True, start_incl_out = start_incl,
                             end_incl_in = False, end_incl_out = end_incl)
    
    def range(self, **conversion_kwargs):
        '''
        Returns genomic range of feature (not adusted for sense)
        '''
        return convert_range([(self.start0, self.end)], **conversion_kwargs)[0]
    
    def subset(self, id):
        return AnnotatedFeature(id = id, annotated_fasta = self.annotated_fasta, gff = self)
    
    def range_from_genomic(self, ranges, **conversion_kwargs):
        '''
        outputs ranges relative to sense strand of feature, where first base on sense strand of feature is 0
        ranges: [(start1, end1), (start2, end2), ...]
          - genomic coordinates
          - 0-indexed; start-inclusive, end-exclusive
        [[output options]]
        index: index
        start_incl: start inclusive
        end_incl: end inclusive
        '''
        if self.plus:
            ranges = [(start - self.start0, end - self.start0) for start, end in ranges]
        else:
            ranges = [(abs(end - self.end), abs(start - self.end)) for start, end in ranges]
        return self.convert_range(sorted(ranges), **conversion_kwargs)
    
    def range_to_genomic(self, ranges, **conversion_kwargs):
        '''
        ranges: [(start1, end1), (start2, end2), ...]
          - relative to sense strand of feature, where first base on sense strand of feature is 0
          - 0-indexed; start-inclusive, end-exclusive
        [[output options]]
        index: index
        start_incl: start inclusive
        end_incl: end inclusive
        '''
        if is_range(ranges):
            start, end = ranges
            if self.plus: r = (start + self.start0, end + self.start0)
            else: r = (abs(end - self.end), abs(start - self.end))
            return self.convert_range([r], **conversion_kwargs)[0]
        else:
            return [self.range_to_genomic(r, **conversion_kwargs) for r in ranges]
    
    def feature_range(self, feature, union = True, genomic = False):
        '''
        Returns ranges of entries of a given feature type (e.g. "CDS")
        '''
        feature_ranges = [(ann.start0, ann.end) for ann in self.get_features(feature, index = False)]
        if not genomic:
            feature_ranges = self.range_from_genomic(feature_ranges)
        if union:
            return ranges_union([feature_ranges])
        return feature_ranges
    
    def subfeature_range(self, feature_id, feature = None, genomic = False, union = True):
        '''
        Returns range of given feature type (e.g. "CDS") of a sub-feature (e.g. an mRNA isoform of gene)
        '''
        subfeature_ann = self.subset(feature_id)
        if feature is None: feature = subfeature_ann.feature
        output = subfeature_ann.feature_range(feature, union = union, genomic = genomic)
        if not genomic:
            offset = (subfeature_ann.start0 - self.start0) if self.plus else (self.end - subfeature_ann.end)
            output = [(start + offset, end + offset) for start, end in output]
        return output
    
    def adj_feature_range(self, feature, ranges, complete = True, **feature_range_kwargs):
        feature_bounds = self.feature_range(feature, **feature_range_kwargs)
        feature_pos = ranges_to_pos(feature_bounds)
        ## coopt adjusted_ranges (meant to account for gaps in alignment)
        feature_dummy_seq = ''.join([('A' if p in feature_pos else '-')
                                     for p in range(max(feature_pos))])
        return adjusted_ranges(feature_dummy_seq, *ranges, subtract_gaps = (not complete))
    
    def adj_subfeature_range(self, feature_id, ranges, complete = True, **subfeature_range_kwargs):
        '''
        Convert range in a subfeature to coords in current feature
        '''
        subfeature_bounds = self.subfeature_range(feature_id, genomic = False, **subfeature_range_kwargs)
        subfeature_pos = ranges_to_pos(subfeature_bounds)
        ## coopt adjusted_ranges (meant to account for gaps in alignment)
        subfeature_dummy_seq = ''.join([('A' if p in subfeature_pos else '-')
                                        for p in range(max(subfeature_pos))])
        return adjusted_ranges(subfeature_dummy_seq, *ranges, subtract_gaps = (not complete))
    
    def range_in_genomic_coords(self, ranges, index=0, start_incl=True, end_incl=False):
        '''
        ranges: [(start1, end1), (start2, end2), ...]
          - relative to sense strand of feature, where first base on sense strand of feature is 0
          - 0-indexed; start-inclusive, end-exclusive
        [[output options]]
        index: index
        start_incl: start inclusive
        end_incl: end inclusive
        '''
        if self.plus: ranges = [(start + self.start0, end + self.start0) for start, end in ranges]
        else: ranges = [(self.length - end, self.length - start) for start, end in ranges]
        if not start_incl: ranges = [(start - 1, end) for start, end in ranges]
        if end_incl: ranges = [(start, end - 1) for start, end in ranges]
        return [(start + index, end + index) for start, end in ranges]
    
    def _get_seq(self, feature=None, ranges=[], adj_dir=True, complete=False, translate=False):
        '''
        feature: GFF feature type (e.g. CDS)
        ranges: 0-indexed ranges, relative to sense strand of AnnotatedFeature obj
        adj_dir: return sense strand
        final range from which sequence will be extracted is the intersection of ranges and feature ranges
        '''
        if feature is None: feature = self.feature
        feature_ranges = self.feature_range(feature, union = True, genomic = True)
        if ranges:
            output_ranges = sorted(ranges_intersect(feature_ranges, self.range_to_genomic(ranges)))
        else:
            output_ranges = sorted(feature_ranges)
        if complete:
            output_ranges = [(min(r[0] for r in output_ranges), max(r[1] for r in output_ranges))]
        ## start extracting seq
        seq = Bio.Seq.Seq('')
        for start, end in output_ranges:
            seq += self.assembly[self.chrom][start:end]
        if adj_dir and not self.plus:
            seq = seq.reverse_complement()
        if translate:
            seq = seq.translate(table = self.annotated_fasta.genetic_code)
        return seq
    
    def get_seq(self,
                ## shared options
                fout=None, adj_dir=True, complete=True, translate=False,
                ## _get_seq options
                feature=None, ranges=[],
                ## _get_domain_seq options
                pssm_id=None, by_gene=False, db=None, rpsblast="rpsblast", rps_hits=None, thread=None,
                mktmp=None, quiet=False, fout_gff=None,
                ## seqid options
                seqid_template = "$source|$gene|$feature|$isoform|$domain|$n|$complete|$revcomp",
                apply_template_to_dict = False):
        '''
        fout: path to output file; if not provided, dictionary of seqs will be returned
        '''
        feature = feature if feature is not None else self.feature
        ## select appropriate get_seq function
        ## seqs fmt: {id: {(<start>, <end>): <seq>}}
        if pssm_id and pssm_id != "gene":
            if isinstance(pssm_id, str): pssm_ids = pssm_id.split(',')
            elif isinstance(pssm_id, int): pssm_ids = [pssm_id]
            else: pssm_ids = pssm_id
            seqs = self._get_domain_seq(*pssm_ids, by_gene = by_gene, feature = feature,
                                        adj_dir = adj_dir, complete = complete, translate = translate,
                                        db = db, rpsblast = rpsblast, rps_hits = rps_hits, thread = thread,
                                        fout_gff = fout_gff, mktmp = mktmp, quiet = quiet)
        else:
            seq = self._get_seq(adj_dir = adj_dir, complete = complete, translate = translate,
                                feature = feature, ranges = ranges)
            seqs = {self.id: {(0, self.length): seq}}
        ## format seqid
        if seqid_template and (fout or apply_template_to_dict):
            ## create function to make seqid
            from string import Template
            def mk_seqid(**kwargs):
                return Template(seqid_template).substitute(source = self.annotated_fasta.name,
                                                           gene = self.id, feature = feature,
                                                           complete = ("complete" if complete else "stitched"),
                                                           domain = (pssm_id if pssm_id else '.'),
                                                           revcomp = ("revcomp" if (adj_dir and not self.plus)
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
    
    def _get_domain_seq(self, *pssm_ids, by_gene=False,
                        ## get_seq options
                        feature=None, adj_dir=True, complete=False, translate=False,
                        ## get_domain_ranges options
                        db=None, rpsblast="rpsblast", rps_hits=None, thread=None,
                        mktmp=None, quiet=False, fout_gff=None):
        '''
        Executes get_seq for each 
        '''
        ## get domain ranges (complete) (relative to each isoform)
        domain_ranges = self._get_domain_ranges(*pssm_ids, db = db, rpsblast = rpsblast,
                                                rps_hits = rps_hits, thread = thread,
                                                by_gene = by_gene, genomic = False,
                                                mktmp = mktmp, quiet = quiet, fout_gff = fout_gff)
        ## get seqs; discontiguous domains are NOT concatenated
        domain_seqs = {isoform: {tuple(hit_range): self.subset(isoform)._get_seq(feature, ranges = [hit_range],
                                                                                 complete = complete,
                                                                                 translate = translate)
                                 for hit_range in hit_ranges}
                       for isoform, hit_ranges in domain_ranges.items()}
        return domain_seqs ## {isoform: {<range of first domain>: <seq>, <range of second domain>: <seq>}...}
    
    def _get_domain_ranges(self, *pssm_ids, db=None, rpsblast="rpsblast", rps_hits=None, thread=None,
                           genomic=False, by_gene=False, mktmp=None, quiet=False, fout_gff=None):
        '''
        If self contains mRNA w/ CDS, CDS of each mRNA will be combined and translated,
            searched for the relevant domain (pssm_id: <Pssm_id>) against a domain database (db: <path to db>),
            and coordinates will be mapped back to the sense strand of self
        If rps_hits not provided AND self.rps_hits is None, rpsblast will be executed to generate rps_hits file
            - rps_hits file is basically blast6 with header
            - required fields: qseqid, sseqid. qseqid must be mRNA ID.
        '''
        ## some functions
        mktmpfname = ((lambda suf: tempfile.mkstemp()[1]) if mktmp is None
                      else (lambda suf: mktmp(f"{self.ann.get_attr('ID',fmt=str)}_{suf}")))
        def printi(msg):
            if not quiet: print(msg)
        ## if rps_hits doesn't alreay exist, execute rps-blast
        if rps_hits is None and self.rps_hits is None:
            ## get isoform sequences
            printi("Extracting peptide sequence(s)")
            isoforms = {entry.get_attr("ID", fmt=list)[0]: self.subset(id=entry.get_attr("ID", fmt=list)[0])
                        for entry in self if entry.feature == "mRNA"}
            isoform_seqs = {isoform: isoform_ann._get_seq(feature = "CDS", adj_dir = True, translate = True)
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
                    db = db, cmd = rpsblast,
                    mt_mode = (1 if thread is not None else 0),
                    num_threads = (1 if thread is None else thread))
            ## store path to output
            self.rps_hits = tmp_rpsblast
        ## filter hits for desired domain
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
                start, end = self.range_to_genomic(r, index = 1, end_incl = True)
                domain_ann = copy.deepcopy(self.ann)
                domain_ann.start = start
                domain_ann.end = end
                domain_ann.feature = "domain"
                domain_ann.source = "minorg"
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
    """
    
    def __init__(self, fasta, gff, name = None, attr_mod = {}, genetic_code = 1,
                 memsave = True):
        '''
        fasta: FASTA file of genome assembly
        gff: GFF3 file of genome annotation
        name [default = None]: name of genome
        attr_mod [default = {}]: attribute modifications
           (fmt: {'<feature>': {'<standard attribute field name>': '<nonstandard attribute field name used>'}})
           (e.g.: {'mRNA': {'Parent': 'Locus_id'}})
        genetic_code [default = 1]: NCBI genetic code
           (see: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
           (1 --> Standard Code; 2 --> Vertebrate Mitochondrial Code; etc.)
        '''
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
        
    def parse_gff(self):
        self.annotation = GFF(self.gff, attr_mod = self.attr_mod, memsave = self.memsave)
    
    def subset_annotation(self, ids, fout = None, memsave = None, sort = True):
        if fout is None:
            fout = tempfile.mkstemp()[1]
        subset = subset_ann(gff_beds = {self.name: self._gff}, ids = ids, fout_fmt = "GFF",
                            mk_tmpf_name = lambda x: fout, attr_mods = {self.name: self.attr_mod},
                            memsave = (memsave if memsave is not None else self.memsave),
                            sort = sort)
        self.gff = fout
        self.parse_gff()
        return
    
    def reduce_annotation(self, ids, fout = None, sort = True):
        return self.subset_annotation(ids, fout = fout, sort = sort)
    
    def feature_range(self, feature_id):
        ann = self.gff.get_id(feature_id, output_list = False)
        return (ann.start - 1, ann.end)
    
    def feature_pos(self, feature_id):
        return ranges_to_pos([self.feature_range(feature_id)])
    
    def annotated_feature(self, feature_id):
        annotated_feature = self._annotated_features.get(feature_id, None)
        if annotated_feature is None:
            annotated_feature = AnnotatedFeature(id = feature_id, annotated_fasta = self)
            self._annotated_features[feature_id] = annotated_feature
        return annotated_feature
    
    def get_seq_domain(self, feature_id, feature, domain, **kwargs):
        pass
    
    def _get_feature_seq(self, feature_id, feature=None, fout=None, pssmid='',
                         ## get_seq options
                         complete=False, adj_dir=False, translate=False, by_gene=False,
                         ## get_domain options
                         db=None, rpsblast="rpsblast", rps_hits=None, thread=None,
                         mktmp=None, fout_gff=None, quiet=True,
                         seqid_template = "$source|$gene|$feature|$isoform|$domain|$n|$complete|$revcomp",
                         apply_template_to_dict = False):
        ## get relevant GFF entries
        def invalid_feature_output():
            if fout: empty_file(fout)
            else: return {}
        ## terminate if feature doesn't exist or feature has no data
        try:
            ann = self.annotated_feature(feature_id)
            if not len(ann):
                return invalid_feature_output()
        except InvalidFeatureID:
            return invalid_feature_output()
        ## get sequences
        seqs = ann.get_seq(adj_dir = adj_dir, complete = complete, translate = translate,
                           feature = feature,
                           pssm_id = pssmid, by_gene = by_gene, db = db, mktmp = mktmp,
                           rpsblast = rpsblast, rps_hits = rps_hits, quiet = quiet,
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
    
    def get_feature_seq(self, *feature_id, feature=None, fout=None, pssmid='',
                        ## get_seq options
                        complete=False, adj_dir=False, translate=False, by_gene=False,
                        ## get_domain options
                        db=None, rpsblast="rpsblast", rps_hits=None, thread=None,
                        mktmp=None, fout_gff=None, quiet=True,
                        seqid_template = "$source|$gene|$feature|$isoform|$domain|$n|$complete|$revcomp",
                        apply_template_to_dict = False):
        ## get all seqs
        seqs = {feature_id: self._get_feature_seq(feature_id, feature = feature, pssmid = pssmid,
                                                  complete = complete, adj_dir = adj_dir, translate = translate,
                                                  by_gene = by_gene, db = db, rpsblast = rpsblast,
                                                  rps_hits = rps_hits, thread = thread,
                                                  mktmp = mktmp, fout_gff = fout_gff,
                                                  quiet = quiet,
                                                  seqid_template = seqid_template,
                                                  apply_template_to_dict = (bool(fout) or
                                                                            apply_template_to_dict))
                for feature_id in feature_id}
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
