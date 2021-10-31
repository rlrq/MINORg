import os
import re
import itertools
from exceptions import InvalidPath
from functions import fasta_to_dict, dict_to_fasta
from mafftcommandline_add import MafftCommandline

def get_recursively(d, default, *keys):
    def helper(d, keys):
        key = keys[0]
        if key in d:
            if len(keys) == 1: return d[key]
            else: return helper(d[key], keys[1:])
        else:
            return default
    return helper(d, keys)

def group_by_gene(fa_cds, fa_genomic, directory, sep = '.', verbose = True, logger = None):
    seqs_cds = {seqid: seq for fname in fa_cds for seqid, seq in fasta_to_dict(fname).items()}
    seqs_genomic = {seqid: seq for fname in fa_genomic for seqid, seq in fasta_to_dict(fname).items()}
    genes = {seqid: [] for seqid in seqs_genomic}
    orphan_cds = []
    ## map CDS to GENOMIC
    for seqid in seqs_cds.keys():
        gene = re.match(f"^.+?(?={re.escape(sep)}[^{sep}]+$)", seqid).group(0)
        if not gene in genes:
            orphan_cds.append(seqid)
        else:
            genes[gene].append(seqid)
    ## print orphans
    orphan_genes = sorted(gene for gene, cds in genes.items() if not cds)
    def log(msg):
        if verbose and logger: logger.plain(msg)
        elif logger: logger.fplain(msg)
        elif verbose: print(msg)
    if orphan_genes:
        log("GENOMIC sequences without CDS:", ','.join(orphan_genes))
    if orphan_cds:
        log("CDS sequences without GENOMIC:", ','.join(orphan_cds))
    ## write
    for gene, cds in genes.items():
        to_write = {gene: seqs_genomic[gene], **{seqid_cds: seqs_cds[seqid_cds] for seqid_cds in cds}}
        dict_to_fasta(to_write, f"{directory}/{gene}.fasta")
    return

def aln_to_annotation(directory, fout, sep = '.', outfmt = "GFF", attr_mod = None, verbose = True):
    fnames = [f"{directory}/{fname}" for fname in os.listdir(directory)]
    entries = []
    attr_mod = attr_mod if isinstance(attr_mod, dict) else {} if not attr_mod else eval(attr_mod)
    for fname in fnames:
        entries.extend(extract_annotation(fname, sep = sep, attr_mod = attr_mod))
    open(fout, 'w+').write('\n'.join([entry.generate_str(fmt = outfmt) for entry in entries]) + '\n')
    return

def extract_annotation(fa_aln, sep = '.', attr_mod = {}):
    seqs = fasta_to_dict(fa_aln)
    seqid_ref = min(seqs.keys(), key = lambda s: s.count(sep))
    seqid_cds = [seqid for seqid in seqs if seqid != seqid_ref]
    ## make gene annotation
    gene_annotation = UserAnnotation(seqid_ref, seqid_ref, "gene", attr_mod = attr_mod)
    gene_annotation.set_range(1, len([c for c in seqs[seqid_ref] if c != '-']))
    features = [gene_annotation]
    ## make mRNA and CDS annotations
    for seqid in seqid_cds:
        cds_annotations = []
        seq = seqs[seqid]
        feature_count = 1
        start = None
        ## make CDS annotations
        for pos, base in enumerate(seq):
            if start == None and base != '-':
                start = pos
            elif start != None and (pos == (len(seq) - 1) or seq[pos+1] == '-'):
                new_annotation = UserAnnotation(f"{seqid}.cds{str(feature_count).zfill(6)}",
                                                seqid_ref, "CDS", parent = seqid)
                new_annotation.set_range(start + 1, pos + 1)
                new_annotation.set_phase(seq[:start])
                cds_annotations.append(new_annotation)
                feature_count += 1
                start = None
            else:
                continue
        ## make mRNA annotation
        mRNA_annotation = UserAnnotation(seqid, seqid_ref, "mRNA", parent = seqid_ref, attr_mod = attr_mod)
        mRNA_annotation.set_range(min(ann._start for ann in cds_annotations),
                                  max(ann._end for ann in cds_annotations))
        features.extend([mRNA_annotation] + cds_annotations)
    return features

class UserAnnotation:
    def __init__(self, ID, molecule, feature, parent = None, attr_mod = {}):
        self._attr_mod = attr_mod
        self._ID = ID
        self._molecule = molecule
        self._feature = feature
        self._parent = parent
        self._phase = '.'
        self._score = '.'
        self._source = "user"
        self._start = None ## 1-index, incl
        self._end = None ## 1-index, incl
        self._strand = '+'
    def set_range(self, start, end):
        self._start = start
        self._end = end
    def set_phase(self, seq):
        self._phase = len([c for c in seq if c != '-']) % 3
    def generate_attr(self):
        def get_mod(field):
            feature_mod = get_recursively(self._attr_mod, field, self._feature, field)
            if feature_mod != field: return feature_mod
            else: return get_recursively(self._attr_mod, field, "all", field)
        return ';'.join([f"{k}={v}" for k, v in
                         [(get_mod("ID"), self._ID), (get_mod("Parent"), self._parent)] if v])
    def generate_str(self, fmt = "GFF"):
        if fmt.upper() in {"GFF", "GFF3"}:
            output = self.generate_gff()
        elif fmt.upper() in {"BED"}:
            output = self.generate_bed()
        return '\t'.join(map(str, output))
    def generate_gff(self):
        return [self._molecule, self._source, self._feature, self._start, self._end,
                self._score, self._strand, self._phase, self.generate_attr()]
    def generate_bed(self):
        return [self._molecule, self._start - 1, self._end, self._ID, self._score,
                self._strand, self._source, self._feature, self._phase, self.generate_attr()]

def extend_genome(args, config):
    if not (args.ext_genome and args.ext_cds):
        return
    ext_genome = [x for x in args.ext_genome if os.path.exists(x)]
    ext_cds = [x for x in args.ext_cds if os.path.exists(x)]
    ## raise Error for inaccessible files
    nonexistent_genome = [x for x in args.ext_genome if x not in ext_genome]
    nonexistent_cds = [x for x in args.ext_cds if x not in ext_cds]
    nonexistent = nonexistent_genome + nonexistent_cds
    if nonexistent:
        InvalidPath(','.join(nonexistent))
    ## separate sequences by gene
    aln_dir = config.mkdir("tmp_extend_aln")
    # config.mkdir(aln_dir)
    group_by_gene(ext_cds, ext_genome, aln_dir, sep = '.', verbose = True, logger = config.logfile)
    ## align
    for fasta in [f for f in os.listdir(aln_dir) if re.search("\.fasta$", f)]:
        gene = re.search("^.+(?=\.fasta$)", os.path.basename(fasta)).group(0)
        with open(os.path.join(aln_dir, f"{gene}_mafft.fa"), "w+") as f:
            stdout, stderr = MafftCommandline(args.mafft, input = os.path.join(aln_dir, fasta),
                                              quiet = True, thread = args.thread)()
            f.write(stdout)
        os.remove(os.path.join(aln_dir, fasta))
    ## convert alignment to gff
    aln_bed = config.mkfname("ext.bed", tmp = True)
    aln_to_annotation(aln_dir, fout = aln_bed, sep = '.', outfmt = "bed", attr_mod = args.attr_mod)
    ## combine genomic files
    ext_genome_fa = config.mkfname("ext.fasta", tmp = True)
    seqs_genome = dict(itertools.chain(*[fasta_to_dict(fa).items() for fa in ext_genome]))
    dict_to_fasta(seqs_genome, ext_genome_fa)
    ## update config
    config.extend_reference("Extended", ext_genome_fa, aln_bed)
    # args.bed = new_bed
    # args.reference = new_ref
    ## remove temporary directories
    import shutil
    shutil.rmtree(aln_dir)
    return
        
