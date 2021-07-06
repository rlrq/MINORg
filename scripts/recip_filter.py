import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from data_manip import make_custom_get
from fasta_manip import fasta_to_dict, dict_to_fasta

def set_overlap(a, b):
    for x in a:
        if x in b:
            return True
    return False

## if relax == True, candidate targets with equal bitscore to a target gene and a non-target gene are retained
## elif relax == False, candidate targets are only retained if all max bitscore genes are subset of target genes
def remove_non_max_bitscore(fasta, bed_f, genes, relax = False,
                            colnames = ["chrom", "start", "end", "candidate", "cstart", "cend", "bitscore",
                                        "bed_chrom", "bed_start", "bed_end", "id", "score", "strand", "source",
                                        "feature", "phase", "attributes", "overlap"]):
    genes = set(genes.split(','))
    get = make_custom_get(colnames)
    with open(bed_f, 'r') as f:
        data = {}
        for entry in (x[:-1].split('\t') for x in f.readlines()
                      if get(x[:-1].split('\t'), "feature") in ("gene", "pseudogene")):
            data[get(entry, "candidate")] = data.get(get(entry, "candidate"), []) + [entry]
    ## get largest bitscore for each candidate target
    max_bitscore = {candidate: max(get(data[candidate], "bitscore")) for candidate in data}
    ## identify sequences to discard and print warnings if candidate has max bitscore with target and non-target
    throw = []
    for candidate in data:
        max_bitscore_genes = set(get(entry, "id") for entry in data[candidate]
                                 if get(entry, "bitscore") == max_bitscore[candidate])
        if max_bitscore_genes.issubset(genes): ## if max score genes are subset of target genes
            continue
        else:
            if max_bitscore_genes.isdisjoint(genes): ## if no target genes have max score
                throw.append(candidate)
            else: ## if overlapping but not subset
                if relax:
                    print(f"Warning: candidate target '{candidate}' has hit(s) with bitscore {max_bitscore[candidate]} that overlap(s) with target gene(s) {genes & max_bitscore_genes} and non-target gene(s) {max_bitscore_genes - genes}. This sequence will be retained as 'relax' has been set to True.")
                else:
                    throw.append(candidate)
                    print(f"Warning: candidate target '{candidate}' has hit(s) with bitscore {max_bitscore[candidate]} that overlap(s) with target gene(s) {genes & max_bitscore_genes} and non-target gene(s) {max_bitscore_genes - genes}. This sequence will be removed from the list of candidate targets as 'relax' has been set to False.")
    ## read original candidate targets, filter, and write
    seqs = fasta_to_dict(fasta)
    dict_to_fasta({seq_id: seq for seq_id, seq in seqs.items() if seq_id not in throw}, fasta)
    return
    
