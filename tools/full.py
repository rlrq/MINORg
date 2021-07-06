import os
import sys
import argparse

ORIGINAL_DIR = os.path.abspath(os.getcwd())
SCRIPT_DIR = os.path.abspath(os.path.realpath(__file__))
SCRIPTS_DIR = os.path.join(os.path.split(SCRIPT_DIR)[0], "scripts")
TOOLS_DIR = os.path.join(SCRIPT_DIR, "tools")

raw_args = sys.argv[1:]

parser = argparse.ArgumentParser(prog = "MINORg full",
                                 description = "MINORg: Identify targets and get minimum gRNA set.")
## housekeeping
parser.add_argument("--original-dir", type = str, default = ORIGINAL_DIR)
parser.add_argument("--scripts-dir", type = str, default = SCRIPTS_DIR)
parser.add_argument("--tools-dir", type = str, default = TOOLS_DIR)
## mutually exclusive vars
annotation_parser = parser.add_mutually_exclusive_group()
annotation_parser.add_argument("--bed", type = str)
annotation_parser.add_argument("--gff", type = str)
accession_parser = parser.add_mutually_exclusive_group()
accession_parser.add_argument("-a", "--acc", type = str)
accession_parser.add_argument("-i", "--input", type = str, dest = "accs_f")
target_parser = parser.add_mutually_exclusive_group()
target_parser.add_argument("-g", "--gene", type = str)
target_parser.add_argument("-c", "--cluster", type = str)
target_parser.add_argument("-t", "--target", type = str)
query_parser = parser.add_mutually_exclusive_group()
query_parser.add_argument("-f", "--fasta", type = str)
query_parser.add_argument("-q", "--query", type = str) ## query preferred
## vars
parser.add_argument("--extend-cds", type = str)
parser.add_argument("--extend-genome", type = str)
parser.add_argument("-d", "--dir", type = str)
parser.add_argument("-r", "--reference", type = str)
parser.add_argument("--prefix", type = str)
parser.add_argument("--domain", type = str)
parser.add_argument("--db", type = str)
parser.add_argument("--cluster-lookup", type = str)
parser.add_argument("--cluster-dir", type = str)
# parser.add_argument("--target-db", type = str)
## flags
parser.add_argument("--merge-redundant", action = "store_true", default = False)
parser.add_argument("--check-id-before-merge", action = "store_true", default = False)
parser.add_argument("--relax", action = "store_true", default = False)
parser.add_argument("--auto", action = "store_true", default = False)
parser.add_argument("--check-recip", "--check-reciprocal", action = "store_true", default = False)
parser.add_argument("--screen-ref", action = "store_true", default = False)
parser.add_argument("--unmask-ref", action = "store_false", default = True, dest = "mask_ref")
parser.add_argument("--skip-bg-check", action = "store_false", default = True, dest = "check_bg")

parser.add_argument("--sep", type = str, default = '.')
parser.add_argument("-v", "--verbose", action = "store_true", default = False)
## arguments to pass to other scripts (that we had to read here because of logging reasons)
parser.add_argument("-p", "--pam", type = str, dest = "PAM")
parser.add_argument("-l", "--length", type = int, dest = "LENGTH")
parser.add_argument("-b", "--background", type = str)
parser.add_argument("-e", "--exclude", type = str)
parser.add_argument("-s", "--sets", type = int)
parser.add_argument("--attr-mod", type = dict)
parser.add_argument("--mismatch", type = int)
parser.add_argument("--gap", type = int)
parser.add_argument("--minlen", type = int)
parser.add_argument("--minid", type = float)
parser.add_argument("--merge-within", type = int, default = 100)
parser.add_argument("--sc-algorithm", type = str)
## things to print
parser.add_argument("-h", "--help", action = "store_true", default = False, dest = "print_man")
parser.add_argument("--readme", action = "store_true", default = False, dest = "print_readme")
parser.add_argument("--aliases", action = "store_true", default = False, dest = "print_aliases")
parser.add_argument("--members", action = "store_true", default = False, dest = "print_members")
parser.add_argument("--accessions", action = "store_true", default = False, dest = "print_accessions")


args, unknown = parser.parse_known_args(sys.argv[1:])

# subprocess.check_output(("python3", os.path.join(TOOLS_DIR, (args.tool + ".py")), *unknown,
#                          "--original-dir", ORIGINAL_DIR, "--scripts-dir", SCRIPTS_DIR,
#                          "--tools-dir", TOOLS_DIR))


HEADER_CDD = "qseqid,sseqid,pident,length,qstart,qend"
HEADER_AL70 = "qseqid,sseqid,pident,length,sstart,send"
# DB_CDD_v3_17 = '/mnt/chaelab/shared/blastdb/RPSdb/Cdd'
DB_CDD_v3_18 = "/mnt/chaelab/shared/blastdb/Cdd.v3.18/Cdd"
DB_CDD = DB_CDD_v3_18
DB_AL70 = "/mnt/chaelab/shared/blastdb/anna_lena70.contigs/anna_lena70.contigs"
RPS_DB_DEFAULT = DB_CDD
TARGET_DB_DEFAULT = DB_AL70
CLUSTER_LOOKUP_DEFAULT = "/mnt/chaelab/rachelle/data/NLR/cluster_aliases.txt"
CLUSTER_DIR_DEFAULT = "/mnt/chaelab/rachelle/data/NLR/clusters/cluster_combined"
DIR_DEFAULT = "/mnt/chaelab/$(whoami)/find_gRNA"
DOMAIN_DEFAULT = "gene"
QUERY_DEFAULT = "/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta"
BACKGROUND_DEFAULT = QUERY_DEFAULT
REFERENCE_DEFAULT = "/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta"
BED_DEFAULT = "/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed"
SC_ALGO_DEFAULT = "LAR"
PREFIX_DEFAULT = "findgRNA"
PAM_DEFAULT = "GG"
MISMATCH_DEFAULT = 0
GAP_DEFAULT = 0
MINID_DEFAULT = 95
MINLEN_DEFAULT = 0
MERGE_WITHIN_DEFAULT = 100
LENGTH_DEFAULT = 20
SEP_DEFAULT = '.'
SETS_DEFAULT = 1
