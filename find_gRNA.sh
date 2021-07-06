#!/bin/bash

## ?? TODO: properly enable --query file to be used as background automatically without overwriting --background ## and --nonreference

## TODO: ensure that using --background doesn't override -a when --target is used (i.e. allow bg check in both user-provided background file and VdW's renseq if -b and -a used along with --target)

## TODO: update manual to reflect new background screening behaviour

ORIGINAL_DIR=$(pwd)
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
TOOLS_DIR="${SCRIPT_DIR}/tools"

## if no arguments, show manual
if [[ $# -eq 0 ]]; then
    man -l ${SCRIPT_DIR}/MANUAL_find_gRNA.2 ## TODO lol
    exit 1
fi

## minimum set
if [[ ${1} == 'minimumset' ]]; then
    ${TOOLS_DIR}/minimumset.sh ${@:2}
    exit 1
fi

## because aliases don't really work with non-interactive shell :(
get_seq=/mnt/chaelab/rachelle/scripts/get_seqs_generic/get_seqs.sh
aln2gff=/mnt/chaelab/rachelle/scripts/aln_to_gff/aln2gff.sh

## some built-in supported PSSMs
declare -A DOMAINS
DOMAINS["TIR"]='366714';DOMAINS["RX-CC_like"]='392282';DOMAINS["CC"]='392282';DOMAINS["RPW8"]='384063';DOMAINS["NB-ARC"]='391514';DOMAINS["NBS"]='391514';DOMAINS["Rx_N"]='375519'

HEADER_CDD="qseqid,sseqid,pident,length,qstart,qend"
HEADER_AL70="qseqid,sseqid,pident,length,sstart,send"
# DB_CDD_v3_17='/mnt/chaelab/shared/blastdb/RPSdb/Cdd'
DB_CDD_v3_18='/mnt/chaelab/shared/blastdb/Cdd.v3.18/Cdd'
DB_CDD=${DB_CDD_v3_18}
DB_AL70='/mnt/chaelab/shared/blastdb/anna_lena70.contigs/anna_lena70.contigs'
RPS_DB_DEFAULT=${DB_CDD}
TARGET_DB_DEFAULT=${DB_AL70}
CLUSTER_LOOKUP_DEFAULT='/mnt/chaelab/rachelle/data/NLR/cluster_aliases.txt'
CLUSTER_DIR_DEFAULT='/mnt/chaelab/rachelle/data/NLR/clusters/cluster_combined'
DIR_DEFAULT="/mnt/chaelab/$(whoami)/find_gRNA"
DOMAIN_DEFAULT="gene"
QUERY_DEFAULT='/mnt/chaelab/shared/anna_lena/anna_lena70.contigs.fasta'
BACKGROUND_DEFAULT="${QUERY_DEFAULT}"
REFERENCE_DEFAULT='/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta'
BED_DEFAULT='/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed'
SC_ALGO_DEFAULT='LAR'
PREFIX_DEFAULT='minorg'
PAM_DEFAULT='GG'
MISMATCH_DEFAULT=0
GAP_DEFAULT=0
MINID_DEFAULT=95
MINLEN_DEFAULT=0
MERGE_WITHIN_DEFAULT=100
LENGTH_DEFAULT=20
SEP_DEFAULT='.'
SETS_DEFAULT=1

params="$@"

while (( "$#" )); do
    case "$1" in
        ### -g|--gene and -c|--cluster can be either gene or cluster names (comma-separated), or paths to file containing newline (\n)-separated gene or cluster names
        -g|--gene) GENE="${2}";; ## mutually exclusive with -c|--cluster and -f|--fasta
        -c|--cluster) CLUSTER="${2}";; ## mutually exclusive with -g|--gene and -f|--fasta
        -t|--target) TARGET="${2}";; ## mutually exclusive with -g|--gene -c|--cluster
        -f|--fasta) FASTA="${2}";; ## -q|--query preferred
        -q|--query) QUERY="${2}";;
        -a|--acc) ACC="${2}";;
        -i|--input) ACCS_F="${2}";;
        -p|--pam) PAM="${2}";;
        -l|--length) LENGTH="${2}";;
        -d|--dir) DIR="${2}";;
        # -n|--nonreference) NONREFERENCE="${2}";; ## -q|--query preferred
        -r|--reference) REFERENCE="${2}";;
        -b|--background) BACKGROUND="${2}";;
        -e|--exclude) EXCLUDE="${2}";;
        -s|--sets) SETS="${2}";;
        --gff) GFF="${2}";;
        --bed) BED="${2}";;
        --prefix) PREFIX="${2}";;
        --domain) DOMAIN="${2}";;
        --db) RPS_DB="${2}";;
        --mismatch) MISMATCH="${2}";;
        --gap) GAP="${2}";;
        --minid) MINID="${2}";;
        --minlen) MINLEN="${2}";;
        --merge-within) MERGE_WITHIN="${2}";;
        --sc-algorithm) SC_ALGO="${2}";;
        --cluster-lookup) CLUSTER_LOOKUP="${2}";;
        --cluster-dir) CLUSTER_DIR="${2}";;
        --target-db) TARGET_DB="${2}";;
        --check-id-before-merge) CHECK_ID_PREMERGE='True';;
        --relax) RELAX='True';;
        --auto) AUTO='True';;
        --check-reciprocal|--check-recip) CHECK_RECIP='True';;
        --relax-reciprocal|--relax-recip) RELAX_RECIP='True'; CHECK_RECIP='True';;
        --screen-ref) SCREEN_REF='True';;
        --unmask-ref) MASK_REF='False'; UNMASK_REF='True';;
        --skip-bg-check) CHECK_BG='False';; ## overrides --screen-ref
        --merge-redundant) MERGE_REDUNDANT='True';;
        # --by-gene) BY_GENE='True';; ## TODO: change this to --merge-redundant, where raising this flag merges identical protein sequences, and should allow domain search time to be shorted significantly if there are multiple redundancies. Bypass --by-gene because this seems to combine any overlapping CDS ranges, which could lead to problems if some CDS are in a different frame from an overlapping one (anyway, --by-gene currently breaks cuz of conversion from domain aa range to domain nt range due to different fasta sequence naming convention of --by-gene sequences
        -h|--help) man -l ${SCRIPT_DIR}/MANUAL_find_gRNA.2; exit 0;;
        --readme) cat ${SCRIPT_DIR}/README; exit 0;;
        --aliases) cat /mnt/chaelab/rachelle/data/NLR/cluster_aliases.txt; exit 0;;
        --accessions) cat /mnt/chaelab/shared/anna_lena/accession_map.txt; exit 0;;
        --members) CHECK_MEMBERS="${2}";;
        --attr-mod) ATTR_MOD="${2}";; ## see get_seqs_generic.sh for explanation
        # --extend-gff) GFF_EXT="${2}";; ## TODO: implement this and --extend-bed
        # --extend-bed) BED_EXT="${2}";;
        --extend-cds) CDS_EXT="${2}";;
        --extend-genome) GENOME_EXT="${2}";;
        --sep) SEP="${2}";;
        --report-bg) REPORT_BG='True';;
        -v|--version) echo "MINORg v2.2"; exit 0;;
    esac
    shift
done

BED="${BED:-${BED_DEFAULT}}"
DIR="${DIR:-${DIR_DEFAULT}}"
CLUSTER_LOOKUP="${CLUSTER_LOOKUP:-${CLUSTER_LOOKUP_DEFAULT}}"
CLUSTER_DIR="${CLUSTER_DIR:-${CLUSTER_DIR_DEFAULT}}"
QUERY="${QUERY:-${FASTA}}"
QUERY="${QUERY:-${QUERY_DEFAULT}}"
# NONREFERENCE="${FASTA:-${QUERY}}"
# NONREFERENCE="${NONREFERENCE:-${NONREFERENCE_DEFAULT}}"
# BACKGROUND="${BACKGROUND:-${BACKGROUND_DEFAULT}}"
REFERENCE="${REFERENCE:-${REFERENCE_DEFAULT}}"
DOMAIN="${DOMAIN:-${DOMAIN_DEFAULT}}"
PREFIX="${PREFIX:-${PREFIX_DEFAULT}}"
SC_ALGO="${SC_ALGO:-${SC_ALGO_DEFAULT}}"
RPS_DB="${RPS_DB:-${RPS_DB_DEFAULT}}"
TARGET_DB="${TARGET_DB:-${TARGET_DB_DEFAULT}}"
PAM="${PAM:-${PAM_DEFAULT}}"
MISMATCH="${MISMATCH:-${MISMATCH_DEFAULT}}"
GAP="${GAP:-${GAP_DEFAULT}}"
MINID="${MINID:-${MINID_DEFAULT}}"
MINLEN="${MINLEN:-${MINLEN_DEFAULT}}"
MERGE_WITHIN="${MERGE_WITHIN:-${MERGE_WITHIN_DEFAULT}}"
LENGTH="${LENGTH:-${LENGTH_DEFAULT}}"
CHECK_ID_PREMERGE="${CHECK_ID_PREMERGE:-$(if [[ ${DOMAIN} == 'gene' ]]; then echo 'True'; else echo 'False'; fi)}"
CHECK_BG="${CHECK_BG:-True}"
AUTO="${AUTO:-False}"
RELAX="${RELAX:-False}"
BY_GENE="${BY_GENE:-False}"
SCREEN_REF="${SCREEN_REF:-False}"
MERGE_REDUNDANT="${MERGE_REDUNDANT:-False}"
MASK_REF="${MASK_REF:-${SCREEN_REF}}"
UNMASK_REF="${UNMASK_REF:-False}"
CHECK_RECIP="${CHECK_RECIP:-False}" ## check if each candidate target has a better bitscore to non-target genes, and if so remove it from the list of candidate targets
SEP="${SEP:-${SEP_DEFAULT}}"
SETS="${SETS:-${SETS_DEFAULT}}"
RELAX_RECIP="${RELAX_RECIP:-False}"
REPORT_BG="${REPORT_BG:-False}"


## only instantiate BACKGROUND if not ( TARGET provided & ACCS not provided ) and --skip-bg-check not raised
if [[ "${CHECK_BG}" == "True" ]] && ! ( ! [ -z "${TARGET}" ] && [ -z "${ACC}" ] && [ -z "${ACCS_F}" ] ); then
    ## if BACKGROUND not provided, default to QUERY
    BACKGROUND="${BACKGROUND:-${QUERY}}"
fi


## if checking members
if ! [ -z "${CHECK_MEMBERS}" ]; then
    readarray -d ',' -t clusters_a <<< "$(tr '\n' ',' < ${CLUSTER_LOOKUP} | sed 's/,\+$//g')"
    ## get cluster gene members (use a lookup file)
    for ref_cluster in "${clusters_a[@]}"; do
        ref_cluster_a=( ${ref_cluster} )
        if [[ ";${ref_cluster_a[2]};" =~ ";${CHECK_MEMBERS};" ]]; then
            cat ${CLUSTER_DIR}/${ref_cluster_a[0]}
            exit 0
        fi
    done
    ## if not found:
    echo "${cluster} is not a recognised cluster name."
    exit 0
## throw error if directory not provided
elif [ -z "${DIR}" ]; then
    echo "Directory required. Please provide directory using '-d <path to directory>'"
    exit 1
## throw error if GENE or CLUSTER or QUERY or TARGET are not provided
elif [ -z "${GENE}" ] && [ -z "${CLUSTER}" ] && [ -z "${TARGET}" ] && \
     [[ "${QUERY}" == "${QUERY_DEFAULT}" ]]; then
    echo "Gene ID(s) or cluster names or fasta file required. Please provide gene ID(s) using '-g <gene ID(s)>', cluster name(s) using '-c <cluster name(s)>', a fasta file to query using '-q <path to file>', or a fasta file of target sequences using '-t <path to file>'."
    exit 1
## else if any combination of 2 or more of GENE or CLUSTER or TARGET (note that QUERY is compatible with CLUSTER or GENE but not with TARGET) are provided for some reason
elif (! [ -z "${GENE}" ] && ! [ -z "${CLUSTER}" ]) || \
         # (! [ -z "${GENE}" ] && (! [ -z "${TARGET}" ])) || \
         # (! [ -z "${CLUSTER}" ] && (! [ -z "${TARGET}" ])) || \
         (! [ -z "${TARGET}" ] && (! [[ "${QUERY}" == "${QUERY_DEFAULT}" ]])); then
    # echo "Please only use either '-g <gene ID(s)>', '-c <cluster name(s)>', '-t <path to file>', and not 2 or more at the same time. These parameters are mutually exclusive. '-f <path to file>' is also mutually exclusive with '-t <path to file>'."
    echo "'-g <gene ID(s)>' and '-c <cluster name(s)>' are mutually exclusive. '-t <path to file>' and '-q <path to file>' are also mutually exclusive."
    exit 1
## else if both ACCS_F, ACC are provided for some reason
elif (! [ -z "${ACCS_F}" ] && ! [ -z "${ACC}" ]); then
    echo "Please only use either '-a <accession number>' or '-i <path to file>', and not both at the same time. These parameters are mutually exclusive."
    exit 1
# ## else if either gene or cluster is provided w/ target, but --skip-bg-check raised or --screen-ref not raised
# elif ( (! [ -z "${GENE}" ] || ! [ -z "${CLUSTER}" ] ) &&
#            (! [ -z "${TARGET}" ] ) &&
#            ( [[ "${SCREEN_REF}" == "False" ]] || [[ "${CHECK_BG}" == "False" ]] ) ); then
#     echo "'-g <gene ID(s)>', '-c <cluster name(s)>' should not be used with '-t <path to file>' unless --screen-ref is raised (and --skip-bg-check is NOT raised)."
#     exit 1
## else if either ACCS_F or ACC is provided w/ target, but --skip-bg-check raised
elif ( (! [ -z "${ACCS_F}" ] || ! [ -z "${ACC}" ] ) &&
           (! [ -z "${TARGET}" ] ) && [[ "${CHECK_BG}" == "False" ]] ); then
    echo "'-a <accession number(s)>', '-i <path to file>' should not be used with '-t <path to file>' unless --skip-bg-check is NOT raised."
    exit 1
## else if either gene or cluster is provided, but accessions are missing
elif ( (! [ -z "${GENE}" ] || (! [ -z "${CLUSTER}" ])) && [ -z "${ACCS_F}" ] &&
           [ -z "${ACC}" ] && [[ "${QUERY}" == "${QUERY_DEFAULT}" ]] ); then
    echo "Fasta file or accession name(s) or number(s) required. Please provide a fasta file in which to query for target(s) using '-q <path to file>', a file of accession number(s) using '-i <path to file>', or a comma-separated list of accessions numbers using '-a <accession number>'"
    exit 1
elif ( (! [ -z "${CDS_EXT}" ] && [ -z "${GENOME_EXT}" ]) ||
           ( [ -z "${CDS_EXT}"] && (! [ -z "${GENOME_EXT}" ]))); then
    echo "'--extend-cds <path to file>' must be used with '--extend-genome <path to file>'."
    exit 1
elif ! [[ "${QUERY}" == "${QUERY_DEFAULT}" ]] && ! [ -f "${QUERY}" ]; then
    echo "Fasta file ${QUERY} does not exist."
    exit 1
elif ! [ -z "${TARGET}" ] && ! [ -f "${TARGET}" ]; then
    echo "Target file ${TARGET} does not exist."
    exit 1
elif ! [ -z "${CDS_EXT}" ] && ! [ -f "${CDS_EXT}" ]; then
    echo "Fasta file ${CDS_EXT} (--extend-cds) does not exist."
    exit 1
elif ! [ -z "${GENOME_EXT}" ] && ! [ -f "${GENOME_EXT}" ]; then
    echo "Fasta file ${GENOME_EXT} (--extend-genome) does not exist."
    exit 1
fi

## if --screen-ref raised and --skip-bg-check not raised, let user know that -g or -c will be used for bg check
if ( (! [ -z "${GENE}" ] || ! [ -z "${CLUSTER}" ] ) &&
         (! [ -z "${TARGET}" ] ) && [[ "${SCREEN_REF}" == "True" ]] && [[ "${CHECK_BG}" == "True" ]] ); then
    echo "--screen-ref raised: ${GENE} will be masked in the reference genome during background screening in reference."
fi

## convert values of variables storing file and directory paths to absolute paths
path_vars_to_convert=( "QUERY" "TARGET" "ACCS_F" "DIR" "REFERENCE" "BACKGROUND" "EXCLUDE" "BED" "RPS_DB" "CLUSTER_LOOKUP" "CLUSTER_DIR" "CDS_EXT" "GENOME_EXT" )
for varname in ${path_vars_to_convert[@]}; do
    if ! [ -z "${!varname}" ] && ([ -f "${!varname}" ] || [ -d "${!varname}" ]); then
        eval ${varname}="$(realpath ${!varname})"
    fi
done

## if acc == 'ref', set background to reference (for background checking reasons)
if [[ "${ACC}" == 'ref' ]]; then
    BACKGROUND="${REFERENCE}"
fi

## move to output dir, create if doesn't exist
mkdir -p ${DIR}
cd ${DIR}
DIR=$(pwd)
echo "Output files will be generated in $(pwd)"
tmp_f=${DIR}/tmp.txt
tmp_f2=${DIR}/tmp2.txt

## write log
to_log=''
short_long=( "g,gene,GENE" "c,cluster,CLUSTER" "a,acc,ACC" "i,input,ACCS_F" "d,dir,DIR" "q,query,QUERY" "b,background,BACKGROUND" "r,reference,REFERENCE" "e,exclude,EXCLUDE" "s,sets,SETS" "p,pam,PAM" "l,length,LENGTH" )
long_only=( "bed,BED" "prefix,PREFIX" "domain,DOMAIN" "db,RPS_DB" "target-db,TARGET_DB" "mismatch,MISMATCH" "gap,GAP" "minlen,MINLEN" "minid,MINID" "merge-within,MERGE_WITHIN" "sc-algorithm,SC_ALGO" "cluster-lookup,CLUSTER_LOOKUP" "cluster-dir,CLUSTER_DIR" "check-id-before-merge,CHECK_ID_PREMERGE" "auto,AUTO" "pam,PAM" "length,LENGTH" "relax,RELAX" "screen-ref,SCREEN_REF" "unmask-ref,UNMASK_REF" "check-recip,CHECK_RECIP" "relax-recip,RELAX_RECIP" "attr-mod,ATTR_MOD" "extend-cds,CDS_EXT" "extend-genome,GENOME_EXT" )
for variable in ${short_long[@]}; do
    readarray -d ',' -t v <<< "${variable}"
    varname=${v[2]}
    eval "var=\$${varname}"
    if ! [ -z ${var} ]; then
        to_log+="-${v[0]}|--${v[1]}:\t${var}"
        varname_default=${varname::-1}_DEFAULT
        eval "var_default=\$${varname_default}"
        if ! [ -z ${var_default} ] && [[ ${var} == ${var_default} ]]; then
            to_log+=" (default)\n"
        else
            to_log+="\n"
        fi
    fi
done
for variable in ${long_only[@]}; do
    readarray -d ',' -t v <<< "${variable}"
    varname=${v[1]}
    eval "var=\$${varname}"
    if ! [ -z ${var} ]; then
        to_log+="--${v[0]}:\t${var}"
        varname_default=${varname::-1}_DEFAULT
        eval "var_default=\$${varname_default}"
        if ! [ -z ${var_default} ] && [[ ${var} == ${var_default} ]]; then
            to_log+=" (default)\n"
        else
            to_log+="\n"
        fi
    fi
done
printf -- "${params}\n\n${to_log}" > ${DIR}/${PREFIX}_findgRNA.log

## generate extended GFF/BED & reference genome
if ! [ -z "${GENOME_EXT}" ] && ! [ -z "${CDS_EXT}" ]; then
    dir_ext=${DIR}/extension
    aln_bed=${dir_ext}/aln2gff_tmp.bed
    new_bed=${dir_ext}/ext_$(date +%s).bed
    new_ref=${dir_ext}/ext_reference.fasta
    ${aln2gff} -c ${CDS_EXT} -g ${GENOME_EXT} -d ${dir_ext} -o ${aln_bed} -s ${SEP} \
               --outfmt bed --attr-mod "${ATTR_MOD}"
    cat ${REFERENCE} <(echo) ${GENOME_EXT} > ${new_ref}
    cat ${BED} <(echo) ${aln_bed} > ${new_bed}
    REFERENCE=${new_ref}
    BED=${new_bed}
fi

## generate common get_seq params
if [[ "${BY_GENE}" == "True" ]]; then
    get_seq_common="${get_seq} --acc ref --reference ${REFERENCE} --attr-mod '${ATTR_MOD:-{}}' --bed ${BED} --no-bed --by-gene"
else
    get_seq_common="${get_seq} --acc ref --reference ${REFERENCE} --attr-mod '${ATTR_MOD:-{}}' --bed ${BED} --no-bed"
fi

## conduct reciprocal blast
recip_blast () {
    local output_dir=${1}
    local output_pref=${2}
    local genes=${3}
    local accs_fa=${4}
    local recip_blast6=${5}
    local recip_bed=${6}
    
    echo "Filtering out candidate targets with higher similarity to non-target genes"
    recip_blast6=${output_dir}/${output_pref}_targets_recip.tsv
    recip_bed=${output_dir}/${output_pref}_targets_recip.bed
    blastn -query ${accs_fa} -subject ${REFERENCE} -out ${recip_blast6} \
           -outfmt '6 sseqid sstart send qseqid qstart qend bitscore'
    python3 -c "f = open('${recip_blast6}', 'r'); data = [x[:-1].split('\t') for x in f.readlines()]; f.close(); output = [x[:1] + (x[1:3] if int(x[1]) < int(x[2]) else x[2:0:-1]) + x[3:] for x in data]; f = open('${recip_blast6}', 'w+'); f.write('\n'.join(['\t'.join(x) for x in output])); f.close()"
    bedtools intersect -wao -a ${recip_blast6} -b ${BED} > ${recip_bed}
    python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); from recip_filter import *; remove_non_max_bitscore('${accs_fa}', '${recip_bed}', '${genes}', relax=(${RELAX_RECIP}))"
    rm ${recip_blast6} ${recip_bed}
}

## get reference fa
get_reference_fa () {
    local genes=${1}
    local out_dir=${2}
    local out_pref=${3}
    local fout=${4}
    local fout_cds=${5}
    local fout_pref=${6}
    if ! [ -z ${7} ] && [[ "${7}" != "gene" ]]; then ## check if domain provided
        ## parse domain
        if [[ ! " ${!DOMAINS[@]} " =~ " ${7} " ]]; then
            echo "${7} is not a supported domain name."
            if [[ ${7} =~ ^[0-9]+$ ]]; then
                echo "Attempting to parse ${7} as CDD PSSM-Id."
                domain=${7}
            else
                echo "Unexpected domain input. Please provide a CDD PSSM-Id (numeric) or one of the following supported domain names: $(echo ${!DOMAINS[@]} | sed 's/ /, /g')"
                exit 1
            fi
        else
            domain=${DOMAINS[${7}]}
        fi
        ## some file names
        local aa_fasta=${out_dir}/ref/${out_pref}_ref_protein.fasta
        local domains_tsv=${out_dir}/ref/${out_pref}_${DOMAIN}.tsv

        ## get protein sequences
        echo "Extracting reference domain range(s)"
        # $get_seq --bed ${BED} --gene ${genes} --acc ref --feature CDS --dir ${out_dir} --out ${aa_fasta} --translate --no-bed --attr-mod ${ATTR_MOD:-'{}'} > /dev/null
        ${get_seq_common} --gene ${genes} --feature CDS --dir ${out_dir} --out ${aa_fasta} \
                          --translate > /dev/null
        ## identify identical protein sequences and collapse identical ones temporarily
        python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); from fasta_manip import *; dat = fasta_to_dict('${aa_fasta}'); identicals = {k: set(seqid for seqid, seq in dat.items() if str(seq) == str(v)) for k, v in dat.items()}; identical_sets = set(map(lambda x: tuple(sorted(x)), identicals.values())); dict_to_fasta({seqids[0]: dat[seqids[0]] for seqids in identical_sets}, '${tmp_f}'); open('${tmp_f2}', 'w+').write('\n'.join(['\t'.join(seqids) for seqids in identical_sets]))"
        ## identify domain positions in protein
        rpsblast+ -db ${RPS_DB} -query ${tmp_f} -out ${domains_tsv} \
                  -outfmt "6 $(echo ${HEADER_CDD} | tr ',' ' ')"
        ## generate file for input as 'DOMAIN_F' to getSeq
        awk -v d="${domain}" '{if ($2 ~ d) print $0 "\t" d}' ${domains_tsv} > ${tmp_f} ## filter (overwrite tmp_f as it's no longer needed)
        ## expand filtered output of rpsblast+ to sequences w/ identical protein
        python3 -c "import sys; sys.path.append('/mnt/chaelab/rachelle/src'); from data_manip import *; dat = [line.split('\t') for line in splitlines('${tmp_f}')]; repr_map = [line.split('\t') for line in splitlines('${tmp_f2}')]; output = [[seqid] + line[1:] for seqids in repr_map for seqid in seqids for line in dat if line[0] == seqids[0]]; open('${tmp_f}', 'w+').write('\n'.join(['\t'.join(x) for x in output]))"
        echo -e "$(echo ${HEADER_CDD},domain | tr ',' '\t')\n$(cat ${tmp_f})" > ${domains_tsv} ## add header
        rm ${tmp_f} ${tmp_f2}
        ## get complete domain sequence
        echo "Extracting reference domain sequence(s)"
        # $get_seq --bed ${BED} --gene ${genes} --acc ref --feature CDS --dir ${out_dir} --out ${fout} --domain-file ${domains_tsv} --complete --minlen ${MINLEN} --domain ${domain} --qname-dname "('qseqid', 'domain')" --qstart-qend "('qstart', 'qend')" --adjust-dir --no-bed --attr-mod ${ATTR_MOD:-'{}'} > /dev/null
        ${get_seq_common} --gene ${genes} --feature CDS --dir ${out_dir} --domain-file ${domains_tsv} \
                          --minlen ${MINLEN} --domain ${domain} --qname-dname "('qseqid', 'domain')" \
                          --qstart-qend "('qstart', 'qend')" --adjust-dir --out ${fout} --complete > /dev/null
        ## get CDS-only
        # $get_seq --bed ${BED} --gene ${genes} --acc ref --feature CDS --dir ${out_dir} --out ${fout_cds} --domain-file ${domains_tsv} --minlen ${MINLEN} --domain ${domain} --qname-dname "('qseqid', 'domain')" --qstart-qend "('qstart', 'qend')" --adjust-dir --no-bed --attr-mod ${ATTR_MOD:-'{}'} > /dev/null
        ${get_seq_common} --gene ${genes} --feature CDS --dir ${out_dir} --domain-file ${domains_tsv} \
                          --minlen ${MINLEN} --domain ${domain} --qname-dname "('qseqid', 'domain')" \
                          --qstart-qend "('qstart', 'qend')" --adjust-dir --out ${fout_cds} > /dev/null
        
        rm ${domains_tsv} ## remove temporary file(s)
    else
        ## just get complete CDS + CDS
        echo "Getting reference gene sequence(s)"
        local fout_bed=${fout_pref}_complete.bed
        # $get_seq --bed ${BED} --gene ${genes} --acc ref --feature CDS --dir ${out_dir} --out ${fout} --complete --bed-out ${fout_bed} --adjust-dir --attr-mod ${ATTR_MOD:-'{}'} > /dev/null
        ${get_seq_common} --gene ${genes} --feature CDS --dir ${out_dir} --adjust-dir --out ${fout} \
                          --bed-out ${fout_bed} --complete > /dev/null
        local fout_bed_cds=${fout_pref}_CDS.bed
        # $get_seq --bed ${BED} --gene ${genes} --acc ref --feature CDS --dir ${out_dir} --out ${fout_cds} --bed-out ${fout_bed_cds} --adjust-dir --attr-mod ${ATTR_MOD:-'{}'} > /dev/null
        ${get_seq_common} --gene ${genes} --feature CDS --dir ${out_dir} --adjust-dir --out ${fout_cds} \
                          --bed-out ${fout_bed_cds} > /dev/null
    fi
}

## parse accession IDs into Python-format tuple
if ! [ -z "${ACCS_F}" ]; then ## if file provided
    sed -i 's/\r//g' ${ACCS_F}
    accs_tuple="('$(awk 'NF > 0' ${ACCS_F} | awk 'NR>1{print PREV} {PREV=$0} END{printf("%s",$0)}' | cat | tr '\n' ',' | sed 's/,/\x27,\x27/g')',)"
elif [[ "${ACC}" == 'ref' ]]; then
    accs_tuple="()"
elif [[ "${ACC}" == 'all' ]]; then
    accs_tuple="('$(cat /mnt/chaelab/shared/anna_lena/accession_map.txt | sed 's/,.\+$//' | sed 's/[^0-9]//' | tr '\n' ',' | sed 's/,$//' | sed 's/,/\x27,\x27/g')')"
elif ! [ -z "${ACC}" ]; then
    accs_tuple="('$(echo ${ACC} | sed 's/,/\x27,\x27/g')',)"
else
    accs_tuple="()"
fi

## if --screen-ref raised and --skip-bg-check not raised, let user know that -a will be used for bg check
if ( ( ! [ -z "${ACCS_F}" ] || ! [ -z "${ACC}" ] ) &&
         (! [ -z "${TARGET}" ] ) && [[ "${QUERY}" == "${QUERY_DEFAULT}" ]] &&
         [[ "${CHECK_BG}" == "True" ]] ); then
    echo "--skip-bg-check not raised: ${accs_tuple} sequences in Van de Weyer et al.'s (2019) RenSeq dataset will be screened for off-targets during background screening."
fi

## if --screen-ref raised and --skip-bg-check not raised, let user know that -g/-c will be used for bg check
if ( ( ! [ -z "${GENE}" ] || ! [ -z "${CLUSTER}" ] ) &&
         (! [ -z "${TARGET}" ] ) && [[ "${QUERY}" == "${QUERY_DEFAULT}" ]] &&
         [[ "${CHECK_BG}" == "True" ]] ); then
    echo "--skip-bg-check not raised: ${CLUSTER}${GENE} in Van de Weyer et al.'s (2019) RenSeq dataset will be masked for off-targets during background screening."
fi

## associative array for <output prefix>:<fasta file> combination
declare -a fasta_a

## if target file provided as query
if ! [ -z "${TARGET}" ] && ( ( [[ "${SCREEN_REF}" == 'False' ]] &&
                                   [ -z "${CLUSTER}" ] && [ -z "${GENE}" ] ) ||
                                 [[ "${CHECK_BG}" == 'False' ]] ); then
    
    mkdir -p ${DIR}/${PREFIX}
    fasta_a+=("${DIR}/${PREFIX} ${PREFIX} ${TARGET}")
    echo "checkpt 1"

## else if cluster or gene specified
elif ! [ -z "${CLUSTER}" ] || ! [ -z "${GENE}" ]; then

    declare -A genes_a
    
    if ! [ -z "${CLUSTER}" ]; then
        
        ## extract cluster members
        if [ -f "${CLUSTER}" ]; then
            readarray -d ',' -t clusters <<< "$(tr '\n' ',' < ${CLUSTER} | sed 's/,\+$//g')"
        else
            readarray -d , -t clusters <<< ${CLUSTER}
        fi
        
        ## read cluster lookup table
        readarray -d ',' -t clusters_a <<< "$(tr '\n' ',' < ${CLUSTER_LOOKUP} | sed 's/,\+$//g')"
        for cluster in ${clusters[@]}; do
            fout_pref=${PREFIX}_${cluster}
            ## get cluster gene members (use a lookup file)
            for ref_cluster in "${clusters_a[@]}"; do
                ref_cluster_a=( ${ref_cluster} )
                if [[ ";${ref_cluster_a[2]};" =~ ";${cluster};" ]]; then
                    ## add prefix (key) and members (value) to associative array
                    cluster_members_f=${CLUSTER_DIR}/${ref_cluster_a[0]}
                    genes_a["${fout_pref}"]="$(tr '\n' ',' < ${cluster_members_f} | sed 's/,\+$//g')"
                    continue 2
                fi
            done
            ## if not found:
            echo "${cluster} is not a recognised cluster name. Please use a different alias or -f <path to fasta file> or -g <comma-separated gene(s)> instead."
        done
        
    elif ! [ -z "${GENE}" ]; then
        ## if gene names provided, get sequences of genes
      if [ -f "${GENE}" ]; then
            genes_a["${PREFIX}"]="$(tr '\n' ',' < ${GENE} | sed 's/,\+$//g')"
        else
            genes_a["${PREFIX}"]=${GENE}
        fi
    fi
    
    ## get fasta_a
    for group in ${!genes_a[@]}; do
        output_dir=${DIR}/${group}
        output_pref=${group}_${DOMAIN}
        ref_fa=${output_dir}/ref/${group}_ref_${DOMAIN}_complete.fasta
        cds_fa=${output_dir}/ref/${group}_ref_${DOMAIN}_CDS.fasta
        fout_pref=${output_dir}/ref/${group}_ref_${DOMAIN}
        accs_blast6=${output_dir}/${output_pref}_targets.tsv
        accs_fa=${output_dir}/${output_pref}_targets.fasta
        align_fa=${output_dir}/${group}_${DOMAIN}_mafft.fa
        almask_fa=${output_dir}/${group}_${DOMAIN}_toMask.fa
        genes=${genes_a["${group}"]} ## comma-separated!
        ## make directories
        mkdir -p ${output_dir}/ref
        ## get reference sequences
        echo "Retrieving reference Col-0 sequence(s)"
        ## get ref domain seqs
        get_reference_fa "${genes}" ${output_dir} ${group} ${ref_fa} ${cds_fa} ${fout_pref} ${DOMAIN}
        ## duplicate reference complete files into accs_fa and rename sequences if accs == ref
        ## if using '-a ref', ${ref_fa} is the target (${accs_fa})
        if [[ "${ACC}" == 'ref' ]]; then
            cp ${ref_fa} ${accs_fa}
        elif [[ "${accs_tuple}" != '()' ]]; then
            ## blast to AL70 database
            echo "Searching for homologues in ${accs_tuple}"
            blastn -query ${ref_fa} -db ${DB_AL70} -outfmt "6 $(echo ${HEADER_AL70} | tr ',' ' ')" \
                   -out ${accs_blast6}
            if [ $(wc -l < ${accs_blast6}) -gt 1 ]; then
                ## add header
                echo -e "$(echo ${HEADER_AL70} | tr ',' '\t')\n$(cat ${accs_blast6})" > ${accs_blast6}
                ## extract sequences from required accessions
                echo "Extracting sequences for homologues in ${accs_tuple}"
                python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); from extract_domain_from_AL70 import *; main(blast6_fname='${accs_blast6}', accIDs=${accs_tuple}, fout='${accs_fa}', fasta='${QUERY_DEFAULT}', min_id=${MINID}, min_len=${MINLEN}, merge_within_range=${MERGE_WITHIN}, check_id_before_merge=(${CHECK_ID_PREMERGE}))"
            else
                ## remove the useless file
                ## - so we can use [ -z {accs_blast6} ] later to check if this step succeeded
                rm ${accs_blast6}
                ## if this step is used for target discovery (i.e. not using -t/-q/-f),
                ## # exit because we can't continue w/ no hits.
                if ( [ -z "${TARGET}" ] && [[ "${QUERY}" == "${QUERY_DEFAULT}" ]] ); then
                    echo "No blast hits, exiting programme."
                    exit 2
                fi
            fi
            ## check if best scoring hits are the same as target genes
            if ( [ -f ${accs_blast6} ] &&
                     [ ${CHECK_RECIP} == 'True' ] && ! [ -z ${REFERENCE} ] && ! [ -z ${BED} ] ) ; then
                recip_blast ${output_dir} ${output_pref} ${genes} ${accs_fa} ${recip_blast6} ${recip_bed}
            fi
        fi
        
        ## if user isn't using -t/-f/-q, (i.e. user is searching AL70 database or ref),
        ## # align these AL70 homologues/reference genes
        if ( [ -z "${TARGET}" ] && [[ "${QUERY}" == "${QUERY_DEFAULT}" ]] ); then
            if [ $(grep '>' ${accs_fa} | wc -l) -lt 1 ]; then
                echo "No targets found, exiting programme."
                exit 2
            fi
            ## align CDS to complete CDS
            echo "Aligning reference CDS and reference complete CDS"
            mafft --quiet ${cds_fa} > ${tmp_f}
            mv ${tmp_f} ${align_fa}
            mafft --quiet --add ${ref_fa} ${align_fa} > ${tmp_f}
            echo "Aligning target sequences to reference sequences"
            mafft --quiet --adjustdirectionaccurately --add ${accs_fa} ${tmp_f} > ${align_fa}
            rm ${tmp_f}
            
            fasta_a+=("${output_dir} ${output_pref} ${accs_fa} ${align_fa} ${ref_fa} ${cds_fa}")
            # for_fasta_a="${output_dir} ${output_pref} ${accs_fa} ${align_fa} ${ref_fa} ${cds_fa}"
        else
            ## if using QUERY or TARGET, then we're only in here to get al_mask
            mv ${accs_fa} ${almask_fa}
        fi
        
        if [ -f "${accs_blast6}" ]; then
            rm ${accs_blast6} ## remove file
        fi
        
        ## if using QUERY, blast to QUERY
        if [[ "${QUERY}" != "${QUERY_DEFAULT}" ]]; then
            echo "Searching for homologues in ${QUERY}"
            blastn -query ${ref_fa} -subject ${QUERY} -outfmt "6 $(echo ${HEADER_AL70} | tr ',' ' ')" \
                   -out ${accs_blast6}
            if [ $(wc -l < ${accs_blast6}) -le 1 ]; then
                echo "No blast hits, exiting programme."
                exit 2
            fi
            ## add header
            echo -e "$(echo ${HEADER_AL70} | tr ',' '\t')\n$(cat ${accs_blast6})" > ${accs_blast6}
            ## extract sequences
            echo "Extracting sequences for homologues in ${QUERY}"
            python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); from extract_domain_from_AL70 import *; main(blast6_fname='${accs_blast6}', accIDs=('.',), fout='${accs_fa}', fasta='${QUERY}', min_id=${MINID}, min_len=${MINLEN}, merge_within_range=${MERGE_WITHIN}, check_id_before_merge=(${CHECK_ID_PREMERGE}), pattern=lambda accID:accID)"
            ## check if best scoring hits are the same as target genes
            if ( [ ${CHECK_RECIP} == 'True' ] && ! [ -z ${REFERENCE} ] && ! [ -z ${BED} ] ) ; then
                recip_blast ${output_dir} ${output_pref} ${genes} ${accs_fa} ${recip_blast6} ${recip_bed}
            fi
            
            if [ $(grep '>' ${accs_fa} | wc -l) -lt 1 ]; then
                echo "No targets found, exiting programme."
                exit 2
            fi
            
            ## align CDS to complete CDS
            echo "Aligning reference CDS and reference complete CDS"
            mafft --quiet ${cds_fa} > ${tmp_f}
            mv ${tmp_f} ${align_fa}
            mafft --quiet --add ${ref_fa} ${align_fa} > ${tmp_f}
            echo "Aligning target sequences to reference sequences"
            mafft --quiet --adjustdirectionaccurately --add ${accs_fa} ${tmp_f} > ${align_fa}
            rm ${tmp_f}
                
            fasta_a+=("${output_dir} ${output_pref} ${accs_fa} ${align_fa} ${ref_fa} ${cds_fa} ${almask_fa}")
            
            if [ -f "${accs_blast6}" ]; then
                rm ${accs_blast6} ## remove file
            fi
        ## if using TARGET and SCREEN_REF, generate relevant entry and move on (don't do homologue discovery)
        elif (! [ -z "${TARGET}" ] ); then
            if [ -f ${align_fa} ]; then
                rm ${align_fa} ## reset align_fa to non-existent
            fi
            if [[ "${SCREEN_REF}" == "True" ]]; then
                ## we're only here because we need the reference sequences
                fasta_a+=("${DIR}/${PREFIX} ${PREFIX} ${TARGET} ${align_fa} ${ref_fa} ${cds_fa} ${almask_fa}")
            else
                ## use non-existent file align_fa as placeholder
                fasta_a+=("${DIR}/${PREFIX} ${PREFIX} ${TARGET} ${align_fa} ${align_fa} ${align_fa} ${almask_fa}")
            fi
        fi
    done
fi

## work on fasta_a
for entry in "${fasta_a[@]}"; do
    echo "Finding gRNA"
    entry_a=( $entry ) ## output_dir, prefix, fasta, alignment (also fasta)
    # echo "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); from find_common_gRNA import *; find_cluster_gRNA_in_acc('${entry_a[2]}', ${accs_tuple}, '${entry_a[0]}', background_usr_fname='${BACKGROUND}', manual_check=(not (${AUTO})), fout_pref='${entry_a[1]}', sc_algorithm='${SC_ALGO}', accs_background_fname='${BACKGROUND}', max_mismatch=${MISMATCH}, max_gap=${GAP}, pam='${PAM}', gRNA_len=int(${LENGTH}), alignment_fname='${entry_a[3]}', exclude='${EXCLUDE}', relax=(${RELAX}), reference_fasta = '${REFERENCE}', mask_reference='${entry_a[4]}', cds_fasta='${entry_a[5]}', complete_fasta='${entry_a[4]}', check_bg=${CHECK_BG})"
    # python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); from find_common_gRNA import *; find_cluster_gRNA_in_acc('${entry_a[2]}', ${accs_tuple}, '${entry_a[0]}', background_usr_fname='${BACKGROUND}', manual_check=(not (${AUTO})), fout_pref='${entry_a[1]}', sc_algorithm='${SC_ALGO}', accs_background_fname='${BACKGROUND}', max_mismatch=${MISMATCH}, max_gap=${GAP}, pam='${PAM}', gRNA_len=int(${LENGTH}), alignment_fname='${entry_a[3]}', exclude='${EXCLUDE}', relax=(${RELAX}), reference_fasta = '${REFERENCE}', mask_reference='${entry_a[4]}', cds_fasta='${entry_a[5]}', complete_fasta='${entry_a[4]}', check_bg=${CHECK_BG})"
    ## changed mask_reference to mask_reference=False
    python3 -c "import sys; sys.path.append('${SCRIPT_DIR}/scripts'); from find_common_gRNA import *; find_cluster_gRNA_in_acc('${entry_a[2]}', ${accs_tuple}, '${entry_a[0]}', background_usr_fname='${BACKGROUND}', manual_check=(not (${AUTO})), fout_pref='${entry_a[1]}', sc_algorithm='${SC_ALGO}', accs_background_fname='${BACKGROUND_DEFAULT}', max_mismatch=${MISMATCH}, max_gap=${GAP}, pam='${PAM}', gRNA_len=int(${LENGTH}), alignment_fname='${entry_a[3]}', exclude='${EXCLUDE}', relax=(${RELAX}), reference_fasta = '${REFERENCE}', mask_reference=${MASK_REF}, screen_reference=${SCREEN_REF}, cds_fasta='${entry_a[5]}', complete_fasta='${entry_a[4]}', check_bg=${CHECK_BG}, num_sets=${SETS}, report_bg=${REPORT_BG}, nonref_mask_fname = '${entry_a[6]}')" ## TODO: implement -a -g bg check even w/ -q or -t
done

if ! [ -z "${dir_ext}" ]; then
    rm -r ${dir_ext}
fi

exit 0
