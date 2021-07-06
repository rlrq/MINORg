#!/bin/bash

ORIGINAL_DIR=$(pwd)
TOOLS_DIR="$(dirname "$(readlink -f "$0")")"
SCRIPTS_DIR="${TOOLS_DIR}/../scripts"

## if no arguments, show manual
if [[ $# -eq 0 ]]; then
    man -l ${TOOLS_DIR}/MANUAL_minimumset.1
    exit 1
fi

while (( "$#" )); do
    case "$1" in
        -m|--mapping) MAPPING="${2}";;
        # -t|--target) TARGET="${2}";;
        --grna) FASTA="${2}";;
        -f|--fasta) FASTA="${2}";;
        -s|--sets) SETS="${2}";;
        -d|--dir) DIR="${2}";;
        -e|--exclude) EXCLUDE="${2}";;
        --prefix) PREFIX="${2}";;
        --out-fasta) FASTA_OUT="${2}";;
        --out-mapping) MAPPING_OUT="${2}";;
        --auto) AUTO='True';;
        --input-ver) INPUT_V="${2}";;
        --output-ver) OUTPUT_V="${2}";;
        --sc-algorithm) SC_ALGO="${2}";;
        --ignore-invalid) IGNORE_INVALID='True';;
        --accept-cds-unknown) ACCEPT_CDS_UNKNOWN='True';;
        -h|--help) man -l ${TOOLS_DIR}/MANUAL_minimumset.1; exit 0;;
        -v|--version) echo "minimumset v1"; exit 0;;
    esac
    shift
done

SETS="${SETS:-1}"
AUTO="${AUTO:-False}"
INPUT_V="${INPUT_V:-None}"
OUTPUT_V="${OUTPUT_V:-2}"
SC_ALGO="${SC_ALGO:-LAR}"
IGNORE_INVALID="${IGNORE_INVALID:-False}"
ACCEPT_CDS_UNKNOWN="${ACCEPT_CDS_UNKNOWN:-False}"

if [ -z "${MAPPING}" ]; then
    echo "Target mapping file is required. Please provide the file using '-m <file>'."
    exit 1
elif [ -z "${FASTA}" ]; then
    echo "FASTA file of gRNA sequences is required. Please provide the file using '--grna <file>'."
    exit 1
elif ! ( ( ! [ -z "${DIR}" ] && ! [ -z "${PREFIX}" ] ) ||
           ( ! [ -z "${FASTA_OUT}" ] && ! [ -z "${MAPPING_OUT}" ] ) ); then
    echo "Please use '-d <directory>' and '--prefix <prefix>' to have the programme automatically generate output file names, or '--out-fasta <file>' and '--out-mapping <file>' to explicitly specify output file names."
    exit 1
fi

python3 -c "import sys; sys.path.append('${SCRIPTS_DIR}'); from get_minimum_set import *; get_minimum_sets_from_files_and_write(num_sets = ${SETS}, mapping = '${MAPPING}', targets = '${TARGET}', fasta = '${FASTA}', directory = '${DIR}', exclude_fname = '${EXCLUDE}', prefix = '${PREFIX}', fout_fasta = '${FASTA_OUT}', fout_mapping = '${MAPPING_OUT}', manual_check = (not ${AUTO}), input_map_ver = (${INPUT_V}), output_map_ver = (${OUTPUT_V}), sc_algorithm = '${SC_ALGO}', ignore_invalid = (${IGNORE_INVALID}), accept_unknown_within_cds_status = (${ACCEPT_CDS_UNKNOWN}))"

exit 0
