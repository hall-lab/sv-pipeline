#!/bin/bash

set -eo pipefail

# SET DIRECTORY HERE:
DIR=`pwd`

# SET LOCATION OF SCRIPTS
BIN_DIR=$( dirname "$0" )

# SET PATH TO LIST OF INPUT FILES FOR BATCH
# Each line should be the path to a VCF file for a sample
BATCH_LIST="$1"

# SET BATCH NAME TO USE IN OUTPUT VCF
# This will be inserted into the sample name column of the lsort output
# A dummy genotyped of 1/1 will be used for all variants.
BATCH="$2"

# NAME OF OUTPUT FILE FROM LSORT. WILL BE BGZIPPED.
SORTFILE="$3"

# NAME OF OUTPUT FILE FROM LMERGE. WILL BE BGZIPPED.
MERGEFILE="$4"

LOGDIR=$DIR/batches/$BATCH/logs
mkdir -p $LOGDIR
NAME=$BATCH.lsort_lmerge.1

bsub \
    -M 12000000 \
    -R 'rusage[mem=12000] select[mem>12000]' \
    -J "$NAME" \
    -oo "$LOGDIR/$NAME.log" \
    bash "$BIN_DIR/lsort_lmerge.1.sh" "$BATCH_LIST" "$BATCH" "$SORTFILE" "$MERGEFILE"
