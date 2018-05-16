#!/bin/bash

set -eo pipefail

# SET DIRECTORY HERE:
DIR=`pwd`

# SET LOCATION OF SCRIPTS
BIN_DIR=$( dirname "$0" )

# SET PATH TO LIST OF INPUT FILES FOR BATCH
# Each line should be the path to a VCF file for a sample
BATCH_LIST="$1"

# NAME OF OUTPUT FILE FROM LSORT. WILL BE BGZIPPED.
SORTFILE="$2"

# NAME OF OUTPUT FILE FROM LMERGE. WILL BE BGZIPPED.
MERGEFILE="$3"

LOGDIR=$DIR/logs
mkdir -p $LOGDIR
NAME=AllBatches.lsort_lmerge.2

bsub \
-M 12000000 \
-R 'rusage[mem=12000] select[mem>12000]' \
-J "$NAME" \
-oo "$LOGDIR/$NAME.log" \
bash "$BIN_DIR/lsort_lmerge.2.sh" "$BATCH_LIST" "$SORTFILE" "$MERGEFILE"
