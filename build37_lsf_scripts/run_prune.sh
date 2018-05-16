#!/bin/bash

set -eo pipefail

# SET DIRECTORY HERE:
DIR=`pwd`

# SET LOCATION OF SCRIPTS
BIN_DIR=$( dirname "$0" )

# SET INPUT VCF.
# Must be gzipped.
INPUT="$1"

# SET OUTPUT VCF.
# Will be gzipped
OUTPUT="$2"

LOGDIR="$DIR/log"
mkdir -p $LOGDIR

bsub \
    -M 12000000 \
    -R 'rusage[mem=12000, tmp=750] select[mem>12000 && tmp>750]' \
    -J "prune" \
    -oo "$LOGDIR/prune.log" \
    bash "$BIN_DIR/prune.sh" "${INPUT}" "${OUTPUT}"

