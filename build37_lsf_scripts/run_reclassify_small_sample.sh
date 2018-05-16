#!/bin/bash

set -eo pipefail

# SET DIRECTORY HERE:
DIR=`pwd`

# SET LOCATION OF SCRIPTS
BIN_DIR=$( dirname "$0" )

# INPUT FILE FOR RECLASSIFICATION
# should be gzipped
INPUT="$1"

# FILE CONTAINING SEX FOR EACH SAMPLE
# Use a 1 for males and 2 for females
SEX_FILE="$2"

# FILE CONTAINING RECENTLY DIVERGED REPEATS
# See https://github.com/hall-lab/svtools/blob/master/Tutorial.md#generate-a-repeat-elements-bed-file
REPEATS="$3"

# OUTPUT FILE
# Will be gzipped
OUTPUT="$4"

# TRAINING VARIANTS
TRAINING_FILE="$5"

LOGDIR=$DIR/log
mkdir -p $LOGDIR
NAME=classify.small_sample

bsub \
-M 12000000 \
-R 'rusage[mem=12000] select[mem>12000]' \
-J "$NAME" \
-oo "$LOGDIR/$NAME.log" \
bash "$BIN_DIR/reclassify.sh" "$INPUT" "$SEX_FILE" "$REPEATS" "$OUTPUT" "naive_bayes" "$TRAINING_FILE"

