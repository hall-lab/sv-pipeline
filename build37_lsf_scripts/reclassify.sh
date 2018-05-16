#!/bin/bash

set -eo pipefail

# DEPENDENCIES
# + svtools
# + bgzip

# SET SVTOOLS PATH HERE:
# Edit this to point to your svtools install if it's not currently on your path
SVTOOLS="svtools"

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

# MODE SHOULD MATCH WHETHER OR NOT YOU PROVIDE A TRAINING VARIANT FILE
MODE="$5"
TRAINING_VARS="$6"

TSTRING=""
if [[ -n "$TRAINING_VARS" ]]
then
    TSTRING="-t $TRAINING_VARS"
fi

zcat "$INPUT" \
     | $SVTOOLS classify \
        -g "$SEX_FILE" \
        -a "$REPEATS" \
        -m "$MODE" $TSTRING \
     | bgzip -c > "$OUTPUT"
