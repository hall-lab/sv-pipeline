#!/bin/bash

set -eo pipefail

# SET DIRECTORY HERE:
DIR=`pwd`

# SET LOCATION OF SCRIPTS
BIN_DIR=$( dirname "$0" )

# SET SAMPLEMAP TO PROCESS
# Each line in this file contains a sample name, whitespace and then the path
# to the BAM file for that sample
#
# Sample names are also used to create directory paths and thus should be unique
SAMPLEMAP="$1"

# SET MASTER VCF CONTAINING ALL SITES ACROSS COHORT
MASTER="$2"

LOGDIR=$DIR/post-merge/gt/logs
mkdir -p $LOGDIR

while read SAMPLE BAM
do
    NAME=$SAMPLE.post_merge_gt
    OUTPUT=$DIR/post-merge/gt/$SAMPLE.gt.vcf

    if [ ! -s $OUTPUT ]
    then
        bsub \
        -M 12000000 \
        -R 'rusage[mem=12000] select[mem>12000]' \
        -J "$NAME" \
        -oo "$LOGDIR/$NAME.log" \
        bash "$BIN_DIR/post_merge_genotype.sh" "$SAMPLE" "$BAM" "$MASTER" "$OUTPUT"
    fi
done < $SAMPLEMAP

