#!/bin/bash

# SET DIRECTORY HERE:
DIR=`pwd`

# DEFINE SAMPLE TO BAM MAPPING:
# Each line in this file contains a sample name, whitespace and then the path
# to the BAM file for that sample
#
# Note that splitters and discordants files are assumed to be present in the same
# location as the main BAM, begin with the same prefix as the main BAM, and
# end in splitters.bam and discordants.bam
#
# Sample names are also used to create directory paths and thus should be unique
# within this file
SAMPLEMAP="$1"

# SET PATH TO REFERENCE FASTA USED FOR ALIGNMENTS:
REF="$2"

# SET PATH TO EXCLUDE FILE FOR LUMPY:
EXCLUDE="$3"

# SET PATH TO CUSTOM CONFIGURATION FOR SPEEDSEQ:
# Typically this is used to obtain a newer version of svtyper than what is included
# within speedseq
# Remove this and the -K option below if not needed
CUSTOM_CONFIG=/path/to/config

# SET SPEEDSEQ PATH HERE:
# Edit this to point to your speedseq install if it's not currently on your path
SPEEDSEQ_LOCATION="speedseq"

while read SAMPLE BAM
do
    echo $SAMPLE
    mkdir -p $DIR/lumpy/$SAMPLE/log
    SPL="${BAM%.*}.splitters.bam"
    DCD="${BAM%.*}.discordants.bam"
    if [[ -z "$BAM" ]]
    then
        echo -e "$SAMPLE\t$BAM"
    fi

    if [[ ! -e "$DIR/lumpy/$SAMPLE/$SAMPLE.sv.vcf.gz" ]]
    then
        bsub \
        -M 30000000 \
        -R 'rusage[mem=30000] select[mem>30000]' \
        -J $SAMPLE.lumpy \
        -oo $DIR/lumpy/$SAMPLE/log/$SAMPLE.lumpy.log \
        -eo $DIR/lumpy/$SAMPLE/log/$SAMPLE.lumpy.log \
        "${SPEEDSEQ_LOCATION} sv \
        -o $DIR/lumpy/$SAMPLE/$SAMPLE \
        -T $DIR/lumpy/$SAMPLE/temp \
        -K $CUSTOM_CONFIG \
        -R $REF \
        -B $BAM \
        -S $SPL \
        -D $DCD \
        -v \
        -x $EXCLUDE \
        -d \
        -P \
        -g \
        -k"
    fi
done < $SAMPLEMAP
