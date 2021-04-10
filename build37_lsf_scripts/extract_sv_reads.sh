#!/bin/bash

set -eu

# SET LOGGING DIRECTORY HERE:
LOG=`pwd`

# SET EXTRACT-SV-READS PATH HERE:
EXTRACT_SV_READS="extract-sv-reads"

# SET SAMTOOLS PATH HERE:
SAMTOOLS="samtools"

# SET NUMBER OF THREADS HERE:
THREADS=4

for file in "$@"
do
    NAME=$(basename "$file" .bam)
    if [ -e "$file" ]
    then
        TMPSPLITTER=$NAME.splitters.tmp.bam
        SPLITTER=$NAME.splitters.bam
        TMPDISCORDANT=$NAME.discordants.tmp.bam
        DISCORDANT=$NAME.discordants.bam
        
        if [[ ! -s $DISCORDANT || ! -s $DISCORDANT.bai || ! -s $SPLITTER || ! -s $SPLITTER.bai ]]
        then
            bsub \
                -oo ${LOG}/extract.$NAME.log \
                -n ${THREADS} \
                -M 4000000 \
                -R 'select[mem>4000 && tmp>100] rusage[mem=4000, tmp=100] span[hosts=1]' \
                bash -c "set -e; ${EXTRACT_SV_READS} --input-threads ${THREADS} -e -r -i $file -s $TMPSPLITTER -d $TMPDISCORDANT && ${SAMTOOLS} index $TMPSPLITTER && ${SAMTOOLS} index $TMPDISCORDANT && mv $TMPDISCORDANT $DISCORDANT && mv $TMPDISCORDANT.bai $DISCORDANT.bai && mv $TMPSPLITTER $SPLITTER && mv $TMPSPLITTER.bai $SPLITTER.bai"
        fi
    fi
done

