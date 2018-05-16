#!/bin/bash

set -eo pipefail
BIN_DIR=$( dirname $0 )

# SET DIRECTORY HERE:
DIR=`pwd`

LUMPY_DIR="$DIR/lumpy"
CREATE_COORDS="create_coordinates"

# DEFINE SAMPLE TO BAM MAPPING:
# Each line in this file contains a sample name, whitespace and then the path
# to the BAM file for that sample
#
# Note that splitters and discordants files are assumed to be present in the same
# location as the main BAM, begin with the same prefix as the main BAM, and
# end in splitters.bam and discordants.bam
#
# Sample names are also used to create directory paths and thus should be unique
SAMPLEMAP="$1"

OUTDIR="$DIR/post-merge/cn"
mkdir -p "$OUTDIR/"
LOGDIR="$OUTDIR/logs"
mkdir -p $LOGDIR

COORDINATES="$OUTDIR/coordinates"

while read SAMPLE BAM REST
do 
    echo $SAMPLE
    ROOT_HIST=$LUMPY_DIR/$SAMPLE/temp/cnvnator-temp/${BAM##*/}.hist.root
    if [ ! -s $ROOT_HIST ]
    then
        echo "$ROOT_HIST not found\n" 1>&2
        exit 1
    fi
    GT_VCF="$DIR/post-merge/gt/$SAMPLE.gt.vcf"

    if [ -s $GT_VCF ]
        then 
            if [ ! -s $COORDINATES ]
            then
                $CREATE_COORDS -i $GT_VCF -o $COORDINATES
            fi
            
            OUTFILE="$OUTDIR/$SAMPLE.cn.vcf.gz"
            if [ ! -s $OUTFILE ]
            then
                bsub \
                    -J $SAMPLE.cn \
                    -oo $LOGDIR/$SAMPLE.cn.log \
                "bash $BIN_DIR/copynumber_sample.sh $SAMPLE $COORDINATES $GT_VCF $ROOT_HIST $OUTFILE"
            fi
     fi
done < $SAMPLEMAP
