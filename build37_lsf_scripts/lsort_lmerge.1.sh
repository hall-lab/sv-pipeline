#!/bin/bash

set -eo pipefail

# DEPENDENCIES
# + svtools
# + bgzip
# + awk
# + sed

# SET SVTOOLS PATH HERE:
# Edit this to point to your svtools install if it's not currently on your path
SVTOOLS="svtools"

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

$SVTOOLS lsort -f $BATCH_LIST | bgzip -c > $SORTFILE

zcat $SORTFILE \
     | $SVTOOLS lmerge -i /dev/stdin -f 20 \
     | awk -v bb="$BATCH" 'BEGIN{OFS="\t"}{
       if ($1~/##/) print $0;
       else if($1=="#CHROM") print $0, "FORMAT", bb;
       else if($0!~/NC_007605/) print $0, "GT", "1/1";
     }' \
     | sed 's/SNAME/SNAME1/g' \
     | bgzip -c > $MERGEFILE
