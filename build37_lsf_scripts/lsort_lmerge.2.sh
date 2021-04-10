#!/bin/bash

set -eo pipefail

# DEPENDENCIES
# + svtools
# + bgzip

# SET SVTOOLS PATH HERE:
# Edit this to point to your svtools install if it's not currently on your path
SVTOOLS="svtools"

# SET PATH TO LIST OF INPUT FILES FOR BATCH
# Each line should be the path to a VCF file for a sample
BATCH_LIST="$1"

# NAME OF OUTPUT FILE FROM LSORT. WILL BE BGZIPPED.
SORTFILE="$2"

# NAME OF OUTPUT FILE FROM LMERGE. WILL BE BGZIPPED.
MERGEFILE="$3"

$SVTOOLS lsort -f $BATCH_LIST | bgzip -c > $SORTFILE

zcat $SORTFILE \
     | $SVTOOLS lmerge -i /dev/stdin -f 20 -w carrier_wt \
     | bgzip -c > $MERGEFILE
