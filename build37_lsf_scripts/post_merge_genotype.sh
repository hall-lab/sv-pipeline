#!/bin/bash

set -eo pipefail

# DEPENDENCIES
# + svtyper
# + awk
# + sed

# SET SVTYPER PATH HERE:
# Edit this to point to your svtyper install if it's not currently on your path
SVTYPER="svtyper"

# SET SAMPLE NAME TO PROCESS
SAMPLE="$1"

# SET BAM or CRAM FILE FOR SAMPLE
BAM="$2"

# SET MASTER VCF CONTAINING ALL SITES ACROSS COHORT
MASTER="$3"

# OUTPUT FILE LOCATION
# Is not bgzipped
OUTPUT="$4"

zcat $MASTER \
  | cut -f-8 \
  | awk 'BEGIN{FS=OFS="\t"} { $6 = "."; print }' \
  | $SVTYPER \
    -B $BAM \
  | sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
  | sed 's/SNAME1=[0-9A-Za-z_:\.e,-/]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
  > $OUTPUT
