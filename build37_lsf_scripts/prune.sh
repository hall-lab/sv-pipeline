#!/bin/bash

set -eo pipefail

# SET SVTOOLS PATH HERE:
# Edit this to point to your svtools install if it's not currently on your path
SVTOOLS="svtools"

# SET INPUT VCF.
# Must be gzipped.
INPUT="$1"

# SET OUTPUT VCF.
# Will be gzipped
OUTPUT="$2"

if [[ ! -s "${OUTPUT}" ]]
then
    zcat ${INPUT} \
        | ${SVTOOLS} afreq \
        | ${SVTOOLS} vcftobedpe \
        | ${SVTOOLS} bedpesort \
        | ${SVTOOLS} prune -s -d 100 -e "AF" \
        | ${SVTOOLS} bedpetovcf \
        | ${SVTOOLS} vcfsort \
        | bgzip -c > "${OUTPUT}.tmp" \
        && mv "${OUTPUT}.tmp" "${OUTPUT}"
fi
