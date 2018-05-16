#!/bin/bash

set -eo pipefail

# DEPENDENCIES
# + svtools
# + bgzip
# + sed
# + split

# SET SVTOOLS PATH HERE:
# Edit this to point to your svtools install if it's not currently on your path
SVTOOLS="svtools"

# SET PATH AND NAME TO OUTPUT FILE FROM LMERGE
MERGE_FILE="$1"

# SET PATH TO A LIST OF FILES TO MERGE
# This should be a file containing one line per post-copynumber annotated file
# it will be split for you
LIST_FILE="$2"

# SET NAME OF THE FINAL OUTPUT FILE
FINAL_OUTPUT="cohort.vcf.gz"

if [[ ! -d "paste" ]]
then
    mkdir -p paste.tmp
    mkdir -p log
    split -a 4 -d -l 200 "${LIST_FILE}" paste.tmp/spl && mv paste.tmp paste
fi

SUBS=final_paste.txt
DEPS=''

for FILE in `ls paste`
do
    INPUT="`pwd`/paste/$FILE"
    OUTPUT="`pwd`/paste/$FILE.vcf.gz"
    if [[ ! -s "${OUTPUT}" ]]
    then
        cmd="${SVTOOLS} vcfpaste -q -m ${MERGE_FILE} -f ${INPUT} | bgzip -c > ${OUTPUT}.tmp"
        bsub \
            -J ${FILE}.paste \
            -oo "`pwd`/log/${FILE}.log" \
            -M 13000000 \
            -R 'select[mem>13000 && tmp>100] rusage[mem=13000, tmp=100] span[hosts=1]' \
            bash -c "$cmd && mv ${OUTPUT}.tmp ${OUTPUT}"
    fi
    echo "$OUTPUT" >> $SUBS
    DEPS="$DEPS $FILE.paste"
done

DEPS=($DEPS)

DEP_STRING="done(${DEPS[0]})"
unset DEPS[0]
for MORE in ${DEPS[*]}
do
    DEP_STRING="$DEP_STRING && done($MORE)"
done

final_cmd="${SVTOOLS} vcfpaste -m ${MERGE_FILE} -q -f $SUBS | svtools vcfsort | bgzip -c > ${FINAL_OUTPUT}.tmp"
bsub \
    -w "${DEP_STRING}" \
    -oo "`pwd`/log/final_paste.log" \
    -M 13000000 \
    -R 'select[mem>13000 && tmp>100] rusage[mem=13000, tmp=100] span[hosts=1]' \
    bash -c "$final_cmd && mv ${FINAL_OUTPUT}.tmp ${FINAL_OUTPUT}"
