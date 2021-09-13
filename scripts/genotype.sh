#!/bin/bash

input_cram=$1  
input_cram_index=$2
basename=$3
ref_cache=$4
input_vcf=$5

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

log "Symlinking ${input_cram} to ${basename}.cram"
ln -s ${input_cram} ${basename}.cram

log "Symlinking ${input_cram_index} to ${basename}.cram.crai"
ln -s ${input_cram_index} ${basename}.cram.crai

log "Setting up the reference sequence cache"
# build the reference sequence cache
tar -zxf ${ref_cache}
export REF_PATH=./cache/%2s/%2s/%s
export REF_CACHE=./cache/%2s/%2s/%s

if [[ -e "${basename}.cram.json" ]]; then
    log "Deleting existing ${basename}.cram.json"
    rm -f ${basename}.cram.json
fi

log "Start running svtyper command"

#(set -exo pipefail; \
# zcat ${input_vcf} \
#    | svtyper \
#      -B ${basename}.cram \
#      -l ${basename}.cram.json \
#    | bgzip -c \
#    > ${basename}.gt.vcf.gz) 2>/dev/null 1>/dev/null

(exit 137;)

rc=$?

log "Finished running svtyper command"
log "original svtyper command return code is: ${rc}"

if [[ ${rc} == "137" ]]
then
  echo "OutOfMemory" >&2 ;
  exit 1;
else
  exit ${rc};
fi
