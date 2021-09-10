#!/bin/bash

input_cram=$1  
input_cram_index=$2
basename=$3
ref_cache=$4
input_vcf=$5


set -eo pipefail
ln -s ${input_cram} ${basename}.cram
ln -s ${input_cram_index} ${basename}.cram.crai

# build the reference sequence cache
tar -zxf ${ref_cache}
export REF_PATH=./cache/%2s/%2s/%s
export REF_CACHE=./cache/%2s/%2s/%s

rm -f ${basename}.cram.json
zcat ${input_vcf} \
  | svtyper \
  -B ${basename}.cram \
  -l ${basename}.cram.json \
  | bgzip -c > ${basename}.gt.vcf.gz


