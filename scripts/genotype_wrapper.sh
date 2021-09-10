#!/bin/bash

input_cram=$1  
input_cram_index=$2
basename=$3
ref_cache=$4
input_vcf=$5
genotype_script=$6

$(/bin/bash ${genotype_script} ${input_cram} ${input_cram_index} ${basename} ${ref_cache} ${input_vcf})
rc=$?
if [[ ${rc} == "137" ]]
then
  echo "OutOfMemory" >&2 ;
  exit 1;
else
  echo "exit code ${rc}";
  exit ${rc};
fi
