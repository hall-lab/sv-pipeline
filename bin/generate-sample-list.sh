#!/bin/bash 

# this is the original data used in the 20200225 call set
# https://github.com/hall-lab/gatk4-cromwell-starter/raw/wustl-gatk4-2020-02-25-call-set/master-gvcf-list.tsv.gz

# curl -L -O https://github.com/hall-lab/gatk4-cromwell-starter/raw/wustl-gatk4-2020-02-25-call-set/master-gvcf-list.tsv.gz

gzcat master-gvcf-list.tsv.gz \
    | ggrep -P '\schr1\s' \
    | while IFS= read -r line; do 
        sample=$(echo ${line}| cut -d\  -f1); 
        gvcf=$(echo ${line}| cut -d\  -f3); 
        cohort=$(basename $(dirname $(dirname ${gvcf})))
        printf "%s\t%s\t%s\n" ${cohort} ${sample}
    done \
    | perl -lane '$_ =~ s/^Sample_\w+/new-1000G/g; print $_' \
    | tee data/20200225-gatk-callset-sample-data.tsv

# rm master-gvcf-list.tsv.gz
