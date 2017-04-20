
# Cohort SV detection pipeline

# Table of contents
1. [Overview](#overview)
2. [WDL scripts](#wdl-scripts)
3. [Dockerfiles](#dockerfiles)

# Overview
This repository contains pipeline scripts for structural variation detection in large cohorts. The pipeline is designed for Illumina paired-end whole genome sequencing data, preferably with at least 30x sequence coverage. Data inputs should be a set of sorted CRAM files, aligned with BWA-MEM.

This pipeline detects structural variation based on breakpoint sequence evidence using the LUMPY algorithm. Structural variant (SV) breakpoints are then unified and merged using the [SVTools](https://github.com/hall-lab/svtools) workflow, followed by re-genotyping with [SVTyper](https://github.com/hall-lab/svtyper) and read-depth annotation with [CNVnator](https://github.com/abyzovlab/CNVnator). Finally, SV types are reclassified based on the concordance between read-depth and breakpoint genotype.

Additional details on the SVTools pipeline are available in the [SVTools tutorial](https://github.com/hall-lab/svtools/blob/master/Tutorial.md).

![Workflow](images/workflow.wdl.v04.low-01.png?raw=true "Workflow")

# WDL scripts

Pipeline scripts (in [WDL format](https://software.broadinstitute.org/wdl/)) are available in the [scripts](scripts) directory.

While the SV pipeline can be run in its entirety via the [SV_Pipeline_Full.wdl](scripts/SV_Pipeline_Full.wdl) script, we recommend running the pipeline in three stages to enable intermediate quality control checkpoints.

## 1. [Pre_Merge_SV.wdl](scripts/Pre_Merge_SV.wdl)

For each sample:
  - Generate LUMPY inputs from an aligned and sorted CRAM file
  - SV discovery with LUMPY
  - Preliminary SV genotyping with SVTyper
  - Generate CNVnator histogram files

At this stage, we recommend performing quality control checks on each sample before merging them into the cohort-level VCF (step 2).

## 2. [Merge_SV.wdl](scripts/Merge_SV.wdl)

Merge

## 3. [Post_Merge_SV.wdl](scripts/Post_Merge_SV.wdl)

Post-merge

# Dockerfiles

Generator scripts for the requisite Docker containers.



