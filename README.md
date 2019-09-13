
# Cohort SV detection pipeline

# Table of contents
1. [Overview](#overview)
2. [WDL scripts](#wdl-scripts)
3. [Docker images](#docker-images)

# Overview
This repository contains pipeline scripts for structural variation detection in large cohorts. The pipeline is designed for Illumina paired-end whole genome sequencing data, preferably with at least 30x sequence coverage. Data inputs should be a set of sorted CRAM files, aligned with BWA-MEM.

This pipeline detects structural variation based on breakpoint sequence evidence using both the LUMPY and Manta algorithms. Structural variant (SV) breakpoints are then unified and merged using the [SVTools](https://github.com/hall-lab/svtools) workflow, followed by re-genotyping with [SVTyper](https://github.com/hall-lab/svtyper) and read-depth annotation with [CNVnator](https://github.com/abyzovlab/CNVnator). Finally, SV types are reclassified based on the concordance between read-depth and breakpoint genotype.

Additional details on the SVTools pipeline are available in the [SVTools tutorial](https://github.com/hall-lab/svtools/blob/master/Tutorial.md).

![Workflow](images/workflow.wdl.v04.low-01.png?raw=true "Workflow")

# WDL scripts

Pipeline scripts (in [WDL format](https://software.broadinstitute.org/wdl/)) are available in the [scripts](scripts) directory. These scripts can be launched using [Cromwell](https://github.com/broadinstitute/cromwell) (version 25 or later).

While the SV pipeline can be run in its entirety via the [SV_Pipeline_Full.wdl](scripts/SV_Pipeline_Full.wdl) script, we recommend running the pipeline in three stages to enable intermediate quality control checkpoints.

## 1. [Pre_Merge_SV.wdl](scripts/Pre_Merge_SV.wdl)

For each sample:
  - SV discovery with LUMPY using the [smoove](https://github.com/brentp/smoove) wrapper
  - Preliminary SV genotyping with SVTyper (also done within the smoove wrapper)
  - SV discovery with [Manta](https://github.com/Illumina/manta), including insertions
  - Generate [CNVnator](https://github.com/abyzovlab/CNVnator) histogram files

After this step, we recommend performing quality control checks on each sample before merging them into the cohort-level VCF (step 2).

## 2. [Merge_SV.wdl](scripts/Merge_SV.wdl)

This step merges the sample-level VCF files from step 1 using the LUMPY breakpoint probability curves to produce a single cohort-level VCF.

## 3. [Post_Merge_SV.wdl](scripts/Post_Merge_SV.wdl)

This step re-genotypes each sample at the sites in the cohort-level VCF file from step 2, and then combines the results into a set of final VCFs, split by variant type for efficiency (deletions, insertions, breakends, and other:duplications+inversions).

For each sample:
  - Re-genotype each SV using SVTyper (note that insertion calls from Manta are taken from the per-sample genotypes and not processed with SVTyper)
  - Annotate the read-depth at each SV using CNVnator
  - Generate a .ped file of sample names and sexes

For the cohort:
  - Combine the re-genotyped VCFs into a single cohort-level VCF
  - Prune overlapping SVs
  - Classify SV type based on the concordance between variant genotypes and read-depths
  - Sort and index the VCF

# Docker images

- Docker images for this pipeline are available at https://hub.docker.com/u/halllab.
- Dockerfiles for these containers are available in the [docker](docker) directory.
- WDL test scripts for each of these Docker containers are available in the [test](test) directory.
