
# Cohort SV detection pipeline

# Table of contents
1. [Overview](#overview)
2. [WDL scripts](#wdl-scripts)
3. [Dockerfiles](#dockerfiles)

# Overview
This repository contains pipeline scripts for structural variation detection in large cohorts. The pipeline is designed for Illumina paired-end whole genome sequencing data, preferably with at least 30x sequence coverage. Data inputs should be a set of sorted CRAM files, aligned with BWA-MEM.

This pipeline detects structural variation based on breakpoint sequence evidence using the LUMPY algorithm. Structural variant (SV) breakpoints are then unified and merged using the [SVTools](https://github.com/hall-lab/svtools) workflow, followed by re-genotyping with [SVTyper](https://github.com/hall-lab/svtyper) and read-depth annotation with [CNVnator](https://github.com/abyzovlab/CNVnator). Finally, SV types are reclassified based on the concordance between read-depth and breakpoint genotype.

![Workflow](images/workflow.wdl.v04.low-01.png?raw=true "Workflow")

# WDL scripts

Pipeline scripts (in [WDL format](https://software.broadinstitute.org/wdl/)) are available in the [scripts](scripts) directory.

# Dockerfiles

Generator scripts for the requisite Docker containers.



