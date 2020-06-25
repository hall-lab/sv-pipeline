# External Data Information

## `cram-list` directory

This is a set of crams identified by @kkanchi to run the pre-merge pipeline on.  They were originally placed here on the MGI cluster:

    /gscmnt/gc2802/halllab/dlarson/jira/BIO-2169_Build38_Realign/Callset_2020/crams_input

### Cohort Breakdowns

      2144 Afib_downsample_crams.tsv
      3363 Brazilian.cram.tsv
       545 Corogene_CAD.crams.tsv
       303 Costarican.cramlist.persamplesv.tsv
       430 Dyslipidemia.crams.tsv
      2956 Finrisk.crams.tsv
       465 Finrisk.downsample_crams.tsv
        30 Indiana.cramlist.persamplesv.tsv
      2036 MEC.cram.tsv
      2817 Metsim.cram.persamplesv.tsv
       517 SCCS.cramlist.persamplesv.tsv
     15606 total


Removed 9 Costarican samples, since they failed GATK QC and also causing error completing the premerge-SV due to lots of 'Abnormal range' in the CNVnator_Histogram.log file and also low coverage
H_WX-1001-1001
H_WX-1087-1087
H_WX-1115-1115
H_WX-1119-1119
H_WX-1152-1152
H_WX-1176-1176
H_WX-1220-1220
H_WX-705-705
H_WX-903-903

Corogene- 5 samples
Dyslipidemia -7
SCCS -1
Brazilian -2
Afib -4
Finrisk -4

### 1000G cram generation list

This was the command used to generate the 1000G parent cram tsv list:

      BUCKET_PATH="gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_13607/Project_CCDG_13607_B01_GRM_WGS.cram.2019-02-06"
      for p in $(gsutil -u washu-genome-inh-dis-analysis ls ${BUCKET_PATH}); do 
          sample=$(basename ${p} | awk -F_ '{print $2;}' );
          echo "${p}analysis/${sample}.final.cram";
      done >data/external/cram-list/1000G.parents.cramlist.tsv

This was the command used to generate the 1000G children cram tsv list:

      BUCKET_PATH="gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_14151/Project_CCDG_14151_B01_GRM_WGS.cram.2020-02-12/"
      for p in $(gsutil -u washu-genome-inh-dis-analysis ls ${BUCKET_PATH}); do
          sample=$(basename ${p} | awk -F_ '{print $2;}' );
          echo "${p}analysis/${sample}.final.cram";
      done >data/external/cram-list/1000G.children.cramlist.tsv

### Metsim downsample cram generation list

This was the command used to generate the Metsim downsampled cram tsv list:

     BUCKET_PATH="gs://wustl-ccdg-cram-data/METSIM_downsample_crams/"
     for p in $(gsutil -u washu-genome-inh-dis-analysis ls ${BUCKET_PATH}); do
         sample=$(basename ${p});
         echo "${p}${sample}.cram";
     done >data/external/cram-list/Metsim.downsample.cram.persamplesv.tsv

The actual downsampling was handled by @tmooney.
