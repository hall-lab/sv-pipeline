import "Pre_Merge_QC_per_sample.wdl" as per_sample

workflow Pre_Merge_QC_Test {
  Array[File] lumpy_vcfs
  Array[File] manta_vcfs
  Array[File] cnvnator_vcfs
  String cohort
  String center


  # system inputs
  Int preemptible_tries

  scatter (i in range(length(lumpy_vcfs))) {
    File lumpy_vcf = lumpy_vcfs[i]
    File manta_vcf = manta_vcfs[i]
    File cnvnator_vcf = cnvnator_vcfs[i]

    call per_sample.Pre_Merge_QC_Per_Sample {
      input:
        lumpy_vcf = lumpy_vcf,
        manta_vcf = manta_vcf,
        cnvnator_vcf = cnvnator_vcf,
        cohort = cohort ,
        center = center ,
        preemptible_tries = preemptible_tries
     }
  }

  output {

    Array[File] lumpy_counts = Pre_Merge_QC_Per_Sample.lumpy_counts
    Array[File] manta_counts = Pre_Merge_QC_Per_Sample.manta_counts
    Array[File] cnvnator_counts = Pre_Merge_QC_Per_Sample.cnvnator_counts

  }
}
