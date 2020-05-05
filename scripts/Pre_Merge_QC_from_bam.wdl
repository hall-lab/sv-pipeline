import "Pre_Merge_QC_per_sample_from_bam.wdl" as per_sample



workflow Pre_Merge_QC {
  Array[File] aligned_crams
  String cohort
  String center


  # system inputs
  Int preemptible_tries

  scatter (i in range(length(aligned_crams))) {
    File aligned_cram = aligned_crams[i]

    call per_sample.Pre_Merge_QC_per_sample_from_bam {
      input:
        aligned_cram = aligned_cram,
        cohort = cohort ,
        center = center ,
        preemptible_tries = preemptible_tries
     }
  } 
    
  output {
    Array[File] lumpy_counts = Pre_Merge_QC_per_sample_from_bam.lumpy_counts
  }
}
