import "SV_Tasks.wdl" as SV

workflow Cnvnator_only {
  # data inputs
  Array[File] aligned_crams
  Array[File] aligned_cram_indices
  String aligned_cram_suffix

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache

  # system inputs
  Int disk_size
  Int preemptible_tries

  scatter (i in range(length(aligned_crams))) {
    File aligned_cram = aligned_crams[i]
    File aligned_cram_index = aligned_cram_indices[i]
    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")
    
    call SV.CNVnator_Histogram {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  output {
    CNVnator_Histogram.output_cn_hist_root
  }
}
