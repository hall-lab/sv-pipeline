import "SV_Tasks.wdl" as SV

workflow Genotype_only {
  # data inputs
  Array[File] aligned_crams
  Array[File] aligned_cram_indices
  Array[File] input_vcfs
  String aligned_cram_suffix

  # reference inputs
  File ref_cache

  # system inputs
  Int disk_size
  Int preemptible_tries

  scatter (i in range(length(aligned_crams))) {
    File aligned_cram = aligned_crams[i]
    File aligned_cram_index = aligned_cram_indices[i]
    File input_vcf = input_vcfs[i]

    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")
    
    call SV.Genotype as Genotype_Unmerged {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = input_vcf,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  output {
    Genotype_Unmerged.output_vcf
  }
}
