import "SV_Tasks.wdl" as SV

workflow Post_Merge_Sx_only_SV {
  # data inputs
  Array[File] cn_hist_roots

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache

  # system inputs
  Int disk_size
  Int preemptible_tries

  # Re-genotype and call copy number for each sample on the merged SV VCF
  scatter (i in range(length(cn_hist_roots))) {
    
    File cn_hist_root = cn_hist_roots[i]

    call SV.Get_Sex {
      input:
      input_cn_hist_root = cn_hist_root,
      ref_fasta_index = ref_fasta_index,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }
}
