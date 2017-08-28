import "SV_Tasks.wdl" as SV

workflow Post_Merge_Cn_only_SV {
  # data inputs
  Array[File] aligned_crams
  String aligned_cram_suffix
  Array[File] cn_hist_roots
  Array[File] genotype_vcf_gzs
  File coordinates

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache

  # system inputs
  Int disk_size
  Int preemptible_tries

  # Re-genotype and call copy number for each sample on the merged SV VCF
  scatter (i in range(length(aligned_crams))) {
    
    File aligned_cram = aligned_crams[i]
    File cn_hist_root = cn_hist_roots[i]
    File genotype_vcf_gz = genotype_vcf_gzs[i]
    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")

    call SV.Get_Sample_Name {
      input:
      input_cram = aligned_cram,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }

    call SV.Copy_Number {
      input:
      basename = basename,
      sample = Get_Sample_Name.sample,
      input_vcf = genotype_vcf_gz,
      coordinates = coordinates,
      input_cn_hist_root = cn_hist_root,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  output {
    Copy_Number.output_vcf_gz
  }
}
