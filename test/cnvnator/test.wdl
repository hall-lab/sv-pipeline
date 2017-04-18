import "../../scripts/SV_Tasks.wdl" as SV

workflow Test_Copy_Number {
  # data inputs
  String basename
  String sample
  File input_cram
  File input_cram_index
  File input_vcf

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache

  # system inputs
  Int disk_size
  Int preemptible_tries

  # -----------------------------------
  # test CNVnator
  call SV.CNVnator_Histogram {
    input:
    basename = basename,
    input_cram = input_cram,
    input_cram_index = input_cram_index,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_cache = ref_cache,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call SV.Copy_Number {
    input:
    basename = basename,
    sample = sample,
    input_vcf = input_vcf,
    input_cn_hist_root = CNVnator_Histogram.output_cn_hist_root,
    ref_cache = ref_cache,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }
  
  # ------------------------------------
  # generate .ped file

  Array[File] aligned_crams
  Array[File] cn_hist_roots
  String cohort_name

  scatter (i in range(length(aligned_crams))) {
    File aligned_cram = aligned_crams[i]
    File cn_hist_root = cn_hist_roots[i]
    
    call SV.Get_Sample_Name {
      input:
      input_cram = aligned_cram,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }

    call SV.Get_Sex {
      input:
      input_cn_hist_root = cn_hist_root,
      ref_fasta_index = ref_fasta_index,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }
  
  call SV.Make_Pedigree_File {
    input:
    sample_array = Get_Sample_Name.sample,
    sex_array = Get_Sex.sex,
    output_ped_basename = cohort_name,
    disk_size = 1
  }
}
