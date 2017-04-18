import "SV_Detect.tasks.wdl" as SV

workflow Merge_Sample_SVs {
  # data inputs
  Array[File] aligned_crams
  String aligned_cram_suffix
  String cohort_name
  String final_vcf_name

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache
  File exclude_regions
  File mei_annotation_bed

  # system inputs
  Int disk_size
  Int preemptible_tries

  call SV.L_Sort_VCF_Variants {
    input:
    input_vcfs = Genotype_Unmerged.output_vcf,
    output_vcf_basename = cohort_name + ".lsort",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call SV.L_Merge_VCF_Variants {
    input:
    input_vcf_gz = L_Sort_VCF_Variants.output_vcf_gz,
    output_vcf_basename = cohort_name + ".lmerge",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }


  # Outputs that will be retained when execution is complete
  output {
  }
}
