import "SV_Tasks.wdl" as SV

workflow Merge_SV {
  # data inputs
  Array[File] input_vcfs
  String cohort_name

  # system inputs
  Int disk_size
  Int preemptible_tries

  call SV.L_Sort_VCF_Variants {
    input:
    input_vcfs = input_vcfs,
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

  output {
    L_Merge_VCF_Variants.output_vcf_gz
  }
}
