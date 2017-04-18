import "../../scripts/SV_Tasks.wdl" as SV

workflow Test_SVTools {
  # data inputs
  Array[File] input_pre_merged_vcfs
  Array[File] input_post_merged_vcfs
  File pedigree_file
  String cohort_name
  String final_vcf_name

  # reference inputs
  File mei_annotation_bed

  # system inputs
  Int disk_size
  Int preemptible_tries

  call SV.L_Sort_VCF_Variants {
    input:
    input_vcfs = input_pre_merged_vcfs,
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

  call SV.Paste_VCF {
    input:
    input_vcfs = input_post_merged_vcfs,
    output_vcf_basename = cohort_name + ".merged.gt.cn",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call SV.Prune_VCF {
    input:
    input_vcf_gz = Paste_VCF.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call SV.Classify {
    input:
    input_vcf_gz = Prune_VCF.output_vcf_gz,
    input_ped = pedigree_file,
    mei_annotation_bed = mei_annotation_bed,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned.class",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call SV.Sort_Index_VCF {
    input:
    input_vcf_gz = Classify.output_vcf_gz,
    output_vcf_name = final_vcf_name,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }
}
