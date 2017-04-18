import "SV_Detect.tasks.wdl" as SV

workflow Post_Merge_SV_Genotype {
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

# Re-genotype and call copy number for each sample on the merged SV VCF
  scatter (i in range(length(Extract_Reads.output_cram))) {
    
    File aligned_cram = Extract_Reads.output_cram[i]
    File aligned_cram_index = Extract_Reads.output_cram_index[i]
    String basename = sub(sub(aligned_cram, "gs://.*/", ""), aligned_cram_suffix + "$", "")

    call Genotype as Genotype_Merged {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = L_Merge_VCF_Variants.output_vcf_gz,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }

    call SV.Copy_Number {
      input:
      basename = basename,
      sample = basename,
      input_vcf = Genotype_Merged.output_vcf,
      input_cn_hist_root = CNVnator_Histogram.output_cn_hist_root[i],
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  call SV.Paste_VCF {
    input:
    input_vcfs = SV_Copy_Number.output_vcf,
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
    input_ped = Make_Pedigree_File.output_ped,
    mei_annotation_bed = mei_annotation_bed,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned.class",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call SV.Sort_Index_VCF {
    input:
    input_vcf_gz = SV_Classify.output_vcf_gz,
    output_vcf_name = final_vcf_name,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    Make_Pedigree_File.*
    SV_Genotype_Unmerged.output_vcf
    SV_Genotype_Unmerged.output_lib
    CNVnator_Histogram.*
    L_Merge_VCF_Variants.*
    Paste_VCF.*
    Prune_VCF.*
    Sort_Index_VCF.*
  }
}
