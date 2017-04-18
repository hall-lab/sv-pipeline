import "SV_Tasks.wdl" as SV

# SV detection workflow
workflow SV_Detect {
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

  scatter (aligned_cram in aligned_crams) {

    String basename = sub(sub(aligned_cram, "gs://.*/", ""), aligned_cram_suffix + "$", "")
    
    call Extract_Reads {
      input:
      input_cram = aligned_cram,
      basename = basename,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call Lumpy {
      input:
      basename = basename,
      input_cram = Extract_Reads.output_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      input_splitters_bam = Extract_Reads.output_splitters_bam,
      input_splitters_bam_index = Extract_Reads.output_splitters_bam_index,
      input_discordants_bam = Extract_Reads.output_discordants_bam,
      input_discordants_bam_index = Extract_Reads.output_discordants_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      exclude_regions = exclude_regions,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call SV_Genotype as SV_Genotype_Unmerged {
      input:
      basename = basename,
      input_cram = Extract_Reads.output_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      input_vcf = Lumpy.output_vcf,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call CNVnator_Histogram {
      input:
      basename = basename,
      input_cram = Extract_Reads.output_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }

    call Get_Sample_Name {
      input:
      input_cram = Extract_Reads.output_cram,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }

    call Get_Sex {
      input:
      input_cn_hist_root = CNVnator_Histogram.output_cn_hist_root,
      ref_fasta_index = ref_fasta_index,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  call Make_Pedigree_File {
    input:
    sample_array = Get_Sample_Name.sample,
    sex_array = Get_Sex.sex,
    output_ped_basename = cohort_name,
    disk_size = 1    
  }

  call L_Sort_VCF_Variants {
    input:
    input_vcfs = SV_Genotype_Unmerged.output_vcf,
    output_vcf_basename = cohort_name + ".lsort",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call L_Merge_VCF_Variants {
    input:
    input_vcf_gz = L_Sort_VCF_Variants.output_vcf_gz,
    output_vcf_basename = cohort_name + ".lmerge",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  # Re-genotype and call copy number for each sample on the merged SV VCF
  scatter (i in range(length(Extract_Reads.output_cram))) {
    
    File aligned_cram = Extract_Reads.output_cram[i]
    File aligned_cram_index = Extract_Reads.output_cram_index[i]
    String basename = sub(sub(aligned_cram, "gs://.*/", ""), aligned_cram_suffix + "$", "")

    call SV_Genotype as SV_Genotype_Merged {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = L_Merge_VCF_Variants.output_vcf_gz,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }

    call SV_Copy_Number {
      input:
      basename = basename,
      sample = basename,
      input_vcf = SV_Genotype_Merged.output_vcf,
      input_cn_hist_root = CNVnator_Histogram.output_cn_hist_root[i],
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  call Paste_VCF {
    input:
    input_vcfs = SV_Copy_Number.output_vcf,
    output_vcf_basename = cohort_name + ".merged.gt.cn",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call Prune_VCF {
    input:
    input_vcf_gz = Paste_VCF.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call SV_Classify {
    input:
    input_vcf_gz = Prune_VCF.output_vcf_gz,
    input_ped = Make_Pedigree_File.output_ped,
    mei_annotation_bed = mei_annotation_bed,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned.class",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call Sort_Index_VCF {
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
