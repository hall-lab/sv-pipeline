version 1.0
import "SV_Tasks.wdl" as SV

workflow Post_Merge_SV {
  # data inputs
  input {
  	File cram_list
  	File index_list
  	File cn_hist_roots_list
  	String aligned_cram_suffix
  	File merged_vcf
  	String cohort_name
  	String final_vcf_name

  	# reference inputs
  	File ref_fasta
  	File ref_fasta_index
  	File ref_cache
  	File mei_annotation_bed

  	# system inputs
  	Int disk_size
  	Int preemptible_tries
  }
  Array[File] aligned_crams = read_lines(cram_list)
  Array[File] aligned_cram_indices = read_lines(index_list)
  Array[File] cn_hist_roots = read_lines(cn_hist_roots_list)

  # Re-genotype and call copy number for each sample on the merged SV VCF
  scatter (i in range(length(aligned_crams))) {
    
    File aligned_cram = aligned_crams[i]
    File aligned_cram_index = aligned_cram_indices[i]
    File cn_hist_root = cn_hist_roots[i]
    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")

    call SV.Get_Sample_Name {
      input:
      input_cram = aligned_cram,
      preemptible_tries = preemptible_tries
    }

    call SV.Get_Sex {
      input:
      input_cn_hist_root = cn_hist_root,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Copy_Number {
      input:
      basename = basename,
      sample = Get_Sample_Name.sample,
      input_vcf = Genotype_Merged.output_vcf,
      input_cn_hist_root = cn_hist_root,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }
  }
  
  call SV.Make_Pedigree_File {
    input:
    sample_array = Get_Sample_Name.sample,
    sex_array = Get_Sex.sex,
    output_ped_basename = cohort_name,
  }

  call SV.Paste_VCF {
    input:
    input_vcfs = Copy_Number.output_vcf,
    output_vcf_basename = cohort_name + ".merged.gt.cn",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call SV.Prune_VCF {
    input:
    input_vcf_gz = Paste_VCF.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned",
    preemptible_tries = preemptible_tries
  }

  call SV.Classify {
    input:
    input_vcf_gz = Prune_VCF.output_vcf_gz,
    input_ped = Make_Pedigree_File.output_ped,
    mei_annotation_bed = mei_annotation_bed,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned.class",
    preemptible_tries = preemptible_tries
  }

  call SV.Sort_Index_VCF {
    input:
    input_vcf_gz = Classify.output_vcf_gz,
    output_vcf_name = final_vcf_name,
    preemptible_tries = preemptible_tries
  }

  output {
    File output_ped = Make_Pedigree_File.output_ped
    File output_vcf = Sort_Index_VCF.output_vcf_gz
    File output_vcf_index = Sort_Index_VCF.output_vcf_gz_index
  }
}
