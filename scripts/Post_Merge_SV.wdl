version 1.0
import "https://raw.githubusercontent.com/hall-lab/sv-pipeline/terra-compatible/scripts/SV_Tasks.wdl" as SV

workflow Post_Merge_SV {
  # data inputs
  input {
  	Array[File] aligned_crams
  	Array[File] aligned_cram_indices
  	Array[File] cn_hist_roots
  	Array[File] manta_vcfs
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
  	Int preemptible_tries
  }

  call SV.Split_By_Type {
    input:
    input_vcf = merged_vcf,
    output_vcf_prefix = cohort_name + ".merged",
    preemptible_tries = preemptible_tries
  }

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

    call SV.Genotype as Genotype_Merged_BND {
      input:
      basename = basename + ".bnd",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = Split_By_Type.bnd_vcf,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged_DEL {
      input:
      basename = basename + ".del",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = Split_By_Type.del_vcf,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Take_Original_Genotypes as Genotype_Merged_INS {
      input:
      sample_name = Get_Sample_Name.sample,
      original_per_sample_vcf = manta_vcfs[i],
      basename = basename + ".ins",
      input_vcf = Split_By_Type.ins_vcf,
      input_variant_to_sname_mapping = Split_By_Type.ins_split,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged_OTHER {
      input:
      basename = basename + ".other",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = Split_By_Type.other_vcf,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Copy_Number as Copy_Number_DEL {
      input:
      basename = basename + ".del",
      sample = Get_Sample_Name.sample,
      input_vcf = Genotype_Merged_DEL.output_vcf,
      input_cn_hist_root = cn_hist_root,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Copy_Number as Copy_Number_OTHER {
      input:
      basename = basename + ".other",
      sample = Get_Sample_Name.sample,
      input_vcf = Genotype_Merged_OTHER.output_vcf,
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

  call SV.Paste_VCF as Paste_VCF_BND {
    input:
    input_vcfs = Genotype_Merged_BND.output_vcf,
    output_vcf_basename = cohort_name + ".merged.gt.bnd",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_DEL {
    input:
    input_vcfs = Copy_Number_DEL.output_vcf,
    output_vcf_basename = cohort_name + ".merged.gt.cn.del",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_INS {
    input:
    input_vcfs = Genotype_Merged_INS.output_vcf,
    output_vcf_basename = cohort_name + ".merged.gt.ins",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_OTHER {
    input:
    input_vcfs = Copy_Number_OTHER.output_vcf,
    output_vcf_basename = cohort_name + ".merged.gt.cn.other",
    preemptible_tries = preemptible_tries
  }

  call SV.Prune_VCF as Prune_VCF_BND{
    input:
    input_vcf_gz = Paste_VCF_BND.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.pruned.bnd",
    preemptible_tries = preemptible_tries
  }

  call SV.Prune_VCF as Prune_VCF_DEL{
    input:
    input_vcf_gz = Paste_VCF_DEL.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned.del",
    preemptible_tries = preemptible_tries
  }

  call SV.Prune_VCF as Prune_VCF_INS{
    input:
    input_vcf_gz = Paste_VCF_INS.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.pruned.ins",
    preemptible_tries = preemptible_tries
  }

  call SV.Prune_VCF as Prune_VCF_OTHER{
    input:
    input_vcf_gz = Paste_VCF_OTHER.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned.other",
    preemptible_tries = preemptible_tries
  }

  call SV.Classify as Classify_DEL{
    input:
    input_vcf_gz = Prune_VCF_DEL.output_vcf_gz,
    input_ped = Make_Pedigree_File.output_ped,
    mei_annotation_bed = mei_annotation_bed,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned.class.del",
    preemptible_tries = preemptible_tries
  }

  call SV.Classify as Classify_OTHER{
    input:
    input_vcf_gz = Prune_VCF_OTHER.output_vcf_gz,
    input_ped = Make_Pedigree_File.output_ped,
    mei_annotation_bed = mei_annotation_bed,
    output_vcf_basename = cohort_name + ".merged.gt.cn.pruned.class.other",
    preemptible_tries = preemptible_tries
  }

  call SV.Sort_Index_VCF as Sort_Index_VCF_BND {
    input:
    input_vcf_gz = Prune_VCF_BND.output_vcf_gz,
    output_vcf_name = final_vcf_name + ".bnd.vcf.gz",
    preemptible_tries = preemptible_tries
  }

  call SV.Sort_Index_VCF as Sort_Index_VCF_DEL {
    input:
    input_vcf_gz = Classify_DEL.output_vcf_gz,
    output_vcf_name = final_vcf_name + ".del.vcf.gz",
    preemptible_tries = preemptible_tries
  }

  call SV.Sort_Index_VCF as Sort_Index_VCF_INS {
    input:
    input_vcf_gz = Prune_VCF_INS.output_vcf_gz,
    output_vcf_name = final_vcf_name + ".ins.vcf.gz",
    preemptible_tries = preemptible_tries
  }

  call SV.Sort_Index_VCF as Sort_Index_VCF_OTHER {
    input:
    input_vcf_gz = Classify_OTHER.output_vcf_gz,
    output_vcf_name = final_vcf_name + ".other.vcf.gz",
    preemptible_tries = preemptible_tries
  }

  call SV.Filter_Index as Filter_Index_BND {
    input:
    input_vcf_gz = Sort_Index_VCF_BND.output_vcf_gz,
    output_vcf_name = final_vcf_name + ".bnd.vcf.gz",
    preemptible_tries = preemptible_tries
  }

  call SV.Filter_Index as Filter_Index_DEL {
    input:
    input_vcf_gz = Sort_Index_VCF_DEL.output_vcf_gz,
    output_vcf_name = final_vcf_name + ".del.vcf.gz",
    preemptible_tries = preemptible_tries
  }

  call SV.Filter_Index as Filter_Index_INS {
    input:
    input_vcf_gz = Sort_Index_VCF_INS.output_vcf_gz,
    output_vcf_name = final_vcf_name + ".ins.vcf.gz",
    preemptible_tries = preemptible_tries
  }

  call SV.Filter_Index as Filter_Index_OTHER {
    input:
    input_vcf_gz = Sort_Index_VCF_OTHER.output_vcf_gz,
    output_vcf_name = final_vcf_name + ".other.vcf.gz",
    preemptible_tries = preemptible_tries
  }

  output {
    File output_ped = Make_Pedigree_File.output_ped
    File output_vcf_bnd = Filter_Index_BND.output_vcf_gz
    File output_vcf_index_bnd = Filter_Index_BND.output_vcf_gz_index
    File output_vcf_del = Filter_Index_DEL.output_vcf_gz
    File output_vcf_ins = Filter_Index_INS.output_vcf_gz
    File output_vcf_index_other = Filter_Index_OTHER.output_vcf_gz_index
    File output_vcf_other = Filter_Index_OTHER.output_vcf_gz
    File output_vcf_index_del = Filter_Index_DEL.output_vcf_gz_index
    File output_vcf_index_ins = Filter_Index_INS.output_vcf_gz_index
  }
}
