version 1.0
import "SV_Tasks.wdl" as SV

workflow Post_Merge_SV_Batched {
  # data inputs
  input {
  	Array[File] aligned_crams
  	Array[File] aligned_cram_indices
  	Array[File] cn_hist_roots
  	Array[File] manta_vcfs
  	String aligned_cram_suffix
  	File merged_vcf_del
        File merged_vcf_bnd
        File merged_vcf_ins
        File merged_vcf_other
        File ins_split
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
      input_vcf = merged_vcf_bnd,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged_DEL {
      input:
      basename = basename + ".del",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf_del,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Take_Original_Genotypes as Genotype_Merged_INS {
      input:
      sample_name = Get_Sample_Name.sample,
      original_per_sample_vcf = manta_vcfs[i],
      basename = basename + ".ins",
      input_vcf = merged_vcf_ins,
      input_variant_to_sname_mapping = ins_split,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged_OTHER {
      input:
      basename = basename + ".other",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf_other,
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

    call SV.Remove_Sname as Remove_Sname_DEL {
      input:
      basename = basename + ".del",
      input_vcf = Copy_Number_DEL.output_vcf,
      preemptible_tries = preemptible_tries
    }

    call SV.Remove_Sname as Remove_Sname_OTHER {
      input:
      basename = basename + ".other",
      input_vcf = Copy_Number_OTHER.output_vcf,
      preemptible_tries = preemptible_tries
    }

    call SV.Remove_Sname as Remove_Sname_INS {
      input:
      basename = basename + ".ins",
      input_vcf = Genotype_Merged_INS.output_vcf,
      preemptible_tries = preemptible_tries
    }

    call SV.Remove_Sname as Remove_Sname_BND {
      input:
      basename = basename + ".bnd",
      input_vcf = Genotype_Merged_BND.output_vcf,
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
    input_vcfs = Remove_Sname_BND.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.bnd",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_DEL {
    input:
    input_vcfs = Remove_Sname_DEL.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.cn.del",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_INS {
    input:
    input_vcfs = Remove_Sname_INS.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.ins",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_OTHER {
    input:
    input_vcfs = Remove_Sname_OTHER.output_vcf_gz,
    output_vcf_basename = cohort_name + ".merged.gt.cn.other",
    preemptible_tries = preemptible_tries
  }

  output {
    Array[File] gt_del = Remove_Sname_DEL.output_vcf_gz
    Array[File] gt_other = Remove_Sname_OTHER.output_vcf_gz
    Array[File] gt_bnd = Remove_Sname_BND.output_vcf_gz
    Array[File] gt_ins = Remove_Sname_INS.output_vcf_gz
    File output_ped = Make_Pedigree_File.output_ped
    File output_vcf_bnd = Paste_VCF_BND.output_vcf_gz
    File output_vcf_del = Paste_VCF_DEL.output_vcf_gz
    File output_vcf_ins = Paste_VCF_INS.output_vcf_gz
    File output_vcf_other = Paste_VCF_OTHER.output_vcf_gz
  }
}
