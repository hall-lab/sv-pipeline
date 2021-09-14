version 1.0
import "SV_Tasks.wdl" as SV

workflow Post_Merge_SV_freeze3 {
  # data inputs
  input {
  	Array[File] aligned_crams
  	Array[File] aligned_cram_indices
	Array[String] sample_names
  	Array[File] cn_hist_roots
  	Array[File] manta_vcfs
  	String aligned_cram_suffix
  	File merged_vcf_ins
	File merged_vcf_1
	File merged_vcf_2
	File merged_vcf_3
	File merged_vcf_4
	File merged_vcf_5
	File merged_vcf_6
	File ins_split
  	String batch
	String cohort
	String center


  	# reference inputs
  	File ref_fasta
  	File ref_fasta_index
  	File ref_cache

  	# system inputs
  	Int preemptible_tries
  }


  # Re-genotype and call copy number for each sample on the merged SV VCF
  scatter (i in range(length(aligned_crams))) {
    
    File aligned_cram = aligned_crams[i]
    File aligned_cram_index = aligned_cram_indices[i]
    String sample_name = sample_names[i]
    File cn_hist_root = cn_hist_roots[i]
    File manta_vcf = manta_vcfs[i]
    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")


    call SV.Get_Sex {
      input:
      input_cn_hist_root = cn_hist_root,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged1 {
      input:
      basename = basename + ".chunk1",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf_1,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged2 {
      input:
      basename = basename + ".chunk2",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf_2,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged3 {
      input:
      basename = basename + ".chunk3",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf_3,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged4 {
      input:
      basename = basename + ".chunk4",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf_4,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged5 {
      input:
      basename = basename + ".chunk5",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf_5,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Genotype as Genotype_Merged6 {
      input:
      basename = basename + ".chunk6",
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = merged_vcf_6,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Take_Original_Genotypes as Genotype_Merged_INS {
      input:
      sample_name = sample_name,
      original_per_sample_vcf = manta_vcf,
      basename = basename + ".chunk7",
      input_vcf =  merged_vcf_ins,
      input_variant_to_sname_mapping = ins_split,
      preemptible_tries = preemptible_tries
    }


    call SV.Copy_Number as Copy_Number3 {
      input:
      basename = basename + ".chunk3",
      sample = sample_name,
      input_vcf = Genotype_Merged3.output_vcf,
      input_cn_hist_root = cn_hist_root,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Copy_Number as Copy_Number4 {
      input:
      basename = basename + ".chunk4",
      sample = sample_name,
      input_vcf = Genotype_Merged4.output_vcf,
      input_cn_hist_root = cn_hist_root,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }
    call SV.Copy_Number as Copy_Number5 {
      input:
      basename = basename + ".chunk5",
      sample = sample_name,
      input_vcf = Genotype_Merged5.output_vcf,
      input_cn_hist_root = cn_hist_root,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

  }
  
  call SV.Make_Pedigree_File {
    input:
    sample_array = sample_names,
    sex_array = Get_Sex.sex,
    output_ped_basename = cohort + "." + batch,
  }

  call SV.Paste_VCF as Paste_VCF1 {
    input:
    input_vcfs = Genotype_Merged1.output_vcf,
    output_vcf_basename =  cohort + "." + batch + ".merged.chunk1",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF2 {
    input:
    input_vcfs = Genotype_Merged2.output_vcf,
    output_vcf_basename =  cohort + "." + batch + ".merged.chunk2",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF6 {
    input:
    input_vcfs = Genotype_Merged6.output_vcf,
    output_vcf_basename =  cohort + "." + batch + ".merged.chunk6",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF3 {
    input:
    input_vcfs = Copy_Number3.output_vcf,
    output_vcf_basename =  cohort + "." + batch + ".merged.chunk3",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF4 {
    input:
    input_vcfs = Copy_Number4.output_vcf,
    output_vcf_basename =  cohort + "." + batch + ".merged.chunk4",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF5 {
    input:
    input_vcfs = Copy_Number5.output_vcf,
    output_vcf_basename =  cohort + "." + batch + ".merged.chunk5",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF7 {
    input:
    input_vcfs = Genotype_Merged_INS.output_vcf,
    output_vcf_basename = cohort + "." + batch + ".merged.chunk7",
    preemptible_tries = preemptible_tries
  }

  output {
    File output_ped = Make_Pedigree_File.output_ped
    File pasted_vcf1 = Paste_VCF1.output_vcf_gz
    File pasted_vcf2 = Paste_VCF2.output_vcf_gz
    File pasted_vcf3 = Paste_VCF3.output_vcf_gz
    File pasted_vcf4 = Paste_VCF4.output_vcf_gz
    File pasted_vcf5 = Paste_VCF5.output_vcf_gz
    File pasted_vcf6 = Paste_VCF6.output_vcf_gz
    File pasted_vcf7 = Paste_VCF7.output_vcf_gz
    Array[File] cn3 = Copy_Number3.output_vcf
    Array[File] cn4 = Copy_Number4.output_vcf
    Array[File] cn5 = Copy_Number5.output_vcf
    Array[File] gt1 = Genotype_Merged1.output_vcf
    Array[File] gt2 = Genotype_Merged2.output_vcf
    Array[File] gt3 = Genotype_Merged3.output_vcf
    Array[File] gt4 = Genotype_Merged4.output_vcf
    Array[File] gt5 = Genotype_Merged5.output_vcf
    Array[File] gt6 = Genotype_Merged6.output_vcf
    Array[File] gt7 = Genotype_Merged_INS.output_vcf
    Array[String] sex = Get_Sex.sex
  }
}
