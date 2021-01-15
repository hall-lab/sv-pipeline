version 1.0
import "SV_Tasks.wdl" as SV

workflow Post_Merge_SV_Batched_No_Paste_No_Name {
  # data inputs
  input {
  	Array[File] aligned_crams
  	Array[File] aligned_cram_indices
  	Array[File] cn_hist_roots
  	Array[File] manta_vcfs
        Array[File] sample_names
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
    #String sample_name = sample_names[i]

    #call SV.Get_Sample_Name {
    #  input:
    #  input_cram = aligned_cram,
    #  preemptible_tries = preemptible_tries
    #}

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
      sample_name = sample_names[i],
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
      sample = sample_names[i],
      input_vcf = Genotype_Merged_DEL.output_vcf,
      input_cn_hist_root = cn_hist_root,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.Copy_Number as Copy_Number_OTHER {
      input:
      basename = basename + ".other",
      sample = sample_names[i],
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
  

  output {
    Array[File] gt_del = Remove_Sname_DEL.output_vcf_gz
    Array[File] gt_other = Remove_Sname_OTHER.output_vcf_gz
    Array[File] gt_bnd = Remove_Sname_BND.output_vcf_gz
    Array[File] gt_ins = Remove_Sname_INS.output_vcf_gz
  }
}
