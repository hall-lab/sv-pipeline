version 1.0
import "SV_Tasks.wdl" as SV

workflow Post_Merge_SV_Paste {
  # data inputs
  input {
  	Array[File] aligned_crams
  	Array[File] aligned_cram_indices
  	Array[File] cn_hist_roots
  	Array[File] vcf_del
        Array[File] vcf_ins
        Array[File] vcf_other
        Array[File] vcf_bnd
        Array[String] sample_name

        File merged_vcf_del
        File merged_vcf_bnd
        File merged_vcf_ins
        File merged_vcf_other

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

    call SV.Get_Sex {
      input:
      input_cn_hist_root = cn_hist_roots[i],
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = preemptible_tries
    }
  }
  
  call SV.Make_Pedigree_File {
    input:
    sample_array = sample_name,
    sex_array = Get_Sex.sex,
    output_ped_basename = cohort_name,
  }

  call SV.Paste_VCF as Paste_VCF_BND {
    input:
    input_vcfs = vcf_bnd,
    output_vcf_basename = cohort_name + ".merged.gt.bnd",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_DEL {
    input:
    input_vcfs = vcf_del,
    output_vcf_basename = cohort_name + ".merged.gt.cn.del",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_INS {
    input:
    input_vcfs = vcf_ins,
    output_vcf_basename = cohort_name + ".merged.gt.ins",
    preemptible_tries = preemptible_tries
  }

  call SV.Paste_VCF as Paste_VCF_OTHER {
    input:
    input_vcfs = vcf_other,
    output_vcf_basename = cohort_name + ".merged.gt.cn.other",
    preemptible_tries = preemptible_tries
  }

  output {
    File output_ped = Make_Pedigree_File.output_ped
    File output_vcf_bnd = Paste_VCF_BND.output_vcf_gz
    File output_vcf_del = Paste_VCF_DEL.output_vcf_gz
    File output_vcf_ins = Paste_VCF_INS.output_vcf_gz
    File output_vcf_other = Paste_VCF_OTHER.output_vcf_gz
  }
}
