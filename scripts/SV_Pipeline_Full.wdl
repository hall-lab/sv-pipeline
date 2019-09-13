version 1.0
import "Pre_Merge_SV.wdl" as premerge
import "Merge_SV.wdl" as merge
import "Post_Merge_SV.wdl" as postmerge

workflow SV_Pipeline_Full {
  input {
    Array[File] aligned_crams
    String aligned_cram_suffix
    File ref_fasta
    File ref_fasta_index
    File ref_cache
    File? call_regions_bed
    File? call_regions_bed_index
    File exclude_regions
    String cohort
    String center
    String final_vcf_name
    Int preemptible_tries
  }

  call premerge.Pre_Merge_SV {
    input:
      aligned_crams = aligned_crams
      aligned_cram_suffix = aligned_cram_suffix
      ref_fasta = ref_fasta
      ref_fasta_index = ref_fasta_index
      ref_cache = ref_cache
      call_regions_bed = call_regions_bed
      call_regions_bed_index = call_regions_bed_index
      exclude_regions = exclude_regions
      cohort = cohort
      center = center
      preemptible_tries = preemptible_tries
  }

  call merge.Merge_SV {
    input {
      manta_input_vcfs = premerge.Pre_Merge_SV.manta_vcfs,
      smoove_inputs_vcfs = premerge.Pre_Merge_SV.smoove_vcfs,
      preemptible_tries = preemptible_tries
    }
  }

  call postmerge.Post_Merge_SV {
    input {
      aligned_crams = aligned_crams,
      aligned_cram_indices = premerge.Pre_Merge_SV.cram_indices,
      cn_hist_roots = premerge.Pre_Merge_SV.cnvnator_cn_hist_roots,
      manta_vcfs = premerge.Pre_Merge_SV.manta_original_vcfs,
      aligned_cram_suffix = aligned_cram_suffix,
      merged_vcf = merge.Merge_SV.output_vcf,
      cohort_name = cohort,
      final_vcf_name = final_vcf_name
    }
  }

  output {
    File output_ped = postmerge.Post_Merge_SV.output_ped
    File output_vcf_bnd = postmerge.Post_Merge_SV.output_vcf_bnd
    File output_vcf_index_bnd = postmerge.Post_Merge_SV.output_vcf_index_bnd
    File output_vcf_del = postmerge.Post_Merge_SV.output_vcf_del
    File output_vcf_ins = postmerge.Post_Merge_SV.output_vcf_ins
    File output_vcf_index_other = postmerge.Post_Merge_SV.output_vcf_index_other
    File output_vcf_other = postmerge.Post_Merge_SV.output_vcf_other
    File output_vcf_index_del = postmerge.Post_Merge_SV.output_vcf_index_del
    File output_vcf_index_ins = postmerge.Post_Merge_SV.output_vcf_index_ins
  }
}
