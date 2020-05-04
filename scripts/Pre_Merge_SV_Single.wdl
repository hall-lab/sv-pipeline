version 1.0
import "https://raw.githubusercontent.com/hall-lab/sv-pipeline/terra-compatible-hja/scripts/Pre_Merge_SV_per_sample.wdl" as per_sample
import "https://raw.githubusercontent.com/hall-lab/sv-pipeline/terra-compatible-hja/scripts/Pre_Merge_QC_per_sample.wdl" as qc
import "https://raw.githubusercontent.com/hall-lab/sv-pipeline/terra-compatible-hja/scripts/SV_Tasks.wdl" as SV

workflow Pre_Merge_SV_Single1 {
  input {
    File aligned_cram
    String aligned_cram_suffix

    # reference inputs
    File ref_fasta
    File ref_fasta_index
    File ref_cache
    File? call_regions_bed
    File? call_regions_bed_index
    File exclude_regions
    String cohort
    String center

    # system inputs
    Int preemptible_tries
  }


  call per_sample.Pre_Merge_SV_Per_Sample {
    input:
      aligned_cram = aligned_cram,
      aligned_cram_suffix = aligned_cram_suffix,
	    ref_fasta = ref_fasta,
	    ref_fasta_index = ref_fasta_index,
      call_regions_bed = call_regions_bed,
      call_regions_bed_index = call_regions_bed_index,
	    ref_cache = ref_cache,
	    exclude_regions = exclude_regions,
	    preemptible_tries = preemptible_tries
  }
  
  call qc.Pre_Merge_QC_Per_Sample {
    input:
      manta_vcf = Pre_Merge_SV_Per_Sample.manta_vcf,
      lumpy_vcf = Pre_Merge_SV_Per_Sample.smoove_vcf,
      cnvnator_vcf = Pre_Merge_SV_Per_Sample.cnvnator_output_cn_txt,
      cohort = cohort,
      center = center,
	preemptible_tries = preemptible_tries
  }


  output {
    File cram_index = Pre_Merge_SV_Per_Sample.cram_index
    File manta_vcf = Pre_Merge_SV_Per_Sample.manta_vcf
    File manta_tbi = Pre_Merge_SV_Per_Sample.manta_tbi
    File manta_original_vcf = Pre_Merge_SV_Per_Sample.manta_original_vcf
    File manta_original_tbi = Pre_Merge_SV_Per_Sample.manta_original_tbi
    File cnvnator_cn_hist_root = Pre_Merge_SV_Per_Sample.cnvnator_cn_hist_root
    File cnvnator_output_cn_txt_file = Pre_Merge_SV_Per_Sample.cnvnator_output_cn_txt
    File cnvnator_cn_bed_file = Pre_Merge_SV_Per_Sample.cnvnator_cn_bed
    File smoove_vcf = Pre_Merge_SV_Per_Sample.smoove_vcf
    File smoove_csi = Pre_Merge_SV_Per_Sample.smoove_csi
    File lumpy_count = Pre_Merge_QC_Per_Sample.lumpy_counts
    File manta_count = Pre_Merge_QC_Per_Sample.manta_counts
  }
}
