import "Pre_Merge_SV_per_sample.wdl" as per_sample
import "Pre_Merge_QC_per_sample.wdl" as qc

workflow Pre_Merge_SV_v2 {
  File cram_list
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
  Array[File] aligned_crams = read_lines(cram_list)

  scatter (i in range(length(aligned_crams))) {
    File aligned_cram = aligned_crams[i]

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
  }

  output {
    Array[File] cram_indices = Pre_Merge_SV_Per_Sample.cram_index
    Array[File] manta_vcfs = Pre_Merge_SV_Per_Sample.manta_vcf
    Array[File] manta_tbis = Pre_Merge_SV_Per_Sample.manta_tbi
    Array[File] cnvnator_cn_hist_roots = Pre_Merge_SV_Per_Sample.cnvnator_cn_hist_root
    Array[File] cnvnator_output_cn_txt_files = Pre_Merge_SV_Per_Sample.cnvnator_output_cn_txt
    Array[File] cnvnator_cn_bed_files = Pre_Merge_SV_Per_Sample.cnvnator_cn_bed
    Array[File] smoove_vcfs = Pre_Merge_SV_Per_Sample.smoove_vcf
    Array[File] smoove_csis = Pre_Merge_SV_Per_Sample.smoove_csi
    Array[File] lumpy_counts = Pre_Merge_QC_Per_Sample.lumpy_counts
    Array[File] manta_counts = Pre_Merge_QC_Per_Sample.cnvnator_counts
    Array[File] cnvnator_counts = Pre_Merge_QC_Per_Sample.cnvnator_counts
  }
}
