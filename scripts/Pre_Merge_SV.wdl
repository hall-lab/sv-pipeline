import "SV_Tasks.wdl" as SV

workflow Pre_Merge_SV {
  # data inputs
  Array[File] aligned_crams
  String aligned_cram_suffix

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache
  File exclude_regions

  # system inputs
  Int disk_size
  Int preemptible_tries

  scatter (aligned_cram in aligned_crams) {
    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")
    
    call SV.Extract_Reads {
      input:
      input_cram = aligned_cram,
      basename = basename,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call SV.Lumpy {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      input_splitters_bam = Extract_Reads.output_splitters_bam,
      input_splitters_bam_index = Extract_Reads.output_splitters_bam_index,
      input_discordants_bam = Extract_Reads.output_discordants_bam,
      input_discordants_bam_index = Extract_Reads.output_discordants_bam_index,
      ref_cache = ref_cache,
      exclude_regions = exclude_regions,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call SV.Genotype as Genotype_Unmerged {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      input_vcf = Lumpy.output_vcf,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call SV.CNVnator_Histogram {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  output {
    Extract_Reads.output_cram_index
    CNVnator_Histogram.output_cn_hist_root
    CNVnator_Histogram.output_cn_txt
    CNVnator_Histogram.output_cn_bed
    Genotype_Unmerged.output_vcf
  }
}
