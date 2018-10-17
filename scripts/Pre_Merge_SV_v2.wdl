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
  Int preemptible_tries

  scatter (aligned_cram in aligned_crams) {
    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")

    call SV.Index_Cram {
      input:
        basename = basename,
        input_cram = aligned_cram,
        ref_cache = ref_cache,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }
    
    call SV.Smoove {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = Index_Cram.output_cram_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      exclude_regions = exclude_regions,
      preemptible_tries = preemptible_tries
    }

    call SV.Manta {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = Index_Cram.output_cram_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }

    call SV.CNVnator_Histogram {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = Index_Cram.output_cram_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }
  }

  output {
    Smoove.output_vcf
    Smoove.output_csi
    Manta.output_vcf
    Manta.output_tbi
    CNVnator_Histogram.output_cn_hist_root
    CNVnator_Histogram.output_cn_txt
    CNVnator_Histogram.output_cn_bed
  }
}
