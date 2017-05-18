import "SV_Tasks.wdl" as SV

workflow Lumpy_and_Genotype {
  # data inputs
  Array[File] aligned_crams
  Array[File] aligned_cram_indices
  Array[File] splitters_bams
  Array[File] splitters_bam_indices
  Array[File] discordants_bams
  Array[File] discordants_bam_indices
  String aligned_cram_suffix

  # reference inputs
  File ref_cache
  File exclude_regions

  # system inputs
  Int disk_size
  Int preemptible_tries

  scatter (i in range(length(aligned_crams))) {
    File aligned_cram = aligned_crams[i]
    File aligned_cram_index = aligned_cram_indices[i]
    File splitters_bam = splitters_bams[i]
    File splitters_bam_index = splitters_bam_indices[i]
    File discordants_bam = discordants_bams[i]
    File discordants_bam_index = discordants_bam_indices[i]

    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")
    
    call SV.Lumpy {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_splitters_bam = splitters_bam,
      input_splitters_bam_index = splitters_bam_index,
      input_discordants_bam = discordants_bam,
      input_discordants_bam_index = discordants_bam_index,
      ref_cache = ref_cache,
      exclude_regions = exclude_regions,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call SV.Genotype as Genotype_Unmerged {
      input:
      basename = basename,
      input_cram = aligned_cram,
      input_cram_index = aligned_cram_index,
      input_vcf = Lumpy.output_vcf,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  output {
    Genotype_Unmerged.output_vcf
  }
}
