import "../../../scripts/SV_Tasks.wdl" as SV

workflow Test_Lumpy {
  # data inputs
  String basename
  File input_cram
  File input_cram_index
  File input_splitters_bam
  File input_splitters_bam_index
  File input_discordants_bam
  File input_discordants_bam_index

  # reference inputs
  File ref_cache
  File exclude_regions

  # system inputs
  Int disk_size
  Int preemptible_tries

  call SV.Lumpy {
    input:
    basename = basename,
    input_cram = input_cram,
    input_cram_index = input_cram_index,
    input_splitters_bam = input_splitters_bam,
    input_splitters_bam_index = input_splitters_bam_index,
    input_discordants_bam = input_discordants_bam,
    input_discordants_bam_index = input_discordants_bam_index,
    ref_cache = ref_cache,
    exclude_regions = exclude_regions,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }
}
