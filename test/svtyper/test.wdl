import "../../scripts/SV_Tasks.wdl" as SV

workflow Test_Genotype {
  # data inputs
  String basename
  File input_cram
  File input_cram_index
  File input_vcf

  # reference inputs
  File ref_cache

  # system inputs
  Int disk_size
  Int preemptible_tries

  call SV.Genotype as Genotype_Merged {
    input:
    basename = basename,
    input_cram = input_cram,
    input_cram_index = input_cram_index,
    input_vcf = input_vcf,
    ref_cache = ref_cache,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }
}
