import "../../../scripts/SV_Tasks.wdl" as SV

workflow Test_Extract_Reads {
  # data inputs
  String basename
  File input_cram

  # reference inputs
  File ref_cache

  # system inputs
  Int disk_size
  Int preemptible_tries

  call SV.Extract_Reads {
    input:
    input_cram = input_cram,
    basename = basename,
    ref_cache = ref_cache,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }
}
