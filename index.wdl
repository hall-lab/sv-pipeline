# index
task Index_Cram {
  File input_cram
  String basename
  Int disk_size
  Int preemptible_tries
  
  command {
    mv ${input_cram} ${basename}.cram
    samtools index ${basename}.cram
  }
  runtime {
    docker: "cc2qe/samtools:v1"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram = "${basename}.cram"
    File output_cram_index = "${basename}.cram.crai"
  }
}

workflow SV_Detect {
  File aligned_cram
  String aligned_cram_suffix
  String basename = sub(sub(aligned_cram, "gs://.*/", ""), aligned_cram_suffix + "$", "")

  Int disk_size
  Int preemptible_tries

  # Because of a wdl/cromwell bug this is not currently valid so we have to sub(sub()) in each task
  # String basename = sub(sub(unmapped_bam, "gs://.*/", ""), unmapped_bam_suffix + "$", "")

  call Index_Cram {
    input:
    input_cram = aligned_cram,
    basename = basename,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }
}
