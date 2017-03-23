# index
task Index_Cram {
  File input_cram
  Int disk_size
  Int preemptible_tries
  
  command {
    samtools index ${input_cram}
  }
  runtime {
    docker: "cc2qe/sv-pipeline:v1"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram_index = "${input_cram}.crai"
  }
}

workflow SV_Detect {
  call Index_Cram
}
