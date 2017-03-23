# index
task Index_Cram {
  File input_cram
  Int disk_size
  Int preemptible_tries
  
  command {
    samtools index ${input_cram}
  }
  runtime {
    docker: "quay.io/cancercollaboratory/dockstore-tool-samtools-index:latest"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram_index = "${input_cram}.crai"
  }
}

workflow SV_Detect {
  File input_cram
  Int disk_size
  Int preemptible_tries

  call Index_Cram {
    input:
    input_cram = input_cram,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }
}
