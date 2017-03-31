task L_Sort_VCF_Variants {
  Array[String] input_vcfs_string_array
  Array[File] input_vcfs_file_array
  Int disk_size
  Int preemptible_tries
  String output_vcf_basename

  command {
    # strip the "gs://" prefix from the file paths
    cat ${write_lines(input_vcfs_string_array)} \
      | sed 's/^gs:\/\//\.\//g' \
      > input_vcfs_file.txt

    svtools lsort \
      -b 200 \
      -f input_vcfs_file.txt \
      > ${output_vcf_basename}.vcf
  }

  runtime {
    docker: "cc2qe/svtools:v1"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${output_vcf_basename}.vcf"
  }
}

# SV detection workflow
workflow SV_Detect {
  # system inputs
  Int disk_size
  Int preemptible_tries

  call L_Sort_VCF_Variants {
    input:
    disk_size = disk_size,
    preemptible_tries = preemptible_tries,
    output_vcf_basename = "cohort"
  }

}
