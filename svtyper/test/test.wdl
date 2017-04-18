task Genotype {
  String basename
  File input_cram
  File input_cram_index
  File input_vcf
  File ref_cache
  Int disk_size
  Int preemptible_tries
  
  command {
    mv ${input_cram} ${basename}.cram
    mv ${input_cram_index} ${basename}.cram.crai

    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    rm -f ${basename}.cram.json
    zless ${input_vcf} \
      | svtyper \
      -q \
      -B ${basename}.cram \
      -l ${basename}.cram.json \
      > ${basename}.gt.vcf
  }
  
  runtime {
    docker: "halllab/svtyper:v0.1.3-6756a0a"
    cpu: "1"
    memory: "6.5 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.gt.vcf"
    File output_lib = "${basename}.cram.json"
  }
}

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

  call Genotype as Genotype_Merged {
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
