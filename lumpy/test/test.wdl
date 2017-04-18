task Lumpy {
  String basename
  File input_cram
  File input_cram_index
  File input_splitters_bam
  File input_splitters_bam_index
  File input_discordants_bam
  File input_discordants_bam_index

  File ref_cache
  File exclude_regions
  Int disk_size
  Int preemptible_tries

  command {
    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    lumpyexpress \
      -P \
      -T ${basename}.temp \
      -o ${basename}.vcf \
      -B ${input_cram} \
      -S ${input_splitters_bam} \
      -D ${input_discordants_bam} \
      -x ${exclude_regions} \
      -k \
      -v
  }

  runtime {
    docker: "halllab/lumpy:v0.2.13-2d611fa"
    cpu: "1"
    memory: "8 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.vcf"
  }
}

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

  call Lumpy {
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
