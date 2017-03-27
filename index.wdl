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

# extract split/discordant reads
task Extract_Reads {
  File input_cram
  String basename
  File ref_fasta
  Int disk_size
  Int preemptible_tries

  command {
    extract-sv-reads \
      -e \
      -r \
      -T ${ref_fasta} \
      -i ${input_cram} \
      -s ${basename}.splitters.bam \
      -d ${basename}.discordants.bam
    samtools index ${basename}.splitters.bam
    samtools index ${basename}.discordants.bam
  }

  runtime {
    docker: "cc2qe/extract-sv-reads:v1"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_splitters_bam = "${basename}.splitters.bam"
    File output_splitters_bam_index = "${basename}.splitters.bam.bai"
    File output_discordants_bam = "${basename}.discordants.bam"
    File output_discordants_bam_index = "${basename}.discordants.bam.bai"
  }
}

# run LUMPY
# change this to lumpyexpress + svtyper --cc
task Lumpy {
  String basename
  File input_cram
  File input_cram_index
  File input_splitters_bam
  File input_spittlers_bam_index
  File input_discordants_bam
  File input_discordants_bam_index
  
  # ref_fasta unnecessary with lumpyexpress
  File ref_fasta
  File exclude_regions
  Int disk_size
  Int preemptible_tries

  command {
    speedseq sv \
      -P \
      -g \
      -T ${basename}.temp \
      -o ${basename} \
      -B ${input_cram} \
      -S ${input_splitters_bam} \
      -D ${input_discordants_bam} \
      -R ${ref_fasta} \
      -x ${exclude_regions} \
      -v \
      -k
  }
  
  runtime {
    docker: "cc2qe/sv-pipeline:v1"
    cpu: "1"
    memory: "8 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  
  # after convert to lumpy express, output the library JSON file as well
  output {
    File output_sv_vcf = "${basename}.sv.vcf.gz"
    File output_sv_vcf_index = "${basename}.sv.vcf.gz.tbi"
  }
}

# SV detection workflow
workflow SV_Detect {
  File aligned_cram
  String aligned_cram_suffix
  String basename = sub(sub(aligned_cram, "gs://.*/", ""), aligned_cram_suffix + "$", "")

  Int disk_size
  Int preemptible_tries

  File ref_fasta
  File ref_fasta_index

  # Because of a wdl/cromwell bug this is not currently valid so we have to sub(sub()) in each task
  # String basename = sub(sub(unmapped_bam, "gs://.*/", ""), unmapped_bam_suffix + "$", "")

  # call Index_Cram {
  #   input:
  #   input_cram = aligned_cram,
  #   basename = basename,
  #   disk_size = disk_size,
  #   preemptible_tries = preemptible_tries
  # }
  
  # call Extract_Reads {
  #   input:
  #   input_cram = Index_Cram.output_cram,
  #   basename = basename,
  #   disk_size = disk_size,
  #   preemptible_tries = preemptible_tries,
  #   ref_fasta = ref_fasta
  # }

  # input_cram = Index_Cram.output_cram,
  # input_cram_index = Index_Cram.output_cram_index,
  # input_splitters_bam = Extract_Reads.output_splitters_bam,
  # input_splitters_bam_index = Extract_Reads.output_splitters_bam_index,
  # input_discordants_bam = Extract_Reads.output_discordants_bam,
  # input_discordants_bam_index = Extract_Reads.output_discordants_bam_index,

  call Lumpy {
    input:
    basename = basename,
    input_cram = input_cram,
    input_cram_index = input_cram_index,
    input_splitters_bam = input_splitters_bam,
    input_splitters_bam_index = input_splitters_bam_index,
    input_discordants_bam = input_discordants_bam,
    input_discordants_bam_index = input_discordants_bam_index,
    ref_fasta = ref_fasta,
    exclude_regions = exclude_regions,
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }
}
