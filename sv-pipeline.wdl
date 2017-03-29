# extract split/discordant reads
task Extract_Reads {
  File input_cram
  String basename
  File ref_cache
  Int disk_size
  Int preemptible_tries

  command {
    mv ${input_cram} ${basename}.cram
    
    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    # index the CRAM
    samtools index ${basename}.cram

    extract-sv-reads \
      -e \
      -r \
      -i ${basename}.cram \
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
    File output_cram = "${basename}.cram"
    File output_cram_index = "${basename}.cram.crai"
    File output_splitters_bam = "${basename}.splitters.bam"
    File output_splitters_bam_index = "${basename}.splitters.bam.bai"
    File output_discordants_bam = "${basename}.discordants.bam"
    File output_discordants_bam_index = "${basename}.discordants.bam.bai"
  }
}

# LUMPY SV discovery
task Lumpy {
  String basename
  File input_cram
  File input_cram_index
  File input_splitters_bam
  File input_splitters_bam_index
  File input_discordants_bam
  File input_discordants_bam_index
  
  File ref_fasta
  File ref_fasta_index
  File exclude_regions
  Int disk_size
  Int preemptible_tries

  command {
    lumpyexpress \
      -P \
      -T ${basename}.temp \
      -o ${basename}.vcf \
      -B ${input_cram} \
      -S ${input_splitters_bam} \
      -D ${input_discordants_bam} \
      -R ${ref_fasta} \
      -x ${exclude_regions} \
      -k \
      -v
  }
  
  runtime {
    docker: "cc2qe/lumpy:v1"
    cpu: "1"
    memory: "8 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  
  output {
    File output_vcf = "${basename}.vcf"
  }
}

task SV_Genotype {
  String basename
  File input_cram
  File input_cram_index
  File input_vcf
  File ref_fasta
  File ref_fasta_index
  Int disk_size
  Int preemptible_tries
  
  command {
    mv ${input_cram} ${basename}.cram
    mv ${input_cram_index} ${basename}.cram.crai
    rm -f ${basename}.cram.json
    svtyper \
      -i ${input_vcf} \
      -B ${basename}.cram \
      -T ${ref_fasta} \
      -l ${basename}.cram.json \
      > ${basename}.gt.vcf
  }
  
  runtime {
    docker: "cc2qe/lumpy:v1"
    cpu: "1"
    memory: "8 GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_vcf = "${basename}.gt.vcf"
    File output_lib = "${basename}.cram.json"
  }
}

task CNVnator_Histogram {
  String basename
  File input_cram
  File input_cram_index
  File ref_fasta
  File ref_fasta_index
  File ref_cache
  String ref_chrom_dir = "cnvnator_chroms"
  Int disk_size
  Int preemptible_tries
  Int threads = 4
  
  command <<<
    mv ${input_cram} ${basename}.cram
    mv ${input_cram_index} ${basename}.cram.crai

    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s
    
    # Create directory of chromosome FASTA files for CNVnator
    mkdir -p ${ref_chrom_dir}
    awk -v CHROM_DIR=${ref_chrom_dir} 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > CHROM_DIR"/"CHROM".fa" }' ${ref_fasta}

    cnvnator_wrapper.py \
      -T cnvnator.out \
      -o ${basename}.cn \
      -t ${threads} \
      -w 100 \
      -b ${basename}.cram \
      -c ${ref_chrom_dir} \
      -g GRCh38 \
      --cnvnator cnvnator

    # infer the sex of the sample
    samtools idxstats ${basename}.cram \
      | awk '$1=="chrX" { print $1":0-"$2 } END { print "exit"}'
      | cnvnator -root cnvnator.out/${basename}.cram.hist.root -genotype 100 \
      | grep -v "^Assuming male" \
      | awk -v SAMPLE=${basename} '{ printf(SAMPLE"\t%.0f\t%f\n", $4,$4); }' \
      > ${basename}.sex.txt
  >>>

  runtime {
    docker: "cc2qe/cnvnator:v1"
    cpu: threads
    memory: "26 GB"
    disks: "local-disk " + disk_size + " HDD" 
  }

  output {
    File output_cn_root = "cnvnator.out/${basename}.cram.root"
    File output_cn_hist_root = "cnvnator.out/${basename}.cram.hist.root"
    File output_cn_txt = "${basename}.cn.txt"
    File output_cn_bed = "${basename}.cn.bed"
    File output_sex = "${basename}.sex.txt"
  }
}

task Merge_Cohort_VCFs {
  Array[File] input_vcfs
  Int disk_size
  Int preemptible_tries

  command {
    svtools lsort \
      basename.vcf \
      | bgzip -c \
      > /data/sorted.vcf.gz
  }
}

# SV detection workflow
workflow SV_Detect {
  # data inputs
  Array[File] aligned_crams
  String aligned_cram_suffix

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache
  File exclude_regions

  # system inputs
  Int disk_size
  Int preemptible_tries

  scatter (aligned_cram in aligned_crams) {
    
    String basename = sub(sub(aligned_cram, "gs://.*/", ""), aligned_cram_suffix + "$", "")
    
    call Extract_Reads {
      input:
      input_cram = aligned_cram,
      basename = basename,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call Lumpy {
      input:
      basename = basename,
      input_cram = Extract_Reads.output_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      input_splitters_bam = Extract_Reads.output_splitters_bam,
      input_splitters_bam_index = Extract_Reads.output_splitters_bam_index,
      input_discordants_bam = Extract_Reads.output_discordants_bam,
      input_discordants_bam_index = Extract_Reads.output_discordants_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      exclude_regions = exclude_regions,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call SV_Genotype {
      input:
      basename = basename,
      input_cram = Extract_Reads.output_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      input_vcf = Lumpy.output_vcf,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
    
    call CNVnator_Histogram {
      input:
      basename = basename,
      input_cram = Extract_Reads.output_cram,
      input_cram_index = Extract_Reads.output_cram_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_cache = ref_cache,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }
  
  # call Merge_Cohort_VCFs {
  #   input:
  # }
}
