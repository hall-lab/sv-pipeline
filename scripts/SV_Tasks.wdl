# get the sample (SM) field from a CRAM file
task Get_Sample_Name {
  File input_cram
  Int disk_size
  Int preemptible_tries

  command {
    samtools view -H ${input_cram} \
      | grep -m 1 '^@RG' | tr '\t' '\n' \
      | grep '^SM:' | sed 's/^SM://g'
  }

  runtime {
    docker: "halllab/extract-sv-reads@sha256:192090f72afaeaaafa104d50890b2fc23935c8dc98988a9b5c80ddf4ec50f70c"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    String sample = read_string(stdout())
  }
}

# infer the sex of a sample based on chrom X copy number
task Get_Sex {
  File input_cn_hist_root
  File ref_fasta_index
  Int disk_size
  Int preemptible_tries
  
  command <<<
    cat ${ref_fasta_index} \
      | awk '$1=="chrX" { print $1":0-"$2 } END { print "exit"}' \
      | cnvnator -root ${input_cn_hist_root} -genotype 100 \
      | grep -v "^Assuming male" \
      | awk '{ printf("%.0f\n",$4); }'
  >>>

  runtime {
    docker: "halllab/cnvnator@sha256:c41e9ce51183fc388ef39484cbb218f7ec2351876e5eda18b709d82b7e8af3a2"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD" 
    preemptible: preemptible_tries
  }

  output {
    String sex = read_string(stdout())
  }
}

# Create pedigree file from samples, with sex inferred from
# CNVnator X chrom copy number
task Make_Pedigree_File {
  Array[String] sample_array
  Array[String] sex_array
  String output_ped_basename
  Int disk_size

  command <<<
    paste ${write_lines(sample_array)} ${write_lines(sex_array)} \
      | awk '{ print $1,$1,-9,-9,$2,-9 }' OFS='\t' \
      > ${output_ped_basename}.ped
  >>>

  runtime {
    docker: "ubuntu@sha256:edf05697d8ea17028a69726b4b450ad48da8b29884cd640fec950c904bfb50ce"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_ped = "${output_ped_basename}.ped"
  }
}

# extract split/discordant reads
task Extract_Reads {
  File input_cram
  String basename
  File ref_cache
  Int disk_size
  Int preemptible_tries

  command {
    ln -s ${input_cram} ${basename}.cram
    
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
    docker: "halllab/extract-sv-reads@sha256:192090f72afaeaaafa104d50890b2fc23935c8dc98988a9b5c80ddf4ec50f70c"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
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

  File ref_cache
  File exclude_regions
  Int disk_size
  Int preemptible_tries

  command {
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

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
    docker: "halllab/lumpy@sha256:59ce7551307a54087e57d5cec89b17511d910d1fe9fa3651c12357f0594dcb07"
    cpu: "1"
    memory: "8 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.vcf"
  }
}

# Smoove wrapper
task Smoove {
  String basename
  File input_cram
  File input_cram_index

  File ref_fasta
  # Can't tell yet if we need the ref cache
  #File ref_cache
  File exclude_regions

  Int disk_size
  Int preemptible_tries
  
  command {
    set -e
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

    # build the reference sequence cache
    # tar -zxf ${ref_cache}
    # export REF_PATH=./cache/%2s/%2s/%s
    # export REF_CACHE=./cache/%2s/%2s/%s
    
    smoove call \
      --name ${basename} \
      --exclude ${exclude_regions} \
      --fasta ${ref_fasta} \
      --genotype \
      ${basename}.cram
  }

  runtime {
    docker: "halllab/smoove@sha256:ee1a3d994b0d6a1eeed91894049ce4d883da40912c35a0aefeb6e26ad0ea34d1",
    cpu: "1"
    memory: "4 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}-smoove.genotyped.vcf.gz"
    File output_csi = "${basename}-smoove.genotyped.vcf.gz.csi"
  }
}

task Genotype {
  String basename
  File input_cram
  File input_cram_index
  File input_vcf
  File ref_cache
  Int disk_size
  Int preemptible_tries
  
  command {
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    rm -f ${basename}.cram.json
    zless ${input_vcf} \
      | svtyper \
      -B ${basename}.cram \
      -l ${basename}.cram.json \
      > ${basename}.gt.vcf
  }
  
  runtime {
    docker: "halllab/svtyper@sha256:21d757e77dfc52fddeab94acd66b09a561771a7803f9581b8cca3467ab7ff94a"
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

task Copy_Number {
  String basename
  String sample
  File input_vcf
  File input_cn_hist_root
  File ref_cache
  Int disk_size
  Int preemptible_tries

  command {
    create_coordinates \
      -i ${input_vcf} \
      -o coordinates.txt

    svtools copynumber \
      -i ${input_vcf} \
      -s ${sample} \
      --cnvnator cnvnator \
      -w 100 \
      -r ${input_cn_hist_root} \
      -c coordinates.txt \
      > ${basename}.cn.vcf
  }
  
  runtime {
    docker: "halllab/cnvnator@sha256:c41e9ce51183fc388ef39484cbb218f7ec2351876e5eda18b709d82b7e8af3a2"
    cpu: "1"
    memory: "4 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.cn.vcf"
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
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

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
  >>>

  runtime {
    docker: "halllab/cnvnator@sha256:c41e9ce51183fc388ef39484cbb218f7ec2351876e5eda18b709d82b7e8af3a2"
    cpu: threads
    memory: "26 GB"
    disks: "local-disk " + disk_size + " HDD" 
    preemptible: preemptible_tries
  }

  output {
    File output_cn_hist_root = "cnvnator.out/${basename}.cram.hist.root"
    File output_cn_txt = "${basename}.cn.txt"
    File output_cn_bed = "${basename}.cn.bed"
  }
}

task L_Sort_VCF_Variants {
  Array[File] input_vcfs
  File input_vcfs_file = write_lines(input_vcfs)
  String output_vcf_basename
  Int disk_size
  Int preemptible_tries

  command {
    # strip the "gs://" prefix from the file paths
    cat ${input_vcfs_file} \
      | sed 's/^gs:\/\//\.\//g' \
      > input_vcfs_file.local_map.txt

    svtools lsort \
      -b 200 \
      -f input_vcfs_file.local_map.txt \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:f2f3f9c788beb613bc26c858f897694cd6eaab450880c370bf0ef81d85bf8d45"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task L_Merge_VCF_Variants {
  File input_vcf_gz
  String output_vcf_basename
  Int disk_size
  Int preemptible_tries

  command {
    zcat ${input_vcf_gz} \
      | svtools lmerge \
      -i /dev/stdin \
      -f 20 \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }
  
  runtime {
    docker: "halllab/svtools@sha256:f2f3f9c788beb613bc26c858f897694cd6eaab450880c370bf0ef81d85bf8d45"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Paste_VCF {
  Array[File] input_vcfs
  File input_vcfs_file = write_lines(input_vcfs)
  String output_vcf_basename
  Int disk_size
  Int preemptible_tries

  command {
    # strip the "gs://" prefix from the file paths
    cat ${input_vcfs_file} \
      | sed 's/^gs:\/\//\.\//g' \
      > input_vcfs_file.local_map.txt
    
    svtools vcfpaste \
      -f input_vcfs_file.local_map.txt \
      -q \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:f2f3f9c788beb613bc26c858f897694cd6eaab450880c370bf0ef81d85bf8d45"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Prune_VCF {
  File input_vcf_gz
  String output_vcf_basename
  Int disk_size
  Int preemptible_tries

  command {
    zcat ${input_vcf_gz} \
      | svtools afreq \
      | svtools vcftobedpe \
      | svtools bedpesort \
      | svtools prune -s -d 100 -e 'AF' \
      | svtools bedpetovcf \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:f2f3f9c788beb613bc26c858f897694cd6eaab450880c370bf0ef81d85bf8d45"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Classify {
  File input_vcf_gz
  File input_ped
  String output_vcf_basename
  File mei_annotation_bed
  Int disk_size
  Int preemptible_tries

  command {
    cat ${input_ped} \
      | cut -f 2,5 \
      > sex.txt

    zcat ${input_vcf_gz} \
      | svtools classify \
      -g sex.txt \
      -a ${mei_annotation_bed} \
      -m large_sample \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:f2f3f9c788beb613bc26c858f897694cd6eaab450880c370bf0ef81d85bf8d45"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Sort_Index_VCF {
  File input_vcf_gz
  String output_vcf_name
  Int disk_size
  Int preemptible_tries

  command {
    zcat ${input_vcf_gz} \
      | svtools vcfsort \
      | bgzip -c \
      > ${output_vcf_name}

    tabix -p vcf -f ${output_vcf_name}
  }

  runtime {
    docker: "halllab/svtools@sha256:f2f3f9c788beb613bc26c858f897694cd6eaab450880c370bf0ef81d85bf8d45"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_name}"
    File output_vcf_gz_index = "${output_vcf_name}.tbi"
  }
}

