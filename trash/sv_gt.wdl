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
    docker: "cc2qe/extract-sv-reads:v1"
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
    docker: "cc2qe/cnvnator:v1"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD" 
    preemptible: preemptible_tries
  }

  output {
    String sex = read_string(stdout())
  }
}

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
    docker: "ubuntu:14.04"
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
    preemptible: preemptible_tries
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
    preemptible: preemptible_tries
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
    memory: "4 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.gt.vcf"
    File output_lib = "${basename}.cram.json"
  }
}

task SV_Copy_Number {
  String basename
  File input_vcf
  File input_cn_hist_root
  File ref_cache
  Int disk_size
  Int preemptible_tries

  command {
    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    cat ${input_vcf} | create_coordinates -o coordinates.txt

    svtools copynumber \
      --cnvnator cnvnator \
      -s ${basename} \
      -w 100 \
      -r ${input_cn_hist_root} \
      -c coordinates.txt \
      -i ${input_vcf} \
      > ${basename}.cn.vcf
  }
  
  runtime {
    docker: "cc2qe/sv-pipeline:v1"
    cpu: "1"
    memory: "1 GB"
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
  >>>

  runtime {
    docker: "cc2qe/cnvnator:v1"
    cpu: threads
    memory: "26 GB"
    disks: "local-disk " + disk_size + " HDD" 
    preemptible: preemptible_tries
  }

  output {
    File output_cn_root = "cnvnator.out/${basename}.cram.root"
    File output_cn_hist_root = "cnvnator.out/${basename}.cram.hist.root"
  }
}

task Sort_VCF_Variants {
  Array[File] input_vcfs
  String output_vcf_basename
  Int disk_size
  Int preemptible_tries

  command {
    svtools lsort \
      -b 200 \
      -f ${write_lines(input_vcfs)} \
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

task Merge_VCF_Variants {
  File input_vcf
  String output_vcf_basename
  Int disk_size
  Int preemptible_tries

  command {
    svtools lmerge \
      -i ${input_vcf} \
      -f 20 \
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

task Paste_VCF {
  Array[File] input_vcfs
  String output_vcf_basename
  Int disk_size
  Int preemptible_tries

  command {
    svtools vcfpaste \
      -f ${write_lines(input_vcfs)} \
      -q \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "cc2qe/svtools:v1"
    cpu: "1"
    memory: "1 GB"
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
    docker: "cc2qe/svtools:v1"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task SV_Classify {
  File input_vcf_gz
  String output_vcf_basename
  File mei_annotation_bed
  Int disk_size
  Int preemptible_tries

  command {
    zcat ${input_vcf_gz} \
      | svtools classify \
      -g ceph.sex.txt \
      -a ${mei_annotation_bed} \
      -m large_sample \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "cc2qe/svtools:v1"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Sort_Index_VCF {
  File input_vcf_gz
  File output_vcf_name
  Int disk_size
  Int preemptible_tries

  command {
    touch ${output_vcf_name}
  }

  runtime {
    docker: "cc2qe/svtools:v1"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_name}"
  }
}

# ============================
# SV detection workflow
workflow SV_Detect {
  # data inputs
  Array[File] aligned_crams
  String aligned_cram_suffix
  String cohort_name
  String final_vcf_name

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache
  File exclude_regions
  File mei_annotation_bed

  # system inputs
  Int disk_size
  Int preemptible_tries

  scatter (aligned_cram in aligned_crams) {

    String basename = sub(sub(aligned_cram, "gs://.*/", ""), aligned_cram_suffix + "$", "")
    

    # input_cram = Extract_Reads.output_cram,
    # input_cram_index = Extract_Reads.output_cram_index,
    # input_vcf = Lumpy.output_vcf,
    
    call SV_Genotype as SV_Genotype_Unmerged {
      input:
      basename = basename,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
 
    # input_cram = Extract_Reads.output_cram,   
    call Get_Sample_Name {
      input:
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }

    # input_cn_hist_root = CNVnator_Histogram.output_cn_hist_root,

    call Get_Sex {
      input:
      ref_fasta_index = ref_fasta_index,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }
  }

  call Make_Pedigree_File {
    input:
    sample_array = Get_Sample_Name.sample,
    sex_array = Get_Sex.sex,
    output_ped_basename = cohort_name,
    disk_size = 1    
  }

  call Sort_VCF_Variants {
    input:
    input_vcfs = SV_Genotype_Unmerged.output_vcf,
    output_vcf_basename = cohort_name + ".lsort.vcf",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

  call Merge_VCF_Variants {
    input:
    input_vcf = Sort_VCF_Variants.output_vcf,
    output_vcf_basename = cohort_name + ".lmerge.vcf",
    disk_size = disk_size,
    preemptible_tries = preemptible_tries
  }

}
