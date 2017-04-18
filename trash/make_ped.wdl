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

# This method is currently very hacky due to
# Cromwell/Google Genomics issues:
# 1) Get_Sex and Get_Sample_Name must output files, not strings
# 2) Bug in parsing Array[File] requires manual replacement
#    of "gs://" prefix
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


# ============================
# SV detection workflow
workflow SV_Detect {
  # data inputs
  Array[File] aligned_crams
  String aligned_cram_suffix
  String cohort_name

  # reference inputs
  File ref_fasta
  File ref_fasta_index
  File ref_cache
  File exclude_regions
  File mei_annotation_bed

  # system inputs
  Int disk_size
  Int preemptible_tries

  # other
  Array[File] histograms

  scatter (i in range(length(aligned_crams))) {
    File aligned_cram = aligned_crams[i]
    
    String basename = sub(sub(aligned_cram, "gs://.*/", ""), aligned_cram_suffix + "$", "")
    
    call Get_Sample_Name {
      input:
      input_cram = aligned_cram,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
    }

    call Get_Sex {
      input:
      input_cn_hist_root = histograms[i],
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

  output {
    Make_Pedigree_File.*
  }
}
