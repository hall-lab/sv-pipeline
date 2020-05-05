import "Pre_Merge_QC_per_sample.wdl" as per_sample

task cat {
  Array[File] lumpy_counts
  Array[File] manta_counts
  Array[File] cnvnator_counts
  String cohort
  String center

  command <<<
   cat ${sep=" "lumpy_counts} | gzip -c  > ${center}.${cohort}.lumpy_counts.txt.gz
   cat ${sep=" "manta_counts} | gzip -c > ${center}.${cohort}.manta_counts.txt.gz
  >>>

  runtime {
    docker: "ubuntu@sha256:edf05697d8ea17028a69726b4b450ad48da8b29884cd640fec950c904bfb50ce"
    cpu: "1"
    memory: "4 GB"
    disks: "local-disk 4 HDD"
  }

  output {
    File lumpy_out = "${center}.${cohort}.lumpy_counts.txt.gz"
    File manta_out = "${center}.${cohort}.manta_counts.txt.gz"
    #File cnvator_out = "${center}.${cohort}.cnvnator_counts.txt.gz"

  }
}


workflow Pre_Merge_QC {
  Array[File] lumpy_vcfs
  Array[File] manta_vcfs
  Array[File] cnvnator_vcfs
  String cohort
  String center


  # system inputs
  Int preemptible_tries

  scatter (i in range(length(lumpy_vcfs))) {
    File lumpy_vcf = lumpy_vcfs[i]
    File manta_vcf = manta_vcfs[i]
    File cnvnator_vcf = cnvnator_vcfs[i]

    call per_sample.Pre_Merge_QC_Per_Sample {
      input:
        lumpy_vcf = lumpy_vcf,
        manta_vcf = manta_vcf,
        cnvnator_vcf = cnvnator_vcf,
        cohort = cohort ,
        center = center ,
        preemptible_tries = preemptible_tries
     }
  }
  
  call cat { input:
     lumpy_counts = Pre_Merge_QC_Per_Sample.lumpy_counts,
     manta_counts = Pre_Merge_QC_Per_Sample.manta_counts,
     cnvnator_counts = Pre_Merge_QC_Per_Sample.cnvnator_counts,
     cohort = cohort,
     center = center
  }
 
    
  output {
    Array[File] lumpy_counts = Pre_Merge_QC_Per_Sample.lumpy_counts
    Array[File] manta_counts = Pre_Merge_QC_Per_Sample.manta_counts
    Array[File] cnvnator_counts = Pre_Merge_QC_Per_Sample.cnvnator_counts
  }
}
