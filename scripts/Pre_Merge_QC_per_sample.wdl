version 1.0
import "SV_Tasks.wdl" as SV

workflow Pre_Merge_QC_Per_Sample {
  input {
    # data inputs
    File manta_vcf
    File lumpy_vcf
    File cnvnator_vcf
    String cohort
    String center

    # system inputs
    Int preemptible_tries
    String basename = sub(sub(lumpy_vcf, "^.*/", ""), ".vcf.gz" + "$", "")
  }

  call SV.Count_Lumpy {
    input:
    cohort = cohort,
    center = center,
    basename = basename,
    input_vcf = lumpy_vcf, 
    preemptible_tries = preemptible_tries
  }

 call SV.Count_Manta {
   input:
   cohort = cohort,
   center = center,
   basename = basename,
   input_vcf = manta_vcf,
   preemptible_tries = preemptible_tries
 }
 
 call SV.Count_Cnvnator {
   input:
   cohort = cohort,
   center = center,
   basename = basename,
   input_vcf = cnvnator_vcf,
   preemptible_tries = preemptible_tries
 }   

  output {
    File lumpy_counts = Count_Lumpy.output_counts
    File manta_counts = Count_Manta.output_counts
    File cnvnator_counts = Count_Cnvnator.output_counts
  }
}
