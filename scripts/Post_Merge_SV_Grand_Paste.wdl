version 1.0
import "SV_Tasks.wdl" as SV

workflow Post_Merge_SV_Grand_Paste {
  # data inputs
  input {

  	Array[File] vcf_del
        Array[File] vcf_ins
        Array[File] vcf_other
        Array[File] vcf_bnd
        Array[File] ped
        Array[String] cohort

        String final_name

  	# system inputs
  	Int preemptible_tries
  }


  # Re-genotype and call copy number for each sample on the merged SV VCF
  scatter (i in range(length(vcf_del))) {

  #  call SV.Append_Cohort as Append_Cohort_DEL {
  #    input:
  #      cohort = cohort[i],
  #      input_vcf = vcf_del[i],
  #      basename = cohort[i] + ".pasted_del",
  #      preemptible_tries = preemptible_tries
  #  }
  #
  #  call SV.Append_Cohort as Append_Cohort_INS {
  #    input:
  #      cohort = cohort[i],
  #      input_vcf = vcf_ins[i],
  #      basename = cohort[i] + ".pasted_ins",
  #      preemptible_tries = preemptible_tries
  #  }
  #
  #  call SV.Append_Cohort as Append_Cohort_BND {
  #    input:
  #      cohort = cohort[i],
  #      input_vcf = vcf_bnd[i],
  #      basename = cohort[i] + ".pasted_bnd",
  #      preemptible_tries = preemptible_tries
  #  }
  #
  #  call SV.Append_Cohort as Append_Cohort_OTHER {
  #    input:
  #      cohort = cohort[i],
  #      input_vcf = vcf_other[i],
  #      basename = cohort[i] + ".pasted_other",
  #      preemptible_tries = preemptible_tries
  #  }

   call SV.Append_Cohort_Pedigree_File as Append_Cohort_Pedigree_File {
     input:
       cohort = cohort[i],
       ped = ped[i],
       basename = cohort[i],
       preemptible_tries = preemptible_tries
   }
 }


  call SV.Cat_Files {
    input:
      input_peds = Append_Cohort_Pedigree_File.reheadered_ped,
      basename = final_name,
      preemptible_tries = preemptible_tries
  }

  #call SV.Paste_VCF as Paste_VCF_BND {
  #  input:
  #  input_vcfs = Append_Cohort_BND.reheadered_vcf,
  #  output_vcf_basename = final_name + ".pasted.gt.bnd",
  #  preemptible_tries = preemptible_tries
  #}

#  call SV.Paste_VCF as Paste_VCF_DEL {
#    input:
#    input_vcfs = Append_Cohort_DEL.reheadered_vcf,
#    output_vcf_basename = final_name + ".pasted.gt.cn.del",
#    preemptible_tries = preemptible_tries
#  }
#
#  call SV.Paste_VCF as Paste_VCF_INS {
#    input:
#    input_vcfs = Append_Cohort_INS.reheadered_vcf,
#    output_vcf_basename = final_name + ".pasted.gt.ins",
#    preemptible_tries = preemptible_tries
#  }
#
#  call SV.Paste_VCF as Paste_VCF_OTHER {
#    input:
#    input_vcfs = Append_Cohort_OTHER.reheadered_vcf,
#    output_vcf_basename = final_name + ".pasted.gt.cn.other",
#    preemptible_tries = preemptible_tries
#  }

  output {
    File output_ped = Cat_Files.output_ped
#    File output_vcf_bnd = Paste_VCF_BND.output_vcf_gz
#    File output_vcf_del = Paste_VCF_DEL.output_vcf_gz
#    File output_vcf_ins = Paste_VCF_INS.output_vcf_gz
#    File output_vcf_other = Paste_VCF_OTHER.output_vcf_gz
  }
}
