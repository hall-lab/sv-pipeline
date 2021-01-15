version 1.0
import "SV_Tasks.wdl" as SV

workflow Index_Crams {
  # data inputs
  input {
  	Array[File] aligned_crams
        String aligned_cram_suffix

  	# reference inputs
  	File ref_fasta
  	File ref_fasta_index
  	File ref_cache

  	# system inputs
  	Int preemptible_tries
  }

  scatter (i in range(length(aligned_crams))) {
    
    File aligned_cram = aligned_crams[i]
    String basename = sub(sub(aligned_cram, "^.*/", ""), aligned_cram_suffix + "$", "")

    call SV.Index_Cram {
      input:
      basename = basename,
      input_cram = aligned_cram,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
    }
  }


  output {
    Array[File] cram_index = Index_Cram.output_cram_index
  }
}
