import "../../scripts/SV_Tasks.wdl" as SV

workflow Test_Simple_Workflow {
    File input_cram

    Int disk_size
    Int preemptible_tries

    call SV.Get_Sample_Name {
        input:
        input_cram = input_cram,
        disk_size = disk_size,
        preemptible_tries = preemptible_tries
    }
}
