task Print_String {
  command {
    echo hello
  }
  runtime {
    docker: "ubuntu:14.04"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + 10 + " HDD"
  }
  output {
    String out_string = read_string(stdout())
  }
}

workflow My_Workflow {
  call Print_String {
  }
}
