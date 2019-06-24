# get the sample (SM) field from a CRAM file
task Get_Sample_Name {
  File input_cram
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
    disks: "local-disk " + ceil( size(input_cram, "GB") + 2.0) + " HDD"
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
  Int preemptible_tries

  command <<<
    cat ${ref_fasta_index} \
      | awk '$1=="chrX" { print $1":0-"$2 } END { print "exit"}' \
      | cnvnator -root ${input_cn_hist_root} -genotype 100 \
      | grep -v "^Assuming male" \
      | awk '{ printf("%.0f\n",$4); }'
  >>>

  runtime {
    docker: "halllab/cnvnator@sha256:8bf4fa64a288c5647a9a6b1ea90d14e76f48a3e16c5bf98c63419bb7d81c8938"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk 4 HDD"
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

  command <<<
    paste ${write_lines(sample_array)} ${write_lines(sex_array)} \
      | awk '{ print $1,$1,-9,-9,$2,-9 }' OFS='\t' \
      > ${output_ped_basename}.ped
  >>>

  runtime {
    docker: "ubuntu@sha256:edf05697d8ea17028a69726b4b450ad48da8b29884cd640fec950c904bfb50ce"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk 4 HDD"
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

# index a CRAM
task Index_Cram {
  File input_cram
  String basename
  File ref_cache
  Int preemptible_tries

  command {
    ln -s ${input_cram} ${basename}.cram

    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    # index the CRAM
    samtools index ${basename}.cram
  }

  runtime {
    docker: "halllab/samtools@sha256:5e6b0430a7ad25f68e5c46a9fa9c0ebba0f9af8ebf5aebe94242954d812a4e68"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + ceil( size(input_cram, "GB") + size(ref_cache, "GB") * 5 + 1.0) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_cram_index = "${basename}.cram.crai"
  }
}

task Count_Lumpy_VCF {
  File input_vcf
  String basename
  String cohort

  Int preemptible_tries

  command <<<
    set -eo pipefail

    echo -e "Cohort\tSample\tGenotype\tType\tLength"
    bcftools query -e 'INFO/SECONDARY=1' -f "[${cohort}\t%SAMPLE\t%GT\t%INFO/SVTYPE\t%INFO/SVLEN\n]" | bgzip -c > ${basename}.counts.txt.gz
  >>>

  runtime {
    docker: "halllab/bcftools@sha256:955cbf93e35e5ee6fdb60e34bb404b7433f816e03a202dfed9ceda542e0d8906"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + ceil( size(input_vcf, "GB") * 2) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_count = "${basename}.counts.txt.gz"
  }
}

task Count_Lumpy {
  String basename
  File input_vcf
  Int preemptible_tries
  String cohort
  String center

  command <<<
    set -eo pipefail
   
     bcftools query  -f "[%CHROM\t${cohort}\t${center}\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SR\t%SAMPLE\t%GT\n]"  ${input_vcf} \
     | awk 'BEGIN{OFS="\t"}{if($1~/chr[1-9]+/ && $1!~/_/) {
       svlen=$6;
       if($6<0) svlen=-1*$6;
       len_bin=">=1kb"
       if(svlen<1000) len_bin="<1kb";
       if($7>0) $7="SR>=1";
       else $7="SR=0";
       print $1, $2, $3, $4, $5, len_bin, $7, $8, $9;}}' \
     | sort -k1,9 \
     | uniq -c \
     | awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $5, $6, $7, $8, $9, $1}'  > ${basename}.lumpy.counts.1.txt
  >>>

  runtime {
    docker: "halllab/bcftools@sha256:955cbf93e35e5ee6fdb60e34bb404b7433f816e03a202dfed9ceda542e0d8906"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + ceil( size(input_vcf, "GB") * 2) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_counts = "${basename}.lumpy.counts.1.txt"
  }
}

task Count_Manta {
  String basename
  File input_vcf
  Int preemptible_tries
  String cohort
  String center

  command <<<
    set -eo pipefail

     bcftools query  -f "[%CHROM\t${cohort}\t${center}\t%FILTER\t%INFO/SVTYPE\t%SAMPLE\t%GT\n]"  ${input_vcf} \
     | awk 'BEGIN{OFS="\t"}{if($1~/chr[1-9]+/ && $1!~/_/ && $4=="PASS") print $0;}' \
     | sort -k1,7 \
     | uniq -c \
     | awk 'BEGIN{OGS="\t"}{print $2, $3, $4, $5, $6, $7, $8,  $1}'  > ${basename}.manta.counts.1.txt
  >>>

  runtime {
    docker: "halllab/bcftools@sha256:955cbf93e35e5ee6fdb60e34bb404b7433f816e03a202dfed9ceda542e0d8906"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + ceil( size(input_vcf, "GB") * 2) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_counts = "${basename}.manta.counts.1.txt"
  }
}

task Count_Cnvnator {
  String basename
  File input_vcf
  Int preemptible_tries
  String cohort
  String center

  command <<<
    set -eo pipefail

     cat ${input_vcf} \
     | awk -v bn="${basename}" -v cohort="${cohort}" -v ctr="${center}" 'BEGIN{OFS="\t"}{
       if($2~/chr[1-9]+/ && $2!~/_/) { 
         len_bin=">2kb";
         if($3<200) len_bin="<200bp";
         else if ($3<500) len_bin="<500bp";
         else if ($3<1000) len_bin="<1kb";
         else if ($3<2000) len_bin="<2kb";
         print cohort, ctr,  bn, $1, int($4*10)/10, len_bin, $3;}
     }' \
     | sort -k1,6 \
     | bedtools groupby -g 1,2,3,4,5,6 -c 7,7 -o count,sum \
     > ${basename}.cnvnator.counts.1.txt
  >>>

  runtime {
    docker: "ernfrid/bedtools@sha256:be4c7bcf70ff0fec467245729251652d92842e4bd255479ae817a1ea2d696168"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + ceil( size(input_vcf, "GB") * 2) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_counts = "${basename}.cnvnator.counts.1.txt"
  }
}

# flagstat a CRAM
task Flagstat {
  File input_cram
  File input_cram_index
  String basename
  Int preemptible_tries

  command {
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

    # index the CRAM
    samtools flagstat ${basename}.cram > ${basename}.flagstat
  }

  runtime {
    docker: "halllab/samtools@sha256:5e6b0430a7ad25f68e5c46a9fa9c0ebba0f9af8ebf5aebe94242954d812a4e68"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + ceil( size(input_cram, "GB") + size(input_cram_index, "GB") + 1.0) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File flagstat = "${basename}.flagstat"
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

task Manta {    
  File input_cram
  File input_cram_index
  File ref_fasta
  File ref_fasta_index
  File ref_cache
  File? call_regions_bed
  File? call_regions_bed_index
  String basename
  Int preemptible_tries

  # Manta requires 2GB per thread for scheduling, but in typical cases uses less than this
  # see https://github.com/Illumina/manta/issues/38
  # Setting below derives CPU count from machine
  # Sets RAM to unlimited to jobs are scheduled only 
  # with respect to cores
  # If a task starts to fail then we can adjust the machine resources to get it
  # to succeed without adjusting the command
  # Note that we are converting to BAM on the fly as CRAM is showing extreme memory usage in some situations. See https://github.com/Illumina/manta/issues/154.
  # Note also that we are specifying an inflation factor of 4, but padding with 20GB of data. This is aimed to get us over 100GB of SSD for better performance on small samples.

  command {
    set -e
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    ${"touch " + call_regions_bed_index} 

    samtools view -hb -@8 ${basename}.cram -o ${basename}.bam 
    samtools index -@8 ${basename}.bam 

    configManta.py \
    --referenceFasta=${ref_fasta} \
    --runDir=MantaWorkflow \
    --bam=${basename}.bam ${"--callRegions=" + call_regions_bed}
    MantaWorkflow/runWorkflow.py -m local -g "unlimited"
    mv MantaWorkflow/results/variants/diploidSV.vcf.gz ${basename}.vcf.gz
    mv MantaWorkflow/results/variants/diploidSV.vcf.gz.tbi ${basename}.vcf.gz.tbi
    zcat ${basename}.vcf.gz | /opt/hall-lab/python-2.7.15/bin/python /opt/hall-lab/doctor_manta.1.py -m 700 | /opt/hall-lab/htslib-1.9/bin/bgzip -c > ${basename}.doctored.vcf.gz
    /opt/hall-lab/htslib-1.9/bin/tabix -p vcf ${basename}.doctored.vcf.gz
    tar -czvf ${basename}.MantaWorkflow.tgz MantaWorkflow
  }
  runtime {
    docker: "apregier/manta_samtools@sha256:cdaf55d96ebb8473383d2190d1a762991c54ca85a851e473a5b1538cb0c3fed1"
    cpu: "8"
    memory: "16 GiB"
    disks: "local-disk " + ceil( size(input_cram, "GB") * 4 + size(input_cram_index, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_cache, "GB") * 5 + 20.0) + " SSD"
    preemptible: preemptible_tries
  }
  output {
    File output_vcf = "${basename}.doctored.vcf.gz"
    File output_tbi = "${basename}.doctored.vcf.gz.tbi"
    File workflow_tgz = "${basename}.MantaWorkflow.tgz"
  }
}

# Smoove wrapper
task Smoove {
  String basename
  File input_cram
  File input_cram_index

  File ref_fasta
  File ref_fasta_index
  File ref_cache
  File exclude_regions

  Int preemptible_tries


  command {
    set -e
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s
    
    export SMOOVE_NO_MAX_CI=TRUE

    smoove call \
      --name ${basename} \
      --exclude ${exclude_regions} \
      --fasta ${ref_fasta} \
      --noextrafilters \
      --genotype \
      ${basename}.cram

    if [ ! -e ${basename}.histo ]; then
      mv *.histo ${basename}.histo
      mv *.split.bam ${basename}.split.bam
      mv *.split.bam.bai ${basename}.split.bam.bai
      mv *.disc.bam ${basename}.disc.bam
      mv *.disc.bam.bai ${basename}.disc.bam.bai
    fi
  }

  runtime {
    docker: "brentp/smoove@sha256:c839ed223462a1c1ae26e7acc27f28f0f67b4581d80a06823895f295ad2bdaf4"
    cpu: "1"
    memory: "2.5 GiB"
    disks: "local-disk " + ceil( size(input_cram, "GB") + size(input_cram_index, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(exclude_regions, "GB") + size(input_cram, "GB") * 0.30 + size(ref_cache, "GB") * 5) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}-smoove.genotyped.vcf.gz"
    File output_csi = "${basename}-smoove.genotyped.vcf.gz.csi"
    File output_histogram = "${basename}.histo"
    File lumpy_script = "${basename}-lumpy-cmd.sh"
    File splitters = "${basename}.split.bam"
    File splitters_index = "${basename}.split.bam.bai"
    File discordants = "${basename}.disc.bam"
    File discordants_index = "${basename}.disc.bam.bai"
  }
}

task Genotype {
  String basename
  File input_cram
  File input_cram_index
  File input_vcf
  File ref_cache
  Int preemptible_tries

  command {
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    rm -f ${basename}.cram.json
    zcat ${input_vcf} \
      | svtyper \
      -B ${basename}.cram \
      -l ${basename}.cram.json \
      > ${basename}.gt.vcf
  }

  runtime {
    docker: "halllab/svtyper@sha256:8839f10f48c254470077ffc2c8f4c6cf78afa8ff6cf34129c1a53c39bfe1987f"
    cpu: "1"
    memory: "6.5 GB"
    disks: "local-disk " + ceil( size(input_cram, "GB") + size(input_vcf, "GB") +  size(ref_cache, "GB") * 5 + 20.0) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.gt.vcf"
    File output_lib = "${basename}.cram.json"
  }
}

task Genotype_Zip {
  String basename
  File input_cram
  File input_cram_index
  File input_vcf
  File ref_cache
  Int preemptible_tries

  command {
    ln -s ${input_cram} ${basename}.cram
    ln -s ${input_cram_index} ${basename}.cram.crai

    # build the reference sequence cache
    tar -zxf ${ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    rm -f ${basename}.cram.json
    zcat ${input_vcf} \
      | svtyper \
      -B ${basename}.cram \
      -l ${basename}.cram.json \
      | bgzip -c \
      > ${basename}.gt.vcf.gz
  }

  runtime {
    docker: "halllab/svtyper@sha256:8839f10f48c254470077ffc2c8f4c6cf78afa8ff6cf34129c1a53c39bfe1987f"
    #docker: "gcr.io/washu-genome-inh-dis-analysis/svtyper:v0.7.0-5fc3076"
    cpu: "1"
    memory: "6.5 GB"
    disks: "local-disk " + ceil( size(input_cram, "GB") + size(ref_cache, "GB") * 5 + 20.0) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.gt.vcf.gz"
    File output_lib = "${basename}.cram.json"
  }
}

task Copy_Number {
  String basename
  String sample
  File input_vcf
  File input_cn_hist_root
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
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "4 GB"
    disks: "local-disk " + 35 + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.cn.vcf"
  }
}

task Zip {
  File input_vcf
  String basename
  Int preemptible_tries

  command {
    cat ${input_vcf} \
    | bgzip -c > ${basename}.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "4 GB"
    disks: "local-disk " +  ceil(2.0*size(input_vcf, "GB")) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.gz"
  }
}

task Copy_Number_Zip {
  String basename
  String sample
  File input_vcf
  File input_cn_hist_root
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
      | bgzip -c \
      > ${basename}.cn.vcf.gz
  }

  runtime {
   
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "6.5 GB"
    disks: "local-disk " + 35 + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.cn.vcf.gz"
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
  Int preemptible_tries
  Int threads = 4
  # Add 7G of pad of the chromosome directory and ~2-3 GB of output files
  Int disk_size = ceil( size(input_cram, "GB") + size(input_cram_index, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_cache, "GB") * 5 + 7.0 )

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
    docker: "halllab/cnvnator@sha256:8bf4fa64a288c5647a9a6b1ea90d14e76f48a3e16c5bf98c63419bb7d81c8938"
    cpu: threads
    memory: "16 GB"
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
   sleep 1

    svtools lsort \
      -b 200 \
      -f input_vcfs_file.local_map.txt \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: 30
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
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Filter_Del {
  File input_vcf_gz
  String output_vcf_basename
  Int disk_size
  Int preemptible_tries
  
  command <<<
    set -eo pipefail

    bcftools view -i '(SVTYPE!="DEL" || SVLEN>1000 || SVLEN<-1000 || INFO/SR>0)' ${input_vcf_gz}  | bgzip -c >  ${output_vcf_basename}.vcf.gz
  >>>

  runtime {
    docker: "halllab/bcftools@sha256:955cbf93e35e5ee6fdb60e34bb404b7433f816e03a202dfed9ceda542e0d8906"
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
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "12 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 0
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Remove_INS {
  File input_vcf_gz
  String output_vcf_basename
  Int preemptible_tries

  command <<<
    zcat ${input_vcf_gz} \
    | awk '{if($5!="<INS>") print $0}' \
    | bgzip -c \
    > ${output_vcf_basename}.vcf.gz
  >>>

  runtime {
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " +  2*ceil( size(input_vcf_gz, "GB")) + " HDD"
    preemptible: preemptible_tries
  }
  
  output {
    File output_vcf_gz =  "${output_vcf_basename}.vcf.gz"
  }
}

task Prune_VCF_Output_Bedpe {
  File input_vcf_gz
  String output_vcf_basename
  Int preemptible_tries

  command {
    zcat ${input_vcf_gz} \
      | svtools afreq \
      | svtools vcftobedpe \
      | svtools bedpesort \
      | svtools prune -s -d 100 -e 'AF' \
      | bgzip -c \ 
      > ${output_vcf_basename}.bedpe.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "8 GB"
    disks: "local-disk " +  4*ceil( size(input_vcf_gz, "GB")) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_bedpe_gz = "${output_vcf_basename}.bedpe.gz"
  }
}

task Prune_VCF {
  File input_vcf_gz
  String output_vcf_basename
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
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " +  3*ceil( size(input_vcf_gz, "GB")) + " HDD"
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
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " +  3*ceil( size(input_vcf_gz, "GB")) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Concat_VCF {
  File input_del_vcf_gz
  File input_dup_vcf_gz
  File input_bnd_vcf_gz
  String output_vcf_basename
  Int preemptible_tries

  command {
    zcat ${input_dup_vcf_gz} \
    | grep "^##"  > header.txt

    zcat ${input_del_vcf_gz} ${input_dup_vcf_gz} ${input_bnd_vcf_gz} \
    | grep -v "^##" \
    | cat header.txt - \
    | bgzip -c > ${output_vcf_basename}.vcf.gz
  }
  runtime {
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " + 3*ceil( size(input_del_vcf_gz, "GB") + size(input_bnd_vcf_gz, "GB") + size(input_dup_vcf_gz, "GB")) + "  HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}
  
task Sort_Index_VCF {
  File input_vcf_gz
  String output_vcf_name
  Int preemptible_tries

  command {
    zcat ${input_vcf_gz} \
      | svtools vcfsort \
      | bgzip -c \
      > ${output_vcf_name}

    tabix -p vcf -f ${output_vcf_name}
  }

  runtime {
    docker: "halllab/svtools@sha256:7571b6e9cbfeba7ebfdefd490e8315ed5742ad034ca075d1f70fc422786cdff3"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " + 3*ceil( size(input_vcf_gz, "GB")) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_name}"
    File output_vcf_gz_index = "${output_vcf_name}.tbi"
  }
}

