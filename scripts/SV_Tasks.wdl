version 1.0
# get the sample (SM) field from a CRAM file
task Split_By_Type {
  input {
    File input_vcf
    String output_vcf_prefix
    Int preemptible_tries
  }
  command <<<
    set -eo pipefail
    zcat ~{input_vcf} | grep -v "random	" | grep -v "alt	" | grep -v "decoy	" | grep -v "EBV	" | grep -v "^chrUn" | grep -v "^HLA"| /opt/hall-lab/vawk/vawk -v svtype=BND --header '{if(I$SVTYPE==svtype) print $0;}' | /opt/hall-lab/htslib-1.9/bin/bgzip -c > ~{output_vcf_prefix}.bnd.vcf.gz
    zcat ~{input_vcf} | grep -v "random	" | grep -v "alt	" | grep -v "decoy	" | grep -v "EBV	" | grep -v "^chrUn" | grep -v "^HLA"| /opt/hall-lab/vawk/vawk -v svtype=DEL --header '{if(I$SVTYPE==svtype) print $0;}' | /opt/hall-lab/htslib-1.9/bin/bgzip -c > ~{output_vcf_prefix}.del.vcf.gz
    zcat ~{input_vcf} | grep -v "random	" | grep -v "alt	" | grep -v "decoy	" | grep -v "EBV	" | grep -v "^chrUn" | grep -v "^HLA"| /opt/hall-lab/vawk/vawk -v svtype=INS --header '{if(I$SVTYPE==svtype) print $0;}' | /opt/hall-lab/htslib-1.9/bin/bgzip -c > ~{output_vcf_prefix}.ins.vcf.gz
    zcat ~{input_vcf} | grep -v "random	" | grep -v "alt	" | grep -v "decoy	" | grep -v "EBV	" | grep -v "^chrUn" | grep -v "^HLA"| /opt/hall-lab/vawk/vawk --header '{if(I$SVTYPE!="DEL" && I$SVTYPE!="BND" && I$SVTYPE!="INS") print $0;}' | /opt/hall-lab/htslib-1.9/bin/bgzip -c > ~{output_vcf_prefix}.other.vcf.gz
    zcat ~{output_vcf_prefix}.ins.vcf.gz | \
    /opt/hall-lab/vawk/vawk '{ct=split(I$SNAME, spl, ","); for(ii=1; ii<=ct; ii++) print $3, spl[ii], $9}' | \
    /opt/hall-lab/htslib-1.9/bin/bgzip -c > ~{output_vcf_prefix}.ins_split.txt.gz
  >>>
  runtime {
    docker: "halllab/vcf_bed_utils@sha256:09c18a5827d67891792ffc110627c7fa05b2262df4b91d6967ad6e544f41e8ec"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + ceil( size(input_vcf, "GB") * 2) + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File bnd_vcf = "${output_vcf_prefix}.bnd.vcf.gz"
    File del_vcf = "${output_vcf_prefix}.del.vcf.gz"
    File ins_vcf = "${output_vcf_prefix}.ins.vcf.gz"
    File other_vcf = "${output_vcf_prefix}.other.vcf.gz"
    File ins_split = "${output_vcf_prefix}.ins_split.txt.gz"
  }
}

task Get_Sample_Name {
  input {
  	File input_cram
  	Int preemptible_tries
  }

  command {
    set -eo pipefail
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
  input {
  	File input_cn_hist_root
  	File ref_fasta_index
  	Int preemptible_tries
  }

  command <<<
    set -eo pipefail
    cat ~{ref_fasta_index} \
      | awk '$1=="chrX" { print $1":0-"$2 } END { print "exit"}' \
      | cnvnator -root ~{input_cn_hist_root} -genotype 100 \
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
  input {
    Array[String] sample_array
    Array[String] sex_array
    String output_ped_basename
    File sample_file = write_lines(sample_array)
    File sex_file = write_lines(sex_array)
  }

  command <<<
    set -eo pipefail
    paste ~{sample_file} ~{sex_file} \
      | awk '{ print $1,$1,-9,-9,$2,-9 }' OFS='\t' \
      > ~{output_ped_basename}.ped
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

# index a CRAM
task Index_Cram {
  input {
    File input_cram
    String basename
    File ref_cache
    Int preemptible_tries
  }

  command {
    set -eo pipefail
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

task Filter_Index {
  input {
    File input_vcf_gz
    String output_vcf_name
    Int preemptible_tries
  }

  command <<<
	set -eo pipefail
	FILTERLINE='##FILTER=<ID=LOW,Description="Test Low quality filter">'
	zcat ~{input_vcf_gz} | \
		/opt/hall-lab/vawk/vawk '{ \
		split(I$STRANDS,x,","); \
		split(x[1],y,":"); \
		split(x[2],z,":"); \
		if (I$SVTYPE=="INS" && I$NSAMP>0) { \
		I$MSQ=QUAL/I$NSAMP; \
		gsub("MSQ=0.00", "MSQ="I$MSQ, $8) \
		} \
		if ((I$SVTYPE=="DEL" || I$SVTYPE=="DUP" || I$SVTYPE=="MEI") && \
		I$MSQ>=100 && sqrt((I$SVLEN)*(I$SVLEN))>=50){ \
		$7="PASS"; print $0; \
		}  else if ( I$SVTYPE=="INV" && $6>=100 && (I$SR/I$SU)>=0.1 && (I$PE/I$SU)>=0.1 && (y[2]/I$SU)>0.1 && (z[2]/I$SU)>0.1 && sqrt((I$SVLEN)*(I$SVLEN))>=50){ \
		$7="PASS"; print $0; \
		} else if ( I$SVTYPE=="BND" && $9 !~ /CN/ && I$MSQ>=500){ \
		$7="PASS"; print $0; \
		} else if ( I$SVTYPE=="BND" && $9 ~ /CN/ && I$MSQ>=250){ \
		$7="PASS"; print $0; \
		} else if ( I$SVTYPE=="INS" && I$MSQ>=100 && I$SVLEN >=50) { \
		$7="PASS"; print $0; \
		} else { \
		$7="LOW"; print $0; \
		} \
	}' |  cat <(zcat ~{input_vcf_gz} | sed -n '/^#[^#]/q;p') <(echo $FILTERLINE) <(zgrep -m 1 '^#CHROM' ~{input_vcf_gz}) - | /opt/hall-lab/htslib-1.9/bin/bgzip -c > ~{output_vcf_name}
	/opt/hall-lab/htslib-1.9/bin/tabix -p vcf -f ~{output_vcf_name}
  >>>

  runtime {
    docker: "halllab/vcf_bed_utils@sha256:09c18a5827d67891792ffc110627c7fa05b2262df4b91d6967ad6e544f41e8ec"
    cpu: "1"
    memory: "1 GB"
    disks: "local-disk " + ceil( size(input_vcf_gz, "GB") * 2) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_name}"
    File output_vcf_gz_index = "${output_vcf_name}.tbi"
  }
  
}

task Count_Lumpy {
  input {
    String basename
    File input_vcf
    Int preemptible_tries
    String cohort
    String center
  }

  command <<<
    set -eo pipefail
   
     bcftools query  -f "[%CHROM\t~{cohort}\t~{center}\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SR\t%SAMPLE\t%GT\n]"  ~{input_vcf} \
     | awk 'BEGIN{OFS="\t"}{if($1~/chr[1-9]+/ && $1!~/_/) {
       svlen=$6;
       if($6<0 && $6!=".") svlen=-1*$6;
       len_bin=">=1kb"
       if(svlen<1000) len_bin="<1kb";
       if($7>0) $7="SR>=1";
       else $7="SR=0";
       print $1, $2, $3, $4, $5, len_bin, $7, $8, $9;}}' \
     | sort -k1,9 \
     | uniq -c \
     | awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $5, $6, $7, $8, $9, $1}'  > ~{basename}.lumpy.counts.1.txt
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
  input {
    String basename
    File input_vcf
    Int preemptible_tries
    String cohort
    String center
  }

  command <<<
    set -eo pipefail

     bcftools query  -f "[%CHROM\t~{cohort}\t~{center}\t%FILTER\t%INFO/SVTYPE\t%SAMPLE\t%GT\n]"  ~{input_vcf} \
     | awk 'BEGIN{OFS="\t"}{if($1~/chr[1-9]+/ && $1!~/_/ && $4=="PASS") print $0;}' \
     | sort -k1,7 \
     | uniq -c \
     | awk 'BEGIN{OGS="\t"}{print $2, $3, $4, $5, $6, $7, $8,  $1}'  > ~{basename}.manta.counts.1.txt
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

task Manta {    
  input {
    File input_cram
    File input_cram_index
    File ref_fasta
    File ref_fasta_index
    File ref_cache
    File? call_regions_bed
    File? call_regions_bed_index
    String basename
    Int preemptible_tries
  }

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
    set -eo pipefail
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
    docker: "halllab/manta_samtools@sha256:d39fac59a2c06f808d115c65b9c191baf5f249769d317263ae3cd19e2c74d20e"
    cpu: "8"
    memory: "16 GiB"
    disks: "local-disk " + ceil( size(input_cram, "GB") * 4 + size(input_cram_index, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_cache, "GB") * 5 + 20.0) + " SSD"
    preemptible: preemptible_tries
  }
  output {
    File output_vcf = "${basename}.doctored.vcf.gz"
    File output_tbi = "${basename}.doctored.vcf.gz.tbi"
    File original_vcf = "${basename}.vcf.gz"
    File original_tbi = "${basename}.vcf.gz.tbi"
    File workflow_tgz = "${basename}.MantaWorkflow.tgz"
  }
}

# Smoove wrapper
task Smoove {
  input {
    String basename
    File input_cram
    File input_cram_index

    File ref_fasta
    File ref_fasta_index
    File ref_cache
    File exclude_regions

    Int preemptible_tries
  }

  command {
    set -eo pipefail
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
  input {
    String basename
    File input_cram
    File input_cram_index
    File input_vcf
    File ref_cache
    Int preemptible_tries
  }

  command {
    set -eo pipefail
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
      | bgzip -c > ${basename}.gt.vcf.gz
  }

  runtime {
    docker: "halllab/svtyper@sha256:8ebb0508bc63a2a32d22b4a3e55453222560daa30b7cc14a4f1189cb311d5922"
    cpu: "1"
    memory: "15 GB"
    disks: "local-disk " + ceil( size(input_cram, "GB") + size(input_vcf, "GB") +  size(ref_cache, "GB") * 5 + 20.0) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.gt.vcf.gz"
    File output_lib = "${basename}.cram.json"
  }
}

task Take_Original_Genotypes {
  input {
    String sample_name
    String basename
    File input_vcf
    File input_variant_to_sname_mapping
    File original_per_sample_vcf
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail
    zcat ~{input_variant_to_sname_mapping} \
      | /opt/hall-lab/vawk/vawk -v sname="~{sample_name}" 'BEGIN{OFS="\t"}{ \
    split($2, spl, ":"); \
    if(spl[1]==sname) { \
        print $1, spl[1], spl[2]":"spl[3]":"spl[4]":"spl[5]":"spl[6]":"spl[7]":"spl[8]; \
        } \
      }' \
      | /opt/hall-lab/io/zjoin -a stdin -b <(paste -d ":" <(zcat ~{original_per_sample_vcf} | grep -v "^#" | cut -f 3,9-) <(zcat ~{original_per_sample_vcf} | grep -v "^#" | cut -f 4,5 | tr "\t" ":") <(zcat ~{original_per_sample_vcf} | /opt/hall-lab/vawk/vawk '{svlen=I$SVLEN; if(svlen==""){svlen="."} print svlen}') | sed 's/:SR/:SR:OREF:OALT:OSVLEN/') -1 3 -2 1 \
      | cut -f 1,5- \
      | awk -v sname="~{sample_name}" 'BEGIN{OFS="\t"; print "ID", "FORMAT", sname;}{ \
           print $0; \
        }' \
      | /opt/hall-lab/htslib-1.9/bin/bgzip -c > temp

    zcat ~{input_vcf} \
      | /opt/hall-lab/io/zjoin -r -p "##" -a stdin -b <(zcat temp | sort -k1,1 | /opt/hall-lab/bin/bedtools groupby -g 1 -c 2,3 -o first,first ) -1 3 -2 1 \
      | cut -f -8,10- \
      | /opt/hall-lab/vawk/vawk --header 'BEGIN{OFS="\t"}{if($9=="NA") {$9="GT:FT:GQ:PL:PR:SR:OREF:OALT:OSVLEN"; $10="0/0:.:.:.:.:.:.:.:.";} print $0;}' \
      | sed 's/^#CHROM/##FORMAT=<ID=OREF,Number=1,Type=String,Description="Original reference sequence">\n##FORMAT=<ID=OALT,Number=1,Type=String,Description="Original alt sequence">\n##FORMAT=<ID=OSVLEN,Number=1,Type=Integer,Description="Original SVLEN">\n#CHROM/' \
      | /opt/hall-lab/htslib-1.9/bin/bgzip -c > ~{basename}.gt.vcf.gz
  >>> 

  runtime {
    docker: "halllab/vcf_bed_utils@sha256:09c18a5827d67891792ffc110627c7fa05b2262df4b91d6967ad6e544f41e8ec"
    cpu: "1"
    memory: "15 GB"
    disks: "local-disk " + ceil( size(original_per_sample_vcf, "GB") + size(input_vcf, "GB") +  size(input_variant_to_sname_mapping, "GB") + 20.0) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.gt.vcf.gz"
  }
}

task Copy_Number {
  input {
    String basename
    String sample
    File input_vcf
    File input_cn_hist_root
    File ref_cache
    Int preemptible_tries
  }

  command {
    set -eo pipefail
    zcat ${input_vcf} \
     | create_coordinates \
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
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
    cpu: "1"
    memory: "4 GB"
    disks: "local-disk " + 35 + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = "${basename}.cn.vcf.gz"
  }
}

task CNVnator_Histogram {
  input {
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
  }

  command <<<
    set -eo pipefail
    ln -s ~{input_cram} ~{basename}.cram
    ln -s ~{input_cram_index} ~{basename}.cram.crai

    # build the reference sequence cache
    tar -zxf ~{ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    # Create directory of chromosome FASTA files for CNVnator
    mkdir -p ~{ref_chrom_dir}
    awk -v CHROM_DIR=~{ref_chrom_dir} 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > CHROM_DIR"/"CHROM".fa" }' ~{ref_fasta}

    cnvnator_wrapper.py \
      -T cnvnator.out \
      -o ~{basename}.cn \
      -t ~{threads} \
      -w 100 \
      -b ~{basename}.cram \
      -c ~{ref_chrom_dir} \
      -g GRCh38 \
      --cnvnator cnvnator
  >>>

  runtime {
    docker: "halllab/cnvnator@sha256:8bf4fa64a288c5647a9a6b1ea90d14e76f48a3e16c5bf98c63419bb7d81c8938"
    cpu: threads
    memory: "16 GB"
    disks: "local-disk " + ceil( size(input_cram, "GB") + size(input_cram_index, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_cache, "GB") * 5 + 7.0 ) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_cn_hist_root = "cnvnator.out/${basename}.cram.hist.root"
    File output_cn_txt = "${basename}.cn.txt"
    File output_cn_bed = "${basename}.cn.bed"
  }
}

task L_Sort_VCF_Variants {
  input {
    Array[File] input_vcfs
    File input_vcfs_file = write_lines(input_vcfs)
    String output_vcf_basename
    Int preemptible_tries
  }

  parameter_meta {
    input_vcfs: {
	description: "vcf files to sort together",
        localization_optional: true
    }
  }

  command {
    set -eo pipefail
    # strip the "gs://" prefix from the file paths
    cat ${input_vcfs_file} \
      | sed 's/^gs:\/\//\.\//g' \
      > ${input_vcfs_file}.local_map.txt
   sleep 1

    svtools lsort \
      -b 200 \
      -f ${input_vcfs_file} \
      -t /cromwell_root/bulk_download \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + 1.5*size(input_vcfs, "GB") + " HDD"
    bootDiskSizeGb: 30
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task L_Merge_VCF_Variants {
  input {
    File input_vcf_gz
    String output_vcf_basename
    Int preemptible_tries
  }

  command {
    set -eo pipefail
    zcat ${input_vcf_gz} \
      | svtools lmerge \
      -i /dev/stdin \
      -f 20 \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + 2*size(input_vcf_gz, "GB")+10 + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task L_Merge_VCF_Variants_weighted {
  input {
    File input_vcf_gz
    String output_vcf_basename
    Int preemptible_tries
  }

  command {
    set -eo pipefail
    zcat ${input_vcf_gz} \
      | svtools lmerge \
      -i /dev/stdin \
      -f 20 \
      -w carrier_wt \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + 2*size(input_vcf_gz, "GB")+10 + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Filter_Del {
  input {
    File input_vcf_gz
    String output_vcf_basename
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail

    bcftools view -i '(SVTYPE!="DEL" || SVLEN>1000 || SVLEN<-1000 || INFO/SR>0)' ~{input_vcf_gz}  | bgzip -c >  ~{output_vcf_basename}.vcf.gz
  >>>

  runtime {
    docker: "halllab/bcftools@sha256:955cbf93e35e5ee6fdb60e34bb404b7433f816e03a202dfed9ceda542e0d8906"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + 2*size(input_vcf_gz, "GB")+10 + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Filter_Pass {
  input {
    File input_vcf_gz
    String output_vcf_basename
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail

    bcftools view -f .,PASS ~{input_vcf_gz}  | bgzip -c >  ~{output_vcf_basename}.vcf.gz
  >>>

  runtime {
    docker: "halllab/bcftools@sha256:955cbf93e35e5ee6fdb60e34bb404b7433f816e03a202dfed9ceda542e0d8906"
    cpu: "1"
    memory: "3.75 GB"
    disks: "local-disk " + 2*size(input_vcf_gz, "GB")+10 + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Paste_VCF {
  input {
    Array[File] input_vcfs
    File input_vcfs_file = write_lines(input_vcfs)
    String output_vcf_basename
    Int preemptible_tries
  }
  parameter_meta {
    input_vcfs: {
	description: "vcf files to paste together",
        localization_optional: true
    }
  }

  command {
    set -eo pipefail
    svtools vcfpaste \
      -f ${input_vcfs_file} \
      -q \
      -t /cromwell_root/bulk_download \
      | bgzip -c \
      > ${output_vcf_basename}.vcf.gz
  }

  runtime {
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
    cpu: "1"
    memory: "12 GB"
    disks: "local-disk " + 1.5*size(input_vcfs, "GB") + " HDD"
    preemptible: 0
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Remove_INS {
  input {
    File input_vcf_gz
    String output_vcf_basename
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail
    zcat ~{input_vcf_gz} \
    | awk '{if($5!="<INS>") print $0}' \
    | bgzip -c \
    > ~{output_vcf_basename}.vcf.gz
  >>>

  runtime {
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " +  2*ceil( size(input_vcf_gz, "GB")) + " HDD"
    preemptible: preemptible_tries
  }
  
  output {
    File output_vcf_gz =  "${output_vcf_basename}.vcf.gz"
  }
}

task Prune_VCF {
  input {
    File input_vcf_gz
    String output_vcf_basename
    Int preemptible_tries
  }

  command {
    set -eo pipefail
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
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
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
  input {
    File input_vcf_gz
    File input_ped
    String output_vcf_basename
    File mei_annotation_bed
    Int preemptible_tries
  }

  command {
    set -eo pipefail
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
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " +  10*ceil( size(input_vcf_gz, "GB")) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_basename}.vcf.gz"
  }
}

task Sort_Index_VCF {
  input {
    File input_vcf_gz
    String output_vcf_name
    Int preemptible_tries
  }

  command {
    set -eo pipefail
    zcat ${input_vcf_gz} \
      | svtools vcfsort \
      | bgzip -c \
      > ${output_vcf_name}

    tabix -p vcf -f ${output_vcf_name}
  }

  runtime {
    docker: "halllab/svtools@sha256:38ac08a8685ff58329b72e2b9c366872086d41ef21da84278676e06ef7f1bfbb"
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " + 20*ceil( size(input_vcf_gz, "GB")) + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf_gz = "${output_vcf_name}"
    File output_vcf_gz_index = "${output_vcf_name}.tbi"
  }
}

