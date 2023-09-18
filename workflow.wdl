version 1.0

workflow lr_sv_calling {
    meta {
    author: "Jean Monlong"
        email: "jean.monlong@gmail.com"
        description: "Calling structural variants from long reads"
    }

    parameter_meta {
        FASTQ: "Raw FASTQ reads. Optional. If provided and no BAM/BAI, alignment will be performed."
        BAM: "Aligned reads in BAM. Optional. If provided, no alignment performed."
        BAI: "Aligned reads index. Optional. If provided, no alignment performed."
        ALIGNER: "Which aligner to use? 'ngmlr' (default)"
        CALLER: "Which caller to use? 'sniffles' (default)"
        REFERENCE_FASTA: "Genome reference FASTA file"
        SAMPLE: "Optional. Sample name. "
    }

    input {
        File? FASTQ
        File? BAM
        File? BAI
        String ALIGNER = 'ngmlr'
        String CALLER = 'sniffles'
        File? VNTR
        File REFERENCE_FASTA
        String SAMPLE = "sample"
    }

    if ( defined(FASTQ) && !defined(BAM) ){
        if ( ALIGNER == 'ngmlr') {
            call alignReadsNGMLR {
                input:
                fastq=select_first([FASTQ]),
                reference_fa=REFERENCE_FASTA,
                sample=SAMPLE
            }
        }

        File alignedBamFile = select_first([alignReadsNGMLR.bam])

        call sortReads {
            input:
            bam=alignedBamFile,
            sample=SAMPLE
        }
        
    }

    File bamFile = select_first([BAM, sortReads.sorted_bam])
    File baiFile = select_first([BAI, sortReads.sorted_bai])

    if (CALLER == 'sniffles') {
        call callSVsniffles {
            input:
            bam=bamFile,
            bai=baiFile,
            sample=SAMPLE,
            vntr=VNTR
        }
    }

    File svVcf = select_first([callSVsniffles.sv_vcf])
    
    output {
        File sv_vcf = svVcf
        File? bam = sortReads.sorted_bam
        File? bai = sortReads.sorted_bai
    }
}

##
#### TASKS
##

task alignReadsNGMLR {
    input {
        File fastq
        File reference_fa
        String sample = ""
        Int threadCount = 32
        Int memoryGB = 16
        String dockerImage="quay.io/jmonlong/ngmlr"
        Int disk_size = 3 * round(size(fastq, 'G') + size(reference_fa, 'G')) + 20
    }

    Int threadSort = if threadCount > 8 then 4 else 1
    Int threadAlign = if threadCount > 8 then threadCount - 4 else threadCount - 1
    
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        ln -s ~{reference_fa} ref.fa

        ngmlr -t ~{threadAlign} -x ont --skip-write -r ref.fa -q ~{fastq} | samtools sort -@ ~{threadSort} -o ~{sample}_ngmlr.bam
        samtools index -@ ~{threadCount} ~{sample}_ngmlr.bam ~{sample}_ngmlr.bai
    >>>
    output {
        File bam = "~{sample}_ngmlr.bam"
        File bai = "~{sample}_ngmlr.bai"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}

task sortReads {
    input {
        File bam
        String sample = ""
        Int threadCount = 8
        Int memoryGB = 8
        String dockerImage="quay.io/biocontainers/samtools@sha256:6f88956b747a67b2a39a3ff72c4de30e665239ee11db610624dd4298e30db1bf"
        Int disk_size = 3 * round(size(bam, 'G')) + 20
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        samtools sort -@ ~{threadCount} ~{bam} -o ~{sample}_sorted.bam
        samtools index -@ ~{threadCount} ~{sample}_sorted.bam ~{sample}_sorted.bai
    >>>
    output {
        File sorted_bam = "~{sample}_sorted.bam"
        File sorted_bai = "~{sample}_sorted.bai"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}

task callSVsniffles {
    input {
        File bam
        File bai
	File vntr = ""
        String sample = ""
	Int minSvLen = 40
        Int threadCount = 8
        Int memoryGB = 32
        String dockerImage = "quay.io/biocontainers/sniffles@sha256:feb1c41eae608ebc2c2cb594144bb3c221b87d9bb691d0af6ad06b49fd54573a"
        Int disk_size = 3 * round(size(bam, 'G')) + 20
  }
  
  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    TRF_ARG=""
    if [ ! -z ~{vntr} ]
    then
       TRF_ARG="--tandem-repeats ~{vntr}"
    fi

    ln -s ~{bam} reads.bam
    ln -s ~{bai} reads.bam.bai
    
    sniffles -i reads.bam -v ~{sample}_sniffles.vcf -t ~{threadCount} ${TRF_ARG} --minsvlen ~{minSvLen} 2>&1 | tee sniffles.log
    gzip ~{sample}_sniffles.vcf
  >>>

  output {
      File sv_vcf = "~{sample}_sniffles.vcf.gz"
      File log = "sniffles.log"
  }

  runtime {
      docker: dockerImage
      cpu: threadCount
      disks: "local-disk " + disk_size + " SSD"
      memory: memoryGB + "GB"
  }
}
