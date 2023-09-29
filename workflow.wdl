version 1.0

workflow lr_sv_calling {
    meta {
    author: "Jean Monlong"
        email: "jean.monlong@gmail.com"
        description: "Calling structural variants from long reads"
    }

    parameter_meta {
        FASTQ: "Raw FASTQ reads. Optional. If provided and no BAM/BAI, alignment will be performed."
        FASTQ_IS_UBAM: "Is the FASTQ file actually a unmapped BAM file? Default is false"
        BAM: "Aligned reads in BAM. Optional. If provided, no alignment performed."
        BAI: "Aligned reads index. Optional. If provided, no alignment performed."
        ALIGNER: "Which aligner to use? 'ngmlr' (default), 'minimap2', or 'winnowmap2'"
        CALLER: "Which caller to use? 'sniffles' (default)"
        REFERENCE_FASTA: "Genome reference FASTA file"
        SAMPLE: "Optional. Sample name. "
        READS_PER_CHUNK: "Number of input reads per chunk (to scale up alignment). Default is 700000"
        VNTR: "Optional. BED annotation of VNTR location for Sniffles."
        REPK_WM: "Optional. K-mer file used by Winnowmap2. Will be recomputed if not provided."
    }

    input {
        File? FASTQ
        Boolean FASTQ_IS_UBAM = false
        File? BAM
        File? BAI
        String ALIGNER = 'ngmlr'
        String CALLER = 'sniffles'
        File? VNTR
        File? REPK_WM
        File REFERENCE_FASTA
        String SAMPLE = "sample"
        Int READS_PER_CHUNK = 700000
    }

    if ( defined(FASTQ) && !defined(BAM) ){

        call splitReads {
	    input:
	    reads = select_first([FASTQ]),
            readsPerChunk = READS_PER_CHUNK,
            ubam=FASTQ_IS_UBAM
	}

        if ( ALIGNER == 'ngmlr') {
            scatter (readChunk in splitReads.readChunks){
                call alignReadsNGMLR {
                    input:
                    fastq=readChunk,
                    reference_fa=REFERENCE_FASTA,
                    sample=SAMPLE
                }
            }
        }

        if ( ALIGNER == 'minimap2') {
            scatter (readChunk in splitReads.readChunks){
                call alignReadsMinimap2 {
                    input:
                    fastq=readChunk,
                    reference_fa=REFERENCE_FASTA,
                    sample=SAMPLE
                }
            }
        }

        if ( ALIGNER == 'winnowmap2') {
            if (!defined(REPK_WM)){
                call prepareWinnowmap2Ref {
                    input:
                    reference_fa=REFERENCE_FASTA
                }
            }

            File repk_txt = select_first([REPK_WM, prepareWinnowmap2Ref.repk])
            
            scatter (readChunk in splitReads.readChunks){
                call alignReadsWinnowmap2 {
                    input:
                    fastq=readChunk,
                    reference_fa=REFERENCE_FASTA,
                    sample=SAMPLE,
                    repk_txt=repk_txt
                }
            }
        }

        Array[File] alignedBamFiles = select_first([alignReadsNGMLR.bam, alignReadsMinimap2.bam, alignReadsWinnowmap2.bam])

        call mergeBAM {
	    input:
	    bams = alignedBamFiles,
            outname = SAMPLE,
	}
        
    }

    File bamFile = select_first([BAM, mergeBAM.bam])
    File baiFile = select_first([BAI, mergeBAM.bai])

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
        File? bam = mergeBAM.bam
        File? bai = mergeBAM.bai
        File? repk = prepareWinnowmap2Ref.repk
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
        Int thread_count = 32
        Int memory_gb = 32
        String dockerImage="quay.io/jmonlong/ngmlr@sha256:ce25d81d1a44f7bcdacef0008ac7542c5e9885d074f6446d6a8549219c91807e"
        Int disk_size = 3 * round(size(fastq, 'G') + size(reference_fa, 'G')) + 20
        Int preemptible = 1
    }

    Int threadSort = if thread_count > 8 then 4 else 1
    Int threadAlign = if thread_count > 8 then thread_count - 4 else thread_count - 1
    
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

        ngmlr -t ~{threadAlign} -x ont --skip-write -r ref.fa -q ~{fastq} --no-progress | samtools sort -@ ~{threadSort} -o ~{sample}_ngmlr.bam
        samtools index -@ ~{thread_count} ~{sample}_ngmlr.bam ~{sample}_ngmlr.bai
    >>>
    output {
        File bam = "~{sample}_ngmlr.bam"
        File bai = "~{sample}_ngmlr.bai"
    }
    runtime {
        preemptible: preemptible
        docker: dockerImage
        cpu: thread_count
        disks: "local-disk " + disk_size + " SSD"
        memory: memory_gb + "GB"
    }
}

task alignReadsMinimap2 {
    input {
        File fastq
        File reference_fa
        String sample = ""
        Int thread_count = 32
        Int memory_gb = 16
        String dockerImage="mkolmogo/card_minimap2:2.23"
        Int disk_size = 3 * round(size(fastq, 'G') + size(reference_fa, 'G')) + 20
    }

    Int threadSort = if thread_count > 8 then 4 else 1
    Int threadAlign = if thread_count > 8 then thread_count - 4 else thread_count - 1
    
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

        minimap2 -a -t ~{threadAlign} -Y -x map-ont ref.fa ~{fastq} | samtools sort -@ ~{threadSort} -o ~{sample}_minimap2.bam
        samtools index -@ ~{thread_count} ~{sample}_minimap2.bam ~{sample}_minimap2.bai
    >>>
    output {
        File bam = "~{sample}_minimap2.bam"
        File bai = "~{sample}_minimap2.bai"
    }
    runtime {
        preemptible: 2
        docker: dockerImage
        cpu: thread_count
        disks: "local-disk " + disk_size + " SSD"
        memory: memory_gb + "GB"
    }
}

task prepareWinnowmap2Ref {
    input {
        File reference_fa
        Int thread_count = 4
        Int memory_gb = 8
        String dockerImage="quay.io/jmonlong/winnowmap@sha256:3c1a1d43451552c30c2d76b326ae83d8db9e51676af69a4116fb212c16132c31"
        Int disk_size = 3 * round(size(reference_fa, 'G')) + 20
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

        ln -s ~{reference_fa} ref.fa

        meryl count k=15 output merylDB ref.fa
        meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
    >>>
    output {
        File repk = "repetitive_k15.txt"
    }
    runtime {
        preemptible: 2
        docker: dockerImage
        cpu: thread_count
        disks: "local-disk " + disk_size + " SSD"
        memory: memory_gb + "GB"
    }
}

task alignReadsWinnowmap2 {
    input {
        File fastq
        File reference_fa
        File repk_txt
        String sample = ""
        Int thread_count = 32
        Int memory_gb = 16
        String dockerImage="quay.io/jmonlong/winnowmap@sha256:3c1a1d43451552c30c2d76b326ae83d8db9e51676af69a4116fb212c16132c31"
        Int disk_size = 3 * round(size(fastq, 'G') + size(reference_fa, 'G') + size(repk_txt, 'G')) + 20
    }

    Int threadSort = if thread_count > 8 then 4 else 1
    Int threadAlign = if thread_count > 8 then thread_count - 4 else thread_count - 1
    
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

        winnowmap -W ~{repk_txt} -ax map-ont -Y -t ~{threadAlign} ref.fa ~{fastq} | samtools sort -@ ~{threadSort} -o ~{sample}_winnowmap2.bam
        samtools index -@ ~{thread_count} ~{sample}_winnowmap2.bam ~{sample}_winnowmap2.bai
    >>>
    output {
        File bam = "~{sample}_winnowmap2.bam"
        File bai = "~{sample}_winnowmap2.bai"
    }
    runtime {
        preemptible: 2
        docker: dockerImage
        cpu: thread_count
        disks: "local-disk " + disk_size + " SSD"
        memory: memory_gb + "GB"
    }
}

task callSVsniffles {
    input {
        File bam
        File bai
	File vntr = ""
        String sample = ""
	Int minSvLen = 40
        Int thread_count = 8
        Int memory_gb = 32
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
    
    sniffles -i reads.bam -v ~{sample}_sniffles.vcf -t ~{thread_count} ${TRF_ARG} --minsvlen ~{minSvLen} --output-rnames --min-alignment-length 500 2>&1 | tee sniffles.log
    gzip ~{sample}_sniffles.vcf
  >>>

  output {
      File sv_vcf = "~{sample}_sniffles.vcf.gz"
      File log = "sniffles.log"
  }

  runtime {
      preemptible: 2
      docker: dockerImage
      cpu: thread_count
      disks: "local-disk " + disk_size + " SSD"
      memory: memory_gb + "GB"
  }
}

task splitReads {
    input {
        File reads
        Int readsPerChunk
        Boolean ubam = false
        Int thread_count = 8
        Int disk_size = 5 * round(size(reads, "G")) + 20
    }

    Int gzThread_count = if thread_count > 1 then thread_count - 1 else 1
    Int chunkLines = readsPerChunk * 4
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        if [ ~{ubam} == "true" ]
        then
            samtools fastq -TMm,Ml,MM,ML ~{reads} | split -l ~{chunkLines} --filter='pigz -p ~{gzThread_count} > ${FILE}.fq.gz' - "fq_chunk.part."
        else
            gzip -cd ~{reads} | split -l ~{chunkLines} --filter='pigz -p ~{gzThread_count} > ${FILE}.fq.gz' - "fq_chunk.part."
        fi
    >>>
    output {
        Array[File] readChunks = glob("fq_chunk.part.*")
    }
    runtime {
        preemptible: 2
        time: 120
        cpu: thread_count
        memory: "4 GB"
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1_pigz"
    }
}

task mergeBAM {
    input {
        Array[File] bams
        String outname = "merged"
        Int thread_count = 8
        Int disk_size = round(5 * size(bams, 'G')) + 20
        Int memory_gb = 8
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        samtools merge -f -p -c --threads ~{thread_count} -b ~{write_lines(bams)} ~{outname}.bam
        samtools index ~{outname}.bam
    >>>
    output {
        File bam = "~{outname}.bam"
        File bai = "~{outname}.bam.bai"
    }
    runtime {
        preemptible: 2
        memory: memory_gb + " GB"
        cpu: thread_count
        disks: "local-disk " + disk_size + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}
