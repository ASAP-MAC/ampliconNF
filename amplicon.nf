#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.raw_dir = "./raw"
params.fastq_files = "${params.raw_dir}/*_{1,2}.fastq.gz"
params.results_dir = "./results"
params.paired = true
params.fwd_primer = "FWDPRIMER"
params.rev_primer = "REVPRIMER"

process fastQC {
    
    tag "FASTQC on ${sample_id}"
    
    publishDir "${params.results_dir}/FastQC"
    
    input:
    tuple val(sample_id), path(fastq_in)
    
    output:
    // path("fastqc_${sample_id}_logs")// into fastqc_ch
    // tuple val(sample_id),
    // path("${sample_id}_{1,2}_fastqc.html"),
    // path("${sample_id}_{1,2}_fastqc.zip")// , 
    // path("${sample_id}_{1,2}.fastqc.gz")
    path("${sample_id}_{1,2}_fastqc.{html,zip}")

    script:
    """
    fastqc ${fastq_in}
    """
}

process multiQC {

    publishDir "${params.results_dir}/MultiQC"

    input:
    path(fastqc_ch)

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc .
    """
}

/*
process cutadapt {

    tag "Cutadapt on ${sample_id}"

    publishDir "${params.results_dir}/Cutadapt"

    input:
    tuple val(sample_id), path(fastq_in)

    output:
    path("${sample_id}*.fastq.gz")

    script:
    """
    if ${params.paired}; then
        cutadapt \
            -a ${params.fwd_primer} \
            -A ${params.rev_primer} \
            -o ${sample_id}.1.fastq.gz \
            -p ${sample_id}.2.fastq.gz \
            ${fastq_in}
    else
        cutadapt \
            -a ${params.fwd_primer} \
            -o ${sample_id}.fastq.gz \
            ${fastq_in}
    fi
    """
}
*/

process DADA2 {

    publishDir "${params.results_dir}/DADA2"

    input:
    tuple val(sample_id), path(fastq_in)

    output:
    path("${sample_id}*")

    script:
    """
    #!/usr/bin/env Rscript

    #print("testing using R")
    #print(paste0("Input:", "${sample_id}", ", ", "${fastq_in}"))

    library(dada2)

    reads <- unlist(strsplit("${fastq_in}", split = " "))
    #print(reads)
    
    if ("${params.paired}" == "true") {
        out <- filterAndTrim(fwd = file.path(reads[1]),
                     filt = paste0("${sample_id}", ".R1.filtered.fastq.gz"),
                     rev = file.path(reads[2]),
                     filt.rev = paste0("${sample_id}", ".R2.filtered.fastq.gz"),
                     minLen = 50,
                     maxN = 0,
                     maxEE = 2,
                     truncQ = 2,
                     trimRight = 0,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE,
                     verbose = TRUE)       
    } else {
        out <- filterAndTrim(fwd = file.path(reads[1]),
                     filt = paste0("${sample_id}", ".R1.filtered.fastq.gz"),
                     minLen = 50,
                     maxN = 0,
                     maxEE = 2,
                     truncQ = 2,
                     trimRight = 0,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE,
                     verbose = TRUE)   
    }

    write.csv(out, paste0("${sample_id}", ".trimmed.txt"))
    
    """
}

workflow {
    // file_channel = channel.fromPath(params.fastq_files, checkIfExists: true)
    //                     .map {it -> [it.simpleName, it]}
    // file_channel.view()

    reads_ch = channel.fromFilePairs(params.fastq_files, checkIfExists: true)

    // Block 1: Initial Quality Control Reports
    fastqc_ch = fastQC(reads_ch)
    multiQC(fastqc_ch.collect())

    // Block 2: Primer Removal
    // cutadapt(reads_ch)

    // Block 3: Quality Trimming and Filtering
    DADA2(reads_ch)

    // Block 4: Final Quality Control Reports

    // Block 5: Phylogenetic Tree Generation

    // Block 6: Packaging of Results

}
