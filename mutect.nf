#! /usr/bin/env nextflow

nextflow.enable_dsl = 2

params.fastq1 = 'foo_1.fastq.gz'
params.fastq2 = 'foo_2.fastq.gz'

container_fastqc = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

process fastqc {
    container container_fastqc
    input:
        path infile
    output:
        path outfile

    """
    fastqc --version
    """
}

workflow {
    def input_ch = Channel.fromPath(params.fastq1)
    fastqc(input_ch)
}