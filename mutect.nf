params.fastq = ''
params.contaminants = ''
params.limits = ''

process fastqc {
    container params.steps.fastqc.container
    input:
        path infile

    script:
    """
    mkdir ${params.outputdir}
    echo "infile is ${infile}"
    echo "container is ${params.steps.fastqc.container}"
    fastqc \
    --outdir ${params.outputdir} \
    ${infile}
    """
}

workflow {
    def fastq_ch = Channel.fromPath(params.fastq)
    fastqc(fastq_ch)
}
