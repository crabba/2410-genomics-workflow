params.fastq = ''
params.contaminants = ''
params.limits = ''

process fastqc {
    container params.steps.fastqc.container
    input:
    tuple val(sample_id), path(reads)

    script:
    """
    mkdir ${params.outputdir}
    echo "sample_id ${sample_id}, reads ${reads}"
    echo "container is ${params.steps.fastqc.container}"
    fastqc \
    --outdir ${params.outputdir} \
    ${reads}
    """
}

process fastq_to_bam {
    container params.steps.fastq_to_bam.container
    input:
    tuple val(sample_id), path(reads)

    script:
    """
    echo "infiles are ${sample_id}, ${reads}"
    java -jar /root/fgbio-2.3.0.jar \
    FastqToBam \
    --input=${reads} \
    --output=out.bam \
    --sample mysample \
    --library mylibrary
    """
}

workflow {
    def fastq_pairs_ch = Channel.fromFilePairs(params.fastq)
    fastqc(fastq_pairs_ch)
    fastq_to_bam(fastq_pairs_ch)
}
