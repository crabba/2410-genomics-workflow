params.fastq = ''
// params.fastq2 = 'foo_2.fastq.gz'
params.outputdir = 'out'
params.contaminants = ''
params.limits = ''

container_fastqc = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

process fastqc {
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    input:
        path infile

    script:
    """
    mkdir ${params.outputdir}
    echo "infile is ${infile}"
    fastqc \
    --outdir ${params.outputdir} \
    ${infile}
    """
}

workflow {
    def fastq_ch = Channel.fromPath(params.fastq)
//    def fastq_ch = Channel.fromPath("/home/ec2-user/H06HDADXX130110.1.ATCACGAT.20k_reads_*.fastq.gz")
    fastqc(fastq_ch)
}
