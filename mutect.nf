params.contaminants = ''
params.fastq = ''
params.intervals = ''
params.limits = ''

fgbio_cmd = "java -jar /root/fgbio-2.3.0.jar"
gatk_cmd = "java -jar /gatk/gatk.jar"

process fastqc {
    container params.steps.fastqc.container
    cpus params.steps.fastqc.cpus
    memory params.steps.fastqc.memory

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
    cpus params.steps.fastq_to_bam.cpus
    memory params.steps.fastq_to_bam.memory

    input:
    tuple val(sample_id), path(reads)

    output:
    path 'fastq_to_bam_out.bam', emit: out_bam

    script:
    """
    echo "infiles are ${sample_id}, ${reads}"
    ${fgbio_cmd} \
    FastqToBam \
    --input=${reads} \
    --output=fastq_to_bam_out.bam \
    --sample mysample \
    --library mylibrary
    """
}

process sort_bam {
    container params.steps.sort_bam.container
    cpus params.steps.sort_bam.cpus
    memory params.steps.sort_bam.memory

    input:
    path in_bam

    output:
    path 'sort_bam_out.bam', emit: out_bam

    script:
    """
    ${fgbio_cmd} \
    SortBam \
    --input=${in_bam} \
    --output=sort_bam_out.bam    
    """
}

process sam_to_fastq {
    container params.steps.sam_to_fastq.container
    cpus params.steps.sam_to_fastq.cpus
    memory params.steps.sam_to_fastq.memory

    input:
    path in_bam

    output:
    path 'sam_to_fastq_out_R?.fastq.gz', emit: out_fastq

    script:
    """
    ${gatk_cmd} \
    SamToFastq \
    --INPUT ${in_bam} \
    --FASTQ sam_to_fastq_out_R1.fastq.gz \
    --SECOND_END_FASTQ sam_to_fastq_out_R2.fastq.gz
    """
}

process bwa_mem2 {
    container params.steps.bwa_mem2.container
    cpus params.steps.bwa_mem2.cpus
    memory params.steps.bwa_mem2.memory

    input:
    path in_fastq

    output:
    path 'bwa_mem2_out.sam', emit: out_sam
    
    script:
    """
    bwa-mem2 \
    index \
    ${in_fastq}

    bwa-mem2 \
    mem \
    -o bwa_mem2_out.sam \
    ${in_fastq}

    """
}

process merge_bam_alignment {
    container params.steps.merge_bam_alignment.container
    cpus params.steps.merge_bam_alignment.cpus
    memory params.steps.merge_bam_alignment.memory

    input:
    path aligned_bam
    path unmapped_bam

    output:
    path 'merge_bam_alignment_out.bam', emit: out_bam
    
    script:
    """
    """
}

process group_reads_by_umi {
    container params.steps.group_reads_by_umi.container
    cpus params.steps.group_reads_by_umi.cpus
    memory params.steps.group_reads_by_umi.memory

    input:
    path in_bam

    output:
    path 'group_reads_by_umi_out.bam', emit: out_bam
    
    script:
    """
    """
}

process call_molecular_consensus_reads {
    container params.steps.call_molecular_consensus_reads.container
    cpus params.steps.call_molecular_consensus_reads.cpus
    memory params.steps.call_molecular_consensus_reads.memory

    input:
    path in_bam

    output:
    path 'call_molecular_consensus_reads_out.bam', emit: out_bam
    
    script:
    """
    """
}

process filter_consensus_reads {
    container params.steps.filter_consensus_reads.container
    cpus params.steps.filter_consensus_reads.cpus
    memory params.steps.filter_consensus_reads.memory

    input:
    path in_bam

    output:
    path 'filter_consensus_reads_out.bam', emit: out_bam
    
    script:
    """
    """
}

process cutadapt {
    container params.steps.cutadapt.container
    cpus params.steps.cutadapt.cpus
    memory params.steps.cutadapt.memory

    input:
    path in_bam

    output:
    path 'cutadapt_out_R?.fastq.gz', emit: out_fastq
    
    script:
    """
    """
}

workflow {
    def fastq_pairs_ch = Channel.fromFilePairs(params.fastq)
    def intervals_ch = Channel.fromFile(params.intervals)
    fastqc(fastq_pairs_ch)
    fastq_to_bam(fastq_pairs_ch) | sam_to_fastq | bwa_mem2
    sort_bam(fastq_to_bam.out_bam)
    merge_bam_alignment(bwa_mem2.out_sam, sort_bam.out_bam) | group_reads_by_umi | call_molecular_consensus_reads \
    | filter_consensus_reads | sam_to_fastq | cutadapt | bwa_mem2 | sort_sam
    mutect2(sort_sam.out_bam, intervals_ch)
}
