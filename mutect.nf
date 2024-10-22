params.contaminants = ''
params.fastq = ''
params.intervals = ''
params.limits = ''


vals_03 = params.steps.fastqc_03.parameters
vals_04 = params.steps.fastq_to_bam_04.parameters
vals_05 = params.steps.fastqc_05.parameters
vals_07 = params.steps.sam_to_fastq_07.parameters
vals_08 = params.steps.bwa_mem2_08.parameters

process fastqc_03 {
    container params.steps.fastqc_03.container
    cpus params.steps.fastqc_03.cpus
    memory params.steps.fastqc_03.memory

    input:
    tuple val(sample_id), path(reads)

    shell:
    """
    mkdir ${params.outputdir}
    echo "sample_id ${sample_id}, reads ${reads}"
    echo "container is ${params.steps.fastqc_03.container}"
    fastqc \
    --kmers !{vals_03['kmers']} \
    --min_length !{vals_03['min_length']} \
    --nogroup !{vals_03['nogroup']} \
    --outdir ${params.outputdir} \
    ${reads}
    """
    // These params currently filed under Too Hard as they require input files
    //    --contaminants !{vals_03['contaminants']} \
    //    --adapters !{vals_03['adapters']} \
    //    --limits !{vals_03['limits']} \

}

process fastq_to_bam_04 {
    container params.steps.fastq_to_bam_04.container
    cpus params.steps.fastq_to_bam_04.cpus
    memory params.steps.fastq_to_bam_04.memory

    input:
    tuple val(sample_id), path(reads)

    output:
    path 'fastq_to_bam_out.bam', emit: out_bam

    shell:
    """
    echo "fastq_to_bam_04: sample_id ${sample_id}, reads ${reads}"
    ${params.cmd.fgbio} \
    FastqToBam \
    --comment=!{vals_04['comment']} \
    --description=!{vals_04['description']} \
    --input=${reads[0]} ${reads[1]} \
    --library=!{vals_04['sample']} \
    --output=fastq_to_bam_out.bam \
    --platform=!{vals_04['platform']} \
    --platform-model=!{vals_04['platform-model']} \
    --platform-unit=!{vals_04['platform-unit']} \
    --predicted-insert-size=!{vals_04['predicted-insert-size']} \
    --read-group-id=!{vals_04['read-group-id']} \
    --run-date=!{vals_04['run-date']} \
    --sample=!{vals_04['sample']} \
    --sequencing-center=!{vals_04['sequencing-center']} \
    --sort=!{vals_04['sort']} \
    --umi-tag=!{vals_04['umi-tag']}
    """
}

process fastqc_05 {
    container params.steps.fastqc_05.container
    cpus params.steps.fastqc_05.cpus
    memory params.steps.fastqc_05.memory

    input:
    tuple val(sample_id), path(reads)

    shell:
    """
    mkdir ${params.outputdir}
    echo "sample_id ${sample_id}, reads ${reads}"
    echo "container is ${params.steps.fastqc.container}"
    fastqc \
    --outdir ${params.outputdir} \
    ${reads}
    """
}

process sort_bam_06 {
    container params.steps.sort_bam_06.container
    cpus params.steps.sort_bam_06.cpus
    memory params.steps.sort_bam_06.memory

    input:
    path in_bam

    output:
    path 'sort_bam_out.bam', emit: out_bam

    shell:
    """
    ${params.cmd.fgbio} \
    SortBam \
    --input=${in_bam} \
    --output=sort_bam_out.bam    
    """
}

process sam_to_fastq_07 {
    container params.steps.sam_to_fastq_07.container
    cpus params.steps.sam_to_fastq_07.cpus
    memory params.steps.sam_to_fastq_07.memory

    input:
    path in_bam

    output:
    path 'sam_to_fastq_out_R?.fastq.gz', emit: out_fastq

    shell:
    """
    ${params.cmd.gatk} \
    SamToFastq \
    --CLIPPING_MIN_LENGTH !{vals_07['CLIPPING_MIN_LENGTH']} \
    --FASTQ !{vals_07['FASTQ']} \
    --INCLUDE_NON_PF_READS !{vals_07['INCLUDE_NON_PF_READS']} \
    --INCLUDE_NON_PRIMARY_ALIGNMENTS !{vals_07['INCLUDE_NON_PRIMARY_ALIGNMENTS']} \
    --INPUT ${in_bam} \
    --INTERLEAVE !{vals_07['INTERLEAVE']} \
    --QUALITY !{vals_07['QUALITY']} \
    --READ1_MAX_BASES_TO_WRITE !{vals_07['READ1_MAX_BASES_TO_WRITE']} \
    --READ1_TRIM !{vals_07['READ1_TRIM']} \
    --READ2_MAX_BASES_TO_WRITE !{vals_07['READ2_MAX_BASES_TO_WRITE']} \
    --READ2_TRIM !{vals_07['READ2_TRIM']} \
    --RE_REVERSE !{vals_07['RE_REVERSE']} \
    --RG_TAG !{vals_07['RG_TAG']} \
    --SECOND_END_FASTQ !{vals_07['SECOND_END_FASTQ']}
    """
    // Parameters requiring valid values
    //    --CLIPPING_ACTION !{vals_07['CLIPPING_ACTION']} \
    //    --CLIPPING_ATTRIBUTE !{vals_07['CLIPPING_ATTRIBUTE']} \
    //    --COMPRESS_OUTPUTS_PER_RG !{vals_07['COMPRESS_OUTPUTS_PER_RG']} \
    //    --OUTPUT_PER_RG !{vals_07['OUTPUT_PER_RG']} \

}

process bwa_mem2_08 {
    container params.steps.bwa_mem2_08.container
    cpus params.steps.bwa_mem2_08.cpus
    memory params.steps.bwa_mem2_08.memory

    input:
    path in_fastq

    output:
    path 'bwa_mem2_out.sam', emit: out_sam
    
    shell:
    """
    echo "bwa_mem2_08: in_fastq ${in_fastq}"
    bwa-mem2 \
    index \
    ${in_fastq[0]} ${in_fastq[1]}

    bwa-mem2 \
    mem \
    -o bwa_mem2_out.sam \
    ${in_fastq[0]} ${in_fastq[1]}

    """
}

process merge_bam_alignment_09 {
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

process group_reads_by_umi_10 {
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

process call_molecular_consensus_reads_11 {
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

process filter_consensus_reads_12 {
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

process sam_to_fastq_13 {
    container params.steps.sam_to_fastq.container
    cpus params.steps.sam_to_fastq.cpus
    memory params.steps.sam_to_fastq.memory

    input:
    path in_bam

    output:
    path 'sam_to_fastq_out_R?.fastq.gz', emit: out_fastq

    script:
    """
    ${params.cmd.gatk} \
    SamToFastq \
    --INPUT ${in_bam} \
    --FASTQ sam_to_fastq_out_R1.fastq.gz \
    --SECOND_END_FASTQ sam_to_fastq_out_R2.fastq.gz
    """
}

process cutadapt_14 {
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
    fastqc_03(fastq_pairs_ch)
    // fastq_to_bam_04(fastq_pairs_ch)
    // fastq_to_bam_04(fastq_pairs_ch) | sam_to_fastq_07
    fastq_to_bam_04(fastq_pairs_ch) | sam_to_fastq_07 | bwa_mem2_08
    // sort_bam_06(fastq_to_bam.out_bam)
    // merge_bam_alignment_09(bwa_mem2.out_sam, sort_bam.out_bam) | group_reads_by_umi_10 \
    // | call_molecular_consensus_reads_11 | filter_consensus_reads_12 | sam_to_fastq_13 | cutadapt_14
    // | bwa_mem2 | sort_sam
    // mutect2(sort_sam.out_bam, intervals_ch)
}
