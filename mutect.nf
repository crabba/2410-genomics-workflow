nextflow.enable.dsl = 2

params.contaminants = ''
params.fastq = ''
params.intervals = ''
params.limits = ''

ecr_prefix = params.onOmics ? params.omics_ecr_prefix : ""
publish_dir = params.onOmics ? "/mnt/workflow/pubdir" : "."

vals_03 = params.steps.fastqc_03.parameters
vals_04 = params.steps.fastq_to_bam_04.parameters
vals_05 = params.steps.fastqc_05.parameters
vals_06 = params.steps.sort_bam_06.parameters
vals_07 = params.steps.sam_to_fastq_07.parameters
vals_08 = params.steps.bwa_mem2_08.parameters
vals_09 = params.steps.merge_bam_alignment_09.parameters
vals_10 = params.steps.group_reads_by_umi_10.parameters

process fastqc_03 {
    container { ecr_prefix + params_steps_fastqc_03_container }
    cpus params_steps_fastqc_03_cpus
    memory params_steps_fastqc_03_memory
    publishDir publish_dir

    input:
    tuple val(sample_id), path(reads)

    shell:
    """
    mkdir ${params.outputdir}
    fastqc \
    --kmers 7 \
    --min_length null \
    --nogroup false \
    --outdir ${params.outputdir} \
    ${reads}
    """
    // These params currently filed under Too Hard as they require input files
    //    --contaminants !{vals_03['contaminants']} \
    //    --adapters !{vals_03['adapters']} \
    //    --limits !{vals_03['limits']} \

}

process fastq_to_bam_04 {
    container { ecr_prefix + params_steps_fastq_to_bam_04_container }
    cpus params_steps_fastq_to_bam_04_cpus
    memory params_steps_fastq_to_bam_04_memory
    publishDir publish_dir

    input:
    tuple val(sample_id), path(reads)

    output:
    path 'fastq_to_bam_out.bam', emit: bam

    shell:
    """
    echo "fastq_to_bam_04: sample_id ${sample_id}, reads ${reads}"
    ${params.cmd.fgbio} \
    FastqToBam \
    --comment=sample-comment \
    --description=sample-description \
    --input=${reads[0]} ${reads[1]} \
    --library=b \
    --output=fastq_to_bam_out.bam \
    --platform=sample-platform \
    --platform-model=sample-platform-model \
    --platform-unit=sample-platform-unit \
    --predicted-insert-size=111 \
    --read-group-id=c \
    --run-date=2024-10-21 \
    --sample=a \
    --sequencing-center=sample-sequencing-center \
    --sort=false \
    --umi-tag=", "
    """
}

process fastqc_05 {
    container { ecr_prefix + params_steps_fastqc_05_container }
    cpus params_steps_fastqc_05_cpus
    memory params_steps_fastqc_05_memory
    publishDir publish_dir

    input:
    tuple val(sample_id), path(reads)

    shell:
    """
    mkdir ${params.outputdir}
    echo "sample_id ${sample_id}, reads ${reads}"
    echo "container is ${params_steps_fastqc_container}"
    fastqc \
    --outdir ${params.outputdir} \
    ${reads}
    """
}

process sort_bam_06 {
    container { ecr_prefix + params_steps_sort_bam_06_container }
    cpus params_steps_sort_bam_06_cpus
    memory params_steps_sort_bam_06_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'sort_bam_out.bam', emit: bam

    shell:
    """
    java -jar /root/fgbio-2.3.0.jar \
    SortBam \
    --input=${in_bam} \
    --max-records-in-ram=${vals_06['max-records-in-ram']} \
    --output=sort_bam_out.bam \
    --sort-order=${vals_06['sort-order']}
    """
}

process sam_to_fastq_07 {
    container { ecr_prefix + params_steps_sam_to_fastq_07_container }
    cpus params_steps_sam_to_fastq_07_cpus
    memory params_steps_sam_to_fastq_07_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'sam_to_fastq_out_R?.fastq.gz', emit: fastq

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
    container { ecr_prefix + params_steps_bwa_mem2_08_container }
    cpus params_steps_bwa_mem2_08_cpus
    memory params_steps_bwa_mem2_08_memory
    publishDir publish_dir

    input:
    path in_fastq
    path reference_sequence
    path reference_index
    path reference_bwt
    path reference_ann
    path reference_amb
    path reference_pac
    path reference_0123

    output:
    path 'bwa_mem2_out.sam', emit: sam
    
    shell:
    '''
    ref_stem=`basename !{reference_sequence}`
    echo "bwa_mem2_08: in_fastq !{in_fastq}, ref_stem $ref_stem"

    bwa-mem2 \
    mem \
    -o bwa_mem2_out.sam \
    !{reference_sequence} \
    !{in_fastq[0]} !{in_fastq[1]}
    '''
    // bwa-mem2 \
    // index \
    // -p ${ref_stem} \
    // !{reference_sequence}
}

process merge_bam_alignment_09 {
    container { ecr_prefix + params_steps_merge_bam_alignment_09_container }
    cpus params_steps_merge_bam_alignment_09_cpus
    memory params_steps_merge_bam_alignment_09_memory
    publishDir publish_dir

    input:
    path aligned_bam
    path unmapped_bam
    path reference_sequence
    path reference_index
    path reference_dict

    output:
    path 'merge_bam_alignment_out.bam', emit: out_bam
    
    shell:
    """
    echo "merge_bam_alignment_09 reference_sequence ${reference_sequence}"
    java -jar /gatk/gatk.jar \
    MergeBamAlignment \
    --ALIGNED_BAM ${aligned_bam} \
    --OUTPUT merge_bam_alignment_out.bam \
    --REFERENCE_SEQUENCE ${reference_sequence} \
    --UNMAPPED_BAM ${unmapped_bam}
    """
}

process group_reads_by_umi_10 {
    container { ecr_prefix + params_steps_group_reads_by_umi_10_container }
    cpus params_steps_group_reads_by_umi_10_cpus
    memory params_steps_group_reads_by_umi_10_memory

    input:
    path in_bam

    output:
    path 'group_reads_by_umi_out.bam', emit: out_bam
    
    shell:
    '''
    !{params.cmd.fgbio} \
    GroupReadsByUmi \
    --input=!{in_bam} \
    --output=group_reads_by_umi_out.bam \
    --strategy=!{vals_10['strategy']} \
    '''
}

process call_molecular_consensus_reads_11 {
    container { ecr_prefix + params_steps_call_molecular_consensus_reads_container }
    cpus params_steps_call_molecular_consensus_reads_cpus
    memory params_steps_call_molecular_consensus_reads_memory

    input:
    path in_bam

    output:
    path 'call_molecular_consensus_reads_out.bam', emit: out_bam
    
    script:
    """
    """
}

process filter_consensus_reads_12 {
    container { ecr_prefix + params_steps_filter_consensus_reads_container }
    cpus params_steps_filter_consensus_reads_cpus
    memory params_steps_filter_consensus_reads_memory

    input:
    path in_bam

    output:
    path 'filter_consensus_reads_out.bam', emit: out_bam
    
    script:
    """
    """
}

process sam_to_fastq_13 {
    container { ecr_prefix + params_steps_sam_to_fastq_container }
    cpus params_steps_sam_to_fastq_cpus
    memory params_steps_sam_to_fastq_memory

    input:
    path in_bam

    output:
    path 'sam_to_fastq_out_R?.fastq.gz', emit: out_fastq

    script:
    """
    java -jar /gatk/gatk.jar \
    SamToFastq \
    --INPUT ${in_bam} \
    --FASTQ sam_to_fastq_out_R1.fastq.gz \
    --SECOND_END_FASTQ sam_to_fastq_out_R2.fastq.gz
    """
}

process cutadapt_14 {
    container { ecr_prefix + params_steps_cutadapt_container }
    cpus params_steps_cutadapt_cpus
    memory params_steps_cutadapt_memory

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
    def ref_ch = Channel.fromPath('Homo_sapiens_assembly38.fasta')
    def ref_idx_ch = Channel.fromPath('Homo_sapiens_assembly38.fasta.fai')
    def ref_idx_bwt = Channel.fromPath('Homo_sapiens_assembly38.fasta.bwt.2bit.64')
    def ref_dict_ch = Channel.fromPath('Homo_sapiens_assembly38.dict')
    def ref_ann = Channel.fromPath('Homo_sapiens_assembly38.fasta.ann')
    def ref_amb = Channel.fromPath('Homo_sapiens_assembly38.fasta.amb')
    def ref_pac = Channel.fromPath('Homo_sapiens_assembly38.fasta.pac')
    def ref_0123 = Channel.fromPath('Homo_sapiens_assembly38.fasta.0123')
    

    fastqc_03(fastq_pairs_ch)
    fastq_to_bam_04(fastq_pairs_ch) | sam_to_fastq_07
    bwa_mem2_08(sam_to_fastq_07.out.fastq, ref_ch, ref_idx_ch, ref_idx_bwt, ref_ann, ref_amb, ref_pac, ref_0123)
    sort_bam_06(fastq_to_bam_04.out.bam)
    merge_bam_alignment_09(bwa_mem2_08.out.sam, sort_bam_06.out.bam, ref_ch, ref_idx_ch, ref_dict_ch)

    // merge_bam_alignment_09(bwa_mem2.out_sam, sort_bam.out_bam) | group_reads_by_umi_10 \
    // | call_molecular_consensus_reads_11 | filter_consensus_reads_12 | sam_to_fastq_13 | cutadapt_14
    // | bwa_mem2 | sort_sam
    // mutect2(sort_sam.out_bam, intervals_ch)
}
