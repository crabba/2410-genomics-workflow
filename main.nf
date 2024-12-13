nextflow.enable.dsl = 2

params.contaminants = ''
params.fastq = ''
params.intervals = ''
params.limits = ''

ecr_prefix = params.onOmics ? params.omics_ecr_prefix : ""
publish_dir = params.onOmics ? "/mnt/workflow/pubdir" : "."

process fastqc_03 {
    container { ecr_prefix + params.steps_03_fastqc_container }
    cpus params.steps_03_fastqc_cpus
    memory params.steps_03_fastqc_memory
    publishDir publish_dir

    input:
    // tuple val(sample_id), path(reads)
    path fastq_1
    path fastq_2

    output:
    

    shell:
    '''
    mkdir !{params.outputdir}
    fastqc \
    --kmers 7 \
    --min_length null \
    --nogroup false \
    --outdir !{params.outputdir} \
    !{fastq_1} !{fastq_2}
    '''
    // These params currently filed under Too Hard as they require input files
    //    --contaminants 
    //    --adapters 
    //    --limits 
    // !{reads}

}

process fastq_to_bam_04 {
    container { ecr_prefix + params.steps_04_fastq_to_bam_container }
    cpus params.steps_04_fastq_to_bam_cpus
    memory params.steps_04_fastq_to_bam_memory
    publishDir publish_dir

    input:
    // tuple val(sample_id), path(reads)
    path fastq_1
    path fastq_2

    output:
    path 'fastq_to_bam_out.bam', emit: bam

    shell:
    '''
    java -jar /root/fgbio-2.3.0.jar \
    FastqToBam \
    --comment="sample-comment" \
    --description="sample-description" \
    --input=!{fastq_1} !{fastq_2} \
    --library=b \
    --output="fastq_to_bam_out.bam" \
    --platform="sample-platform" \
    --platform-model="sample-platform-model" \
    --platform-unit="sample-platform-unit" \
    --predicted-insert-size=111 \
    --read-group-id="c" \
    --read-structures=12M11S+T +T \
    --run-date=2024-10-21 \
    --sample="a" \
    --sequencing-center="sample-sequencing-center" \
    --sort=false \
    --umi-tag="RX"
    '''
    //     --input=!{reads[0]} !{reads[1]} \
}

process fastqc_05 {
    container { ecr_prefix + params.steps_05_fastqc_container }
    cpus params.steps_05_fastqc_cpus
    memory params.steps_05_fastqc_memory
    publishDir publish_dir

    input:
    tuple val(sample_id), path(reads)

    shell:
    '''
    mkdir !{params.outputdir}
    fastqc \
    --outdir !{params.outputdir} \
    !{reads}
    '''
}

process sort_bam_06 {
    container { ecr_prefix + params.steps_06_sort_bam_container }
    cpus params.steps_06_sort_bam_cpus
    memory params.steps_06_sort_bam_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'sort_bam_out.bam', emit: bam

    shell:
    '''
    java -jar /root/fgbio-2.3.0.jar \
    SortBam \
    --input=!{in_bam} \
    --max-records-in-ram=111111 \
    --output="sort_bam_out.bam" \
    --sort-order="Queryname"
    '''
}

process sam_to_fastq_07 {
    container { ecr_prefix + params.steps_07_sam_to_fastq_container }
    cpus params.steps_07_sam_to_fastq_cpus
    memory params.steps_07_sam_to_fastq_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'sam_to_fastq_out_R?.fastq.gz', emit: fastq

    shell:
    '''
    java -jar /gatk/gatk.jar \
    SamToFastq \
    --CLIPPING_MIN_LENGTH "0" \
    --FASTQ "sam_to_fastq_out_R1.fastq.gz" \
    --INCLUDE_NON_PF_READS false \
    --INCLUDE_NON_PRIMARY_ALIGNMENTS false \
    --INPUT !{in_bam} \
    --INTERLEAVE false \
    --QUALITY null \
    --READ1_MAX_BASES_TO_WRITE null \
    --READ1_TRIM "0" \
    --READ2_MAX_BASES_TO_WRITE null \
    --READ2_TRIM "0" \
    --RE_REVERSE true \
    --RG_TAG "PU" \
    --SECOND_END_FASTQ "sam_to_fastq_out_R2.fastq.gz"
    '''
    // Parameters requiring valid values
    //    --CLIPPING_ACTION 
    //    --CLIPPING_ATTRIBUTE 
    //    --COMPRESS_OUTPUTS_PER_RG 
    //    --OUTPUT_PER_RG 

}

process bwa_mem2_08 {
    container { ecr_prefix + params.steps_08_bwa_mem2_container }
    cpus params.steps_08_bwa_mem2_cpus
    memory params.steps_08_bwa_mem2_memory
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
    -t !{params.steps_08_bwa_mem2_threads} \
    !{reference_sequence} \
    !{in_fastq[0]} !{in_fastq[1]}
    '''
    // bwa-mem2 \
    // index \
    // -p ${ref_stem} \
    // !{reference_sequence}
}

process merge_bam_alignment_09 {
    container { ecr_prefix + params.steps_09_merge_bam_alignment_container }
    cpus params.steps_09_merge_bam_alignment_cpus
    memory params.steps_09_merge_bam_alignment_memory
    publishDir publish_dir

    input:
    path aligned_bam
    path unmapped_bam
    path reference_sequence
    path reference_index
    path reference_dict

    output:
    path 'merge_bam_alignment_out.bam', emit: bam
    
    shell:
    '''
    java -jar /gatk/gatk.jar \
    MergeBamAlignment \
    --ADD_MATE_CIGAR !{params.steps_09_merge_bam_alignment_add_mate_cigar} \
    --ALIGNED_READS_ONLY !{params.steps_09_merge_bam_alignment_aligned_reads_only} \
    --ALIGNER_PROPER_PAIR_FLAGS !{params.steps_09_merge_bam_alignment_aligner_proper_pair_flags} \
    --PRIMARY_ALIGNMENT_STRATEGY !{params.steps_09_merge_bam_alignment_primary_alignment_strategy} \
    --SORT_ORDER !{params.steps_09_merge_bam_alignment_sort_order} \
    --ALIGNED_BAM !{aligned_bam} \
    --OUTPUT merge_bam_alignment_out.bam \
    --REFERENCE_SEQUENCE !{reference_sequence} \
    --UNMAPPED_BAM !{unmapped_bam}
    '''
}

process group_reads_by_umi_10 {
    container { ecr_prefix + params.steps_10_group_reads_by_umi_container }
    cpus params.steps_10_group_reads_by_umi_cpus
    memory params.steps_10_group_reads_by_umi_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'group_reads_by_umi_out.bam', emit: bam
    
    shell:
    '''
    java -jar /root/fgbio-2.3.0.jar \
    GroupReadsByUmi \
    --edits=1 \
    --include-non-pf-reads="false" \
    --include-secondary="false" \
    --include-supplementary="false" \
    --input=!{in_bam} \
    --mark-duplicates="false" \
    --output=group_reads_by_umi_out.bam \
    --strategy="adjacency"
    '''
    // --family_size_histogram_on="false" \
    // --assign-tag="" \
    // --raw-tag="" \
}

process call_molecular_consensus_reads_11 {
    container { ecr_prefix + params.steps_11_call_molecular_consensus_reads_container }
    cpus params.steps_11_call_molecular_consensus_reads_cpus
    memory params.steps_11_call_molecular_consensus_reads_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'call_molecular_consensus_reads_out.bam', emit: bam
    
    shell:
    '''
    java -jar /root/fgbio-2.3.0.jar \
    CallMolecularConsensusReads \
    --consensus-call-overlapping-bases="false" \
    --debug="false" \
    --input=!{in_bam} \
    --min-reads=6 \
    --output=call_molecular_consensus_reads_out.bam \
    --output-per-base-tags="true"
    '''
    // "rejects_output_on": "false" - Could this be '-r PathToBam, --rejects=PathToBam'?
    // Galaxy tool ver 2.2.1, my ver 2.3.0
    // --sort-order="" \
    // --error-rate-pre-umi=null \
    // --max-reads=null \
    // --min-input-base-quality=null \
    // --error-rate-post-umi=null \
    // --read-group-id=null \
    // --read-name-prefix=null \
    // --tag=null
}

process filter_consensus_reads_12 {
    container { ecr_prefix + params.steps_12_filter_consensus_reads_container }
    cpus params.steps_12_filter_consensus_reads_cpus
    memory params.steps_12_filter_consensus_reads_memory

    input:
    path in_bam
    path reference_sequence

    output:
    path 'filter_consensus_reads_out.bam', emit: bam
    
    shell:
    '''
    java -Xmx4096m -jar /root/fgbio-2.3.0.jar \
    FilterConsensusReads \
    --input=!{in_bam} \
    --max-base-error-rate="1.0" \
    --max-read-error-rate="1.0" \
    --min-base-quality=10 \
    --min-reads=6 \
    --output=filter_consensus_reads_out.bam \
    --ref=!{reference_sequence} \
    --require-single-strand-agreement="false" \
    --reverse-per-base-tags="false"
    '''
    // 'ref_cond' section of the GA file has "ref": "hg19"
    // --max-no-call-fraction=null \
    // --min-mean-base-quality=null \
    // --sort-order=null
}

process sam_to_fastq_13 {
    container { ecr_prefix + params.steps_13_sam_to_fastq_container }
    cpus params.steps_13_sam_to_fastq_cpus
    memory params.steps_13_sam_to_fastq_memory

    input:
    path in_bam
    path reference_sequence

    output:
    path 'sam_to_fastq_out_R?.fastq.gz', emit: fastq

    shell:
    '''
    java -jar /gatk/gatk.jar \
    SamToFastq \
    --CLIPPING_MIN_LENGTH "0" \
    --COMPRESSION_LEVEL 5 \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
    --FASTQ "sam_to_fastq_out_R1.fastq.gz" \
    --INCLUDE_NON_PF_READS false \
    --INCLUDE_NON_PRIMARY_ALIGNMENTS false \
    --INPUT !{in_bam} \
    --INTERLEAVE false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --QUALITY null \
    --QUIET false \
    --READ1_MAX_BASES_TO_WRITE null \
    --READ1_TRIM "0" \
    --READ2_MAX_BASES_TO_WRITE null \
    --READ2_TRIM "0" \
    --RE_REVERSE true \
    --REFERENCE_SEQUENCE !{reference_sequence} \
    --RG_TAG "PU" \
    --SECOND_END_FASTQ "sam_to_fastq_out_R2.fastq.gz" \
    --USE_JDK_DEFLATER false \
    --USE_JDK_INFLATER false \
    --VALIDATION_STRINGENCY "STRICT" \
    --VERBOSITY "INFO"
    '''
    //         "GA4GH_CLIENT_SECRETS": "client_secrets.json",
    //         "arguments_file": {
    //        "CLIPPING_ACTION": "",
    //        "CLIPPING_ATTRIBUTE": "",
    //        "UNPAIRED_FASTQ_sel": "false"

    // Argument 'FASTQ' cannot be used in conjunction with argument(s) OUTPUT_PER_RG COMPRESS_OUTPUTS_PER_RG
    // --OUTPUT_PER_RG false \
    // --COMPRESS_OUTPUTS_PER_RG false \

}

process cutadapt_14 {
    container { ecr_prefix + params.steps_cutadapt_container }
    cpus params.steps_cutadapt_cpus
    memory params.steps_cutadapt_memory

    input:
    path in_bam

    output:
    path 'cutadapt_out_R?.fastq.gz', emit: out_fastq
    
    shell:
    '''
    '''
}

workflow {
    // def fastq_pairs_ch = Channel.fromFilePairs(params.fastq)
    def fastq_1_ch = Channel.fromPath(params.fastq_1)
    def fastq_2_ch = Channel.fromPath(params.fastq_2)

    // Hard coded for now
    // def ref_prefix = "s3://crabba-raas-data/genomics/ref/hg38/"
    def ref_prefix = params.ref_prefix
    def ref_ch      = Channel.fromPath(ref_prefix + 'Homo_sapiens_assembly38.fasta')
    def ref_idx_ch  = Channel.fromPath(ref_prefix + 'Homo_sapiens_assembly38.fasta.fai')
    def ref_idx_bwt = Channel.fromPath(ref_prefix + 'Homo_sapiens_assembly38.fasta.bwt.2bit.64')
    def ref_dict_ch = Channel.fromPath(ref_prefix + 'Homo_sapiens_assembly38.dict')
    def ref_ann     = Channel.fromPath(ref_prefix + 'Homo_sapiens_assembly38.fasta.ann')
    def ref_amb     = Channel.fromPath(ref_prefix + 'Homo_sapiens_assembly38.fasta.amb')
    def ref_pac     = Channel.fromPath(ref_prefix + 'Homo_sapiens_assembly38.fasta.pac')
    def ref_0123    = Channel.fromPath(ref_prefix + 'Homo_sapiens_assembly38.fasta.0123')

    fastqc_03(fastq_1_ch, fastq_2_ch)
    fastq_to_bam_04(fastq_1_ch, fastq_2_ch) | sam_to_fastq_07
    bwa_mem2_08(sam_to_fastq_07.out.fastq, ref_ch, ref_idx_ch, ref_idx_bwt, ref_ann, ref_amb, ref_pac, ref_0123)
    sort_bam_06(fastq_to_bam_04.out.bam)
    merge_bam_alignment_09(bwa_mem2_08.out.sam, sort_bam_06.out.bam, ref_ch, ref_idx_ch, ref_dict_ch) | group_reads_by_umi_10 | call_molecular_consensus_reads_11
    filter_consensus_reads_12(call_molecular_consensus_reads_11.out.bam, ref_ch)
    sam_to_fastq_13(filter_consensus_reads_12.out.bam, ref_ch) | cutadapt_14 \


    // merge_bam_alignment_09(bwa_mem2.out_sam, sort_bam.out_bam) | group_reads_by_umi_10 \
    // | call_molecular_consensus_reads_11 | filter_consensus_reads_12 | sam_to_fastq_13 | cutadapt_14
    // | bwa_mem2 | sort_sam
    // mutect2(sort_sam.out_bam, intervals_ch)
}
