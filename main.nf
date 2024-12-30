nextflow.enable.dsl = 2

params.contaminants = ''
params.fastq = ''
params.intervals = ''
params.limits = ''

ecr_prefix = params.onOmics ? params.omics_ecr_prefix : ""
publish_dir = params.onOmics ? "/mnt/workflow/pubdir" : "."

process step03_fastqc {
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

process step04_fastq_to_bam {
    container { ecr_prefix + params.steps_04_fastq_to_bam_container }
    cpus params.steps_04_fastq_to_bam_cpus
    memory params.steps_04_fastq_to_bam_memory
    publishDir publish_dir

    input:
    // tuple val(sample_id), path(reads)
    path fastq_1
    path fastq_2

    output:
    path 'step04_out_fastq_to_bam.bam', emit: bam

    shell:
    '''
    java -jar /root/fgbio-2.3.0.jar \
    FastqToBam \
    --comment="sample-comment" \
    --description="sample-description" \
    --input=!{fastq_1} !{fastq_2} \
    --library=b \
    --output="step04_out_fastq_to_bam.bam" \
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

process step05_fastqc {
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

process step06_sort_bam {
    container { ecr_prefix + params.steps_06_sort_bam_container }
    cpus params.steps_06_sort_bam_cpus
    memory params.steps_06_sort_bam_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'step06_out_sort_bam.bam', emit: bam

    shell:
    '''
    java -jar /root/fgbio-2.3.0.jar \
    SortBam \
    --input=!{in_bam} \
    --max-records-in-ram=111111 \
    --output="step06_out_sort_bam.bam" \
    --sort-order="Queryname"
    '''
}

process step07_sam_to_fastq {
    container { ecr_prefix + params.steps_07_sam_to_fastq_container }
    cpus params.steps_07_sam_to_fastq_cpus
    memory params.steps_07_sam_to_fastq_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'step07_out_sam_to_fastq_R?.fastq.gz', emit: fastq

    shell:
    '''
    java -jar /gatk/gatk.jar \
    SamToFastq \
    --CLIPPING_MIN_LENGTH "0" \
    --FASTQ "step07_out_sam_to_fastq_R1.fastq.gz" \
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
    --SECOND_END_FASTQ "step07_out_sam_to_fastq_R2.fastq.gz"
    '''
    // Parameters requiring valid values
    //    --CLIPPING_ACTION 
    //    --CLIPPING_ATTRIBUTE 
    //    --COMPRESS_OUTPUTS_PER_RG 
    //    --OUTPUT_PER_RG 

}

process step08_bwa_mem2 {
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
    path 'step08_out_bwa_mem2.sam', emit: sam
    
    shell:
    '''
    bwa-mem2 \
    mem \
    -o step08_out_bwa_mem2.sam \
    -t !{params.steps_08_bwa_mem2_threads} \
    !{reference_sequence} \
    !{in_fastq[0]} !{in_fastq[1]}
    '''
    // bwa-mem2 \
    // index \
    // -p ${ref_stem} \
    // !{reference_sequence}
}

process step09_merge_bam_alignment {
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
    path 'step09_out_merge_bam_alignment.bam', emit: bam
    
    shell:
    '''
    java -jar /gatk/gatk.jar \
    MergeBamAlignment \
    --ADD_MATE_CIGAR true \
    --ALIGNED_READS_ONLY false \
    --ALIGNER_PROPER_PAIR_FLAGS false \
    --PRIMARY_ALIGNMENT_STRATEGY BestMapq \
    --SORT_ORDER coordinate \
    --ALIGNED_BAM !{aligned_bam} \
    --OUTPUT step09_out_merge_bam_alignment.bam \
    --REFERENCE_SEQUENCE !{reference_sequence} \
    --UNMAPPED_BAM !{unmapped_bam}
    '''
}

process step10_group_reads_by_umi {
    container { ecr_prefix + params.steps_10_group_reads_by_umi_container }
    cpus params.steps_10_group_reads_by_umi_cpus
    memory params.steps_10_group_reads_by_umi_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'step10_out_group_reads_by_umi.bam', emit: bam
    
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
    --output=step10_out_group_reads_by_umi.bam \
    --strategy="adjacency"
    '''
    // --family_size_histogram_on="false" \
    // --assign-tag="" \
    // --raw-tag="" \
}

process step11_call_molecular_consensus_reads {
    container { ecr_prefix + params.steps_11_call_molecular_consensus_reads_container }
    cpus params.steps_11_call_molecular_consensus_reads_cpus
    memory params.steps_11_call_molecular_consensus_reads_memory
    publishDir publish_dir

    input:
    path in_bam

    output:
    path 'step11_out_call_molecular_consensus_reads.bam', emit: bam
    
    shell:
    '''
    java -jar /root/fgbio-2.3.0.jar \
    CallMolecularConsensusReads \
    --consensus-call-overlapping-bases="false" \
    --debug="false" \
    --input=!{in_bam} \
    --min-reads=2 \
    --output=step11_out_call_molecular_consensus_reads.bam \
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

process step12_filter_consensus_reads {
    container { ecr_prefix + params.steps_12_filter_consensus_reads_container }
    cpus params.steps_12_filter_consensus_reads_cpus
    memory params.steps_12_filter_consensus_reads_memory

    input:
    path in_bam
    path reference_sequence

    output:
    path 'step12_out_filter_consensus_reads.bam', emit: bam
    
    shell:
    '''
    java -Xmx4096m -jar /root/fgbio-2.3.0.jar \
    FilterConsensusReads \
    --input=!{in_bam} \
    --max-base-error-rate="1.0" \
    --max-read-error-rate="1.0" \
    --min-base-quality=10 \
    --min-reads=1 \
    --output=step12_out_filter_consensus_reads.bam \
    --ref=!{reference_sequence} \
    --require-single-strand-agreement="false" \
    --reverse-per-base-tags="false"
    '''
    // 'ref_cond' section of the GA file has "ref": "hg19"
    // --max-no-call-fraction=null \
    // --min-mean-base-quality=null \
    // --sort-order=null
    // --min-reads=6 - set to 1 for testing with smaller input files
}

process step13_sam_to_fastq {
    container { ecr_prefix + params.steps_13_sam_to_fastq_container }
    cpus params.steps_13_sam_to_fastq_cpus
    memory params.steps_13_sam_to_fastq_memory

    input:
    path in_bam
    path reference_sequence

    output:
    path 'step13_out_sam_to_fastq_R?.fastq.gz', emit: fastq

    shell:
    '''
    java -jar /gatk/gatk.jar \
    SamToFastq \
    --CLIPPING_MIN_LENGTH "0" \
    --COMPRESSION_LEVEL 5 \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
    --FASTQ "step13_out_sam_to_fastq_R1.fastq.gz" \
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
    --SECOND_END_FASTQ "step13_out_sam_to_fastq_R2.fastq.gz" \
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

process step14_cutadapt {
    container { ecr_prefix + params.steps_14_cutadapt_container }
    cpus params.steps_14_cutadapt_cpus
    memory params.steps_14_cutadapt_memory

    input:
    path in_bam

    output:
    path 'step14_cutadapt_out_R?.fastq.gz', emit: fastq
    
    shell:
    '''
    cutadapt \
    --error-rate="0.1" \
    --times=1 \
    --overlap=7 \
    --match-read-wildcards \
    --minimum-length=20 \
    --maximum-length=0 \
    --pair-filter="any" \
    -a 21nt=CAAAACGCAATACTGTACTGG \
    -g 11nt=ATGACTCGAAT \
    --output step14_cutadapt_out_R1.fastq.gz \
    --paired-output step14_cutadapt_out_R2.fastq.gz \
    !{in_bam}

    '''
    // Unimplemented:            "front_adapters2": [
    // Need correct values: -a, -g, --cut, --output, --paired-output, --match-read-wildcards
    // --no-indels="false" \
    // --no-trim="false" \
    // --mask-adapter="false" \
    // --discard-trimmed="false" \
    // --discard-untrimmed="false" \

    // --cut 11 \
    // --cut -21 \

}

process step15_bwa_mem2 {
    container { ecr_prefix + params.steps_15_bwa_mem2_container }
    cpus params.steps_15_bwa_mem2_cpus
    memory params.steps_15_bwa_mem2_memory
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
    path 'step15_out_bwa_mem2.sam', emit: sam
    
    shell:
    '''
    bwa-mem2 \
    mem \
    -w 100 \
    -L 5 \
    -W 10000 \
    -E 1 \
    -O 6 \
    -M \
    -A 1 \
    -T 30 \
    -k 19 \
    -B 4 \
    -t 16 \
    -U 9 \
    -d 100 \
    -o step15_out_bwa_mem2.sam \
    -t !{params.steps_15_bwa_mem2_threads} \
    !{reference_sequence} \
    !{in_fastq[0]} !{in_fastq[1]}
    '''
}

process step16_sortsam {
    container { ecr_prefix + params.steps_16_sortsam_container }
    cpus params.steps_16_sortsam_cpus
    memory params.steps_16_sortsam_memory

    input:
    path in_bam

    output:
    path 'step16_out_sortsam.bam', emit: bam
    path 'step16_out_sortsam.bam.bai', emit: bai

    shell:
    '''
    echo "Before SortSam"
    pwd
    ls -l
    
    java -jar /gatk/gatk.jar SortSam \
    --COMPRESSION_LEVEL 5 \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
    --INPUT !{in_bam} \
    --MAX_RECORDS_IN_RAM 2000000 \
    --OUTPUT step16_intermed_sortsam.sam \
    --QUIET false \
    --SORT_ORDER "coordinate" \
    --USE_JDK_DEFLATER false \
    --USE_JDK_INFLATER false \
    --VALIDATION_STRINGENCY "LENIENT" \
    --VERBOSITY "INFO"

    echo "Before samtools addreplacerg"
    ls -l

    samtools addreplacerg \
    -r '@RG\tID:samplename\tSM:samplename' \
    step16_intermed_sortsam.sam \
    -o step16_out_sortsam.sam

    echo "Before samtools sort"
    ls -l

    samtools sort \
    -@ !{params.steps_16_sortsam_cpus} \
    -o step16_out_sortsam.bam \
    step16_out_sortsam.sam

    echo "Before samtools index"
    ls -l

    samtools index \
    step16_out_sortsam.bam

    echo "After samtools index"
    ls -l

    '''
    // mutect2 without addreplacerg gives 'java.lang.IllegalArgumentException: samples cannot be empty'
    // Probably missing RG insertion in a previous step
    // https://www.biostars.org/p/424817/

    // mutect2_17 without samtools index gives: 
    // 'A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.'
    // 'Please index all input files'

    //         "GA4GH_CLIENT_SECRETS": "client_secrets.json",
}

process step17_mutect2 {
    container { ecr_prefix + params.steps_17_mutect2_container }
    cpus params.steps_17_mutect2_cpus
    memory params.steps_17_mutect2_memory

    input:
    path in_bam
    path in_bai
    path reference_sequence
    path reference_index
    path reference_bwt
    path reference_ann
    path reference_amb
    path reference_pac
    path reference_0123
    path reference_dict

    output:
    path 'step17_out_mutect2.vcf.gz', emit: vcf

    shell:
    '''
    echo "Before Mutect2"
    pwd
    ls -l

    java -Xmx16384M -jar /gatk/gatk.jar \
    Mutect2 \
    --input !{in_bam} \
    --output step17_out_mutect2.vcf.gz \
    --reference !{reference_sequence}
    '''
    // Filed under Too Hard
    //         "gvcf_lod_band_rpt": [
    // 	    --kmer_size: null, \
    //      --num-pruning-samples "1" \
    // --active-probability-threshold "0.002" \
    // --adaptive-pruning-initial-error-rate "0.001" \
    // --allele-informative-reads-overlap-margin "2" \
    // --allow-non-unique-kmers-in-ref "false" \
    // --bam-writer-type "CALLED_HAPLOTYPES" \
    // --disable-adaptive-pruning "false" \
    // --disable-tool-default-annotations "false" \
    // --dont-increase-kmer-sizes-for-cycles "false" \
    // --dont-use-soft-clipped-bases "false" \
    // --emit-ref-confidence "NONE" \
    // --enable-all-annotations "false" \
    // --force-active "false" \
    // --force-call-filtered-alleles "false" \
    // --independent-mates "false" \
    // --max-mnp-distance "1" \
    // --max-num-haplotypes-in-population "128" \
    // --max-prob-propagation-distance "50" \
    // --max-suspicious-reads-per-alignment-start "0" \
    // --max-unpruned-variants "100" \
    // --min-dangling-branch-length "4" \
    // --min-pruning "2" \
    // --minimum-allele-fraction "0.0" \
    // --pair-hmm-gap-continuation-penalty "10" \
    // --pair-hmm-implementation "FASTEST_AVAILABLE" \
    // --pcr-indel-model "CONSERVATIVE" \
    // --phred-scaled-global-read-mismapping-rate "45" \
    // --pruning-lod-threshold "2.302585092994046" \
    // --recover-all-dangling-branches "false" \
    // --smith-waterman "JAVA" \
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

    // Put this back in - only omitted to speed up cached test runs (this step does not cache on AHO)
    // step03_fastqc(fastq_1_ch, fastq_2_ch)
    step04_fastq_to_bam(fastq_1_ch, fastq_2_ch) | step07_sam_to_fastq
    step08_bwa_mem2(step07_sam_to_fastq.out.fastq, ref_ch, ref_idx_ch, ref_idx_bwt, ref_ann, ref_amb, ref_pac, ref_0123)
    step06_sort_bam(step04_fastq_to_bam.out.bam)
    step09_merge_bam_alignment(step08_bwa_mem2.out.sam, step06_sort_bam.out.bam, ref_ch, ref_idx_ch, ref_dict_ch) | step10_group_reads_by_umi | step11_call_molecular_consensus_reads
    step12_filter_consensus_reads(step11_call_molecular_consensus_reads.out.bam, ref_ch)
    step13_sam_to_fastq(step12_filter_consensus_reads.out.bam, ref_ch) | step14_cutadapt
    step15_bwa_mem2(step14_cutadapt.out.fastq, ref_ch, ref_idx_ch, ref_idx_bwt, ref_ann, ref_amb, ref_pac, ref_0123) | step16_sortsam
    step17_mutect2(step16_sortsam.out.bam, step16_sortsam.out.bai, ref_ch, ref_idx_ch, ref_idx_bwt, ref_ann, ref_amb, ref_pac, ref_0123, ref_dict_ch)
}
