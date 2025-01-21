nextflow.enable.dsl=2

process filterReads {
    scratch true
    label 'process_medium'
    publishDir("${params.unmapped_bam_dir}", mode: 'copy')
    conda "${params.fastp_env}"
    errorStrategy 'retry'
    maxRetries 3

    input:
    path r1_fastq
    path r2_fastq

    output:
    path("${r1_fastq.baseName.split('\\.')[0]}.R1.UNMAPPED.FASTP.FILTERED.fastq.gz"), emit: r1_fastq
    path("${r2_fastq.baseName.split('\\.')[0]}.R2.UNMAPPED.FASTP.FILTERED.fastq.gz"), emit: r2_fastq

    script:
    """
    fastp -l 45 --adapter_fasta ${params.adapters} --cut_tail -i ${r1_fastq} -w 16 -o ${r1_fastq.baseName.split('\\.')[0]}.R1.UNMAPPED.FASTP.FILTERED.fastq.gz
    fastp -l 45 --adapter_fasta ${params.adapters} --cut_tail -i ${r2_fastq} -w 16 -o ${r2_fastq.baseName.split('\\.')[0]}.R2.UNMAPPED.FASTP.FILTERED.fastq.gz
    """
}



