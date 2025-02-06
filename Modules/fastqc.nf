nextflow.enable.dsl=2

process FASTQC {
    scratch true
    label 'process_medium'
    publishDir("${params.fastqc_dir}", mode: 'copy')
    conda "${params.samtools_env}"
    errorStrategy 'retry'
    maxRetries 3
    input:
    path(r1_fastq)
    path(r2_fastq)

    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip

    script:
    """
    # Running FastQC on the extracted R1 and R2 fastq.gz files
    fastqc -t 8 -o ./  ${r1_fastq} ${r2_fastq}
    """
}

