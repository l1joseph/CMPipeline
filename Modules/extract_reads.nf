nextflow.enable.dsl=2

process extractReads {
    scratch true
    label 'extract_reads'
    publishDir("${params.unmapped_bam_dir}", mode: 'copy')
    conda "${params.samtools_env}"
    
    input:
    val(meta)  // Assuming 'meta' contains the necessary information like 'bam' and 'patient'

    output:
    path("${meta.patient}.R1.UNMAPPED.fastq.gz"), emit: r1_fastq
    path("${meta.patient}.R2.UNMAPPED.fastq.gz"), emit: r2_fastq


    script:
    """
    # Extract reads and split into R1 and R2 fastq.gz files
    samtools view -f 4 -O BAM ${meta.bam} | samtools bam2fq \
        -1 ${meta.patient}.R1.UNMAPPED.fastq.gz \
        -2 ${meta.patient}.R2.UNMAPPED.fastq.gz

    """
}