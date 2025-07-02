nextflow.enable.dsl=2

process filterReads {
    scratch true
    label 'filter_reads'
    publishDir("${params.unmapped_bam_dir}", mode: 'copy')
    conda "${params.fastp_env}"
    errorStrategy 'retry'
    maxRetries 3

    input:
    // input as tuple
    tuple val(sampleID), path(r1_fastq), path(r2_fastq)

    // Defining output as tuple (sampleID, r1_filtered, r2_filtered)
    output:
    tuple val(sampleID), path("${sampleID}.R1.UNMAPPED.FASTP.FILTERED.fastq.gz"), path("${sampleID}.R2.UNMAPPED.FASTP.FILTERED.fastq.gz")

    script:
    """

    R1="${r1_fastq.baseName.split('\\.')[0].toString()}.R1.UNMAPPED.FASTP.FILTERED.fastq.gz"
    R2="${r2_fastq.baseName.split('\\.')[0].toString()}.R2.UNMAPPED.FASTP.FILTERED.fastq.gz"

    if [[ -f "${params.unmapped_bam_dir}/\$R1" && -f "${params.unmapped_bam_dir}/\$R2" ]]; then
        echo "Skipping filterReads: Found existing \$R1, \$R2"
        if [[ ! -f "\$R1" ]]; then
            ln -s "${params.unmapped_bam_dir}/\$R1" . 2>/dev/null || cp "${params.unmapped_bam_dir}/\$R1" .
        fi
        if [[ ! -f "\$R2" ]]; then
            ln -s "${params.unmapped_bam_dir}/\$R2" . 2>/dev/null || cp "${params.unmapped_bam_dir}/\$R2" .
        fi
        exit 0
    fi

    fastp -l 45 --adapter_fasta ${params.adapters} --cut_tail \
        -i ${r1_fastq.toString()} -w 16 -o "\$R1"

    fastp -l 45 --adapter_fasta ${params.adapters} --cut_tail \
        -i ${r2_fastq.toString()} -w 16 -o "\$R2"
    """
}

