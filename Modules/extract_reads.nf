nextflow.enable.dsl=2

process extractReads {
    scratch true
    label 'extract_reads'
    publishDir("${params.unmapped_bam_dir}", mode: 'copy')
    conda "${params.samtools_env}"

    input:
    val(meta)  // Assuming 'meta' contains the necessary information like 'bam' and 'patient'

    output:
    tuple val(meta.patient), path("${meta.patient}.R1.UNMAPPED.fastq.gz"), path("${meta.patient}.R2.UNMAPPED.fastq.gz")
    
    script:
    """
    R1="${meta.patient}.R1.UNMAPPED.fastq.gz"
    R2="${meta.patient}.R2.UNMAPPED.fastq.gz"

    # Skip condition: if R1, R2 files already exist, link and exit
    if [[ -f "${params.unmapped_bam_dir}/\$R1" && -f "${params.unmapped_bam_dir}/\$R2" ]]; then
        echo "Skipping extractReads: Found existing \$R1, \$R2 in publishDir"

        # If not present in the current work directory, link or copy
        if [[ ! -f "\$R1" ]]; then
            ln -s "${params.unmapped_bam_dir}/\$R1" . 2>/dev/null || cp "${params.unmapped_bam_dir}/\$R1" .
        fi
        if [[ ! -f "\$R2" ]]; then
            ln -s "${params.unmapped_bam_dir}/\$R2" . 2>/dev/null || cp "${params.unmapped_bam_dir}/\$R2" .
        fi
        exit 0
    fi

    # Extract reads and split into R1 and R2 fastq.gz files
    samtools view -f 4 -O BAM ${meta.bam} | samtools bam2fq \
        -1 "\$R1" \
        -2 "\$R2"

    """
}