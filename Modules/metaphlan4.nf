nextflow.enable.dsl=2

process metaphlan4 {
    label 'process_high'
    scratch true
    publishDir "${params.metaphlan4_dir}", mode: 'copy'
    conda "${params.metaphlan4_env}"  // Set conda environment

    input:
    tuple path(r1_fastq), path(r2_fastq)
 
    output:
    path "*.metagenome.bowtie2.bz2"
    path "*.sam.bz2"
    path "*.profiled_metagenome.txt"

    script:

    // Define the sample name from the input file name
    def SAMPLE_NAME = r1_fastq.baseName.split('\\.')[0]

    """
    echo "[DEBUG] PATH is: \$PATH"
    echo "[DEBUG] which metaphlan:"
    which metaphlan || echo "metaphlan not found"
    metaphlan --version || true

    out1="${SAMPLE_NAME}.metagenome.bowtie2.bz2"
    out2="${SAMPLE_NAME}.sam.bz2"
    out3="${SAMPLE_NAME}.profiled_metagenome.txt"

    # Skip condition
    if [[ -f "${params.metaphlan4_dir}/\$out1" && -f "${params.metaphlan4_dir}/\$out2" && -f "${params.metaphlan4_dir}/\$out3" ]]; then
        echo "Skipping metaphlan4: Found \$out1, \$out2, \$out3"
        for f in "\$out1" "\$out2" "\$out3"; do
            if [[ ! -f "\$f" ]]; then
                ln -s "${params.metaphlan4_dir}/\$f" . 2>/dev/null || cp "${params.metaphlan4_dir}/\$f" .
            fi
        done
        exit 0
    fi

    # Actual execution
    # Combine R1 and R2 fastq files into a single file
    cat "${r1_fastq}" "${r2_fastq}" > "${SAMPLE_NAME}.trimmed.fastq.gz"

    gunzip "${SAMPLE_NAME}.trimmed.fastq.gz"

    # Run MetaPhlAn4
    metaphlan "${SAMPLE_NAME}.trimmed.fastq" \\
        --bowtie2db "${params.metaphlan_db}" \\
        --index mpa_vJun23_CHOCOPhlAnSGB_202307 \\
        --bowtie2out "${SAMPLE_NAME}.metagenome.bowtie2.bz2" \\
        -s "${SAMPLE_NAME}.sam.bz2" \\
        --nproc 4 \\
        --input_type fastq \\
        -o "${SAMPLE_NAME}.profiled_metagenome.txt"
    """
}
