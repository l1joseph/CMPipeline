nextflow.enable.dsl=2

process humann3 {
    label 'process_high'
    scratch true
    publishDir "${params.humann3_dir}", mode: 'copy'
    conda "${params.humann3_env}"  // Set conda environment

    input:
    tuple path(r1_fastq), path(r2_fastq)
 
    output:
    path "*.log"
    path "*.metaphlan_profile.tsv"
    path "*.genefamilies.tsv"
    path "*.pathabundance.tsv"


    script:

    // Define the sample name from the input file name
    def SAMPLE_NAME = r1_fastq.baseName.split('\\.')[0]

    """
    echo "[DEBUG] PATH is: \$PATH"
    echo "[DEBUG] which humann3:"
    which metaphlan || echo "humann3 not found"
    humann --version || true

    out1="${SAMPLE_NAME}.genefamilies.tsv"
    out2="${SAMPLE_NAME}.pathabundance.tsv"

    # Skip condition
    if [[ -f "${params.humann3_dir}/\$out1" && -f "${params.humann3_dir}/\$out2" ]]; then
        echo "Skipping humann3: Found \$out1, \$out2"
        for f in "\$out1" "\$out2"; do
            if [[ ! -f "\$f" ]]; then
                ln -s "${params.humann3_dir}/\$f" . 2>/dev/null || cp "${params.humann3_dir}/\$f" .
            fi
        done
        exit 0
    fi

    # Actual execution
    # Combine R1 and R2 fastq files into a single file
    cat "${r1_fastq}" "${r2_fastq}" > "${SAMPLE_NAME}.trimmed.fastq.gz"

    gunzip "${SAMPLE_NAME}.trimmed.fastq.gz"

    # Run humann3
    humann --input "${SAMPLE_NAME}.trimmed.fastq" \\
        --search-mode 'uniref90' \\
        --nucleotide-database "${params.humann3_nucleotide_db}" \\
        --protein-database "${params.humann3_protein_db}" \\
         –metaphlan-options "–bowtie2db /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/metaphlan -x mpa_vJun23_CHOCOPhlAnSGB_202307" \\
        --o-log ${SAMPLE_NAME}.log \\
        --output-basename "${SAMPLE_NAME}" \\
        --output "${params.humann3_dir}" \\
        --threads 4 \\
        --input-format fastq \\
        --taxonomic-profile "${SAMPLE_NAME}.profiled_metagenome.txt"
    """
}