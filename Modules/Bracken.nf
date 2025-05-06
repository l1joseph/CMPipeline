nextflow.enable.dsl=2

process Bracken {

    scratch true
    label 'process_high_disk'
    publishDir "${params.krakenuniq_bracken_dir}", mode: 'copy'
    conda "${params.krakenuniq_bracken_env}"

    input:
    tuple path(r1_fastq), path(r2_fastq)

    output:
    path "*.krakenuniq.report.txt"
    path "*.classified.fasta"
    path "*.unclassified.fasta"
    path "*.bracken.*.report.txt"
    path "*.bracken.*.krakenreport.txt"


    script:

    // Define the sample name from the input file name
    def SAMPLE_NAME = r1_fastq.baseName.split('\\.')[0]

    """

    KRAKEN_REPORT="${SAMPLE_NAME}.krakenuniq.report.txt"
    CLASSIFIED="${SAMPLE_NAME}.classified.fasta"
    UNCLASSIFIED="${SAMPLE_NAME}.unclassified.fasta"
    OUTPUT="${SAMPLE_NAME}.krakenuniq.output.txt"
    BOUT_G="${SAMPLE_NAME}.bracken.G.report.txt"
    BOUT_S="${SAMPLE_NAME}.bracken.S.report.txt"

    # Skip condition
    if [[ -f "${params.krakenuniq_bracken_dir}/\$KRAKEN_REPORT" && -f "${params.krakenuniq_bracken_dir}/\$CLASSIFIED" && -f "${params.krakenuniq_bracken_dir}/\$UNCLASSIFIED" ]]; then
        # Also check if .bracken.G.report.txt and .bracken.S.report.txt exist
        if [[ -f "${params.krakenuniq_bracken_dir}/\$BOUT_G" && -f "${params.krakenuniq_bracken_dir}/\$BOUT_S" ]]; then
            echo "Skipping Bracken: Found existing files for ${SAMPLE_NAME}"

            # Link/Copy
            for f in "\$KRAKEN_REPORT" "\$CLASSIFIED" "\$UNCLASSIFIED" "${SAMPLE_NAME}.bracken.G.report.txt" "${SAMPLE_NAME}.bracken.S.report.txt" "${SAMPLE_NAME}.bracken.G.krakenreport.txt" "${SAMPLE_NAME}.bracken.S.krakenreport.txt"; do
                [ -f "${params.krakenuniq_bracken_dir}/\$f" ] || continue
                if [[ ! -f "\$f" ]]; then
                    ln -s "${params.krakenuniq_bracken_dir}/\$f" . 2>/dev/null || cp "${params.krakenuniq_bracken_dir}/\$f" .
                fi
            done
            exit 0
        fi
    fi

    # Actual execution

    echo "SAMPLE_NAME is: $SAMPLE_NAME"
    cat "${r1_fastq}" "${r2_fastq}" > "${SAMPLE_NAME}.trimmed.fastq.gz"
    gunzip "${SAMPLE_NAME}.trimmed.fastq.gz"

    CLASSIFIED_FASTA="${SAMPLE_NAME}.classified.fasta"
    UNCLASSIFIED_FASTA="${SAMPLE_NAME}.unclassified.fasta"
    REPORT="${SAMPLE_NAME}.krakenuniq.report.txt"
    OUTPUT="${SAMPLE_NAME}.krakenuniq.output.txt"

    krakenuniq --db "${params.kraken_db}" --threads 4 --report-file \${REPORT} --output \${OUTPUT} \
        --classified-out \${CLASSIFIED_FASTA} --unclassified-out \${UNCLASSIFIED_FASTA} "${SAMPLE_NAME}.trimmed.fastq"

    for lvl in G S; do
        bracken_output="${SAMPLE_NAME}.bracken.\${lvl}.report.txt"
        bracken_kraken_report="${SAMPLE_NAME}.bracken.\${lvl}.krakenreport.txt"
        bracken -d "${params.kraken_db}" -i \${REPORT} -o \${bracken_output} -w \${bracken_kraken_report} -r 50 -l \${lvl} -t 2
    done

    echo "Done processing ${SAMPLE_NAME}"

    """
}