nextflow.enable.dsl=2

process Bracken {
    
    scratch true
    label 'process_high'
    publishDir "${params.krakenuniq_bracken_dir}", mode: 'copy'
    conda "${params.krakenuniq_bracken_env}"

    input:
    val(r1_fastq)
    val(r2_fastq)

    output:
    path "*.krakenuniq.report.txt"
    path "*.classified.fasta"
    path "*.unclassified.fasta"
    path "*.bracken.*.report.txt"
    path "*.bracken.*.krakenreport.txt", emit: krakenreport 

    script:

    // Define the sample name from the input file name
    def SAMPLE_NAME = r1_fastq.baseName.split('\\.')[0]
    """
    
    echo "SAMPLE_NAME is: $SAMPLE_NAME"
    cat "${r1_fastq}" "${r2_fastq}" > "${SAMPLE_NAME}.trimmed.fastq.gz"
    gunzip "${SAMPLE_NAME}.trimmed.fastq.gz"

    CLASSIFIED_FASTA="${SAMPLE_NAME}.classified.fasta"
    UNCLASSIFIED_FASTA="${SAMPLE_NAME}.unclassified.fasta"
    REPORT="${SAMPLE_NAME}.krakenuniq.report.txt"
    OUTPUT="${SAMPLE_NAME}.krakenuniq.output.txt"

    krakenuniq --db "${params.kraken_db}" --threads 32 --report-file \${REPORT} --output \${OUTPUT} \
        --classified-out \${CLASSIFIED_FASTA} --unclassified-out \${UNCLASSIFIED_FASTA} "${SAMPLE_NAME}.trimmed.fastq"

    for lvl in G S; do
        bracken_output="${SAMPLE_NAME}.bracken.\${lvl}.report.txt"
        bracken_kraken_report="${SAMPLE_NAME}.bracken.\${lvl}.krakenreport.txt"
        bracken_kraken_mpa_report="${SAMPLE_NAME}.bracken.\${lvl}.mpa.report.txt"
        bracken -d "${params.kraken_db}" -i \${REPORT} -o \${bracken_output} -w \${bracken_kraken_report} -r 50 -l \${lvl} -t 2
        python ${params.krakentools_pack}/kreport2mpa.py -r \${bracken_kraken_report} -o \${bracken_kraken_mpa_report} --display-header
    done

    echo "Done processing ${SAMPLE_NAME}"

    """
}