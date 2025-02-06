nextflow.enable.dsl=2

process metaphlan4 {

    label 'process_medium'
    scratch true
    publishDir "${params.metaphlan4_dir}", mode: 'copy'
    conda "${params.metaphlan4_env}"  // Set conda environment

    input:
    val(r1_fastq)
    val(r2_fastq)
 
    output:
    path "*.metagenome.bowtie2.bz2"
    path "*.sam.bz2"
    path "*.profiled_metagenome.txt", emit: metagenome_file 

    script:

    // Define the sample name from the input file name
    def SAMPLE_NAME = r1_fastq.baseName.split('\\.')[0]
    """

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
