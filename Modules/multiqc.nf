nextflow.enable.dsl=2

process multiqc {

    scratch true
    label 'process_medium'
    publishDir("${params.multiqc_dir}", mode: 'copy')
    conda "${params.multiqc_env}"  // Add the conda environment if required for multiqc


    input:
    val(input_files)  // The list of input files

    output:
    path("multiqc_report.html")  // Output report for MultiQC
    path("multiqc_data.tsv")     // Output data file for MultiQC

    script:
    """
    # Run MultiQC on all input files
    files_list=\$(echo ${input_files.join(' ')})
    multiqc -v -p -ip -f --data-dir --data-format tsv --no-report --cl-config "max_table_rows: 10000" --outdir ./ \${files_list} --force
    """
}



