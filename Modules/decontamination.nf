// Modules/decontamination.nf

nextflow.enable.dsl=2

process Decontamination {
    
    scratch true
    label 'process_medium'
    conda "${params.krakenuniq_bracken_env}"
    
    label 'process_high_disk'
    publishDir "${params.krakenuniq_bracken_dir}", mode: 'copy'
    conda "${params.decontam_env}"

    publishDir "${params.decontam_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(prefix), path(otu_table), path(metadata)

    output:
    path "${prefix}.*.csv"

    script:
    def decontam_script = "${projectDir}/scripts/decontamination.R"
    """
    ${decontam_script} \\
        --otu_table ${otu_table} \\
        --metadata ${metadata} \\
        --prefix ${prefix}
    """
}