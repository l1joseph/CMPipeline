process decontamination {
    tag "$dataset_type"
    publishDir "${params.decontam_dir}/${dataset_type}", mode: 'copy'
    
    conda params.decontam_env
    
    input:
    tuple val(dataset_type), path(otu_table), path(metadata)
    
    output:
    tuple val(dataset_type), path("${dataset_type}.decontam_pkg_decontaminated.csv"), emit: decontam_table
    path("decontamination_plots/*"), emit: plots
    path("${dataset_type}_decontam.log"), emit: log
    
    script:
    """
    # Create output directory for plots
    mkdir -p decontamination_plots
    
    # Run decontamination
    Rscript ${params.scripts}/decontamination.R \\
        --otu_table ${otu_table} \\
        --metadata ${metadata} \\
        --prefix ${dataset_type} \\
        --threshold ${params.decontam_threshold} \\
        --min_prevalence ${params.decontam_min_prevalence} \\
        --var_batch ${params.decontam_batch_var} \\
        > ${dataset_type}_decontam.log 2>&1
    
    # Check if the script completed successfully
    if [ ! -f "${dataset_type}.decontam_pkg_decontaminated.csv" ]; then
        echo "ERROR: Decontamination failed to produce output file"
        exit 1
    fi
    """
    
    stub:
    """
    # Create dummy output for testing
    touch ${dataset_type}.decontam_pkg_decontaminated.csv
    mkdir -p decontamination_plots
    touch decontamination_plots/test_plot.png
    echo "Stub mode: decontamination" > ${dataset_type}_decontam.log
    """
}
