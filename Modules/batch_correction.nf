process batch_correction {
    tag "$dataset_type"
    publishDir "${params.batch_correction_dir}/${dataset_type}", mode: 'copy'
    
    conda params.batch_correction_env
    
    input:
    tuple val(dataset_type), path(otu_table), path(metadata)
    
    output:
    tuple val(dataset_type), path("corrected/ConQuR_default.tsv"), emit: corrected_table
    path("normalized/*"), emit: normalized_tables
    path("pcoa_plots/*"), emit: plots
    path("permanova_summary.tsv"), emit: permanova_summary
    path("${dataset_type}_batch_correction.log"), emit: log
    
    script:
    """
    # Create output directories
    mkdir -p corrected normalized pcoa_plots
    
    # Run batch correction and normalization
    Rscript ${params.scripts}/batch_correction_normalization.R \\
        --otu ${otu_table} \\
        --meta ${metadata} \\
        --prefix ${dataset_type} \\
        --batch_var ${params.batch_var} \\
        --covariates "${params.covariates}" \\
        ${params.tumor_only ? '--tumor_only' : ''} \\
        > ${dataset_type}_batch_correction.log 2>&1
    
    # Check if the script completed successfully
    if [ ! -f "permanova_summary.tsv" ]; then
        echo "ERROR: Batch correction failed to produce permanova summary"
        exit 1
    fi
    """
    
    stub:
    """
    # Create dummy outputs for testing
    mkdir -p corrected normalized pcoa_plots
    touch corrected/ConQuR_default.tsv
    touch normalized/test_clr_pseudo.tsv
    touch pcoa_plots/test.png
    touch permanova_summary.tsv
    echo "Stub mode: batch_correction" > ${dataset_type}_batch_correction.log
    """
}
