nextflow.enable.dsl=2

process process_metaphlan4 {

    scratch true
    publishDir "${params.metaphlan4_dir}", mode: 'copy'
    conda "${params.metaphlan4_env}"  // Use the conda environment for MetaPhlAn
    errorStrategy 'retry'

    input:
    path profiled_files  // Input directory containing *.profiled_metagenome.txt files

    output:
    path "merged_abundance_table.txt"
    path "merged_abundance_table_genus.txt"
    path "merged_abundance_table_species.txt"
    path "merged_abundance_table_SGB.txt"

    script:
    """
    # Define the directory where profiled files are located
    DIR=${profiled_files}

    # Run the MetaPhlAn script to merge tables
    python3 "${params.metaphlan_pack}/metaphlan/utils/merge_metaphlan_tables.py" \${DIR}/*.profiled_metagenome.txt > merged_abundance_table.txt

    # Genus level
    grep -E "g__|PD" merged_abundance_table.txt \\
    | grep -v "s__" \\
    | sed "s/^.*|//g" \\
    | sed "s/.profiled_metagenome//g" \\
    > merged_abundance_table_genus.txt

    # Species level
    grep -E "s__|PD" merged_abundance_table.txt \\
    | grep -v "t__" \\
    | sed "s/^.*|//g" \\
    | sed "s/.profiled_metagenome//g" \\
    > merged_abundance_table_species.txt

    # SGB level
    grep -E "t__|PD" merged_abundance_table.txt \\
    | sed "s/^.*|//g" \\
    | sed "s/.profiled_metagenome//g" \\
    > merged_abundance_table_SGB.txt
    """
}

process process_bracken {

    // Specify the resources needed
    label 'process_medium'
    conda "${params.kraken_env}" // Use the appropriate Conda environment
    errorStrategy 'terminate'

    input:
    path bracken_species_reports, type: 'file' // List of species-level Kraken reports
    path bracken_genus_reports, type: 'file'   // List of genus-level Kraken reports

    output:
    path "${params.result_path}/crc.bracken.species.mpa.report.txt", emit: species_combined_report
    path "${params.result_path}/crc.bracken.genus.mpa.report.txt", emit: genus_combined_report

    script:
    """
    # Set up result directories
    RESULT_PATH="${params.result_path}/mpa.reports"
    mkdir -p \$RESULT_PATH

    # Species-level processing
    for file in ${bracken_species_reports}; do
        echo "Processing species report: \$file"
        SAMPLE=\$(basename \$file .bracken.species.krakenreport.txt)
        python ${params.kraken_tools}/kreport2mpa.py -r \$file -o \$RESULT_PATH/\$SAMPLE.bracken.species.mpa.report.txt --display-header
    done

    # Genus-level processing
    for file in ${bracken_genus_reports}; do
        echo "Processing genus report: \$file"
        SAMPLE=\$(basename \$file .bracken.genus.krakenreport.txt)
        python ${params.kraken_tools}/kreport2mpa.py -r \$file -o \$RESULT_PATH/\$SAMPLE.bracken.genus.mpa.report.txt --display-header
    done

    # Combine genus-level MPA reports
    GENUS_FILES=\$(ls \$RESULT_PATH/*.bracken.genus.mpa.report.txt)
    python ${params.kraken_tools}/combine_mpa.py --input \${GENUS_FILES[@]} --output ${params.result_path}/crc.bracken.genus.mpa.report.txt

    # Combine species-level MPA reports
    SPECIES_FILES=\$(ls \$RESULT_PATH/*.bracken.species.mpa.report.txt)
    python ${params.kraken_tools}/combine_mpa.py --input \${SPECIES_FILES[@]} --output ${params.result_path}/crc.bracken.species.mpa.report.txt
    """
}


process consensus_taxa {


}

process remove_contaminants {

}


