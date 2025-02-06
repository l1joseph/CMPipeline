nextflow.enable.dsl=2

process process_metaphlan {

    scratch true
    publishDir "${params.metaphlan4_dir}", mode: 'copy'
    conda "${params.metaphlan4_env}"  // Use the conda environment for MetaPhlAn
    errorStrategy 'retry'

    input:
    path profiled_files  // Input directory containing *.profiled_metagenome.txt files

    output:
    path "merged_abundance_table.txt", emit: metaphlan_file
    path "merged_abundance_table_genus.txt"
    path "merged_abundance_table_species.txt"
    path "merged_abundance_table_SGB.txt"

    script:
    """
    # Define the directory where profiled files are located

    # Run the MetaPhlAn script to merge tables
    python3 "${params.metaphlan4_pack}/metaphlan/utils/merge_metaphlan_tables.py" *.profiled_metagenome.txt > merged_abundance_table.txt

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

    scratch true
    publishDir "${params.krakenuniq_bracken_dir}", mode: 'copy'

    input:
    path bracken_files  // Input directory containing *.profiled_metagenome.txt files

    output:
    path "bracken.genus.mpa.report.txt", emit: bracken_genus_file
    path "bracken.species.mpa.report.txt", emit: bracken_species_file

    script:
    def bracken_files_list = bracken_files as List
    def GENUS_FILES = bracken_files_list.findAll { it.name.endsWith('.G.mpa.report.txt') }*.toString().join(' ')
    def SPECIES_FILES = bracken_files_list.findAll { it.name.endsWith('.S.mpa.report.txt') }*.toString().join(' ')

    """
    # Combine genus-level files
    if [ -n "$GENUS_FILES" ]; then
        python ${params.krakentools_pack}/combine_mpa.py --input $GENUS_FILES \\
                                                         --output bracken.genus.mpa.report.txt
    else
        echo "No genus files found." > bracken.genus.mpa.report.txt
    fi

    # Combine species-level files
    if [ -n "$SPECIES_FILES" ]; then
        python ${params.krakentools_pack}/combine_mpa.py --input $SPECIES_FILES \\
                                                         --output bracken.species.mpa.report.txt
    else
        echo "No species files found." > bracken.species.mpa.report.txt
    fi
    """   
}


process consensus_taxa {

    scratch true
    publishDir "${params.consensus_taxa_dir}", mode: 'copy'
    conda "downstream_CMP_env.yml"  // Nextflow will use this Conda environment

    input:
    path metaphlan_file
    path bracken_genus_file
    path bracken_species_file

    output:
    path "bracken.metaphlan.taxa.prop.tumor.pdf"
    path "crc.bracken.metaphlan.common.genus.mpa.report.txt"
    path "crc.bracken.metaphlan.common.species.mpa.report.txt"

    script:
    """
    python "${params.python_scripts}"/compute_consensus_taxa.py --metaphlan ${metaphlan_file} \\
                                     --bracken_genus ${bracken_genus_file} \\
                                     --bracken_species ${bracken_species_file} \\
                                     --output bracken.metaphlan.taxa.prop.tumor.pdf \\
                                     --output_common_genus crc.bracken.metaphlan.common.genus.mpa.report.txt \\
                                     --output_common_species crc.bracken.metaphlan.common.species.mpa.report.txt
    """
}
