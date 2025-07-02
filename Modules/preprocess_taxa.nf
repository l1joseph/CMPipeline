nextflow.enable.dsl=2

process process_metaphlan {

    scratch true
    publishDir "${params.metaphlan4_dir}", mode: 'copy'

    input:
    path profiled_files  // Input directory containing *.profiled_metagenome.txt files

    output:
    tuple path("merged_abundance_table.txt"), 
    path("merged_abundance_table_genus.txt"), 
    path("merged_abundance_table_species.txt"),
    path("merged_abundance_table_SGB.txt")

    script:

    // Prepare space-separated file name lists
    def metaphlan_files = profiled_files.findAll { it.name.endsWith('.profiled_metagenome.txt') }*.getName()
    def metaphlan_str = metaphlan_files.join(' ')


    """
    
    # Run the MetaPhlAn script to merge tables
    python3 "${params.scripts}/merge_metaphlan_tables.py" ${metaphlan_str} > merged_abundance_table.txt
    sed -i "s/.profiled_metagenome//g" merged_abundance_table.txt


    # Genus level
    (head -n 2 merged_abundance_table.txt && \
     grep -E "g__" merged_abundance_table.txt | grep -v "s__") | \
     sed "s/^.*|//g" | sed "s/.profiled_metagenome//g" > merged_abundance_table_genus.txt

    # Species level
    (head -n 2 merged_abundance_table.txt && \
     grep -E "s__" merged_abundance_table.txt | grep -v "t__") | \
     sed "s/^.*|//g" | sed "s/.profiled_metagenome//g" > merged_abundance_table_species.txt

    # SGB level
    (head -n 2 merged_abundance_table.txt && \
     grep -E "t__" merged_abundance_table.txt) | \
     sed "s/^.*|//g" | sed "s/.profiled_metagenome//g" > merged_abundance_table_SGB.txt

    """
}


process process_bracken {

    scratch true
    publishDir "${params.krakenuniq_bracken_dir}", mode: 'copy'

    input:
    path bracken_files

    output:
    tuple path("bracken.genus.mpa.report.txt"),
          path("bracken.species.mpa.report.txt")

    script:
    // Prepare space-separated file name lists for genus and species
    def genus_files = bracken_files.findAll { it.name.endsWith('.G.mpa.krakenreport.txt') }*.getName()
    def species_files = bracken_files.findAll { it.name.endsWith('.S.mpa.krakenreport.txt') }*.getName()

    def genus_str = genus_files.join(' ')
    def species_str = species_files.join(' ')

    """

    # Combine genus-level files
    if [ -n "${genus_str}" ]; then
        python "${params.scripts}"/combine_mpa.py --input ${genus_str} --output bracken.genus.mpa.report.txt
    else
        echo "No genus files found." > bracken.genus.mpa.report.txt
    fi

    # Combine species-level files
    if [ -n "${species_str}" ]; then
        python "${params.scripts}"/combine_mpa.py --input ${species_str} --output bracken.species.mpa.report.txt
    else
        echo "No species files found." > bracken.species.mpa.report.txt
    fi
    """    
}


process consensus_taxa {

    scratch true
    publishDir "${params.consensus_taxa_dir}", mode: 'copy'
    conda "${params.consensus_taxa_env}"

    input:
    tuple path(metaphlan_file), path(bracken_genus_file), path(bracken_species_file)

    output:
    path "bracken.metaphlan.taxa.prop.pdf"
    path "bracken.metaphlan.common.genus.mpa.report.txt"
    path "bracken.metaphlan.common.species.mpa.report.txt"

    script:
    """
    python "${params.scripts}"/compute_consensus_taxa.py --metaphlan ${metaphlan_file} \\
                                     --bracken_genus ${bracken_genus_file} \\
                                     --bracken_species ${bracken_species_file} \\
                                     --output bracken.metaphlan.taxa.prop.pdf \\
                                     --output_common_genus bracken.metaphlan.common.genus.mpa.report.txt \\
                                     --output_common_species bracken.metaphlan.common.species.mpa.report.txt
    """
}
