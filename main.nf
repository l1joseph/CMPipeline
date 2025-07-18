nextflow.enable.dsl=2

// Edit this with your sample.csv path and metadata
params.sample = "/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/sample.csv"
params.meta = ""

// Parameteres for preprocessing. Edit this according to analytical purpose
/// Decontamination
params.decontam_threshold = 0.1
params.decontam_min_prevalence = 0.05
params.decontam_batch_var = "shipment_batch"  // Default batch variable, can be overridden
/// Batch correction & Normalization
params.batch_var = "shipment_batch"  // Batch variable for correction
params.covariates = "age_diag,sex"   // Comma-separated list of covariates
params.tumor_only = true              // Default: process tumor samples only

// Output directories
params.unmapped_bam_dir = "${projectDir}/RESULTS/UNMAPPED_BAM"
params.mapped_reads_dir = "${projectDir}/RESULTS/MAPPED_READS"
params.fastqc_dir = "${projectDir}/RESULTS/FASTQC"
params.multiqc_dir = "${projectDir}/RESULTS/MULTIQC"
params.krakenuniq_bracken_dir = "${projectDir}/RESULTS/BRACKEN"
params.metaphlan4_dir = "${projectDir}/RESULTS/METAPHLAN4"
params.humann3_dir = "${projectDir}/RESULTS/HUMANN3"
params.consensus_taxa_dir = "${projectDir}/RESULTS/CONSENSUS_TAXA"
params.decontam_dir = "${projectDir}/RESULTS/DECONTAM"
params.batch_correction_dir = "${projectDir}/RESULTS/BATCH_CORRECTION"

// Databases and ref files [CHANGE THIS]
params.hg38_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi"
params.t2t_phix_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi"
params.pangenome_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/pangenome_mmi"
params.kraken_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023"
params.metaphlan_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/metaphlan"
params.humann3_nucleotide_db='/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/databases/humann3/chocophlan/'
params.humann3_protein_db='/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/microbiome/databases/humann3/uniref/'
params.adapters="${projectDir}/ref/known_adapters.fna"

// Enviroment paths
params.samtools_env = "./conda_envs/samtools_env.yml"
params.fastp_env = "./conda_envs/fastp_env.yml"
params.minimap2_env = "./conda_envs/minimap2_env.yml"
params.multiqc_env = "./conda_envs/multiqc_env.yml"
params.krakenuniq_bracken_env = "./conda_envs/krakenUniq_bracken_env.yml"
params.consensus_taxa_env = "./conda_envs/consensus_taxa_env.yml"
params.metaphlan4_env = "./conda_envs/metaphlan4_env.yml"
params.humann3_env = "./conda_envs/humann3_env.yml"
params.krakentools_pack ="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/KrakenTools"
params.metaphlan4_pack ="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/MetaPhlAn-4.1.1"
params.decontam_env = "./conda_envs/decontam_env.yml"
params.batch_correction_env = "./conda_envs/batch_correction_env.yml"

// Package and script paths
params.scripts ="${projectDir}/scripts"


// Include the external processes
include { extractReads } from './Modules/extract_reads.nf'
include { FASTQC as FASTQC1 } from './Modules/fastqc.nf'
include { FASTQC as FASTQC2 } from './Modules/fastqc.nf'
include { FASTQC as FASTQCHG38 } from './Modules/fastqc.nf'
include { FASTQC as FASTQCT2T } from './Modules/fastqc.nf'
include { FASTQC as FASTQCPANGENOME } from './Modules/fastqc.nf'
include { filterReads } from './Modules/filter_reads.nf'
include { mapReads as mapReads } from './Modules/map_reads.nf'
include { multiqc as MCR1} from './Modules/multiqc.nf'
include { multiqc as MCR2} from './Modules/multiqc.nf'
include { multiqc as MCHG38} from './Modules/multiqc.nf'
include { multiqc as MCT2T} from './Modules/multiqc.nf'
include { Bracken } from './Modules/Bracken.nf'
include { metaphlan4 } from './Modules/metaphlan4.nf'
include { process_metaphlan; process_bracken; consensus_taxa } from './Modules/preprocess_taxa.nf'
include { humann3 } from './Modules/humann3.nf'
include { decontamination } from './Modules/decontamination.nf'
include { batch_correction } from './Modules/batch_correction.nf'


// Define the workflow
workflow {

    // ------------------- STEP1: HOST DEPLETION ---------------------- //

    // Read and parse the sample sheet
    sample_sheet = nextflow.Channel.fromPath(params.sample)
        .splitCsv(header: true)
        .map { row ->
            row.subMap('patient', 'bam') // Extract relevant metadata
        }

    // Extract reads from BAM files
    extractReads(sample_sheet).set { UNMAPPED_READS }

    UNMAPPED_READS.multiMap { sampleID, r1, r2 -> 
        path_only: tuple(r1, r2)
        whole: tuple(sampleID, r1, r2)
    }
    .set { UNMAPPED_READS_MULTI }

    // Perform FastQC on the extracted fastq files
    FASTQC1(UNMAPPED_READS_MULTI.path_only)    
    
    // Filter poor quality reads using fastp
    filterReads(UNMAPPED_READS_MULTI.whole).set { FILTERED_UNMAPPED_READS }

    FILTERED_UNMAPPED_READS.multiMap { sampleID, r1, r2 ->
        path_only: tuple(r1, r2)
        whole: tuple(sampleID, r1, r2)
    }
    .set { FILTERED_UNMAPPED_READS_MULTI }

   // Perform FastQC on the filtered fastq files 
    FASTQC2(FILTERED_UNMAPPED_READS_MULTI.path_only)

    // gather the list of pangenome .mmi files
    def mmiFiles = []
    def dir = new File("${params.pangenome_db}")
    dir.eachFileRecurse (groovy.io.FileType.FILES) { file ->
        if (file.name.endsWith('.mmi')) {
            mmiFiles << file
        }
    }

    // Processing for both READSs
    mapReads(FILTERED_UNMAPPED_READS_MULTI.whole, mmiFiles).set { MAPPED_READS }

    MAPPED_READS.multiMap { sampleID, r1Hg38, r1T2T, r1Pan, r2Hg38, r2T2T, r2Pan ->
        Hg38: tuple(r1Hg38, r2Hg38) 
        T2T: tuple(r1T2T, r2T2T) 
        PAN: tuple(r1Pan, r2Pan) 
    }
    .set { MAPPED_READS_MULTI }

    // Perforn FastQC on the mapped Reads
    FASTQCHG38(MAPPED_READS_MULTI.Hg38)
    FASTQCT2T(MAPPED_READS_MULTI.T2T)
    FASTQCPANGENOME(MAPPED_READS_MULTI.PAN)

    // ------------------- STEP2: TAXONOMIC CLASSIFICATION -- BRACKEN ---------------------- //

    Bracken(MAPPED_READS_MULTI.PAN).set { BRACKEN_OUT }

    // Collect all kraken reports once all samples are done
    BRACKEN_OUT.map { kraken_report, classified_fasta, unclassified_fasta, bracken_reports, bracken_krakenreports, bracken_mpa_reports -> 
        tuple(bracken_mpa_reports)
    }
    .flatten()   
    .collect()
    .set {Bracken_mpa_files}

    // Run process_bracken once all reports are available
    process_bracken(Bracken_mpa_files).set { BRACKEN_FILES }

    // ------------------- STEP2: TAXONOMIC CLASSIFICATION -- METAPHLAN ---------------------- //

    metaphlan4(MAPPED_READS_MULTI.PAN).set { METAPHLAN_OUT }

    // Collect all kraken reports once all samples are done
    METAPHLAN_OUT.map { bowtie2_files, sam_files, profiled_metagenomes -> 
        tuple(profiled_metagenomes)
    }
    .flatten()   
    .collect()
    .set {metaphlan4_files}

    process_metaphlan(metaphlan4_files).set { METAPHLAN_FILES }


    // -------------------  OPTIONAL STEP3: CONSENSUS TAXA ---------------------- //
    
    // Step 1: Get merged table
    METAPHLAN_FILES
        .map { merged_table, merged_genus, merged_species, merged_SGB -> merged_genus }
        .set { metaphlan4_merged_table }

    // Step 2: Get bracken tuple
    BRACKEN_FILES
        .map { bracken_genus_file, bracken_species_file -> tuple(bracken_genus_file, bracken_species_file) }
        .set { bracken_file_tuple }

    // Step 3: Combine for consensus
    metaphlan4_merged_table
        .combine(bracken_file_tuple)
        .set { consensus_input }

    consensus_taxa(consensus_input).set { CONSENSUS_OUTPUT }  # set CONSENSUS_OUTPUT at here


    // ------------------- OPTIONAL STEP4: DECONTAMINATION ---------------------- //

    // Step 0 : Check if metadata file is provided
    if (!params.meta) {
        exit 1, "Decontamination step requires valid path for --meta parameter"
    }

    // Step 1 : Prepare metadata channel
    metadata_ch = Channel.fromPath(params.meta)

    // Step 2: Extract consensus output files
    CONSENSUS_OUTPUT
        .map { prop_pdf, common_genus, common_species -> 
            [
                tuple("consensus_genus", common_genus),
                tuple("consensus_species", common_species)
            ]
        }
        .flatten()
        .combine(metadata_ch)
        .set { decontam_input }

    // Step 3 : Run decontamination on consensus output
    decontamination(decontam_input).set { DECONTAM_RESULTS }

    // ------------------- OPTIONAL STEP5: BATCH CORRECTION ---------------------- //

    // Step 1 : Prepare input for batch correction
    DECONTAM_RESULTS
        .map { dataset_type, decontam_table, plots, log -> 
            tuple(dataset_type, decontam_table)
        }
        .combine(metadata_ch)
        .set { batch_correction_input }

    // Step 2 : Run batch correction and normalization
    batch_correction(batch_correction_input).set { BATCH_CORRECTED_RESULTS }

}



