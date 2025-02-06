nextflow.enable.dsl=2

// Edit this with your sample.csv path
params.sample = "/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/sample.csv"

// Output directories
params.unmapped_bam_dir = "${projectDir}/RESULTS/UNMAPPED_BAM"
params.mapped_reads_dir = "${projectDir}/RESULTS/MAPPED_READS"
params.fastqc_dir = "${projectDir}/RESULTS/FASTQC"
params.multiqc_dir = "${projectDir}/RESULTS/MULTIQC"
params.krakenuniq_bracken_dir = "${projectDir}/RESULTS/BRACKEN"
params.metaphlan4_dir = "${projectDir}/RESULTS/METAPHLAN4"
params.consensus_taxa_dir = "${projectDir}/RESULTS/CONSENSUS_TAXA"

// Databases and ref files [CHANGE THIS]
params.hg38_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi"
params.t2t_phix_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi"
params.pangenome_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/pangenome_mmi"
params.kraken_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023"
params.metaphlan_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/metaphlan"
params.adapters="./ref/known_adapters.fna"

// Enviroment paths
params.samtools_env = "./conda_envs/samtools_env.yml"
params.fastp_env = "./conda_envs/fastp_env.yml"
params.minimap2_env = "./conda_envs/minimap2_env.yml"
params.multiqc_env = "./conda_envs/multiqc_env.yml"
params.krakenuniq_bracken_env = "./conda_envs/krakenUniq_bracken_env.yml"
params.downstream_CMP_env = "./conda_envs/downstream_CMP_env.yml"
params.metaphlan4_env = "./conda_envs/metaphlan4_env.yml"

// Package and script paths
params.scripts ="./scripts"


// Include Processes
include { extractReads } from './Modules/extract_reads.nf'
include { FASTQC as FASTQC1; FASTQC as FASTQC2; FASTQC as FASTQCHG38;FASTQC as FASTQCT2T;FASTQC as FASTQCPANGENOME; } from './Modules/fastqc.nf'
include { filterReads } from './Modules/filter_reads.nf'
include { mapReads as mapReadsR1;mapReads as mapReadsR2 } from './Modules/map_reads.nf'
include { multiqc as MCR1;multiqc as MCR2;multiqc as MCHG38; multiqc as MCT2T;} from './Modules/multiqc.nf'
include { Bracken } from './Modules/Bracken.nf'
include { metaphlan4 } from './Modules/metaphlan4.nf'
include { process_metaphlan; process_bracken;consensus_taxa; } from './Modules/consensus_taxa.nf'

// Define the workflow
workflow {

    // Read and parse the sample sheet
    sample_sheet = Channel.fromPath(params.sample)
        .splitCsv(header: true)
        .map { row ->
            row.subMap('patient', 'bam') // Extract relevant metadata
        }

    // Extract reads from BAM files
    extractReads(sample_sheet).set { UNMAPPED_READS }

    // Perform FastQC on the extracted fastq files
    FASTQC1(UNMAPPED_READS.r1_fastq, UNMAPPED_READS.r2_fastq)

    // Filter poor quality reads using fastp
    filterReads(UNMAPPED_READS.r1_fastq, UNMAPPED_READS.r2_fastq).set { FILTERED_UNMAPPED_READS }

    // Perform FastQC on the filtered fastq files
    FASTQC2(FILTERED_UNMAPPED_READS.r1_fastq, FILTERED_UNMAPPED_READS.r2_fastq)

    // gather the list of pangenome .mmi files
    def mmiFiles = []
    def dir = new File("${params.pangenome_db}")
    dir.eachFileRecurse (groovy.io.FileType.FILES) { file ->
        if (file.name.endsWith('.mmi')) {
            mmiFiles << file
        }
    }

    // R1 Processing
    mapReadsR1(FILTERED_UNMAPPED_READS.r1_fastq, mmiFiles, 'R1').set { R1_MAPPED }

    // R2 Processing
    mapReadsR2(FILTERED_UNMAPPED_READS.r2_fastq,  mmiFiles, 'R2').set { R2_MAPPED }

    // Perform FastQC on the filtered fastq files
    FASTQCHG38(R1_MAPPED.hg38_fastq, R2_MAPPED.hg38_fastq)

    // Perform FastQC on the filtered fastq files
    FASTQCT2T(R1_MAPPED.t2t_fastq, R2_MAPPED.t2t_fastq)

    // Perform FastQC on the filtered fastq files
    FASTQCPANGENOME(R1_MAPPED.pangenome_fastq, R2_MAPPED.pangenome_fastq)

    // Metaphlan Taxonomic Classification
    metaphlan4(R1_MAPPED.pangenome_fastq, R2_MAPPED.pangenome_fastq)

    // Kracken Taxonomic Classification
    Bracken(R1_MAPPED.pangenome_fastq, R2_MAPPED.pangenome_fastq)

    // PROCESS TAXA ANNOTATIONS
    process_metaphlan(metaphlan4.out.metagenome_file.collect()) 
    process_bracken(Bracken.out.krakenreport.collect()) 

}



