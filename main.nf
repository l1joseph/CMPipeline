nextflow.enable.dsl=2

// Define parameters
params.sample = "/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/sample.csv"
params.unmapped_bam_dir = "${projectDir}/RESULTS/UNMAPPED_BAM"
params.mapped_reads_dir = "${projectDir}/RESULTS/MAPPED_READS"
params.adapters="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/ref/known_adapters.fna"
params.hg38_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi"
params.t2t_phix_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi"
params.pangenome_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/pangenome_mmi"

params.fastqc_dir = "${projectDir}/RESULTS/FASTQC"
params.samtools_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/samtools_env.yml"
params.fastp_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/fastp_env.yml"
params.minimap2_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/minimap2_env.yml"
params.multiqc_dir = "${projectDir}/RESULTS/MULTIQC"
params.multiqc_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/multiqc_env.yml"

// Include the external processes
include { extractReads } from './Modules/extract_reads.nf'
include { FASTQC as FASTQC1 } from './Modules/fastqc.nf'
include { FASTQC as FASTQC2 } from './Modules/fastqc.nf'
include { FASTQC as FASTQCHG38 } from './Modules/fastqc.nf'
include { FASTQC as FASTQCT2T } from './Modules/fastqc.nf'
include { FASTQC as FASTQCPANGENOME } from './Modules/fastqc.nf'
include { filterReads } from './Modules/filter_reads.nf'
include { mapReads as mapReadsR1 } from './Modules/map_reads.nf'
include { mapReads as mapReadsR2 } from './Modules/map_reads.nf'
include { multiqc } from './Modules/multiqc.nf'


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

    // Perform FastQC on the filtered fastq files
    mapReadsR1(FILTERED_UNMAPPED_READS.r1_fastq,mmiFiles).set { R1_MAPPED}

    // Perform FastQC on the filtered fastq files
    mapReadsR2(FILTERED_UNMAPPED_READS.r2_fastq,mmiFiles).set { R2_MAPPED}

    // Perform FastQC on the filtered fastq files
    FASTQCHG38(R1_MAPPED.hg38_fastq, R2_MAPPED.hg38_fastq)

    // Perform FastQC on the filtered fastq files
    FASTQCT2T(R1_MAPPED.t2t_fastq, R2_MAPPED.t2t_fastq)

    // Perform FastQC on the filtered fastq files
    FASTQCPANGENOME(R1_MAPPED.pangenome_fastq, R2_MAPPED.pangenome_fastq)

    // gather the list of pangenome .mmi files
    def input_files = []
    def zip_dir = new File("${params.fastqc_dir}")
    zip_dir.eachFileRecurse (groovy.io.FileType.FILES) { file ->
        if (file.name.endsWith('fastqc.zip')) {
            input_files << file
        }
    }

    // multiqc
    multiqc(input_files)

}



