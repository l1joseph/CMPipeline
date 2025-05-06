nextflow.enable.dsl=2

// Define parameters
// Edit params.sample with real sample list
// params.sample = "/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/99.Raw_data/test_data/250127_sample_test.csv"
params.sample = "/tscc/lustre/restricted/alexandrov-ddn/users/kohjy2000/IAG_microbiome/99.Raw_data/RCC/sample.csv"
params.unmapped_bam_dir = "${projectDir}/RESULTS/UNMAPPED_BAM"
params.mapped_reads_dir = "${projectDir}/RESULTS/MAPPED_READS"
params.adapters="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/ref/known_adapters.fna"
params.hg38_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi"
params.t2t_phix_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi"
params.pangenome_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/pangenome_mmi"
params.kraken_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023"
params.metaphlan_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/metaphlan"

params.fastqc_dir = "${projectDir}/RESULTS/FASTQC"
params.samtools_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/samtools_env.yml"
params.fastp_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/fastp_env.yml"
params.minimap2_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/minimap2_env.yml"
params.multiqc_dir = "${projectDir}/RESULTS/MULTIQC"
params.multiqc_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/multiqc_env.yml"
params.krakenuniq_bracken_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/krakenUniq_bracken_env.yml"
params.krakenuniq_bracken_dir = "${projectDir}/RESULTS/BRACKEN"
params.krakentools_pack ="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/KrakenTools"
params.metaphlan4_pack ="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/MetaPhlAn-4.1.1"
params.metaphlan4_env = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/yml/metaphlan4_env.yml"
params.metaphlan4_dir = "${projectDir}/RESULTS/METAPHLAN4"

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
include { mapBothReads as mapBothReads } from './Modules/map_both_reads.nf'
include { multiqc as MCR1} from './Modules/multiqc.nf'
include { multiqc as MCR2} from './Modules/multiqc.nf'
include { multiqc as MCHG38} from './Modules/multiqc.nf'
include { multiqc as MCT2T} from './Modules/multiqc.nf'
include { Bracken } from './Modules/Bracken.nf'
include { metaphlan4 } from './Modules/metaphlan4.nf'

// Define the workflow
workflow {

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
    mapBothReads(FILTERED_UNMAPPED_READS_MULTI.whole, mmiFiles).set { MAPPED_READS }

    MAPPED_READS.multiMap { sampleID, r1Hg38, r1T2T, r1Pan, r2Hg38, r2T2T, r2Pan ->
        Hg38: tuple(r1Hg38, r2Hg38) 
        T2T: tuple(r1T2T, r2T2T) 
        PAN: tuple(r1Pan, r2Pan) 
    }
    .set { MAPPED_READS_MULTI }

    // Perforn FastQC on the mapped Reads
    // Perform FastQC on the filtered fastq files

    FASTQCHG38(MAPPED_READS_MULTI.Hg38)
    FASTQCT2T(MAPPED_READS_MULTI.T2T)
    FASTQCPANGENOME(MAPPED_READS_MULTI.PAN)

    // metaphlan4(MAPPED_READS_MULTI.PAN)  <- metaphlan pause

    // Kracken Taxonomic Classification
    Bracken(MAPPED_READS_MULTI.PAN)
   // CONSENSUS TAXA ANNOTATIONS      

}



