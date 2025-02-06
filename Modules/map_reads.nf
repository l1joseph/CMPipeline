process mapReads {

    scratch true
    label 'process_medium'
    publishDir("${params.mapped_reads_dir}", mode: 'copy')
    conda "${params.minimap2_env}"

    input:
    path(fastq_file)
    val(mmi_files)  // pass the list of .mmi files
    val(read_type)  // accept read_type as an input (e.g., R1 or R2)

    output:
    path "${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz", emit: hg38_fastq
    path "${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz", emit: t2t_fastq
    path "${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz", emit: pangenome_fastq

    script:
    """

    # Run minimap2 on hg38 reference
    minimap2 -2 -ax sr -t 16 ${params.hg38_db} ${fastq_file} -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz

    # Run minimap2 on t2t_phix reference
    minimap2 -2 -ax sr -t 16 ${params.t2t_phix_db} ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz
    
    # Create temporary file for processing
    cp ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz
    
    # Convert the mmi_files list to a space-separated string
    mmi_files_list=\$(echo ${mmi_files.join(' ')})
    # Iterate over pangenome database files (list passed from main workflow)
    for mmi in \${mmi_files_list}; do
        echo "Running minimap2 on \${mmi}"
        minimap2 -2 -ax sr -t 16 \${mmi} ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.new.fastq.gz

        # Move output to temporary file for next iteration
        mv ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.new.fastq.gz ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz
    done

    # Final move to save the processed output
    mv ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.tmp.fastq.gz ${fastq_file.baseName.split('\\.')[0]}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz

    """    
}


