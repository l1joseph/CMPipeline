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
    // Define the sample name from the input file name
    def SAMPLE_NAME = fastq_file.baseName.split('\\.')[0]

    """
    out1="${SAMPLE_NAME}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz"
    out2="${SAMPLE_NAME}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"
    out3="${SAMPLE_NAME}.${read_type}.UNMAPPED.FASTP.FILTERED.hg38.t2t.pangenome.fastq.gz"



    # Skip condition
    if [[ -f "${params.mapped_reads_dir}/\$out1" && -f "${params.mapped_reads_dir}/\$out2" && -f "${params.mapped_reads_dir}/\$out3" ]]; then
        echo "Skipping mapReads: Found \$out1, \$out2, \$out3 in publishDir"

        for f in "\$out1" "\$out2" "\$out3"; do
            if [[ ! -f "\$f" ]]; then
                ln -s "${params.mapped_reads_dir}/\$f" . 2>/dev/null || cp "${params.mapped_reads_dir}/\$f" .
            fi
        done
        exit 0
    fi
    
    # ========= Actual execution =========

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


