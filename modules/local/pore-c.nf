process digest_align_annotate {
    errorStrategy = 'retry'
    maxRetries 3
    cpus params.threads
    memory "32 GB"
    label 'wfporec'
    input:
        tuple val(meta), path("concatemers.bam"),
            path("concatemers.bam.bci"),
            val(ref), path("reference.fasta.mmi"),
            val(minimap2_settings)
    output:
        tuple val(meta),
            path("${meta.alias}_out.ns.bam"),
            emit: ns_bam
        tuple val(meta),
            path("${meta.alias}.cs.bam"),
            path("${meta.alias}.cs.bam.csi"),
            emit: cs_bam
        tuple val(meta),
            path("${meta.alias}.chromunity.parquet"),
            emit: chromunity_pq, optional: true
        tuple val(meta),
            path("${meta.alias}.pe.bam"),
            emit: paired_end_bam, optional: true
    script:
        args = task.ext.args ?: " "
        if (params.chromunity) {
            args += "--chromunity "
            if (params.chromunity_merge_distance != null) {
                args += "--chromunity_merge_distance ${params.chromunity_merge_distance} "
            }
        }
        if (params.paired_end) {
            args += "--paired_end "
            if (params.filter_pairs) {
                args += "--filter_pairs "
                if (params.paired_end_minimum_distance != null) {
                    args += "--paired_end_minimum_distance ${params.paired_end_minimum_distance} "
                }
                if (params.paired_end_maximum_distance != null) {
                    args += "--paired_end_maximum_distance ${params.paired_end_maximum_distance} "
                }
            }
        }
        if (params.summary_json) {
            args  += "--summary "
        }
        def chunk = task.index - 1
        // 2 threads are recommended for each the pore-c-py processes
        def digest_annotate_threads = params.threads >= 8 ? 2 : 1 
        // if possible use 3 for samtools (--threads 2 + 1)
        def samtools_threads = params.threads >= 8 ? 2 : 1
        // calculate the left over threads for mapping and leave one as samtools will require 3
        def ubam_map_threads = params.threads - (digest_annotate_threads * 2) - samtools_threads - 1
        if (params.chunk_size > 0){
            """
            echo "${ref}"
            bamindex fetch --chunk=${chunk} "concatemers.bam" |
            pore-c-py digest "${meta.cutter}" --header "concatemers.bam" \
            --threads ${digest_annotate_threads} |
            samtools fastq --threads 1 -T '*' |
            minimap2 -ay -t ${ubam_map_threads} ${minimap2_settings} \
            "reference.fasta.mmi" - |
            pore-c-py annotate - "${meta.alias}" --monomers \
            --threads ${digest_annotate_threads}  --stdout true ${args} | \
            tee "${meta.alias}_out.ns.bam" |
            samtools sort -m 1G --threads ${samtools_threads}  -u --write-index -o "${meta.alias}.cs.bam" -  
            """  
        }else{
            """
            pore-c-py digest "concatemers.bam" "${meta.cutter}" --header "concatemers.bam" \
            --threads ${digest_annotate_threads} | 
            samtools fastq --threads 1 -T '*' |
            minimap2 -ay -t ${ubam_map_threads} ${minimap2_settings} \
            "reference.fasta.mmi" - |
            pore-c-py annotate - "${meta.alias}" --monomers \
            --threads ${digest_annotate_threads}  --stdout true ${args} | \
            tee "${meta.alias}_out.ns.bam" |
            samtools sort -m 1G --threads ${samtools_threads}  -u --write-index -o "${meta.alias}.cs.bam" -  
            """  
        }
        
}

process haplotagReads {
    label 'wfporec'
    cpus 2
    memory "16 GB"
    input:
        tuple val(meta),
            path("concatemers.cs.bam"),
            path("concatemers.cs.bam.csi"),
            path("reference.fasta"),
            path("reference.fasta.fai"),
            path(phased_vcf),
            path(phased_vcf_tbi)
    output:
        tuple val(meta),
            path("${meta.sample_id}.ht.bam"),
            path("${meta.sample_id}.ht.bam.csi"),
            emit: "cs_bam"
        tuple val(meta),
            path("${meta.sample_id}.ht.txt.gz"),
            emit: "haplotagged_monomers"
    shell:
        args = task.ext.args ?: "--ignore-read-groups --skip-missing-contigs "
    """
    whatshap haplotag --reference "reference.fasta"  -o "${meta.sample_id}.ht.bam" \
    --output-haplotag-list "${meta.sample_id}.ht.txt.gz" $args "$phased_vcf" "concatemers.cs.bam"
    samtools index -c "${meta.sample_id}.ht.bam"
    """
}

/// gather individual parquets into a single directory
process merge_parquets_to_dataset {
    label 'wfporec'
    cpus 2
    memory "2 GB"
    input:
    tuple val(meta),
          path("to_merge/part?????.parquet")
    output:
        tuple val(meta),
            path("$prefix"),
            emit: "parquets"
    shell:
        prefix = task.ext.prefix ?: "${meta.sample_id}.chromunity.parquet"
    """
    mkdir $prefix
    cp to_merge/part*.parquet  $prefix/
    """
}
