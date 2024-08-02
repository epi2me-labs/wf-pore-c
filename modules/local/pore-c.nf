process digest_align_annotate {
    label 'pore_c_py'
    errorStrategy = 'retry'
    memory { 15.GB * task.attempt }
    maxRetries 1
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    cpus params.threads
    input:
        tuple val(meta), 
            path("concatemers.bam"),
            path("concatemers.bam.bci"),
            val(chunk_index), val(chunk_ref), path("reference.fasta.mmi"),
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
        tuple val(meta),
            path("filtered_reads.txt"),
            emit: filtered_read_ids, optional: true
    script:
        args = task.ext.args ?: " "
        if (params.chromunity) {
            args += "--chromunity "
            if (params.chromunity_merge_distance != null) {
                args += "--chromunity_merge_distance ${params.chromunity_merge_distance} "
            }
        }
        if (params.paired_end | params.bed) {
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
            echo "${chunk_ref}"
            bamindex fetch --chunk=${chunk_index} "concatemers.bam" |
                pore-c-py digest "${meta.cutter}" --max_monomers ${params.max_monomers} --excluded_list "filtered_reads.txt" \
                --header "concatemers.bam" \
                --threads ${digest_annotate_threads} |
            samtools fastq --threads 1 -T '*' |
            minimap2 -ay -t ${ubam_map_threads} ${minimap2_settings} --cap-kalloc 100m --cap-sw-mem 50m \
                "reference.fasta.mmi" - |
            pore-c-py annotate - "${meta.alias}" --monomers \
                --threads ${digest_annotate_threads}  --stdout ${args} | \
            tee "${meta.alias}_out.ns.bam" |
            samtools sort -m 1G --threads ${samtools_threads}  -u --write-index -o "${meta.alias}.cs.bam" -  
            """  
        }else{
            """
            pore-c-py digest "concatemers.bam" "${meta.cutter}" --max_monomers ${params.max_monomers} --excluded_list "filtered_reads.txt" \
                --header "concatemers.bam" \
                --threads ${digest_annotate_threads} | 
            samtools fastq --threads 1 -T '*' |
            minimap2 -ay -t ${ubam_map_threads} ${minimap2_settings} --cap-kalloc 100m --cap-sw-mem 50m \
                "reference.fasta.mmi" - |
            pore-c-py annotate - "${meta.alias}" --monomers \
                --threads ${digest_annotate_threads}  --stdout ${args} | \
            tee "${meta.alias}_out.ns.bam" |
            samtools sort -m 1G --threads ${samtools_threads}  -u --write-index -o "${meta.alias}.cs.bam" -  
            """  
        }
        
}

process haplotagReads {
    label 'wfporec'
    cpus 2
    memory "15 GB"
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
            path("${meta.alias}.ht.bam"),
            path("${meta.alias}.ht.bam.csi"),
            emit: "cs_bam"
        tuple val(meta),
            path("${meta.alias}.ht.txt.gz"),
            emit: "haplotagged_monomers"
    shell:
        args = task.ext.args ?: "--ignore-read-groups --skip-missing-contigs "
    """
    whatshap haplotag --reference "reference.fasta"  -o "${meta.alias}.ht.bam" \
    --output-haplotag-list "${meta.alias}.ht.txt.gz" $args "$phased_vcf" "concatemers.cs.bam"
    samtools index -c "${meta.alias}.ht.bam"
    """
}

/// gather individual parquets into a single directory
process merge_parquets_to_dataset {
    label 'wfporec'
    cpus 2
    memory "4 GB"
    input:
    tuple val(meta),
          path("to_merge/part?????.parquet")
    output:
        tuple val(meta),
            path("$prefix"),
            emit: "parquets"
    shell:
        prefix = task.ext.prefix ?: "${meta.alias}.chromunity.parquet"
    """
    mkdir $prefix
    cp to_merge/part*.parquet  $prefix/
    """
}
