process digest_align_annotate {
    errorStrategy = 'retry'
    maxRetries 3
    maxForks 10
    label 'wfporec'
    input:
        tuple val(meta), path("concatemers.bam"), path("reference.fasta.mmi"),
              val(minimap2_settings)
    output:
        tuple val(meta),
            path("${meta.sample_id}_out.ns.bam"),
            emit: ns_bam
        tuple val(meta),
            path("${meta.sample_id}.cs.bam"),
            path("${meta.sample_id}.cs.bam.csi"),
            emit: cs_bam
        tuple val(meta),
            path("${meta.sample_id}.chromunity.parquet"),
            emit: chromunity_pq, optional: true
        tuple val(meta),
            path("${meta.sample_id}.pe.bam"),
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
        if (params.summary) {
            args  += "--summary "
        }
        test_task = task.ext.suffix
    """
    pore-c-py digest "concatemers.bam" "${meta.enzyme}" \
    --threads ${params.digest_annotate_threads} |
    samtools fastq --threads 1 -T '*'  |
    minimap2 -ay -t ${params.ubam_map_threads} ${minimap2_settings} \
    "reference.fasta.mmi" - |
    pore-c-py annotate - "${meta.sample_id}" --monomers \
    --threads ${params.digest_annotate_threads}  --stdout true ${args}|
    tee "${meta.sample_id}_out.ns.bam" |
    samtools sort -m 1G --threads ${params.digest_annotate_threads}  -u --write-index -o "${meta.sample_id}.cs.bam" - 
    
    """
}

process haplotagReads {
    label 'wfporec'
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
