process digest_concatemers {
    label 'wfporec'
    input:
        tuple val(meta), path("concatemers.bam")
    output:
        tuple val(meta), path("${meta.sample_id}.monomers.bam")
    script:
        args = task.ext.args ?: " "
    """
    pore-c-py digest "concatemers.bam" "${meta.cutter}" "${meta.sample_id}.monomers.bam" $args
    """
}


process minimap2_ubam_namesort {
    label 'wfporec'
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    input:
        tuple val(meta),
              path("concatermers.monomers.bam"),
              path("reference.fasta.mmi"),
              val(minimap2_settings)
    output:
        tuple val(meta), path("${meta.sample_id}.mm2.bam")
    script:
    """
    samtools fastq --threads ${params.ubam_bam2fq_threads} -T '*' "concatermers.monomers.bam" |
    minimap2 -ay -t ${params.ubam_map_threads}  ${minimap2_settings} "reference.fasta.mmi" - |
    samtools sort --threads ${params.ubam_sort_threads} -n -o "${meta.sample_id}.mm2.bam"
    """
}

process annotate_monomers {
    label 'wfporec'
    input:
        tuple val(meta), path("concatemers.mm2.bam")
    output:
        tuple val(meta),
            path("annotated/${meta.sample_id}.ns.bam"),
            emit: ns_bam
        tuple val(meta),
            path("annotated/${meta.sample_id}.cs.bam"),
            path("annotated/${meta.sample_id}.cs.bam.csi"),
            emit: cs_bam
        tuple val(meta),
            path("annotated/${meta.sample_id}.chromunity.parquet"),
            emit: chromunity_pq, optional: true
        tuple val(meta),
            path("annotated/${meta.sample_id}.pe.bam"),
            emit: paired_end_bam, optional: true
    script:
        args = task.ext.args ?: " "
        if (params.chromunity) {
            args += "--chromunity "
            if (params.chromunity_merge_distance != null) {
                args += "--chromunity-merge-distance ${params.chromunity_merge_distance} "
            }
        }
        if (params.paired_end) {
            args += "--paired-end "
            if (params.paired_end_minimum_distance != null) {
                args += "--paired-end-minimum-distance ${params.paired_end_minimum_distance} "
            }
            if (params.paired_end_maximum_distance != null) {
                args += "--paired-end-maximum-distance ${params.paired_end_maximum_distance} "
            }
        }
        if (params.summary) {
            args  += "--summary "
        }
        test_task = task.ext.suffix
        //test_task_prefix = task.ext.prefix
    """
    mkdir annotated
    pore-c-py parse-bam "concatemers.mm2.bam" "annotated/${meta.sample_id}" $args
    samtools sort --write-index -o "annotated/${meta.sample_id}.cs.bam" "annotated/${meta.sample_id}.ns.bam"
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
