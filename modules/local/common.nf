process index_ref_fai {
    cpus 1
    input:
        file reference
    output:
        path "${reference}.fai", emit: reference_index
    """
    samtools faidx ${reference}
    """
}

process index_ref_mmi {
    input:
    path(ref)
    val(minimap_settings)

    output:
    path("${ref}.mmi")

    script:
    """
    minimap2 ${minimap_settings} -d ${ref}.mmi $ref
    """
}



// NOTE -f required to compress symlink
process decompress_ref {
    cpus 1
    input:
        file compressed_ref
    output:
        path "${compressed_ref.baseName}", emit: decompressed_ref
    """
    gzip -df ${compressed_ref}
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish_artifact {
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}

process chunk_ubam {
    input:
    tuple val(meta), path(ubam)
    val chunk_size
    output:
    tuple val(meta), path "batches/*.bam"
    shell:
    """
    mkdir batches
    picard SplitSamByNumberOfReads I=$ubam OUTPUT=batches SPLIT_TO_N_READS=$chunk_size
    """
}


// TODO rewrite as single merge process
process merge_namesorted_bams {
    input:
    tuple val(meta), path('to_merge/src*.bam')
    output:
    tuple val(meta), path("${prefix}.bam")
    shell:
    prefix = task.ext.prefix ?: "${meta.sample_id}.ns"
    """
    samtools merge -n --threads $task.cpus -o ${prefix}.bam --no-PG to_merge/src*.bam
    """
}

process merge_coordsorted_bams {
    input:
    tuple val(meta), path('to_merge/src*.bam')
    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.csi")
    shell:
    prefix = task.ext.prefix ?: "${meta.sample_id}.cs"
    """
    samtools merge --threads $task.cpus -o ${prefix}.bam --write-index --no-PG to_merge/src*.bam
    """
}


process mosdepth_coverage {
    input:
    tuple val(meta),
          path(bam),
          path(bai),
          path(bed)
    output:
    tuple val(meta),
          path("${prefix}.per-base.d4"),
          emit: d4
    tuple val(meta),
          path("${prefix}.regions.bed.gz"),
          path("${prefix}.regions.bed.gz.csi"),
          emit: regions
    tuple val(meta),
          path("${prefix}.thresholds.bed.gz"),
          path("${prefix}.thresholds.bed.gz.csi"),
          emit: thresholds
    tuple val(meta),
          path("${prefix}.mosdepth.*"),
          emit: summaries
    shell:
    prefix = task.ext.prefix ?: "${meta.sample_id}"
    args = task.ext.args ?: "--thresholds 1,10,30,60,100"
    """
    mosdepth --threads $task.cpus --d4 --by $bed $args $prefix $bam
    """
}
