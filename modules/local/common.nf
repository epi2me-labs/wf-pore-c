process index_ref_fai {
    label 'wfporec'
    cpus 1
    input:
        path "reference.fasta"
    output:
        path "reference.fasta.fai", emit: reference_index
    """
    samtools faidx "reference.fasta"
    """
}

process index_ref_mmi {
    label 'wfporec'
    input:
        path "reference.fasta"
        val(minimap_settings)
    output:
        path "reference.fasta.mmi"
    """
    minimap2 ${minimap_settings} -d "reference.fasta.mmi" "reference.fasta"
    """
}

// NOTE -f required to compress symlink
process decompress_ref {
    label 'wfporec'
    cpus 1
    input:
        path compressed_ref
    output:
        path "${compressed_ref.baseName}", emit: decompressed_ref
    """
    gzip -df "${compressed_ref}"
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish_artifact {
    label 'wfporec'
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}

// what is better, this or the one in input_check that uses porecpy? - that is actually being used.
process chunk_bam {
    label 'wfporec'
    input:
        tuple val(meta), path("concatemers.bam")
        val chunk_size
    output:
        tuple val(meta), path "batches/*.bam"
    shell:
    """
    mkdir batches
    picard SplitSamByNumberOfReads I="concatemers.bam" OUTPUT=batches SPLIT_TO_N_READS=${chunk_size}
    """
}

// TODO rewrite as single merge process
process merge_namesorted_bams {
    label 'wfporec'
    cpus params.threads
    input:
        tuple val(meta), path('to_merge/src*.bam')
    output:
        tuple val(meta), path("${prefix}.${suffix}.bam")
    shell:
        suffix = task.ext.suffix ?: "ns"
        prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    samtools cat --threads $task.cpus -o "${prefix}.${suffix}.bam" --no-PG to_merge/src*.bam
    """
}

process merge_coordsorted_bams {
    label 'wfporec'
    input:
        tuple val(meta), path('to_merge/src*.bam')
    output:
        tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.csi")
    shell:
        prefix = task.ext.prefix ?: "${meta.sample_id}.cs"
    """
    samtools merge --threads $task.cpus -o "${prefix}.bam" -p --write-index --no-PG to_merge/src*.bam
    """
}

process mosdepth_coverage {
    label 'wfporec'
    cpus params.threads
    input:
        tuple val(meta),
              path("concatemers.cs.bam"),
              path("concatemers.cs.bam.csi"),
              path("fragments.bed")
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
    mosdepth --threads $task.cpus --d4 --by "fragments.bed" $args $prefix "concatemers.cs.bam"
    """
}
