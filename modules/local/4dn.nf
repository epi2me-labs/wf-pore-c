#!usr/bin/env nextflow
nextflow.enable.dsl = 2

process to_pairs_file {

    input:
    tuple val(meta), path(bam), path(fai), path(bed)

    output:
    tuple val(meta), path(fai), path("${bam.baseName}.pairs.gz"), emit: "pairs"
    tuple val(meta), path("${bam.baseName}.stats.txt"), emit: "stats"

    shell:
    def args = task.ext.args ?: "--drop-sam --drop-seq --expand --add-pair-index --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches"
    """
    pairtools parse2  \
    --output-stats ${bam.baseName}.stats.txt \
    -c $fai  --single-end --readid-transform 'readID.split(":")[0]'  \
    $args $bam | pairtools restrict  -f $bed -o ${bam.baseName}.pairs.gz
    """
}

process merge_pairs {
    input:
    tuple val(meta), path('to_merge/src*.pairs.gz')
    output:
        output:
    tuple val(meta), path("${prefix}.pairs.gz")
    shell:
    prefix = task.ext.prefix ?: "${meta.sample_id}"
    def args = task.ext.args ?: "--concatenate"
    """
    pairtools merge -o ${prefix}.pairs.gz $args to_merge/src*.pairs.gz
    """
}

process merge_pairs_stats {
    input: tuple val(meta), path('to_merge/src*.stats.txt')
    output: tuple val(meta), path("${prefix}.pairs.stats.txt")
    shell:
    prefix = task.ext.prefix ?: "${meta.sample_id}"
    def args = task.ext.args ?: "--merge "
    """
    pairtools stats -o ${prefix}.pairs.stats.txt $args to_merge/src*.stats.txt
    """
}


process create_restriction_bed {
    input: tuple val(enzyme), path(fasta), path(fai)
    output: tuple val(enzyme), path(fai), path("${fasta.baseName}.${enzyme}.fragments.bed")
    shell:
    def args = task.ext.args ?: " "
    """
    cooler digest -o ${fasta.baseName}.${enzyme}.fragments.bed $args $fai $fasta $enzyme
    """
}



process pairsToCooler {
    input:
    tuple val(meta), path(fai), path(pairs), val(min_bin_width)

    output:
    tuple val(meta), path("${pairs.baseName}.cool")

    shell:
    """
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $fai:${min_bin_width} $pairs ${pairs.baseName}.cool
    """
}


process merge_mcools {
    input: tuple val(meta), path('to_merge/src*.cool'), val(resolutions)
    output: tuple val(meta), path("${prefix}.mcool")
    shell:
    prefix = task.ext.prefix ?: "${meta.sample_id}"
    def args = task.ext.args ?: " "
    """
    cooler merge  ${prefix}.cool $args to_merge/src*.cool
    cooler zoomify -r ${resolutions} -o ${prefix}.mcool  ${prefix}.cool
    """
}

