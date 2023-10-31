#!usr/bin/env nextflow
nextflow.enable.dsl = 2


process to_pairs_file {
    label 'wfporec'
    input:
        tuple val(meta), path("monomers.mm2.ns.bam"), path("fasta.fai"), path("fragments.bed")
    output:
        tuple val(meta), path("fasta.fai"), path("${meta.sample_id}.pairs.gz"), emit: "pairs"
        tuple val(meta), path("${meta.sample_id}.stats.txt"), emit: "stats"
    shell:
        def args = task.ext.args ?: "--drop-sam --drop-seq --expand --add-pair-index --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches"
    """
    pairtools parse2  \
    --output-stats "${meta.sample_id}.stats.txt" \
    -c "fasta.fai" --single-end --readid-transform 'readID.split(":")[0]' \
    $args "monomers.mm2.ns.bam" > extract_pairs.tmp
    pairtools restrict  -f "fragments.bed" -o "${meta.sample_id}.pairs.gz" extract_pairs.tmp
    rm -rf extract_pairs.tmp
    """
}


process prepare_hic {
    label 'wfporec'
    input:
        tuple val(meta), path("input.pairs.gz"), path("fasta.fai")
    output:
        path "${meta.sample_id}.hic", emit: hic
    """
    cut -f1,2 fasta.fai > sizes.genome
    pairtools flip input.pairs.gz -c sizes.genome  > flipped.pairs.tmp
    pairtools sort flipped.pairs.tmp > sorted.pairs.tmp
    pairtools dedup sorted.pairs.tmp > dedup.pairs.tmp
    java -jar /home/epi2melabs/juicer_tools_1.22.01.jar pre dedup.pairs.tmp "${meta.sample_id}.hic" sizes.genome
    rm -rf "*.pairs.tmp"
    """
}

process merge_pairs {
    label 'wfporec'
    input:
        tuple val(meta), path('to_merge/{?}.gz')
    output:
        tuple val(meta), path("${prefix}.pairs.gz"), emit: merged_pairs
    shell:
        prefix = task.ext.prefix ?: "${meta.sample_id}"
        def args = task.ext.args ?: "--concatenate"
    """
    # pass a quoted glob, pairtools will do its own globbing
    pairtools merge -o "${prefix}.pairs.gz" $args 'to_merge/*'
    """
}

process merge_pairs_stats {
    label 'wfporec'
    input: 
        tuple val(meta), path('to_merge/src*.stats.txt')
    output: 
        tuple val(meta), path("${prefix}.pairs.stats.txt")
    shell:
        prefix = task.ext.prefix ?: "${meta.sample_id}"
        def args = task.ext.args ?: "--merge "
    """
    pairtools stats -o "${prefix}.pairs.stats.txt" $args to_merge/src*.stats.txt
    """
}

process pair_stats_report {
    label 'wfporec'
    input: 
        tuple val(meta), path("pairs.stats.txt")
    output:
        tuple val(meta), path("${prefix}.pairs.stats.html")
    shell:
        prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    create_pairs_report.py "pairs.stats.txt" "${prefix}.pairs.stats.html"
    """
}

process create_restriction_bed {
    label 'wfporec'
    input:
        tuple val(enzyme), path("reference.fasta"), path("reference.fasta.fai")
    output:
        tuple val(enzyme), path("reference.fasta.fai"), path("fragments.bed")
    shell:
        def args = task.ext.args ?: " "
    """
    cooler digest -o "fragments.bed" $args "reference.fasta.fai" "reference.fasta"  $enzyme
    """
}

process pairsToCooler {
    label 'wfporec'
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
    label 'wfporec'
    input:
        tuple val(meta), path('to_merge/src*.cool'), val(resolutions)
    output:
        tuple val(meta), path("${prefix}.mcool")
    shell:
        prefix = task.ext.prefix ?: "${meta.sample_id}"
        def args = task.ext.args ?: " "
    """
    cooler merge  ${prefix}.cool $args to_merge/src*.cool
    cooler zoomify -r ${resolutions} -o ${prefix}.mcool  ${prefix}.cool
    """
}
