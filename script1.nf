nextflow.enable.dsl = 2

log.info("UBAM: ${params.ubam}")

process splitConcatemers {
    label 'wftemplate'
    input:
    path concatemer
    val enzyme

    output:
    path 'monomers.bam'

    script:
    """
    pore-c2 utils digest-concatemers $concatemer $enzyme monomers.bam
    """
}

process indexGenome {
    label 'wftemplate'
    input:
    tuple val(ref_id), path(ref)
    val minimap_settings

    output:
    path "${ref_id}.mmi"
    path "${ref_id}.fasta.fai"

    script:
    """
    samtools faidx --fai-idx ${ref_id}.fasta.fai $ref
    minimap2 $minimap_settings -d ${ref_id}.mmi $ref
    """
}

process mapMonomers {
    label 'wftemplate'

    input:
    path monomers
    path mmi
    val minimap_settings

    output:
    path 'aligned_monomers.bam'

    script:
    """
    samtools fastq -T '*' $monomers |
    minimap2 -ay  $minimap_settings $mmi - |
    samtools sort -t MI -n -o aligned_monomers.bam
    """
}

process annotateMonomers {
    label 'wftemplate'

    input:
    path aligned_monomers

    output:
    path 'annotated_monomers.ns.bam'
    path 'annotated_monomers.cs.bam'
    path 'annotated_monomers.cs.bam.csi'

    script:
    """
    pore-c2 utils process-monomer-alignments $aligned_monomers annotated_monomers
    samtools sort --write-index -o annotated_monomers.cs.bam annotated_monomers.ns.bam
    """
}

process toPairs {
    label 'wftemplate'

    input:
    path annotated_monomers
    path fai

    output:
    path 'monomers.pairs.gz'
    path 'monomers.stats.txt'

    shell:
    """
    pairtools parse2  \
    -o monomers.pairs.gz \
    --output-stats monomers.stats.txt \
    -c $fai  \
    --single-end  \
    --readid-transform 'readID.split(":")[0]'  \
    --drop-sam  \
    --drop-seq \
    --expand   \
    $annotated_monomers
    """
}

process pairsToCooler {
    label 'wftemplate'

    input:
    path pairs
    path fai
    val resolutions

    output:
    path 'contacts.mcool'

    shell:
    """
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $fai:${resolutions[0]} $pairs contacts.cooler
    cooler zoomify -r ${resolutions.join(',')} -o contacts.mcool  contacts.cooler
    """
}

process haplotagReads {
    label 'wftemplate'

    input:
    path coordsorted_bam
    path coordsorted_bai
    tuple val(ref_id), path(ref)
    path fai
    path phased_vcf
    path phased_vcf_tbi

    output:
    path 'haplotagged.bam'
    path 'haplotagged.reads.txt'
    path 'haplotagged.bam.csi'

    shell:
    """
    whatshap haplotag --reference $ref \
    -o haplotagged.bam \
    --output-haplotag-list haplotagged.reads.txt \
    --ignore-read-groups \
    $phased_vcf $coordsorted_bam

    samtools index -c haplotagged.bam
    """
}

workflow {
    resolutions = Channel.of(params.mcool_resolutions).splitCsv().take(1)
    resolutions.view()

    concatemer_ch = Channel.fromPath(params.ubam)
    monomers = splitConcatemers(concatemer_ch, 'NlaIII')
    monomers.view()

    minimap_settings = Channel.from('-x map-ont')
    ref_ch = Channel.fromPath(params.ref).map({ [it.getSimpleName(), it] })
    (mmi_ch, fai_ch) = indexGenome(ref_ch, minimap_settings)

    aligned_monomers = mapMonomers(monomers, mmi_ch, minimap_settings)
    aligned_monomers.view()

    (
    namesorted_annotated_monomers,
    coordsorted_annotated_monomers,
    coordsorted_annotated_monomers_idx
    ) = annotateMonomers(aligned_monomers)
    namesorted_annotated_monomers.view()

    (pairs, pair_stats) = toPairs(namesorted_annotated_monomers, fai_ch)
    pairs.view()
    pair_stats.view()

    mcool = pairsToCooler(pairs, fai_ch, resolutions)
    mcool.view()

    phased_vcf = Channel.fromPath(params.phased_vcf)
    phased_vcf_tbi = Channel.fromPath(params.phased_vcf + '.tbi')
    (haplotagged_bam, haplotagged_reads) = haplotagReads(
            coordsorted_annotated_monomers,
            coordsorted_annotated_monomers_idx,
            ref_ch,
            fai_ch,
            phased_vcf,
            phased_vcf_tbi)
    haplotagged_bam.view()
}
