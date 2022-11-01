process digest_concatemers {
    input:
    tuple val(meta), path(ubam)

    output:
    tuple val(meta), path("${ubam.baseName}.monomers.bam")

    script:
    """
    pore-c2 utils digest-concatemers $ubam ${meta.cutter} ${ubam.baseName}.monomers.bam
    """
}

//TODO remove
process index_ref_mmi {
    input:
    path ref

    output:
    path "${ref}.mmi"

    script:
    """
    minimap2 ${params.minimap_settings} -d ${ref}.mmi $ref
    """
}

process minimap2_ubam_namesort {
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    input:
        tuple val(meta),
              path(ubam),
              path(mmi),
              val(minimap2_settings)
    output:
        tuple val(meta), path("${ubam.baseName}.mm2.bam")

    script:
    """
    samtools fastq --threads ${params.ubam_bam2fq_threads} -T '*' $ubam |
    minimap2 -ay -t ${params.ubam_map_threads}  ${minimap2_settings} $mmi - |
    samtools sort --threads ${params.ubam_sort_threads} -n -o ${ubam.baseName}.mm2.bam
    """
}

process annotate_monomers {
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta),
          path("annotated/${bam.baseName}.ns.bam"),
          emit: ns_bam
    tuple val(meta),
          path("annotated/${bam.baseName}.cs.bam"),
          path("annotated/${bam.baseName}.cs.bam.csi"),
          emit: cs_bam
    script:
    """
    mkdir annotated
    pore-c2 utils process-monomer-alignments $bam annotated/${bam.baseName}
    samtools sort --write-index -o annotated/${bam.baseName}.cs.bam annotated/${bam.baseName}.ns.bam
    """
}

process haplotagReads {
    label 'wftemplate'

    input:
    tuple val(meta),
          path(bam),
          path(bai),
          path(ref),
          path(fai),
          path(phased_vcf),
          path(phased_vcf_tbi)
    output:
    tuple val(meta),
          path("${bam.baseName}.ht.bam"),
          path("${bam.baseName}.ht.bam.csi"),
          emit: "cs_bam"
    tuple val(meta),
          path("${bam.baseName}.ht.txt.gz"),
          emit: "haplotagged_monomers"
    shell:
    """
    whatshap haplotag \
      --reference $ref \
      -o ${bam.baseName}.ht.bam \
      --output-haplotag-list ${bam.baseName}.ht.txt.gz \
      --ignore-read-groups \
      $phased_vcf $bam
    samtools index -c ${bam.baseName}.ht.bam
    """
}


