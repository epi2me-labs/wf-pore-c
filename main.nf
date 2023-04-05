#!usr/bin/env nextflow
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2
include { ingress } from "./lib/ingress"
include {
    index_ref_fai
    decompress_ref
    publish_artifact
    merge_namesorted_bams
    merge_namesorted_bams as merge_paired_end_bams
    merge_coordsorted_bams
    mosdepth_coverage
} from './modules/local/common'

include {
    digest_align_annotate
    haplotagReads as haplotag_alignments
    merge_parquets_to_dataset
} from './modules/local/pore-c'
include {
   to_pairs_file
   pairsToCooler
   merge_mcools
   merge_pairs
   merge_pairs_stats
   create_restriction_bed
   pair_stats_report
} from './modules/local/4dn'


include { prepare_genome } from "./subworkflows/local/prepare_genome"

process chunk_ubam {
    label "wfporec"
    input:
        tuple val(meta), path(bam)
        val chunk_size
    output:
        tuple val(meta), path("batches/shard*.bam")
    shell:
        args = task.ext.args ?: " "
    """
    mkdir batches
    pore-c-py chunk-bam $bam/*.bam "batches/shard" --chunk_size $chunk_size $args
    """
}

// entrypointworkflow
WorkflowMain.initialise(workflow, params, log)

workflow POREC {
    main:
        if (params.disable_ping == false) {
            Pinguscript.ping_post(workflow, 'start', 'none', params.out_dir, params)
        }
        /// PREPARE INPUTS  ///
        // make sure that one of `--fastq` or `--ubam` was given
        def input_type = ['fastq', 'ubam'].findAll { params[it] }
        if (input_type.size() != 1) {
            error "Only provide one of '--fastq' or '--ubam'."
        }
        input_type = input_type[0]

        sample_data = ingress([
        "input":params[input_type],
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "fastcat_stats": true,
        "fastcat_extra_args": "",
        "input_type": input_type])

        // create channel of input concatemers
        reads = sample_data.map{meta,bam,reads -> [meta, bam]}
        if (params.chunk_size > 0) {
            reads = chunk_ubam(sample_data, channel.value(params.chunk_size)).transpose()
        }
        if (params.params_sheet == null) {
            ch_chunks = reads.map{meta, bam ->
                vcf_file = params.vcf == null ? null : file(params.vcf, checkExists:true)
                tbi_file = vcf_file == null ? null : file(params.vcf + '.tbi')
                [ meta + [enzyme: params.cutter, vcf:vcf_file, tbi:tbi_file], bam]}
        } else {
            // create tuples from params sheet
            per_sample_file = Channel.fromPath(params.params_sheet).splitCsv(header:true)
            per_sample = per_sample_file.map { it ->  
                vcf_file = (it["vcf"] == null || it["vcf"] == '' )? null : file(it["vcf"], checkExists: true)
                tbi_file = vcf_file == null ? null : file(it["vcf"] + '.tbi')
                [it["alias"], it["cutter"], vcf_file, tbi_file]}
            // combine with output of ingress
            combined_samples = reads
            .map { [it[0]["alias"], *it] }
            .combine(per_sample, by: 0)
            .map { it[1..-1] }
            // add tuple values to meta data
            ch_chunks = combined_samples.map{meta, bam, cutter, vcf_file, tbi_file ->
             [meta + [enzyme: cutter, vcf:vcf_file, tbi:tbi_file], bam]}
        }
     
        ref = prepare_genome(params.ref, params.minimap2_settings)
        
        /// RUN PORE-C TOOLS ///
        chunks_refs = ch_chunks.map{it -> tuple(it[0], it[1])}.combine(ref.mmi).combine(ref.minimap2_settings)
        ch_annotated_monomers = digest_align_annotate(chunks_refs)

        // create a fork for samples that have phase info available
        ch_annotated_monomers.cs_bam
            .branch{
                to_haplotag: it[0].vcf != null
                no_haplotag: it[0].vcf == null
            }
            .set { haplotag_fork }
        // haplotag bams when we have VCF available
        (haplotag_fork
            .to_haplotag  // [meta, bam bai]
            .combine(ref.fasta)
            .combine(ref.fai)
            .map(i -> {
                [
                i[0], // meta
                i[1], // bam
                i[2], // bai
                i[3], // fasta
                i[4], // fai
                i[0].vcf, // vcf
                i[0].tbi, // tbi
                ]
            })) | haplotag_alignments | set {haplotagged_monomers}

        // merge haplotagged and non-haplotagged coord-sorted bam chunks
        // back to single channel
        haplotag_fork
            .no_haplotag
            .mix(haplotagged_monomers.cs_bam)
            .set { cs_bam_chunks }

        /// MERGE PORE-C BAMS ///

        // merge coord-sorted bams by sample_id
        cs_bam = merge_coordsorted_bams(
            cs_bam_chunks.map(i -> [i[0], i[1]])
            .groupTuple()
        )
        // merge namesorted bams by sample_id
        ns_bam = merge_namesorted_bams(
            ch_annotated_monomers
            .ns_bam
            .map(i -> [i[0],  i[1]])
            .groupTuple()
        )

        if (params.coverage || params.pairs || params.mcool ) {
            // for each enzyme a bed file of the fragments
            digest_ch = create_restriction_bed(
                ch_chunks.map{meta, bam -> meta.enzyme}
                .unique()
                .combine(ref.fasta)
                .combine(ref.fai)
            )
        }

        /// COVERAGE CALCULATIONS
        if (params.coverage) {
            // calculate coverage on the merged BAM
            digest_ch
                .cross(
                    cs_bam
                    .map(i -> [i[0].cutter, i[0], i[1], i[2]]) // [key, meta, bam, bai]
                )
                .map(i -> [
                    i[1][1], // meta
                    i[1][2], // bam
                    i[1][3], // bai
                    i[0][2], // bed
                ]) | mosdepth_coverage | set{ coverage }
        }
        /// 4DN file formats
        if (params.pairs || params.mcool ) {
            (digest_ch
                .cross(
                    ch_annotated_monomers
                    .ns_bam
                    .map(i -> [i[0].enzyme, i[0], i[1]]) // [key, meta, bam]
                )
                .map(i -> [
                        i[1][1], // meta
                        i[1][2], // bam
                        i[0][1], // fai
                        i[0][2], // bed
                    ])
                ) | to_pairs_file | set {pair_chunks}
            if (params.mcool) {
                mcool_chunks = pairsToCooler(
                    pair_chunks
                    .pairs
                    .combine(Channel.of(params.cool_bin_size))
                )
                mcool = merge_mcools(
                    mcool_chunks
                    .groupTuple()
                    .combine(Channel.of(params.mcool_resolutions))
                )
            }
            if (params.pairs) {
                unsorted_pairs = merge_pairs(
                    pair_chunks.pairs.map(i -> [i[0], i[2]]).groupTuple()
                )
                pairs_stats = merge_pairs_stats(
                    pair_chunks.stats.groupTuple()
                )
                pairs_report = pair_stats_report(
                    pairs_stats
                )
          
            }
        }
        /// CHROMUNITY
        if (params.chromunity) {
            chromunity_pq = merge_parquets_to_dataset(
                ch_annotated_monomers
                .chromunity_pq
                .groupTuple()
            )
        }

        /// Paired end bams
        if (params.paired_end) {
            pe_bam = merge_paired_end_bams(
                ch_annotated_monomers
                .paired_end_bam
                .map(i -> [i[0],  i[1]])
                .groupTuple()
            )
        }

    emit:
        name_sorted_bam = ns_bam
        coord_sorted_bam = cs_bam
}

workflow {
    POREC()
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, 'end', 'none', params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, 'error', "$workflow.errorMessage", params.out_dir, params)
    }
}
