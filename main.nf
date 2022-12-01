#!usr/bin/env nextflow
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include {
    index_ref_fai
    decompress_ref
    publish_artifact
    chunk_ubam
    merge_namesorted_bams
    merge_namesorted_bams as merge_paired_end_bams
    merge_coordsorted_bams
    mosdepth_coverage
} from './modules/local/common'

include {
    digest_concatemers
    minimap2_ubam_namesort as map_monomers
    haplotagReads as haplotag_alignments
    annotate_monomers
    merge_parquets_to_dataset
} from './modules/local/pore-c'
include {
   to_pairs_file
   pairsToCooler
   merge_mcools
   merge_pairs
   merge_pairs_stats
   create_restriction_bed
} from './modules/local/4dn'


include { check_input } from "./subworkflows/local/input_check"
include { prepare_genome } from "./subworkflows/local/prepare_genome"

// entrypointworkflow
WorkflowMain.initialise(workflow, params, log)

workflow POREC {
    main:
        if (params.disable_ping == false) {
            Pinguscript.ping_post(workflow, 'start', 'none', params.out_dir, params)
        }
        /// PREPARE INPUTS  ///
        // create channel of input concatemers
        check_input(
            params.chunk_size
        ).set { ch_chunks }  // [meta, ubam_chunk]

        // scalar channel of genome files.
        ref = prepare_genome(params.ref, params.minimap2_settings)

       /// RUN PORE-C TOOLS ///
        // digest concatemers into monomers
        ch_monomers = digest_concatemers(ch_chunks)
        ch_monomers
            .combine(ref.mmi)
            .combine(ref.minimap2_settings)
            .map{     // meta, ubam, mmi, minimap2_settings
                it -> [it[0], it[1], it[2], it[3]]
            }
            .set{ch_mapping}
        // map monomers against ref genome and add pore-c specific tags
        ch_annotated_monomers  = map_monomers(ch_mapping) | annotate_monomers

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
                ch_chunks
                .map(i -> [i[0].cutter])
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
                    .map(i -> [i[0].cutter, i[0], i[1]]) // [key, meta, bam]
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
