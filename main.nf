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
   prepare_hic
} from './modules/local/4dn'


include { prepare_genome } from "./subworkflows/local/prepare_genome"

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

// bamindex will work with bam or fastq format file as input
process index_ubam {

    label "wfporec"
    input:
        tuple val(meta), path("concatemers.bam")
        val chunk_size
    output:
        tuple val(meta), path("concatemers.bam"), path("concatemers.bam.bci"), path("chunks.csv")
    shell:
        args = task.ext.args ?: " "
    """
    bamindex build -c ${params.chunk_size} -t ${task.cpus} concatemers.bam
    bamindex dump concatemers.bam.bci > chunks.csv
    """
}


process getVersions {
    label "wfporec"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    pairtools --version | sed 's/\\<version\\>//g' >> versions.txt
    whatshap --version | sed 's/^/whatshap,/' >> versions.txt
    pore-c-py --version | sed 's/ /,/' >> versions.txt
    samtools --version | (head -n 1 && exit 0) | sed 's/ /,/' >> versions.txt
    """
}


process getParams {
    label "wfporec"
    cpus 1
    output:
        path "params.json"
    script:
        String paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process makeReport {
    label "wfporec"
    input:
        val metadata
        path "per_read_stats.tsv"
        path "versions/*"
        path "params.json"
    output:
        path "wf-pore-c-report.html"
    script:
        String report_name = "wf-pore-c-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --metadata metadata.json \
        --stats per_read_stats.tsv \
        --versions versions \
        --params params.json

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
        input_reads = sample_data.map{meta,bam,reads -> [meta, bam]}
        if (params.chunk_size > 0) {
            chunks = index_ubam(input_reads, channel.value(params.chunk_size))
            // create tuple for each region
            reads = chunks.map{it -> 
                tuple(
                    it[0], it[1], it[2], 
                    it[3].splitCsv(
                        header: ['region', 'ref'], skip: 1 , sep:' ').ref)}.transpose()
        } else {
            reads = input_reads.combine(Channel.of(tuple(OPTIONAL_FILE, 'NA')))
        }
        if (!params.sample_sheet) {
            ch_chunks = reads.map{meta, bam, index, ref ->
                vcf_file = params.vcf == null ? null : file(params.vcf, checkExists:true)
                tbi_file = vcf_file == null ? null : file(params.vcf + '.tbi')
                [meta + [cutter: params.cutter, vcf:vcf_file, tbi:tbi_file], bam, index, ref]}
        } else {
            // check meta vcf exists, add tbi and convert to files
            per_sample = sample_data.map{meta,bam,reads -> meta
                vcf_file = meta["vcf"] ? file(meta["vcf"], checkExists: true) : null
                tbi_file = vcf_file ? file(meta["vcf"] + '.tbi') : null
                [meta["alias"], vcf_file, tbi_file]}
            // combine with output of ingress
            combined_samples = reads
            .map { [it[0]["alias"], *it] }
            .combine(per_sample, by: 0)
            .map { it[1..-1] }
            // add tuple values to meta data
            ch_chunks = combined_samples.map{meta, bam, index, ref, vcf_file, tbi_file ->
            [meta + [vcf:vcf_file, tbi:tbi_file], bam, index, ref]}
        }
        ref = prepare_genome(params.ref, params.minimap2_settings)
        
        /// RUN PORE-C TOOLS ///
        chunks_refs = ch_chunks.map{it[0..3]}.combine(ref.mmi).combine(ref.minimap2_settings)

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

        if (params.coverage || params.pairs || params.mcool || params.hi_c) {
            // for each cutter a bed file of the fragments
            digest_ch = create_restriction_bed(
                ch_chunks.map{meta, bam, index, ref -> meta.cutter}
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
        if (params.pairs || params.mcool || params.hi_c) {
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
            if (params.pairs || params.hi_c) {
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

        // Make a report
        software_versions = getVersions()
        workflow_params = getParams()
        metadata = input_reads.map { it[0] }.toList()
        per_read_stats = sample_data.map {
            it[2] ? it[2].resolve('per-read-stats.tsv') : null
        }
        | collectFile ( keepHeader: true )

        report = makeReport(
            metadata, per_read_stats, software_versions, workflow_params
        )

        if (params.hi_c){
            hi_c = prepare_hic(merge_pairs.out.merged_pairs.combine(ref.fai))
        }


    emit:
        name_sorted_bam = ns_bam
        coord_sorted_bam = cs_bam
        report = report
}

workflow {
    if (params.containsKey("params_sheet")) {
        error = "`--params_sheet` parameter is deprecated. Use parameter `--sample_sheet` instead."
    }
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
