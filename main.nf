#!usr/bin/env nextflow
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2
include {
    fastq_ingress
    xam_ingress
} from "./lib/ingress"
include {
    index_ref_fai
    decompress_ref
    publish_artifact
    merge_namesorted_bams
    merge_namesorted_bams as merge_paired_end_bams
    merge_coordsorted_bams
    mosdepth_coverage
    get_filtered_out_bam
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
   createBed
   mergeBed
} from './modules/local/4dn'


include { prepare_genome } from "./subworkflows/local/prepare_genome"

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

// bamindex will work with bam or fastq format file as input
process index_bam {
    label "wfporec"
    cpus 4
    memory "8 GB"
    input:
        tuple val(meta), path("concatemers.bam")
        val chunk_size
    output:
        tuple val(meta), path("concatemers.bam"), path("concatemers.bam.bci"), path("indexed_chunks.csv")
    shell:
        args = task.ext.args ?: " "
    """
    bamindex build -c ${params.chunk_size} -t ${task.cpus} concatemers.bam
    bamindex dump concatemers.bam.bci > chunks.csv
    awk -F' ' -v OFS=' ' 'NR == 1 {print "ID", \$0; next} {print (NR-2), \$0}' chunks.csv > indexed_chunks.csv
    """
}


process getVersions {
    label "wfporec"
    cpus 4
    memory "4 GB"
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
    memory "4 GB"
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
    label "wf_common"
    cpus 4
    memory "15 GB"
    input:
        val metadata
        path(stats, stageAs: "stats_*")
        path "versions/*"
        path "params.json"
        val wf_version
    output:
        path "wf-pore-c-report.html"
    script:
        String report_name = "wf-pore-c-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --metadata metadata.json \
        --stats $stats \
        --versions versions \
        --params params.json \
        --wf_version $wf_version
    """
}

// Creates a new directory named after the sample alias and moves the ingress results
// into it. So output folder will contain alias named folders with stats.
process collectIngressResultsInDir {
    label "wf_common"
    input:
        // inputs might be `OPTIONAL_FILE` --> stage in different sub-directories
        // to avoid name collisions
        tuple val(meta),
            path(stats, stageAs: "stats/*")
    output:
        // use sub-dir to avoid name clashes (in the unlikely event of a sample alias
        // being `reads` or `stats`)
        tuple path("out/*"), val("ingress_results")
    script:
    String outdir = "out/${meta["alias"]}"
    String metaJson = new JsonBuilder(meta).toPrettyString()
    String stats = stats.fileName.name == OPTIONAL_FILE.name ? "" : stats
    """
    mkdir -p $outdir
    echo '$metaJson' > metamap.json
    mv metamap.json $stats $outdir
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process publish {
    // publish inputs to output directory
    label "wfporec"
    cpus 1
    memory "4 GB"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}

// entrypointworkflow
WorkflowMain.initialise(workflow, params, log)

workflow POREC {
    main:
        Pinguscript.ping_start(nextflow, workflow, params)
        /// PREPARE INPUTS  ///

        if (params.fastq) {
            sample_data = fastq_ingress([
                "input":params.fastq,
                "sample":params.sample,
                "sample_sheet":params.sample_sheet,
                "analyse_unclassified":params.analyse_unclassified,
                "stats": true,
                "fastcat_extra_args": "",
            ])
            // fastq_ingress doesn't have the index; add one extra null for compatibility.
            // We do not use variable name as assigning variable name with a tuple
            // not matching (e.g. meta, bam, bai, stats <- [meta, bam, stats]) causes
            // the workflow to crash.
            sample_data = sample_data
                .map{
                    it.size() == 4 ? it : [it[0], it[1], null, it[2]]
                }
        } else {
        // if we didn't get a `--fastq`, there must have been a `--bam` (as is codified
        // by the schema)
            sample_data = xam_ingress([
                "input":params.bam,
                "sample":params.sample,
                "sample_sheet":params.sample_sheet,
                "analyse_unclassified":params.analyse_unclassified,
                "keep_unaligned": true,
                "stats": true,
            ])
        }

        // create channel of input chimeric reads
        input_reads = sample_data.map{meta, path, index, stats -> [meta, path]}

        if (params.chunk_size > 0) {
            chunks = index_bam(input_reads, channel.value(params.chunk_size))
            // create tuple for each region
            reads = chunks
                .map{meta, bam, bai, chunk_csv ->
                    tuple(meta, bam, bai,chunk_csv.splitCsv(header: ['index','region', 'ref'], skip: 1 , sep:' '))}
                .transpose()
                .map{ meta, bam, bai, chunk_index ->
                    [meta, bam, bai, chunk_index.index, chunk_index.ref]}
        } else {
            // Add optional file and nulls to satisfy channel structure.
            // These values are ignored in digest_align_annotate
            reads = input_reads.combine(Channel.of(tuple(OPTIONAL_FILE, null, null)))
        }
        if (!params.sample_sheet) {
            ch_chunks = reads.map{meta, bam, index, chunk_index, chunk_ref ->
                vcf_file = params.vcf == null ? null : file(params.vcf, checkExists:true)
                tbi_file = vcf_file == null ? null : file(params.vcf + '.tbi')
                [meta + [cutter: params.cutter, vcf:vcf_file, tbi:tbi_file], bam, index, chunk_index, chunk_ref]}
        } else {
            // check meta vcf exists, add tbi and convert to files
            per_sample = sample_data.map{meta, path, index, stats ->
                vcf_file = meta["vcf"] ? file(meta["vcf"], checkExists: true) : null
                tbi_file = vcf_file ? file(meta["vcf"] + '.tbi') : null
                [meta["alias"], vcf_file, tbi_file]}
            // combine with output of ingress
            combined_samples = reads
            .map { [it[0]["alias"], *it] }
            .combine(per_sample, by: 0)
            .map { it[1..-1] }
            // add tuple values to meta data
            pre_chunks = combined_samples.map{meta, bam, index, chunk_index, chunk_ref, vcf_file, tbi_file ->
            [meta + [vcf:vcf_file, tbi:tbi_file], bam, index, chunk_index, chunk_ref]}
            // use params.cutter if it was missing from user provided sample_sheet
            ch_chunks = pre_chunks.map{ meta, bam, index, chunk_index, chunk_ref -> 
                if (meta.cutter && params.cutter){
                    log.warn("Using cutter: ${meta.cutter} from sample sheet column for ${meta.alias}")
                }
                cutter = meta.cutter ?: params.cutter
                return [ meta + ["cutter": cutter], bam, index, chunk_index, chunk_ref]
            }   
        }
        ref = prepare_genome(params.ref, params.minimap2_settings)
        
        /// RUN PORE-C TOOLS ///
        chunks_refs = ch_chunks.combine(ref.mmi).combine(ref.minimap2_settings)

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

        // merge coord-sorted bams by alias
        cs_bam = merge_coordsorted_bams(
            cs_bam_chunks.map(i -> [i[0], i[1]])
            .groupTuple()
        )
        // merge namesorted bams by alias
        ns_bam = merge_namesorted_bams(
            ch_annotated_monomers
            .ns_bam
            .map(i -> [i[0],  i[1]])
            .groupTuple()
        )

        if (params.coverage || params.pairs || params.mcool || params.hi_c) {
            // for each cutter a bed file of the fragments
            digest_ch = create_restriction_bed(
                ch_chunks.map{meta, bam, index, chunk_index, chunk_ref -> meta.cutter}
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
    
        // get metadata and stats files, keeping them ordered (could do with transpose I suppose)
        sample_data.multiMap{ meta, path, index, stats ->
            meta: meta
            stats: stats
        }.set { for_report }
        metadata = for_report.meta.collect()
        // create a file list of the stats, and signal if its empty or not
        stats = for_report.stats | collect
        report = makeReport(
            metadata, stats, software_versions, workflow_params, workflow.manifest.version
        )

        if (params.hi_c){
            hi_c = prepare_hic(merge_pairs.out.merged_pairs.combine(ref.fai))
        }

        if (params.bed){
            bed_chunks = createBed(ch_annotated_monomers.paired_end_bam)
            mergeBed(bed_chunks.groupTuple())

        }
      

        sample_data
        | map {
            meta, path, index, stats ->
            if (stats) [ meta, stats ]
        }
        | collectIngressResultsInDir


        // Group together lists of filtered reads from all the processed chunks
        named_filtered_read_ids = ch_annotated_monomers.filtered_read_ids.groupTuple().map{ meta, read_ids -> tuple(meta.alias, read_ids)}
        named_reads = input_reads.map{ meta, reads -> tuple(meta.alias, reads)}
        // Combine with input reads
        filtered_reads = named_filtered_read_ids.join(named_reads, remainder:false)
        // Retrieve filtered out BAM from list of filtered reads per sample
        filtered_out = get_filtered_out_bam(filtered_reads)


    emit:
        name_sorted_bam = ns_bam
        coord_sorted_bam = cs_bam
        report = report
        ingress_results = collectIngressResultsInDir.out
}

workflow {
    if (params.containsKey("params_sheet")) {
        error = "`--params_sheet` parameter is deprecated. Use parameter `--sample_sheet` instead."
    }
    POREC()
    publish(POREC.out.ingress_results)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
