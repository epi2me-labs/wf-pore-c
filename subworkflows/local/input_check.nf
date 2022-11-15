nextflow.enable.dsl = 2

process chunk_ubam {
    input:
    tuple val(meta), path(ubam)
    val chunk_size
    output:
    tuple val(meta), path("batches/shard_*.bam")
    shell:
    """
    mkdir batches
    picard SplitSamByNumberOfReads I=$ubam OUTPUT=batches SPLIT_TO_N_READS=$chunk_size
    """
}

workflow check_input {
    take:
    chunk_size
    // path sample_sheet TODO
    main:
        // TODO: parse this from sample sheet to enable multiple
        if (params.sample_sheet == null) {
            // samples for now just pass on command line
            ch_samples = Channel.of(
                [
                sample_id: params.sample_id,
                ubam: params.ubam,
                cutter: params.cutter,
                vcf: params.vcf]
                )
        } else {
            ch_samples = Channel.fromPath(params.sample_sheet).splitCsv(header:true)
        }
        ch_samples \
            .map{ create_concatemer_channel(it) }
            .set{ reads }
        if (chunk_size > 0) {
            reads = chunk_ubam(reads, channel.value(chunk_size)).transpose()
        }
    emit:
       // [meta, ubam]  per chunk
       reads

}

def create_concatemer_channel(row) {
    def meta = [:]
    // TODO: do something like fastqingress for minknow directories
    ubam = file(row.ubam, checkExists: true )
    meta.sample_id = row.sample_id ?: ubam.baseName
    meta.cutter = row.cutter
    // TODO: better check for empty string
    meta.vcf = (row.vcf == null || row.vcf == '') ? null : file(row.vcf, checkExists: true)
    meta.tbi = meta.vcf == null ? null : file(row.vcf + '.tbi')
    if (meta.tbi != null) {
        if (!meta.tbi.exists()) { // TODO: create tbi
            exit 1, "ERROR: VCF must be indexed ${row.vcf}"
    }
    }
    [meta, ubam]
}