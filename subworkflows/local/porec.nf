
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



workflow porec {
    take:
        input_ch   // [[meta, ubam], [meta, ubam],...]
        ref        // path to fasta
        chunk_size // integer channel
    main:
        if (chunk_size > 0) {
            reads = chunk_ubam(reads, chunk_size).transpose()
        } else {
            reads = input_ch
        }


    emit:



}
