#!/usr/bin/env nextflow
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include {
    index_ref_fai
    index_ref_mmi
    decompress_ref
} from '../../modules/local/common'



workflow prepare_genome {
    take:
        ref_param
        minimap2_settings
    main:
        // taken from wf-human-variation
        // Check ref and decompress if needed
        ref = null
        ref_index_fp = null
        if (ref_param.toLowerCase().endsWith('gz')) {
            // gzipped ref not supported by some downstream tools (pyfaidx, cram_cache)
            // easier to just decompress and pass it around rather than confusing the user
            decompress_ref(file(ref_param))
            ref = decompress_ref.out.decompressed_ref
        }
        else {
            ref = Channel.fromPath(ref_param, checkIfExists: true)
            ref_index_fp = file(ref_param + '.fai')
        }
        // Create ref index if required
        if (!ref_index_fp || !ref_index_fp.exists()) {
            index_ref = index_ref_fai(ref)
            ref_index = index_ref.reference_index
        }
        else {
            ref_index = Channel.of(ref_index_fp)
        }
        ref_channel = ref.concat(ref_index).buffer(size: 2)
        // create a minimap2 index, not strictly necessary
        mmi = index_ref_mmi(ref, minimap2_settings)

    emit:
        fasta = ref
        fai = ref_index
        mmi = mmi
        minimap2_settings = minimap2_settings
}
