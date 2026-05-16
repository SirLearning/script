nextflow.enable.dsl=2

include { GENETICS_PIPELINE } from './subworkflows/local/genetics_pipeline.nf'

workflow {
    GENETICS_PIPELINE()
}
