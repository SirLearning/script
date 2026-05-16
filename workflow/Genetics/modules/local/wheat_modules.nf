nextflow.enable.dsl=2

include { integrated_wheat as INTEGRATED_WHEAT } from '../../integrated/integrated_wheat.nf'

workflow RUN_WHEAT_INTEGRATED {
    main:
    INTEGRATED_WHEAT()
}
