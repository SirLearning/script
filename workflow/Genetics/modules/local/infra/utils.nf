nextflow.enable.dsl=2

// Index of shared DSL helpers (job config, vmap4 ref maps, tool paths, PLINK reuse).
// Consumers must include from infra_*.nf directly — Nextflow does not re-export Groovy def via this barrel.

include {
    getJavaSetupScript
    getSoftwareConfig
    getWorkshopJarConfig
    getSamtoolsConfig
    validatePaths
} from './infra_tools.nf'

include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
    hasPlinkBasicInfoForMergedTests
    listMergedSubgenomeSnpQcPlotTuples
    listMergedSubgenomeTestBfileTuples
} from './infra_plink_reuse.nf'

include {
    getJobConfig
    getCallingJobConfig
    getPopulationConfig
    getChrConfig
    getVcfIdFromPath
} from './infra_job_config.nf'

include {
    getServerPopulations
    getTigerJarConfig
} from './infra_tiger.nf'

include {
    getTaxaBamMapFile_v1
    getPopDepTaxaBamFile_v1
    getRefV1SubChr
    getRefV1NameChr
    getRefV1ChrName
    getRefV1ChrOffset
    getRefV1ChrLength
} from './infra_ref_v1.nf'
