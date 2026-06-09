nextflow.enable.dsl=2

/*
 * Router-gap analysis modules (included in plink_genotype_modes historically).
 * Not wired from main.nf yet — import here for future params.mod branches.
 * See doc/TODO.md §2 "Router gap".
 */

include { database as DATABASE } from '../../../modules/local/genotype/database/database.nf'
include { kinship as KINSHIP } from '../../../modules/local/dynamic/kinship.nf'
include { population_structure as POPULATION_STRUCTURE } from '../../../modules/local/dynamic/ps.nf'
include { GWAS } from '../../../modules/local/static/gwas/gwas.nf'
include { HAIL } from '../../../modules/local/genotype/hail/hail.nf'
