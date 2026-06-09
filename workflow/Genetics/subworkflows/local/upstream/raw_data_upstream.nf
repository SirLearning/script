nextflow.enable.dsl=2

/*
 * Upstream raw-data tracks (not routed from main.nf yet).
 * Align: USB BAM transfer + MD5 verification.
 * Calling: FastCall3 disc / blib / scan on BAM + reference.
 *
 * Entry scripts remain standalone:
 *   modules/local/genotype/align/align.nf
 *   modules/local/genotype/calling/caller.nf
 */

include { RUN_ALIGN_USB_TRANSFER } from '../../../modules/local/genotype/align/align.nf'
include { run_FastCall3 } from '../../../modules/local/genotype/calling/caller.nf'
