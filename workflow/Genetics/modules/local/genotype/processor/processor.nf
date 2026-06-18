nextflow.enable.dsl=2

include { getTigerJarConfig } from '../../infra/infra_tiger.nf'
include {
    hasMergedSubgenomeTestPfiles
    hasPlinkBasicInfoForMergedTests
    listMergedSubgenomeTestPfileTuples
    listMergedSubgenomeTestBfileTuples
} from '../../infra/infra_plink_reuse.nf'
include { hasVariantMqForMergedTests } from '../../infra/infra_mq_reuse.nf'
include { hasVariantPopdepForMergedTests } from '../../infra/infra_popdep_reuse.nf'

include {
    format_vcf_bgzip
    format_vcf_bgzip_idx
    arrange_sample_info
    arrange_vcf_wheat_chr_by_awk
    merge_arranged_vcf
    format_vcf_plink
} from './processor_vcf.nf'

include {
    subsampling_pfile_for_test
    subsampling_common_variant_pfile_for_test
    merge_subgenome_test_pfile
} from './processor_test.nf'

include {
    mk_plink_basic_info_camp_pop_with_filter
    mk_plink_basic_info_camp_pop
    plink2_pca
    plink2_tagsnp_prune
    awk_depth_cnv_call
    mk_plink_basic_info
    calc_plink_ld_unphased
    calc_plink_ld_crosschr_random
} from './processor_plink2.nf'

include { mk_vcftools_basic_info } from './processor_legacy.nf'

include { calc_population_depth; popdep_tiger_gz_to_bgzip_tabix; annotate_subgenome_variant_popdep } from './processor_depth.nf'

include { annotate_subgenome_variant_mq } from './processor_mq.nf'

include {
    filter_sample_plink
    filter_sample_after_variant_plink
    filter_variant_plink
    filter_variant_after_sample_plink
    filter_vcf_v0_vcftools
    filter_vcf_v0_bcftools
} from './processor_filter.nf'

include {
    plink2_assess_debug_slice
    quick_count
    bcftools_qc_assess
} from './processor_assess.nf'

// --- Genotype processor workflows ---

workflow test_plink_processor {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    def reuseMerged = hasMergedSubgenomeTestPfiles(params.process_dir)
    def merge_pfile
    def merge_bfile

    if (reuseMerged) {
        log.info "${params.c_green}Reuse merged test pfiles from:${params.c_reset} ${params.process_dir} (skip per-chr thin + merge)"
        merge_pfile = Channel.from(listMergedSubgenomeTestPfileTuples(params.process_dir))
        merge_bfile = Channel.from(listMergedSubgenomeTestBfileTuples(params.process_dir))
        if (hasPlinkBasicInfoForMergedTests(params.process_dir)) {
            log.info "${params.c_green}Reuse PLINK2 basic info from ${params.process_dir}/variant and sample/${params.c_reset}"
            basic_info_out = merge_pfile.multiMap { id, chr, prefix, pgen, psam, pvar ->
                def d = params.process_dir
                smiss: [ id, chr, file("${d}/sample/${id}.info.smiss") ]
                vmiss: [ id, chr, file("${d}/variant/${id}.info.vmiss") ]
                scount: [ id, chr, file("${d}/sample/${id}.info.scount") ]
                gcount: [ id, chr, file("${d}/variant/${id}.info.gcount") ]
                afreq: [ id, chr, file("${d}/variant/${id}.info.afreq") ]
                hardy: [ id, chr, file("${d}/variant/${id}.info.hardy") ]
            }
        } else {
            log.info "${params.c_green}Making basic info from reused merged pfiles.${params.c_reset}"
            basic_info_out = mk_plink_basic_info(merge_pfile)
        }
        ld_out = calc_plink_ld_unphased(merge_pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_pfile)
    } else if (params.process_dir) {
        log.info "${params.c_green}Using main_raw pfiles in dir:${params.c_reset} ${params.process_dir}"
        source_pfiles = vcf_in.map { id, chr, _vcf ->
            def prefix = "${params.process_dir}/${id}.plink2"
            def pgen = file("${prefix}.pgen")
            def psam = file("${prefix}.psam")
            def pvar = file("${prefix}.pvar")
            return [ id, chr, prefix, pgen, psam, pvar ]
        }
        subsampling_out = subsampling_pfile_for_test(source_pfiles, params.thin_rate)

        // Map pfiles to include subgenome info
        def grouped_pfile = subsampling_out.pfile.map { id, chr, prefix, pgen, psam, pvar ->
            def subgenome = "Others"
            def c = chr.toString()
            if (["1","2","7","8","13","14","19","20","25","26","31","32","37","38"].contains(c)) subgenome = "A"
            else if (["3","4","9","10","15","16","21","22","27","28","33","34","39","40"].contains(c)) subgenome = "B"
            else if (["5","6","11","12","17","18","23","24","29","30","35","36","41","42"].contains(c)) subgenome = "D"
            return [ subgenome, id, chr, prefix, pgen, psam, pvar ]
        }.groupTuple(by: 0)

        def merge_out = merge_subgenome_test_pfile(grouped_pfile)
        merge_pfile = merge_out.pfile
        merge_bfile = merge_out.bfile

        log.info "${params.c_green}Making basic info normally.${params.c_reset}"
        basic_info_out = mk_plink_basic_info(merge_out.pfile)
        ld_out = calc_plink_ld_unphased(merge_out.pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_out.pfile)
    } else {
        // preprocess VCF to PLINK format
        preprocess_out = plink_preprocess(vcf_in)
        log.info "${params.c_green}No process_dir specified, subsampling VCFs normally.${params.c_reset}"
        // subsample pfile for testing
        subsampling_out = subsampling_pfile_for_test(preprocess_out.pfile, params.thin_rate)

        // Map pfiles to include subgenome info
        def grouped_pfile = subsampling_out.pfile.map { id, chr, prefix, pgen, psam, pvar ->
            def subgenome = "Others"
            def c = chr.toString()
            if (["1","2","7","8","13","14","19","20","25","26","31","32","37","38"].contains(c)) subgenome = "A"
            else if (["3","4","9","10","15","16","21","22","27","28","33","34","39","40"].contains(c)) subgenome = "B"
            else if (["5","6","11","12","17","18","23","24","29","30","35","36","41","42"].contains(c)) subgenome = "D"
            return [ subgenome, id, chr, prefix, pgen, psam, pvar ]
        }.groupTuple(by: 0)

        def merge_out = merge_subgenome_test_pfile(grouped_pfile)

        log.info "${params.c_green}Making basic info normally.${params.c_reset}"
        basic_info_out = mk_plink_basic_info(merge_out.pfile)
        ld_out = calc_plink_ld_unphased(merge_out.pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_out.pfile)
        merge_pfile = merge_out.pfile
        merge_bfile = merge_out.bfile
    }

    if (reuseMerged && hasVariantMqForMergedTests(params.process_dir)) {
        log.info "${params.c_green}Reuse variant MQ from ${params.process_dir}/variant${params.c_reset}"
        mq_out = merge_pfile.multiMap { id, chr, prefix, pgen, psam, pvar ->
            def d = params.process_dir
            mq: [ id, chr, file("${d}/variant/${id}.mq.info.tsv") ]
        }
    } else {
        log.info "${params.c_green}Annotating variant MQ from frozen site_mq.ref (${params.mq_dir})${params.c_reset}"
        mq_out = annotate_subgenome_variant_mq(merge_pfile)
    }

    if (reuseMerged && hasVariantPopdepForMergedTests(params.process_dir)) {
        log.info "${params.c_green}Reuse variant popdep from ${params.process_dir}/variant${params.c_reset}"
        popdep_out = merge_pfile.multiMap { id, chr, prefix, pgen, psam, pvar ->
            def d = params.process_dir
            popdep: [ id, chr, file("${d}/variant/${id}.popdep.info.tsv") ]
        }
    } else {
        log.info "${params.c_green}Annotating variant popdep from frozen main_raw grids (${params.popdep_dir})${params.c_reset}"
        popdep_out = annotate_subgenome_variant_popdep(merge_pfile)
    }

    emit:
    vcf = params.process_dir ? channel.empty() : preprocess_out.gz_vcf
    merged_bfile = merge_bfile
    merged_pfile = merge_pfile
    smiss = basic_info_out.smiss
    vmiss = basic_info_out.vmiss
    scount = basic_info_out.scount
    gcount = basic_info_out.gcount
    afreq = basic_info_out.afreq
    hardy = basic_info_out.hardy
    ld = ld_out.ld
    ld_cross = ld_cross_out.ld
    mq = mq_out.mq
    popdep = popdep_out.popdep
}

workflow test_plink_camp {
   take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in
    camp_vmap4_map_tsv

    main:
    def reuseMerged = hasMergedSubgenomeTestPfiles(params.process_dir)
    def merge_pfile
    def merge_bfile

    if (reuseMerged) {
        log.info "${params.c_green}Reuse merged test pfiles from:${params.c_reset} ${params.process_dir} (skip per-chr thin + merge)"
        merge_pfile = Channel.from(listMergedSubgenomeTestPfileTuples(params.process_dir))
        merge_bfile = Channel.from(listMergedSubgenomeTestBfileTuples(params.process_dir))
        log.info "${params.c_green}Making basic info (CAMP) from reused merged pfiles.${params.c_reset}"
        basic_info_out = mk_plink_basic_info_camp_pop_with_filter(merge_pfile, camp_vmap4_map_tsv)
        ld_out = calc_plink_ld_unphased(merge_pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_pfile)
    } else if (params.process_dir) {
        log.info "${params.c_green}Using main_raw pfiles in dir:${params.c_reset} ${params.process_dir}"
        source_pfiles = vcf_in.map { id, chr, _vcf ->
            def prefix = "${params.process_dir}/${id}.plink2"
            def pgen = file("${prefix}.pgen")
            def psam = file("${prefix}.psam")
            def pvar = file("${prefix}.pvar")
            return [ id, chr, prefix, pgen, psam, pvar ]
        }
        subsampling_out = subsampling_pfile_for_test(source_pfiles, params.thin_rate)

        def grouped_pfile = subsampling_out.pfile.map { id, chr, prefix, pgen, psam, pvar ->
            def subgenome = "Others"
            def c = chr.toString()
            if (["1","2","7","8","13","14","19","20","25","26","31","32","37","38"].contains(c)) subgenome = "A"
            else if (["3","4","9","10","15","16","21","22","27","28","33","34","39","40"].contains(c)) subgenome = "B"
            else if (["5","6","11","12","17","18","23","24","29","30","35","36","41","42"].contains(c)) subgenome = "D"
            return [ subgenome, id, chr, prefix, pgen, psam, pvar ]
        }.groupTuple(by: 0)

        def merge_out = merge_subgenome_test_pfile(grouped_pfile)
        merge_pfile = merge_out.pfile
        merge_bfile = merge_out.bfile

        basic_info_out = mk_plink_basic_info_camp_pop_with_filter(merge_out.pfile, camp_vmap4_map_tsv)
        ld_out = calc_plink_ld_unphased(merge_out.pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_out.pfile)
    } else {
        preprocess_out = plink_preprocess(vcf_in)
        log.info "${params.c_green}No process_dir specified, subsampling VCFs normally.${params.c_reset}"
        subsampling_out = subsampling_pfile_for_test(preprocess_out.pfile, params.thin_rate)

        def grouped_pfile = subsampling_out.pfile.map { id, chr, prefix, pgen, psam, pvar ->
            def subgenome = "Others"
            def c = chr.toString()
            if (["1","2","7","8","13","14","19","20","25","26","31","32","37","38"].contains(c)) subgenome = "A"
            else if (["3","4","9","10","15","16","21","22","27","28","33","34","39","40"].contains(c)) subgenome = "B"
            else if (["5","6","11","12","17","18","23","24","29","30","35","36","41","42"].contains(c)) subgenome = "D"
            return [ subgenome, id, chr, prefix, pgen, psam, pvar ]
        }.groupTuple(by: 0)

        def merge_out = merge_subgenome_test_pfile(grouped_pfile)

        basic_info_out = mk_plink_basic_info_camp_pop_with_filter(merge_out.pfile, camp_vmap4_map_tsv)
        ld_out = calc_plink_ld_unphased(merge_out.pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_out.pfile)
        merge_pfile = merge_out.pfile
        merge_bfile = merge_out.bfile
    }

    if (reuseMerged && hasVariantMqForMergedTests(params.process_dir)) {
        log.info "${params.c_green}Reuse variant MQ from ${params.process_dir}/variant${params.c_reset}"
        mq_out = merge_pfile.multiMap { id, chr, prefix, pgen, psam, pvar ->
            def d = params.process_dir
            mq: [ id, chr, file("${d}/variant/${id}.mq.info.tsv") ]
        }
    } else {
        log.info "${params.c_green}Annotating variant MQ from frozen site_mq.ref (${params.mq_dir})${params.c_reset}"
        mq_out = annotate_subgenome_variant_mq(merge_pfile)
    }

    if (reuseMerged && hasVariantPopdepForMergedTests(params.process_dir)) {
        log.info "${params.c_green}Reuse variant popdep from ${params.process_dir}/variant${params.c_reset}"
        popdep_out = merge_pfile.multiMap { id, chr, prefix, pgen, psam, pvar ->
            def d = params.process_dir
            popdep: [ id, chr, file("${d}/variant/${id}.popdep.info.tsv") ]
        }
    } else {
        log.info "${params.c_green}Annotating variant popdep from frozen main_raw grids (${params.popdep_dir})${params.c_reset}"
        popdep_out = annotate_subgenome_variant_popdep(merge_pfile)
    }

    emit:
    vcf = params.process_dir ? channel.empty() : preprocess_out.gz_vcf
    merged_bfile = merge_bfile
    merged_pfile = merge_pfile
    smiss = basic_info_out.smiss
    vmiss = basic_info_out.vmiss
    scount = basic_info_out.scount
    gcount = basic_info_out.gcount
    afreq = basic_info_out.afreq
    hardy = basic_info_out.hardy
    ld = ld_out.ld
    ld_cross = ld_cross_out.ld
    mq = mq_out.mq
    popdep = popdep_out.popdep
} 

workflow test_common_thin_processor {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    def reuseMerged = hasMergedSubgenomeTestPfiles(params.process_dir)
    def merge_pfile
    def merge_bfile

    if (reuseMerged) {
        log.info "${params.c_green}Reuse merged test pfiles from:${params.c_reset} ${params.process_dir} (skip common-thin + merge)"
        merge_pfile = Channel.from(listMergedSubgenomeTestPfileTuples(params.process_dir))
        merge_bfile = Channel.from(listMergedSubgenomeTestBfileTuples(params.process_dir))
        if (hasPlinkBasicInfoForMergedTests(params.process_dir)) {
            log.info "${params.c_green}Reuse PLINK2 basic info from ${params.process_dir}/variant and sample/${params.c_reset}"
            basic_info_out = merge_pfile.multiMap { id, chr, prefix, pgen, psam, pvar ->
                def d = params.process_dir
                smiss: [ id, chr, file("${d}/sample/${id}.info.smiss") ]
                vmiss: [ id, chr, file("${d}/variant/${id}.info.vmiss") ]
                scount: [ id, chr, file("${d}/sample/${id}.info.scount") ]
                gcount: [ id, chr, file("${d}/variant/${id}.info.gcount") ]
                afreq: [ id, chr, file("${d}/variant/${id}.info.afreq") ]
                hardy: [ id, chr, file("${d}/variant/${id}.info.hardy") ]
            }
        } else {
            log.info "${params.c_green}Making basic info from reused merged pfiles.${params.c_reset}"
            basic_info_out = mk_plink_basic_info(merge_pfile)
        }
        ld_out = calc_plink_ld_unphased(merge_pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_pfile)
    } else if (params.process_dir) {
        log.info "${params.c_green}Using main_raw pfiles in dir:${params.c_reset} ${params.process_dir}"
        source_pfiles = vcf_in.map { id, chr, _vcf ->
            def prefix = "${params.process_dir}/${id}.plink2"
            def pgen = file("${prefix}.pgen")
            def psam = file("${prefix}.psam")
            def pvar = file("${prefix}.pvar")
            return [ id, chr, prefix, pgen, psam, pvar ]
        }
        log.info "${params.c_green}Building common-thin test pfiles with --geno ${params.hf_geno} --maf ${params.hf_maf} --thin ${params.hf_thin_rate}${params.c_reset}"
        subsampling_out = subsampling_common_variant_pfile_for_test(source_pfiles)

        def grouped_pfile = subsampling_out.pfile.map { id, chr, prefix, pgen, psam, pvar ->
            def subgenome = "Others"
            def c = chr.toString()
            if (["1","2","7","8","13","14","19","20","25","26","31","32","37","38"].contains(c)) subgenome = "A"
            else if (["3","4","9","10","15","16","21","22","27","28","33","34","39","40"].contains(c)) subgenome = "B"
            else if (["5","6","11","12","17","18","23","24","29","30","35","36","41","42"].contains(c)) subgenome = "D"
            return [ subgenome, id, chr, prefix, pgen, psam, pvar ]
        }.groupTuple(by: 0)

        def merge_out = merge_subgenome_test_pfile(grouped_pfile)
        merge_pfile = merge_out.pfile
        merge_bfile = merge_out.bfile

        basic_info_out = mk_plink_basic_info(merge_out.pfile)
        ld_out = calc_plink_ld_unphased(merge_out.pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_out.pfile)
    } else {
        preprocess_out = plink_preprocess(vcf_in)
        log.info "${params.c_green}Building common-thin test pfiles with --geno ${params.hf_geno} --maf ${params.hf_maf} --thin ${params.hf_thin_rate}${params.c_reset}"
        subsampling_out = subsampling_common_variant_pfile_for_test(preprocess_out.pfile)

        def grouped_pfile = subsampling_out.pfile.map { id, chr, prefix, pgen, psam, pvar ->
            def subgenome = "Others"
            def c = chr.toString()
            if (["1","2","7","8","13","14","19","20","25","26","31","32","37","38"].contains(c)) subgenome = "A"
            else if (["3","4","9","10","15","16","21","22","27","28","33","34","39","40"].contains(c)) subgenome = "B"
            else if (["5","6","11","12","17","18","23","24","29","30","35","36","41","42"].contains(c)) subgenome = "D"
            return [ subgenome, id, chr, prefix, pgen, psam, pvar ]
        }.groupTuple(by: 0)

        def merge_out = merge_subgenome_test_pfile(grouped_pfile)

        basic_info_out = mk_plink_basic_info(merge_out.pfile)
        ld_out = calc_plink_ld_unphased(merge_out.pfile)
        ld_cross_out = calc_plink_ld_crosschr_random(merge_out.pfile)
        merge_pfile = merge_out.pfile
        merge_bfile = merge_out.bfile
    }

    if (reuseMerged && hasVariantMqForMergedTests(params.process_dir)) {
        log.info "${params.c_green}Reuse variant MQ from ${params.process_dir}/variant${params.c_reset}"
        mq_out = merge_pfile.multiMap { id, chr, prefix, pgen, psam, pvar ->
            def d = params.process_dir
            mq: [ id, chr, file("${d}/variant/${id}.mq.info.tsv") ]
        }
    } else {
        log.info "${params.c_green}Annotating variant MQ from frozen site_mq.ref (${params.mq_dir})${params.c_reset}"
        mq_out = annotate_subgenome_variant_mq(merge_pfile)
    }

    if (reuseMerged && hasVariantPopdepForMergedTests(params.process_dir)) {
        log.info "${params.c_green}Reuse variant popdep from ${params.process_dir}/variant${params.c_reset}"
        popdep_out = merge_pfile.multiMap { id, chr, prefix, pgen, psam, pvar ->
            def d = params.process_dir
            popdep: [ id, chr, file("${d}/variant/${id}.popdep.info.tsv") ]
        }
    } else {
        log.info "${params.c_green}Annotating variant popdep from frozen main_raw grids (${params.popdep_dir})${params.c_reset}"
        popdep_out = annotate_subgenome_variant_popdep(merge_pfile)
    }

    emit:
    vcf = params.process_dir ? channel.empty() : preprocess_out.gz_vcf
    merged_bfile = merge_bfile
    merged_pfile = merge_pfile
    smiss = basic_info_out.smiss
    vmiss = basic_info_out.vmiss
    scount = basic_info_out.scount
    gcount = basic_info_out.gcount
    afreq = basic_info_out.afreq
    hardy = basic_info_out.hardy
    ld = ld_out.ld
    ld_cross = ld_cross_out.ld
    mq = mq_out.mq
    popdep = popdep_out.popdep
}

workflow plink_processor {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    // 1. preprocess VCF to PLINK format
    preprocess_out = plink_preprocess(vcf_in)

    // 2. make sample and variant basic information
    if (params.process_dir) {
        log.info "${params.c_green}Making basic info from pfiles in dir:${params.c_reset} ${params.process_dir}"
        basic_info_out = preprocess_out.pfile.multiMap { id, chr, _prefix, _pgen, _psam, _pvar ->
            def smiss = file("${params.process_dir}/sample/${id}.info.smiss")
            def vmiss = file("${params.process_dir}/variant/${id}.info.vmiss")
            def scount = file("${params.process_dir}/sample/${id}.info.scount")
            def gcount = file("${params.process_dir}/variant/${id}.info.gcount")
            def afreq = file("${params.process_dir}/variant/${id}.info.afreq")
            def hardy = file("${params.process_dir}/variant/${id}.info.hardy")
            smiss: [ id, chr, smiss ]
            vmiss: [ id, chr, vmiss ]
            scount: [ id, chr, scount ]
            gcount: [ id, chr, gcount ]
            afreq: [ id, chr, afreq ]
            hardy: [ id, chr, hardy ]
        }
        popdep_out = preprocess_out.gz_vcf.multiMap { id, chr, _vcf ->
            def root = params.popdep_dir ?: params.process_dir
            def popdep = file("${root}/variant/${id}.popdep.txt")
            popdep: [ id, chr, popdep ]
        }
    } else {
        log.info "${params.c_green}Making basic info normally.${params.c_reset}"
        basic_info_out = mk_plink_basic_info(preprocess_out.pfile)

        // 3. calculate population depth using TIGER
        def pd_config = getTigerJarConfig(params.popdep_tiger_jar, params.home_dir, params.popdep_tiger_app)
        ch_tiger_config = channel.value(tuple(pd_config.path, pd_config.app_name, pd_config.java_version))
        popdep_gz_out = calc_population_depth(preprocess_out.gz_vcf, ch_tiger_config)
        popdep_out = popdep_tiger_gz_to_bgzip_tabix(popdep_gz_out.tiger_gz)
    }


    emit:
    vcf = preprocess_out.gz_vcf
    plink_bfile = preprocess_out.bfile
    plink_pfile  = preprocess_out.pfile
    smiss = basic_info_out.smiss
    vmiss = basic_info_out.vmiss
    scount = basic_info_out.scount
    gcount = basic_info_out.gcount
    afreq = basic_info_out.afreq
    hardy = basic_info_out.hardy
    popdep = popdep_out.popdep
}

workflow plink_preprocess {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    if (params.process_dir) {
        log.info "${params.c_green}Processed VCFs in dir:${params.c_reset} ${params.process_dir}"
        vcf = vcf_in.map { id, chr, _vcf ->
            def out_vcf = file("${params.process_dir}/${id}.vcf.gz")
            return [ id, chr, out_vcf ]
        }
        pfiles = vcf_in.map { id, chr, _vcf ->
            def prefix = "${params.process_dir}/${id}.plink2"
            def pgen = file("${prefix}.pgen")
            def psam = file("${prefix}.psam")
            def pvar = file("${prefix}.pvar")
            return [ id, chr, prefix, pgen, psam, pvar ]
        }
        bfiles = vcf_in.map { id, chr, _vcf ->
            def prefix = "${params.process_dir}/${id}.plink"
            def bed = file("${prefix}.bed")
            def bim = file("${prefix}.bim")
            def fam = file("${prefix}.fam")
            return [ id, chr, prefix, bed, bim, fam ]
        }
        gz_vcf = [ vcf: vcf ]
        plink_out = [
            bfile: bfiles,
            pfile: pfiles
        ]
    } else {
        gz_vcf = format_vcf_bgzip(vcf_in)
        plink_out = format_vcf_plink(gz_vcf.vcf)
    }

    emit:
    gz_vcf = gz_vcf.vcf
    bfile = plink_out.bfile
    pfile = plink_out.pfile
}

workflow vcf_arrange_merge {
    take:
    // Expect a channel: [ val(id), val(chr), path(vcf) ]
    vcf_in

    main:
    arrange_out = arrange_vcf_wheat_chr_by_awk(vcf_in)
    arrange_ch = arrange_out.vcf.groupTuple(by: 2)
    merge_out = merge_arranged_vcf(arrange_ch)

    emit:
    vcf = merge_out.vcf
}

workflow vcf_bcftools_filter {
    take:
    vcf_in

    main:
    filter_out = filter_vcf_v0_bcftools(vcf_in)

    emit:
    vcf = filter_out.vcf
}
