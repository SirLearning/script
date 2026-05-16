nextflow.enable.dsl=2

/*
 * Table/matrix wheat analytics (formerly standalone entry scripts under integrated/<topic>/).
 * Invoked from main.nf when params.mod starts with wheat_.
 * Python: import run_* from genetics.wheat.* inside process scripts (no python -m CLI).
 */

process WHEAT_SNP_QC {
    tag "wheat_snp_qc"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path snp_table
    val output_prefix

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.wheat.snp_qc import run_snp_qc
    run_snp_qc(
        "${snp_table}",
        "${output_prefix}",
        maf=${params.wheat_snp_qc_maf},
        max_missing=${params.wheat_snp_qc_max_missing},
        min_qual=${params.wheat_snp_qc_min_qual},
    )
    """
}

process WHEAT_POPULATION_STRUCTURE {
    tag "wheat_pca_tsne"
    label 'cpus_2'
    cpus 2
    memory '16.GB'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path geno_matrix
    val output_prefix

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.wheat.population_structure import run_population_structure
    run_population_structure(
        "${geno_matrix}",
        "${output_prefix}",
        n_pcs=${params.wheat_pca_n_pcs},
        tsne_perplexity=${params.wheat_pca_tsne_perplexity},
    )
    """
}

process WHEAT_TAGSNP {
    tag "wheat_tagsnp"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"

    input:
    path genotype
    val output_prefix

    output:
    path "*.tsv"

    script:
    """
    #!/usr/bin/env python
    from genetics.wheat.tagsnp import run_tagsnp_selection
    run_tagsnp_selection(
        "${genotype}",
        "${output_prefix}",
        max_tags=${params.wheat_tagsnp_max_tags},
        ld_threshold=${params.wheat_tagsnp_ld_threshold},
    )
    """
}

process WHEAT_HAPMAP {
    tag "wheat_hapmap"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"

    input:
    path genotype
    val output_prefix

    output:
    path "*.tsv"

    script:
    """
    #!/usr/bin/env python
    from genetics.wheat.hapmap import run_hapmap_build
    run_hapmap_build(
        "${genotype}",
        "${output_prefix}",
        window_size=${params.wheat_hapmap_window_size},
    )
    """
}

process WHEAT_CNV {
    tag "wheat_cnv"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path depth_matrix
    val output_prefix

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.wheat.cnv import run_cnv_calling
    run_cnv_calling(
        "${depth_matrix}",
        "${output_prefix}",
        del_z=${params.wheat_cnv_del_z},
        dup_z=${params.wheat_cnv_dup_z},
    )
    """
}

process WHEAT_GENETIC_MAP {
    tag "wheat_genetic_map"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path marker_table
    val output_prefix

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.wheat.genetic_map import run_genetic_map
    run_genetic_map("${marker_table}", "${output_prefix}")
    """
}

process WHEAT_GWAS {
    tag "wheat_gwas"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path genotype
    path phenotype
    val output_prefix

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.wheat.wheat_gwas import run_gwas
    run_gwas(
        "${genotype}",
        "${phenotype}",
        "${output_prefix}",
        trait="${params.wheat_gwas_trait}",
    )
    """
}

process WHEAT_KGWAS {
    tag "wheat_kgwas"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    path kmer_matrix
    path phenotype
    val output_prefix

    output:
    path "*.tsv"
    path "*.png"

    script:
    """
    #!/usr/bin/env python
    from genetics.wheat.kgwas import run_kgwas
    run_kgwas(
        "${kmer_matrix}",
        "${phenotype}",
        "${output_prefix}",
        trait="${params.wheat_kgwas_trait}",
    )
    """
}

workflow integrated_wheat {
    main:
    if (!params.output_dir) {
        error "integrated_wheat: params.output_dir is required."
    }
    if (!params.job) {
        error "integrated_wheat: params.job is required."
    }
    def prefix = "${params.job}.${params.mod}"

    if (params.mod == 'wheat_snp_qc') {
        if (!params.wheat_table_input) {
            error "wheat_snp_qc requires params.wheat_table_input (path to SNP summary table)."
        }
        WHEAT_SNP_QC(file(params.wheat_table_input), prefix)
    } else if (params.mod == 'wheat_pca_tsne') {
        if (!params.wheat_table_input) {
            error "wheat_pca_tsne requires params.wheat_table_input (genotype matrix with Sample column)."
        }
        WHEAT_POPULATION_STRUCTURE(file(params.wheat_table_input), prefix)
    } else if (params.mod == 'wheat_tagsnp') {
        if (!params.wheat_table_input) {
            error "wheat_tagsnp requires params.wheat_table_input."
        }
        WHEAT_TAGSNP(file(params.wheat_table_input), prefix)
    } else if (params.mod == 'wheat_hapmap') {
        if (!params.wheat_table_input) {
            error "wheat_hapmap requires params.wheat_table_input."
        }
        WHEAT_HAPMAP(file(params.wheat_table_input), prefix)
    } else if (params.mod == 'wheat_cnv') {
        if (!params.wheat_table_input) {
            error "wheat_cnv requires params.wheat_table_input (depth matrix)."
        }
        WHEAT_CNV(file(params.wheat_table_input), prefix)
    } else if (params.mod == 'wheat_genetic_map') {
        if (!params.wheat_table_input) {
            error "wheat_genetic_map requires params.wheat_table_input (Marker, CHR, POS, sample columns)."
        }
        WHEAT_GENETIC_MAP(file(params.wheat_table_input), prefix)
    } else if (params.mod == 'wheat_gwas') {
        if (!params.wheat_gwas_genotype || !params.wheat_gwas_phenotype) {
            error "wheat_gwas requires params.wheat_gwas_genotype and params.wheat_gwas_phenotype."
        }
        WHEAT_GWAS(file(params.wheat_gwas_genotype), file(params.wheat_gwas_phenotype), prefix)
    } else if (params.mod == 'wheat_kgwas') {
        if (!params.wheat_kgwas_kmer_matrix || !params.wheat_kgwas_phenotype) {
            error "wheat_kgwas requires params.wheat_kgwas_kmer_matrix and params.wheat_kgwas_phenotype."
        }
        WHEAT_KGWAS(file(params.wheat_kgwas_kmer_matrix), file(params.wheat_kgwas_phenotype), prefix)
    } else {
        error "Unknown wheat mod '${params.mod}'. Expected one of: wheat_snp_qc, wheat_pca_tsne, wheat_tagsnp, wheat_hapmap, wheat_cnv, wheat_genetic_map, wheat_gwas, wheat_kgwas."
    }
}
