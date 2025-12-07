nextflow.enable.dsl=2

include { HAIL_PCA } from './hail_pca.nf'

workflow population_structure {
    take:
    vcf_in
    config

    main:
    
    plink_bed = Channel.empty()
    pca_results = Channel.empty()

    if (params.tool == 'hail') {
        // Hail PCA
        HAIL_PCA(vcf_in)
        pca_results = HAIL_PCA.out.scores
        
        // Hail Export to PLINK (with LD pruning) for Admixture
        plink_bed = hail_export_plink(vcf_in).bed
        
    } else {
        // Run PCA for Population Structure Analysis
        pca_results = run_pca(vcf_in).pca_results
        
        // Run Admixture Analysis
        // 1. Convert VCF to PLINK BED
        plink_bed = vcf_to_plink(vcf_in).bed
    }
    
    // Common Admixture Pipeline
    // 2. Run Admixture for K=2 to K=10 (or params.max_k)
    // We create a channel of Ks
    def max_k = params.max_k ?: 5
    k_range = Channel.from(2..max_k)
    
    // Combine bed files with K range
    admixture_input = plink_bed.combine(k_range)
    
    // Run Admixture
    admixture_results = run_admixture(admixture_input)
    
    // Collect all Q files for a sample to plot
    // Group by meta (sample id)
    q_files = admixture_results.q_file.groupTuple()
    
    // 3. Plot Structure
    plot_structure(q_files)
    
    // 4. Additional Visualizations (if metadata provided)
    if (params.metadata) {
        metadata_ch = Channel.fromPath(params.metadata)
        
        // Allele Frequency Trend by Latitude
        run_allele_trend(vcf_in, metadata_ch)
        
        // Geographic Map of Samples
        run_geo_map(metadata_ch)
    }
    
    // 5. XP-CLR Gene Analysis (Optional)
    // Requires params.xpclr_genes file: GeneName, Chr, Pos, Threshold(optional), Centromere(optional)
    if (params.xpclr_genes && params.xpclr_results_dir) {
        gene_info = Channel.fromPath(params.xpclr_genes)
            .splitCsv(header: true, sep: '\t')
        
        // We need to find the corresponding XP-CLR output file for each gene
        // This assumes a naming convention or that we pass the file path in the gene info
        
        run_gene_xpclr_plot(gene_info)
    }
    
    // 6. Sample Clustering based on Environmental Data (Optional)
    if (params.env_data) {
        env_ch = Channel.fromPath(params.env_data)
        run_sample_cluster(env_ch)
    }

    emit:
    pca = pca_results
}

process run_sample_cluster {
    tag "SampleCluster"
    publishDir "${params.outdir}/genotype/population_structure/cluster", mode: 'copy'

    input:
    path env_data

    output:
    path "*.pdf"
    path "*.txt"

    script:
    def k = params.cluster_k ?: 5
    """
    Rscript ${params.src_dir}/r/genetics/dynamic/cluster/sample_cluster_env.R \\
        --input ${env_data} \\
        --k ${k} \\
        --out sample_cluster_K${k}
    """
}

process run_gene_xpclr_plot {
    tag "XPCLR_Plot: ${gene.GeneName}"
    publishDir "${params.outdir}/genotype/selection/xpclr_plots", mode: 'copy'

    input:
    val gene

    output:
    path "*.pdf"

    script:
    def xpclr_file = "${params.xpclr_results_dir}/${gene.FileName}"
    def threshold_arg = gene.Threshold ? "--threshold ${gene.Threshold}" : ""
    def centromere_arg = gene.Centromere ? "--centromere ${gene.Centromere}" : ""
    """
    Rscript ${params.src_dir}/r/genetics/dynamic/xpclr/plot_gene_xpclr.R \\
        --xpclr_file ${xpclr_file} \\
        --gene_pos ${gene.Pos} \\
        --gene_name "${gene.GeneName}" \\
        --chr_name "${gene.Chr}" \\
        ${threshold_arg} \\
        ${centromere_arg} \\
        --out ${gene.GeneName}_xpclr.pdf
    """
}

process run_geo_map {
    tag "GeoMap"
    publishDir "${params.outdir}/genotype/population_structure/maps", mode: 'copy'

    input:
    path metadata

    output:
    path "*.pdf"

    script:
    """
    Rscript ${params.src_dir}/r/genetics/dynamic/structure/plot_geo_map.R \\
        --meta ${metadata} \\
        --out sample_distribution_map \\
        --group_col Region
    """
}

process hail_export_plink {
    tag "HailExportPLINK: ${meta.id}"
    publishDir "${params.outdir}/genotype/process", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.bed"), path("${meta.id}.bim"), path("${meta.id}.fam"), emit: bed

    script:
    def prefix = meta.id
    def ref = params.reference_genome ?: ""
    def ref_arg = ref ? "--reference ${ref}" : ""
    """
    python ${params.src_dir}/python/genetics/hail/export_plink.py \\
        --vcf ${vcf} \\
        --out ${prefix} \\
        ${ref_arg}
    """
}

process vcf_to_plink {
    tag "VCF2PLINK: ${meta.id}"
    
    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.bed"), path("${meta.id}.bim"), path("${meta.id}.fam"), emit: bed

    script:
    def prefix = meta.id
    """
    plink2 --vcf ${vcf} --make-bed --out ${prefix} --allow-extra-chr
    """
}

process run_admixture {
    tag "Admixture: ${meta.id} K=${k}"
    publishDir "${params.outdir}/genotype/population_structure/admixture/${meta.id}/K${k}", mode: 'copy'

    input:
    tuple val(meta), path(bed), path(bim), path(fam), val(k)

    output:
    tuple val(meta), path("*.Q"), emit: q_file
    path "*.P"
    path "*.log"

    script:
    def prefix = meta.id
    """
    # Admixture requires .bed .bim .fam
    # It expects the input filename without extension
    # We need to make sure files are named consistently
    
    # Admixture output is named {prefix}.{K}.Q
    admixture --cv ${bed} ${k} > ${prefix}.K${k}.log
    """
}

process plot_structure {
    tag "PlotStructure: ${meta.id}"
    publishDir "${params.outdir}/genotype/population_structure/plots", mode: 'copy'

    input:
    tuple val(meta), path(q_files)

    output:
    path "*.pdf"
    path "*.txt"

    script:
    def prefix = meta.id
    """
    # Create a directory for Q files to pass to R script
    mkdir q_files
    cp ${q_files} q_files/
    
    Rscript ${params.src_dir}/r/genetics/dynamic/structure/plot_structure.R \\
        --input_dir q_files \\
        --pattern "*.Q" \\
        --output_dir . \\
        --prefix ${prefix}_structure
    """
}

process run_allele_trend {
    tag "AlleleTrend: ${meta.id}"
    publishDir "${params.outdir}/genotype/population_structure/allele_trend", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    path metadata

    output:
    path "*.pdf"

    script:
    def prefix = meta.id
    """
    # Convert VCF to 012 format
    vcftools --gzvcf ${vcf} --012 --out ${prefix}
    
    # Run R script
    Rscript ${params.src_dir}/r/genetics/dynamic/structure/plot_allele_latitude.R \\
        --geno ${prefix}.012 \\
        --indv ${prefix}.012.indv \\
        --pos ${prefix}.012.pos \\
        --meta ${metadata} \\
        --out ${prefix}_allele_trend
    """
}

process run_pca {
    tag "PCA: ${meta.id}"
    publishDir "${params.outdir}/genotype/population_structure", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.eigenvec"), emit: pca_results
    path "*.log"

    script:
    def prefix = meta.id
    """
    echo "Running PCA for population structure analysis on ${vcf}" > ${prefix}.log

    # Use plink2 for PCA
    plink2 \\
        --vcf ${vcf} \\
        --pca ${params.pca_k} \\
        --out ${prefix}

    echo "PCA complete. Results are in ${prefix}.eigenvec" >> ${prefix}.log
    """
}
