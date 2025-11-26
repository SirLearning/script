

process xpclr_smooth {
    tag { "${meta.id}:${chrom}" }
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/smooth_raw", mode: 'copy'

    input:
    tuple val(meta), val(chrom), path(raw_xpclr)

    output:
    tuple val(meta), path("${meta.id}_smooth${chrom}.txt"), emit: smooth_txt
    path("*.pdf"), optional: true, emit: smooth_pdf

    script:
    """
    Rscript ${params.src_dir}/r/genetics/xpclr_smooth.r \
        --input ${raw_xpclr} \
        --output ${meta.id}_smooth${chrom}.txt \
        --plot ${meta.id}_smooth${chrom}.pdf \
        --smoothness 2000 \
        --method 2
    """
}

process merge_xp_clr_results {
    tag { meta.id }
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/smooth", mode: 'copy'

    input:
    tuple val(meta), path(xpclr_files)

    output:
    tuple val(meta), path("${meta.id}_smooth_A.txt"), path("${meta.id}_smooth_B.txt"), path("${meta.id}_smooth_D.txt"), emit: lineage_smooth

    script:
    """
    # A Lineage: Chr 1,2,7,8,13,14,19,20,25,26,31,32,37,38
    for i in 1 2 7 8 13 14 19 20 25 26 31 32 37 38; do
        # Find file matching pattern (assuming standard naming or passed explicitly)
        # The input is a list of files. We need to find the one for chr \${i}
        # Assuming file name contains "smooth\${i}.txt" or similar.
        f=\$(ls *smooth\${i}.txt 2>/dev/null || true)
        if [ -n "\$f" ]; then
            sed '1d' "\$f" | awk -v chr=\$i '{print chr"\t"\$0}'
        fi
    done | sed '/NA/d' | sort -k6,6g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > ${meta.id}_smooth_A.txt

    # B Lineage: Chr 3,4,9,10,15,16,21,22,27,28,33,34,39,40
    for i in 3 4 9 10 15 16 21 22 27 28 33 34 39 40; do
        f=\$(ls *smooth\${i}.txt 2>/dev/null || true)
        if [ -n "\$f" ]; then
            sed '1d' "\$f" | awk -v chr=\$i '{print chr"\t"\$0}'
        fi
    done | sed '/NA/d' | sort -k6,6g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > ${meta.id}_smooth_B.txt

    # D Lineage: Chr 5,6,11,12,17,18,23,24,29,30,35,36,41,42
    for i in 5 6 11 12 17 18 23 24 29 30 35 36 41 42; do
        f=\$(ls *smooth\${i}.txt 2>/dev/null || true)
        if [ -n "\$f" ]; then
            sed '1d' "\$f" | awk -v chr=\$i '{print chr"\t"\$0}'
        fi
    done | sed '/NA/d' | sort -k6,6g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > ${meta.id}_smooth_D.txt
    """
}

process select_top_windows {
    tag { meta.id }
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/top5", mode: 'copy'

    input:
    tuple val(meta), path(smooth_A), path(smooth_B), path(smooth_D)

    output:
    tuple val(meta), path("${meta.id}_A.top5.bed"), path("${meta.id}_B.top5.bed"), path("${meta.id}_D.top5.bed"), emit: top5_bed

    script:
    """
    # Function to extract top 5%
    get_top5() {
        input=\$1
        output=\$2
        lines=\$(wc -l < "\$input")
        # Subtract header
        lines=\$((lines - 1))
        if [ "\$lines" -gt 0 ]; then
            top_n=\$(echo "\$lines * 0.05" | bc | awk '{print int(\$1+0.5)}')
            if [ "\$top_n" -lt 1 ]; then top_n=1; fi
            # Sort by Wstat (col 6) descending? The shell script used `sort -k6,6g` (general numeric) which is ascending?
            # Wait, `tail -n` implies the file is sorted ascending, so tail gets the highest values.
            # The previous process sorted `sort -k6,6g`.
            tail -n "\$top_n" "\$input" > "\$output"
        else
            touch "\$output"
        fi
    }

    get_top5 ${smooth_A} ${meta.id}_A.top5.bed
    get_top5 ${smooth_B} ${meta.id}_B.top5.bed
    get_top5 ${smooth_D} ${meta.id}_D.top5.bed
    """
}

process annotate_and_filter {
    tag { meta.id }
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/genes", mode: 'copy'

    input:
    tuple val(meta), path(bed_A), path(bed_B), path(bed_D)
    path gff
    path known_genes

    output:
    path("${meta.id}_*.cloned.gene"), emit: cloned_genes

    script:
    """
    # Annotate and filter for each lineage
    for bed in ${bed_A} ${bed_B} ${bed_D}; do
        lineage=\$(echo "\$bed" | sed 's/.*_\\([ABD]\\)\\.top5.bed/\\1/')
        
        # Intersect with GFF
        # GFF format assumed: col 1=chr, col 4=start, col 5=end, col 9=attributes
        # BED format: chr, start, end, ...
        
        bedtools intersect -a ${gff} -b "\$bed" -wa | \
            awk '{print \$1"\t"\$4"\t"\$5"\t"\$9}' | \
            awk -F";" '{print \$1}' | \
            awk -F"ID=" '{print \$2}' | \
            sort | uniq > "\$bed.gene"
        
        # Filter known genes
        if [ -s "${known_genes}" ]; then
            grep -w -f "${known_genes}" "\$bed.gene" > "${meta.id}_\${lineage}.top5.cloned.gene" || true
        else
            cp "\$bed.gene" "${meta.id}_\${lineage}.top5.cloned.gene"
        fi
    done
    """
}

process plot_xpclr_heatmap {
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/plot", mode: 'copy'

    input:
    path(cloned_gene_files)

    output:
    path("heatmap_format2.txt")

    script:
    """
    # Prepare input for R script
    # We need to combine all files into a format: File Lineage Gene
    
    for f in ${cloned_gene_files}; do
        # Filename format: {id}_{Lineage}.top5.cloned.gene
        # Extract ID and Lineage
        # Example: EU_South_A.top5.cloned.gene -> ID=EU_South, Lineage=A
        
        base=\$(basename "\$f" .top5.cloned.gene)
        lineage=\${base##*_}
        id=\${base%_*}
        
        awk -v id="\$id" -v lin="\$lineage" '{print id"\t"lin"\t"\$1}' "\$f" >> combined_genes.txt
    done

    Rscript ${params.src_dir}/r/genetics/xpclr_heatmap_prep.r \
        --input combined_genes.txt \
        --output heatmap_format2.txt
    """
}

process extract_gene_vcf {
    tag { gene_id }
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/haplotype_vcf", mode: 'copy'

    input:
    tuple val(gene_id), val(chrom), val(start), val(end)
    path vcf_file
    path vcf_index

    output:
    path("5k.${gene_id}.${chrom}.${start}-${end}.recode.vcf"), emit: gene_vcf

    script:
    """
    # Extract region for gene with 5k buffer (already included in start/end if passed correctly)
    # User script: vcftools --gzvcf \$2 --chr \$chr --from-bp \$from --to-bp \$to ...
    
    vcftools --gzvcf ${vcf_file} \
        --chr ${chrom} \
        --from-bp ${start} \
        --to-bp ${end} \
        --recode \
        --recode-INFO-all \
        --out 5k.${gene_id}.${chrom}.${start}-${end}
    """
}

process overlap_xpclr_bayenv {
    tag { meta.id }
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/bayenv_overlap", mode: 'copy'

    input:
    tuple val(meta), path(xpclr_bed_A), path(xpclr_bed_B), path(xpclr_bed_D)
    tuple path(bayenv_bed_A), path(bayenv_bed_B), path(bayenv_bed_D)

    output:
    tuple val(meta), path("${meta.id}_A.overlap.bed"), path("${meta.id}_B.overlap.bed"), path("${meta.id}_D.overlap.bed"), emit: overlapped_beds

    script:
    """
    # Overlap XP-CLR Top5 with BayEnv Top5
    # User script logic: bedtools intersect -b xpclr -a bayenv -wa
    # Note: User script had 'bayenv_xpclr_site' and 'bayenv_xpclr_region'. 
    # Here we implement one intersection (Region overlap).
    
    bedtools intersect -a ${xpclr_bed_A} -b ${bayenv_bed_A} -wa | sort | uniq > ${meta.id}_A.overlap.bed
    bedtools intersect -a ${xpclr_bed_B} -b ${bayenv_bed_B} -wa | sort | uniq > ${meta.id}_B.overlap.bed
    bedtools intersect -a ${xpclr_bed_D} -b ${bayenv_bed_D} -wa | sort | uniq > ${meta.id}_D.overlap.bed
    """
}

process annotate_overlaps {
    tag { meta.id }
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/bayenv_overlap/genes", mode: 'copy'

    input:
    tuple val(meta), path(bed_A), path(bed_B), path(bed_D)
    path gff
    path known_genes
    path nlr_genes

    output:
    path("${meta.id}_*.cloned.gene"), emit: cloned_genes
    path("${meta.id}_*.nlr.gene"), emit: nlr_genes

    script:
    """
    for bed in ${bed_A} ${bed_B} ${bed_D}; do
        lineage=\$(echo "\$bed" | sed 's/.*_\\([ABD]\\)\\.overlap.bed/\\1/')
        
        # Annotate
        bedtools intersect -a ${gff} -b "\$bed" -wa | \
            awk '{print \$1"\t"\$4"\t"\$5"\t"\$9}' | \
            awk -F";" '{print \$1}' | \
            awk -F"ID=" '{print \$2}' | \
            sort | uniq > "\$bed.gene"
        
        # Filter Cloned Genes
        if [ -s "${known_genes}" ]; then
            grep -w -f "${known_genes}" "\$bed.gene" > "${meta.id}_\${lineage}.top5.cloned.gene" || true
        else
            touch "${meta.id}_\${lineage}.top5.cloned.gene"
        fi

        # Filter NLR Genes
        if [ -s "${nlr_genes}" ]; then
            grep -w -f "${nlr_genes}" "\$bed.gene" > "${meta.id}_\${lineage}.top5.nlr.gene" || true
        else
            touch "${meta.id}_\${lineage}.top5.nlr.gene"
        fi
    done
    """
}

process plot_overlap_heatmap {
    publishDir "${params.output_dir}/${params.job}/population_genetics/xp_clr/bayenv_overlap/plot", mode: 'copy'

    input:
    path(cloned_gene_files)

    output:
    path("heatmap_format2.txt"), emit: matrix
    path("heatmap.pdf"), emit: plot

    script:
    """
    # Prepare input for R script (Reuse logic)
    for f in ${cloned_gene_files}; do
        if [ -s "\$f" ]; then
            base=\$(basename "\$f" .top5.cloned.gene)
            lineage=\${base##*_}
            id=\${base%_*}
            awk -v id="\$id" -v lin="\$lineage" '{print id"\t"lin"\t"\$1}' "\$f" >> combined_genes.txt
        fi
    done

    if [ -f combined_genes.txt ]; then
        Rscript ${params.src_dir}/r/genetics/dynamic/xpclr/xpclr_heatmap_prep.r \
            --input combined_genes.txt \
            --output heatmap_format2.txt
        
        Rscript ${params.src_dir}/r/genetics/dynamic/xpclr/plot_heatmap.r \
            --input heatmap_format2.txt \
            --output heatmap.pdf
    else
        touch heatmap_format2.txt heatmap.pdf
    fi
    """
} 
