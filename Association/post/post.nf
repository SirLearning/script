nextflow.enable.dsl = 2

process COLLECT_RESULTS {
    tag "collect results"
    publishDir params.outdir + "/post/${params.trait}", mode: 'copy'

    input:
    val(trigger)

    output:
    path 'summary.tsv'
    path 'plot.manhattan.png', optional: true
    path 'plot.qq.png', optional: true
    path 'annotated.tsv', optional: true

    script:
    def gff = params.gff ?: ''
    """
    set -euo pipefail
    # Build unified summary
    : > summary.tsv
    if ls ${params.outdir}/results/${params.trait}/**/RESULTS.tsv >/dev/null 2>&1; then
      head -n1 $(ls ${params.outdir}/results/${params.trait}/**/RESULTS.tsv | head -n1) > summary.tsv
      awk 'FNR>1' $(ls ${params.outdir}/results/${params.trait}/**/RESULTS.tsv) >> summary.tsv || true
  elif ls ${params.outdir}/results/${params.trait}/**/plink.glm.linear >/dev/null 2>&1; then
      awk 'NR==1 || $0!~"^#"{print $3"\t"$1"\t"$4"\t"$12}' OFS='\t' $(ls ${params.outdir}/results/${params.trait}/**/plink.glm.linear) > summary.tsv || true
    elif ls ${params.outdir}/results/${params.trait}/**/plink.glm.logistic >/dev/null 2>&1; then
      awk 'NR==1 || $0!~"^#"{print $3"\t"$1"\t"$4"\t"$9}' OFS='\t' $(ls ${params.outdir}/results/${params.trait}/**/plink.glm.logistic) > summary.tsv || true
    fi

    # Plots
    if [ -s summary.tsv ]; then
      Rscript ${projectDir}/src/r/gwas/plot_manhattan_qq.R --input summary.tsv --outprefix plot || true
    fi

    # Optional: annotate with GFF gene features if provided
    if [ -n "${gff}" ] && [ -s "${gff}" ] && [ -s summary.tsv ]; then
      awk 'BEGIN{OFS="\t"} NR>1{print $2, $3-1, $3, $1, $4}' summary.tsv > snps.bed
      awk 'BEGIN{OFS="\t"} $3=="gene"{split($9,a,";"); id="."; for(i in a){if(a[i]~/(ID|Name)=/){sub(/^(ID|Name)=/,"",a[i]); id=a[i]; break}}; print $1,$4-1,$5,id}' ${gff} > genes.bed
      if command -v bedtools >/dev/null 2>&1; then
        bedtools intersect -wao -a snps.bed -b genes.bed | awk 'BEGIN{OFS="\t"}{print $4,$1,$3,$5}' > annotated.tsv || true
      else
        awk 'BEGIN{OFS="\t"} FNR==NR{g[$1]=g[$1] sprintf("%s\t%s\t%s\t%s\n", $1,$2,$3,$4); next} {chr=$1; p=$3; split(g[chr],rows,"\n"); hit="."; for(i in rows){split(rows[i],x,"\t"); if(p>=x[2] && p<=x[3]){hit=x[4]; break}}; print $4,chr,p,hit}' genes.bed snps.bed > annotated.tsv || true
      fi
    fi
    """
}

workflow POST {
    take:
    ch_results

    main:
    COLLECT_RESULTS(ch_results.collect())

    emit:
    COLLECT_RESULTS.out
}
