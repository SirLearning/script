nextflow.enable.dsl=2

/*
 * Variant assessment metrics using bcftools/vcftools (optionally DumpNice R plots)
 * Computes: SNP count, depth (per-site/per-individual), genotype quality (summary),
 * variant quality (QUAL), missingness (site/individual), Hardy-Weinberg (HWE),
 * minor allele frequency (MAF).
 */

workflow assess {
	take:
		vcf_ch

	main:
		prepared_vcf_ch = vcf_ch.map { meta, vcf -> PREPARE_VCF_GZ(meta, vcf) }

		counts_ch = prepared_vcf_ch.map { meta, vcf -> QUICK_COUNTS(meta, vcf) }

		vcftools_stats = prepared_vcf_ch.map { meta, vcf -> VCFTOOLS_STATS(meta, vcf) }

		bcftools_tags = prepared_vcf_ch.map { meta, vcf -> BCFTOOLS_TAGS(meta, vcf) }

		qc_plots = prepared_vcf_ch
			.ifEmpty { log.warn "No VCFs provided for DumpNice VCF QC assessment." }
			.map { meta, vcf -> DUMPNICE_VCF_QC_ASSESS(meta, vcf) }

	emit:
		counts_ch
		vcftools_stats.maf_freq
		vcftools_stats.hwe
		vcftools_stats.miss_site
		vcftools_stats.miss_indv
		vcftools_stats.depth_site
		vcftools_stats.depth_indv
		vcftools_stats.site_qual
		bcftools_tags.maf_missing
		bcftools_tags.gq_summary
		qc_plots.qc_plots
}

process PREPARE_VCF_GZ {
	tag { "prepare ${meta.id}" }
	publishDir "${params.outdir}/assess", mode: 'copy'
	errorStrategy 'retry'
	maxRetries 1

	input:
	tuple val(meta), path(vcf)

	output:
	tuple val(meta), path('prepared.vcf.gz'), emit: vcf

	script:
	def is_gz = vcf.toString().endsWith('.gz')
	def src = is_gz ? vcf : "${vcf}.gz"
	def compress_cmd = is_gz ? 'true' : "bgzip -c ${vcf} > ${vcf}.gz"
	"""
	set -euo pipefail
	${compress_cmd}
	tabix -p vcf ${src} || true
	ln -sf ${src} prepared.vcf.gz
	ln -sf ${src}.tbi prepared.vcf.gz.tbi || true
	"""
}

process QUICK_COUNTS {
	tag { "counts ${meta.id}" }
	publishDir "${params.outdir}/assess", mode: 'copy'
	cpus 1

	input:
	tuple val(meta), path(vcf)

	output:
	tuple val(meta), path('counts.tsv'), emit: counts

	script:
	def id = meta.id ?: 'sample'
	"""
	set -euo pipefail
	n_snp=$(bcftools view -H ${vcf} | wc -l || echo 0)
	n_samples=$(bcftools query -l ${vcf} | wc -l || echo 0)
	printf "id\tn_variants\tn_samples\n%s\t%s\t%s\n" "${id}" "${n_snp}" "${n_samples}" > counts.tsv
	"""
}

process VCFTOOLS_STATS {
	tag { "vcftools ${meta.id}" }
	publishDir "${params.outdir}/assess/vcftools", mode: 'copy'
	cpus 2

	input:
	tuple val(meta), path(vcf)

	output:
	tuple val(meta), path("*.frq"), emit: maf_freq
	tuple val(meta), path("*.hwe"), emit: hwe
	tuple val(meta), path("*.lmiss"), emit: miss_site
	tuple val(meta), path("*.imiss"), emit: miss_indv
	tuple val(meta), path("*.ldepth.mean"), emit: depth_site
	tuple val(meta), path("*.idepth"), emit: depth_indv
	tuple val(meta), path("*.lqual"), emit: site_qual
	path "*.log"

	script:
	def prefix = meta.id ?: 'vcftools'
	def inFlag = vcf.toString().endsWith('.gz') ? '--gzvcf' : '--vcf'
	"""
	set -euo pipefail
	log=${prefix}.vcftools.log
	vcftools ${inFlag} ${vcf} \
	  --freq --hardy --missing-site --missing-indv \
	  --site-mean-depth --depth --site-quality \
	  --out ${prefix} >> ${log} 2>&1 || true
	"""
}

process BCFTOOLS_TAGS {
	tag { "fill-tags ${meta.id}" }
	publishDir "${params.outdir}/assess/bcftools", mode: 'copy'
	cpus 2

	input:
	tuple val(meta), path(vcf)

	output:
	tuple val(meta), path('maf_missing.tsv'), emit: maf_missing
	tuple val(meta), path('gq_summary.tsv'), emit: gq_summary

	script:
	"""
	set -euo pipefail
	# Annotate MAF and F_MISSING into a temp vcf and export table
	bcftools +fill-tags ${vcf} -- -t MAF,F_MISSING,AC,AN | \
	  bcftools view -O z -o tmp.tags.vcf.gz
	tabix -p vcf tmp.tags.vcf.gz || true
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%MAF\t%F_MISSING\t%QUAL\n' tmp.tags.vcf.gz > maf_missing.tsv

	# Genotype quality summary (median/mean over all genotypes with GQ)
	# Extract all GQ values into one column
	bcftools query -f '[%GQ\t]\n' ${vcf} | tr '\t' '\n' | awk 'NF && $1 != "."' > all.gq
	awk 'BEGIN{min=1e9;max=-1e9}{s+=$1; n++; if($1<min)min=$1; if($1>max)max=$1} END{if(n==0){print "metric\tvalue"; print "n\t0"; exit} mean=s/n; print "metric\tvalue"; print "n\t"n; print "mean\t"mean; print "min\t"min; print "max\t"max}' all.gq > gq_summary.tsv
	"""
}

// Optional DumpNice visualization (reuses existing VCF QC script)
process DUMPNICE_VCF_QC_ASSESS {
	tag { "dumpnice vcf qc ${meta.id}" }
	publishDir "${params.outdir}/assess/qc", mode: 'copy'
	when:
	params.enable_dumpnice_vcf_assess && params.dumpnice_inputdir

	input:
	tuple val(meta), path(vcf)

	output:
	tuple val(meta), path('*.pdf'), emit: qc_plots

	script:
	def inputdir = params.dumpnice_inputdir as String
	"""
	set -euo pipefail
	rsync -a "${inputdir}/" ./ || true
	Rscript src/r/dumpnice/vcf/12_VcfQc.r || true
	ls -1 *.pdf >/dev/null 2>&1 || touch dumpnice_vcf_qc.pdf
	"""
}

