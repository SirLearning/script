# Association Pipeline (Nextflow DSL2)

This pipeline performs VCF filtering, GWAS association with multiple models (PLINK GLM/LOGISTIC, GAPIT or rMVP for GLM/MLM/FarmCPU/BLINK), and post-analysis (Manhattan/QQ plots, summary).

## Inputs
- --vcf: Input VCF(.gz)
- --pheno: Phenotype table (CSV/TSV) containing an ID column and at least one trait
- --trait: Trait/column name to analyze
- --covariates: Optional covariates table (CSV/TSV) with matching ID column

## Key params
- --pheno_id_col: ID column name in phenotype/covariates (default IID)
- --filter_tool: bcftools (default) or gatk
- --qual, --maf, --max_miss, --region: filtering thresholds
- --models: comma list of models [GLM,MLM,FarmCPU,Blink]
- --gwas_engine: rmvp (default) or gapit
- --outdir: Output directory (default Association/output)
- --gff: Optional GFF for post-annotation (future extension)
- --pheno_process: true/false. If true, preprocess phenotype and compute BLUE/BLUP and stats.
- --pheno_effect: raw|BLUE|BLUP. Select which phenotype to feed into GWAS when pheno_process=true (default raw).
- --envcol/--repcol: Optional columns in phenotype for environment/replicate.
- --fixed/--random: Additional fixed/random effects in mixed model formula for BLUE/BLUP, e.g. "env+rep" and "IID".

## Software requirements
- bcftools, tabix
- plink (1.9+)
- R >= 4.0 with packages: data.table, qqman, optparse, plink2R and one of:
  - rMVP
  - GAPIT3
- GATK (optional when --filter_tool gatk)

## Run
```bash
nextflow run Association/main.nf \
  --vcf path/to/input.vcf.gz \
  --pheno path/to/pheno.tsv \
  --trait Height \
  --pheno_id_col IID \
  --pheno_process true \
  --pheno_effect BLUE \
  --models GLM,MLM,FarmCPU,Blink \
  --gwas_engine rmvp \
  --outdir results/assoc
```

## Outputs
- outdir/filter: filtered VCF
- outdir/run: PLINK files, PCA, prepared pheno/covar
- outdir/results/<trait>/: per-model results
- outdir/post/<trait>/: summary.tsv, Manhattan and QQ plots

## Notes
- Phenotype and covariate files can be CSV or TSV; the pipeline autodetects by splitting on tab or comma.
- If no covariates provided, top 10 PCs from PLINK are used as covariates automatically.
- For binary traits, PLINK logistic will be attempted if GLM fails.

## Troubleshooting
- Ensure sample IDs in phenotype/covariates match PLINK FID (first column written by --double-id).
- Consider pre-harmonizing IDs to a single column IID across all inputs.
```