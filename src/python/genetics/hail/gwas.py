import hail as hl
import argparse
import sys

def run_gwas(vcf_path, pheno_path, output_prefix, response_col, id_col='IID', covar_path=None, covar_cols=None, reference=None):
    hl.init(default_reference=reference if reference else 'GRCh37', log=f'{output_prefix}.hail.log')
    
    print(f"Importing VCF: {vcf_path}")
    if reference:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=reference)
    else:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=None)
    
    print(f"Importing Phenotypes from {pheno_path}")
    # Import phenotype table
    # We assume the file is delimited (tsv/csv). import_table defaults to tsv.
    pheno = hl.import_table(pheno_path, impute=True, key=id_col)
    
    # Annotate cols (samples) with phenotype
    mt = mt.annotate_cols(pheno = pheno[mt.s])
    
    # Prepare covariates
    covariates = [1.0] # Intercept
    
    if covar_path:
        print(f"Importing Covariates from {covar_path}")
        cov = hl.import_table(covar_path, impute=True, key=id_col)
        mt = mt.annotate_cols(cov = cov[mt.s])
        
        if covar_cols:
            for c in covar_cols.split(','):
                c = c.strip()
                if c:
                    print(f"Adding covariate: {c}")
                    covariates.append(mt.cov[c])
    
    print(f"Running Linear Regression for trait: {response_col}")
    # Check if trait exists
    # We can't easily check existence before running in lazy mode without triggering computation, 
    # but Hail will error if column missing.
    
    # Run GWAS
    # y: response variable
    # x: genotype dosage (n_alt_alleles)
    gwas = hl.linear_regression_rows(
        y=mt.pheno[response_col],
        x=mt.GT.n_alt_alleles(),
        covariates=covariates
    )
    
    print(f"Exporting results to {output_prefix}.gwas.tsv")
    gwas.export(f"{output_prefix}.gwas.tsv")
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GWAS using Hail")
    parser.add_argument("--vcf", required=True, help="Input VCF file (bgzipped)")
    parser.add_argument("--pheno", required=True, help="Phenotype file path")
    parser.add_argument("--out", required=True, help="Output prefix")
    parser.add_argument("--response", required=True, help="Response column name (trait)")
    parser.add_argument("--id_col", default="IID", help="Sample ID column name in pheno/covar files")
    parser.add_argument("--covar", help="Covariate file path (optional)")
    parser.add_argument("--covar_cols", help="Comma-separated list of covariate column names")
    parser.add_argument("--reference", default=None, help="Reference genome")
    
    args = parser.parse_args()
    
    run_gwas(
        vcf_path=args.vcf,
        pheno_path=args.pheno,
        output_prefix=args.out,
        response_col=args.response,
        id_col=args.id_col,
        covar_path=args.covar,
        covar_cols=args.covar_cols,
        reference=args.reference
    )
