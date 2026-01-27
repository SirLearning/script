import hail as hl
import argparse
import os

def run_vcf_to_mt(vcf_path, output_path, reference=None):
    hl.init(default_reference=reference if reference else 'GRCh37', log=f'{output_path}.hail.log')
    
    print(f"Importing VCF: {vcf_path}")
    if reference:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=reference)
    else:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=None)
    
    print(f"Writing MatrixTable to {output_path}")
    mt.write(output_path, overwrite=True)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--reference")
    args = parser.parse_args()
    
    run_vcf_to_mt(args.vcf, args.out, args.reference)
