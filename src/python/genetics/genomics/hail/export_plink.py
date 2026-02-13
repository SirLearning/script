import hail as hl
import argparse
import os
import sys
import inspect

def run_export_plink(vcf_path, output_prefix, prune=True, r2=0.1, reference=None):
    hl.init(default_reference=reference if reference else 'GRCh37', log=f'{output_prefix}.hail.log')
    
    print(f"Importing VCF: {vcf_path}")
    if reference:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=reference)
    else:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=None)
    
    if prune:
        print(f"Running LD pruning (r2={r2})...")
        # LD pruning requires variant QC to be run first usually? No, ld_prune takes mt.
        # But usually we want to filter to common variants first.
        # Assuming input is already filtered for MAF/Missingness as per previous steps.
        pruned_variant_table = hl.ld_prune(mt.GT, r2=r2)
        mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))
        print(f"Retained {mt.count_rows()} variants after pruning.")

    print(f"Exporting to PLINK: {output_prefix}")
    hl.export_plink(mt, output_prefix, ind_id=mt.s, fam_id=mt.s) 
    # Note: Hail export_plink creates .bed, .bim, .fam. 
    # output_prefix should not have extension.

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--no_prune", action="store_true")
    parser.add_argument("--r2", type=float, default=0.1)
    parser.add_argument("--reference")
    args = parser.parse_args()
    
    run_export_plink(args.vcf, args.out, prune=not args.no_prune, r2=args.r2, reference=args.reference)
