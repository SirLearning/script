import hail as hl
import argparse
import os

def run_kinship(vcf_path, output_prefix, method='pc_relate', k=10, min_kinship=None, reference=None):
    hl.init(default_reference=reference if reference else 'GRCh37', log=f'{output_prefix}.hail.log')
    
    print(f"Importing VCF: {vcf_path}")
    if reference:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=reference)
    else:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=None)
    
    if method == 'pc_relate':
        print(f"Running PC-Relate with k={k}...")
        # PC-Relate requires PCA first to control for population structure
        # But pc_relate function in Hail takes the number of PCs to use (k) and computes them internally 
        # OR takes pre-computed scores. 
        # Actually, hl.pc_relate takes `k` (int) or `scores` (Table).
        # If k is provided, it runs PCA internally.
        
        # Note: pc_relate requires min_individual_maf to be set if not using pre-computed scores?
        # Let's check docs or assume defaults.
        # It returns a Table with keys (i, j) and field 'kinship'.
        
        rel = hl.pc_relate(mt.GT, min_individual_maf=0.01, k=k, statistics='kinship')
        
        if min_kinship is not None:
            rel = rel.filter(rel.kinship >= min_kinship)
            
        print(f"Exporting Kinship to {output_prefix}.kinship.tsv")
        rel.export(f"{output_prefix}.kinship.tsv")
        
    elif method == 'ibd':
        print("Running Identity By Descent (IBD)...")
        # IBD requires LD pruning usually for best results, but we'll run on input
        # It returns a Table with keys (i, j) and fields Z0, Z1, Z2, PI_HAT
        
        # We might want to filter to common variants first
        mt_common = hl.variant_qc(mt)
        mt_common = mt_common.filter_rows(mt_common.variant_qc.AF[1] > 0.01)
        
        ibd = hl.identity_by_descent(mt_common)
        
        print(f"Exporting IBD to {output_prefix}.ibd.tsv")
        ibd.export(f"{output_prefix}.ibd.tsv")

    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--method", default='pc_relate', choices=['pc_relate', 'ibd'])
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--min_kinship", type=float)
    parser.add_argument("--reference")
    args = parser.parse_args()
    
    run_kinship(args.vcf, args.out, args.method, args.k, args.min_kinship, args.reference)
