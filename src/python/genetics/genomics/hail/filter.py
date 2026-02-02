import hail as hl
import argparse
import os

def run_filter(vcf_path, output_prefix, maf, mac, min_alleles, max_alleles, max_missing, reference=None):
    hl.init(default_reference=reference if reference else 'GRCh37', log=f'{output_prefix}.hail.log')
    
    print(f"Importing VCF: {vcf_path}")
    if reference:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=reference)
    else:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=None)
    
    initial_count = mt.count()
    print(f"Initial count: {initial_count[0]} variants, {initial_count[1]} samples")

    # Annotate with variant QC metrics
    mt = hl.variant_qc(mt)
    
    filters = []
    
    # Max Missing (fraction of missing genotypes allowed)
    # call_rate >= 1 - max_missing
    if max_missing is not None:
        filters.append(mt.variant_qc.call_rate >= (1.0 - max_missing))
        
    # MAF (Minor Allele Frequency)
    # We filter based on the frequency of the minor allele.
    # For biallelic, this is min(AF[0], AF[1]).
    # For multiallelic, it's more complex, but usually we want the 2nd most common allele to be >= MAF?
    # Or just that the site is polymorphic enough.
    # VCFtools --maf: "Include only sites with a Minor Allele Frequency greater than or equal to the --maf value."
    # We will use min(AF) across all alleles (including ref if it's minor? No, usually AF refers to Alt alleles in some contexts, but Hail AF includes Ref).
    # Hail AF is array of frequencies summing to 1.
    # We want the frequency of the second most common allele?
    # Or simply `hl.min(mt.variant_qc.AF) >= maf`?
    # If Ref=0.99, Alt=0.01. Min=0.01. If MAF=0.05, 0.01 < 0.05 -> Filter out. Correct.
    # If Ref=0.01, Alt=0.99. Min=0.01. Filter out. Correct.
    if maf is not None:
        filters.append(hl.min(mt.variant_qc.AF) >= maf)

    # MAC (Minor Allele Count)
    if mac is not None:
        filters.append(hl.min(mt.variant_qc.AC) >= mac)

    # Alleles
    if min_alleles is not None:
        filters.append(hl.len(mt.alleles) >= min_alleles)
    if max_alleles is not None:
        filters.append(hl.len(mt.alleles) <= max_alleles)

    # Apply filters
    if filters:
        combined_filter = filters[0]
        for f in filters[1:]:
            combined_filter = combined_filter & f
        mt = mt.filter_rows(combined_filter)

    final_count = mt.count()
    print(f"Final count: {final_count[0]} variants, {final_count[1]} samples")

    output_file = f"{output_prefix}.vcf.bgz"
    print(f"Exporting to {output_file}")
    hl.export_vcf(mt, output_file)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--maf", type=float)
    parser.add_argument("--mac", type=int)
    parser.add_argument("--min_alleles", type=int)
    parser.add_argument("--max_alleles", type=int)
    parser.add_argument("--max_missing", type=float)
    parser.add_argument("--reference")
    args = parser.parse_args()
    
    run_filter(args.vcf, args.out, args.maf, args.mac, args.min_alleles, args.max_alleles, args.max_missing, args.reference)
