import hail as hl
import argparse
import os

def run_qc(vcf_path, output_prefix, reference=None):
    hl.init(default_reference=reference if reference else 'GRCh37', log=f'{output_prefix}.hail.log')
    
    print(f"Importing VCF: {vcf_path}")
    if reference:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=reference)
    else:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=None)
    
    print("Running Variant QC...")
    mt = hl.variant_qc(mt)
    
    print("Running Sample QC...")
    mt = hl.sample_qc(mt)
    
    # Export Variant QC metrics
    # Fields: dp_stats, gq_stats, call_rate, n_called, n_not_called, n_filtered, n_het, n_non_ref, n_hom_var, n_hom_ref, p_value_hwe, p_value_fit, het_freq_hwe, AF, AC, AN, homozygote_count
    # We select some common ones to export to TSV
    rows = mt.rows()
    rows = rows.flatten()
    
    # Handle locus struct if no reference
    if reference is None:
        rows = rows.transmute(
            contig = rows.locus.contig,
            position = rows.locus.position,
            alleles = hl.str(rows.alleles)
        )
        key_cols = ['contig', 'position', 'alleles']
    else:
        rows = rows.transmute(alleles = hl.str(rows.alleles))
        key_cols = ['locus', 'alleles']

    # Select relevant variant QC columns
    # Note: AF is an array (one per allele). We take AF[1] for biallelic assumption or export as string.
    # For simplicity, let's export call_rate, AF, AC, AN, p_value_hwe, dp_stats.mean, gq_stats.mean
    rows = rows.select(
        *key_cols,
        call_rate = rows.variant_qc.call_rate,
        AF = hl.str(rows.variant_qc.AF),
        AC = hl.str(rows.variant_qc.AC),
        AN = rows.variant_qc.AN,
        p_value_hwe = rows.variant_qc.p_value_hwe,
        mean_dp = rows.variant_qc.dp_stats.mean,
        mean_gq = rows.variant_qc.gq_stats.mean,
        n_het = rows.variant_qc.n_het,
        n_hom_var = rows.variant_qc.n_hom_var
    )
    
    print(f"Exporting Variant QC to {output_prefix}.variant_qc.tsv")
    rows.export(f"{output_prefix}.variant_qc.tsv")
    
    # Export Sample QC metrics
    cols = mt.cols()
    cols = cols.flatten()
    
    # Select relevant sample QC columns
    cols = cols.select(
        s = cols.s,
        call_rate = cols.sample_qc.call_rate,
        n_called = cols.sample_qc.n_called,
        n_not_called = cols.sample_qc.n_not_called,
        n_hom_ref = cols.sample_qc.n_hom_ref,
        n_het = cols.sample_qc.n_het,
        n_hom_var = cols.sample_qc.n_hom_var,
        n_snp = cols.sample_qc.n_snp,
        n_insertion = cols.sample_qc.n_insertion,
        n_deletion = cols.sample_qc.n_deletion,
        n_singleton = cols.sample_qc.n_singleton,
        n_transition = cols.sample_qc.n_transition,
        n_transversion = cols.sample_qc.n_transversion,
        r_ti_tv = cols.sample_qc.r_ti_tv,
        dp_mean = cols.sample_qc.dp_stats.mean,
        gq_mean = cols.sample_qc.gq_stats.mean
    )
    
    print(f"Exporting Sample QC to {output_prefix}.sample_qc.tsv")
    cols.export(f"{output_prefix}.sample_qc.tsv")
    
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--reference")
    args = parser.parse_args()
    
    run_qc(args.vcf, args.out, args.reference)
