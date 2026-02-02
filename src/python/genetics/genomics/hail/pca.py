import hail as hl
import argparse
import sys
import os

def run_pca(vcf_path, output_prefix, k=10, reference=None):
    # Initialize Hail
    # Use a local tmp dir to avoid permission issues on shared clusters if needed
    hl.init(default_reference=reference if reference else 'GRCh37', log=f'{output_prefix}.hail.log')

    print(f"Importing VCF: {vcf_path}")
    if reference:
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=reference)
    else:
        # If no reference is provided, we try to infer or use None. 
        # For non-human organisms, it's often safer to not specify a reference 
        # or create one from the VCF header if Hail supports it, 
        # but import_vcf requires a reference_genome argument usually unless we define one.
        # If we pass None, Hail treats locus as a struct {contig, position} rather than a Locus object.
        mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome=None)

    print("Running PCA...")
    # hwe_normalized_pca works on the GT entry
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=k, compute_loadings=True)

    print(f"Exporting results to {output_prefix}...")
    # Export Scores (Sample PCs)
    # scores is a Table with key 's' (sample) and field 'scores' (array of float64)
    scores = scores.transmute(**{f'PC{i+1}': scores.scores[i] for i in range(k)})
    scores.export(f"{output_prefix}.scores.tsv")

    # Export Loadings (Variant PCs)
    # loadings is a Table with key 'locus', 'alleles' and field 'loadings'
    # If reference was None, locus is a struct. We flatten it.
    if reference is None:
        loadings = loadings.key_by() # Unkey to access fields
        # Assuming locus has contig and position
        loadings = loadings.transmute(
            contig = loadings.locus.contig,
            position = loadings.locus.position,
            alleles = hl.str(loadings.alleles),
            **{f'PC{i+1}': loadings.loadings[i] for i in range(k)}
        )
    else:
        loadings = loadings.transmute(**{f'PC{i+1}': loadings.loadings[i] for i in range(k)})
    
    loadings.export(f"{output_prefix}.loadings.tsv")

    # Export Eigenvalues
    with open(f"{output_prefix}.eigenvalues.txt", 'w') as f:
        for val in eigenvalues:
            f.write(f"{val}\n")
    
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run PCA using Hail")
    parser.add_argument("--vcf", required=True, help="Input VCF file (bgzipped)")
    parser.add_argument("--out", required=True, help="Output prefix")
    parser.add_argument("--k", type=int, default=10, help="Number of PCs")
    parser.add_argument("--reference", default=None, help="Reference genome (e.g., GRCh37, GRCh38). Leave empty for non-standard.")
    
    args = parser.parse_args()
    
    run_pca(args.vcf, args.out, args.k, args.reference)
