"""SnpEff ANN VCF maps onto the VEP-shaped columns coding_mask already reads."""

from pathlib import Path

from genetics.genomics.variant.snpeff import run_snpeff_vcf_to_tab, snpeff_vcf_to_vep_rows


def test_snpeff_vcf_to_vep_rows_takes_first_ann(tmp_path: Path):
    vcf = tmp_path / "sites.ann.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t1-100\tA\tG\t.\t.\t"
        "ANN=G|missense_variant|MODERATE|TraesCS1A02G000100|gene|transcript|"
        "TraesCS1A02G000100.1|protein_coding|1/3|c.10A>G|p.Thr4Ala|10/300|4/100|10/300,"
        "G|synonymous_variant|LOW|TraesCS1A02G000100|gene|transcript|"
        "TraesCS1A02G000100.1|protein_coding|1/3|c.12A>G|p.Thr4Thr|12/300|4/100|12/300\n",
        encoding="utf-8",
    )
    table = snpeff_vcf_to_vep_rows(vcf)
    assert list(table["Uploaded_variation"]) == ["1-100"]
    assert list(table["Gene"]) == ["TraesCS1A02G000100"]
    assert list(table["Consequence"]) == ["missense_variant"]
    written = run_snpeff_vcf_to_tab(vcf, tmp_path / "A")
    assert (tmp_path / "A.snpeff.tsv").is_file()
    assert len(written) == 1
