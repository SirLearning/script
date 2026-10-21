nextflow.enable.dsl=2

process subsampling_pfile_for_test {
    tag "subsample pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    val thin_rate

    output:
    tuple val(id), val(chr), val("${id}.thin"), path("${id}.thin.pgen"), path("${id}.thin.psam"), path("${id}.thin.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > subsample_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --thin ${thin_rate} \\
        --seed \$(date +%s) \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.thin
    """
}

process compute_abd_global_thin_rate {
    tag "ABD global thin rate to ${params.test_set_target_n} sites"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/info", mode: 'copy', pattern: "global_thin_rate.tsv"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path pvars

    output:
    path "thin_rate.txt", emit: rate_file
    path "global_thin_rate.tsv", emit: table

    script:
    """
    python3 <<'PY'
    from pathlib import Path
    import sys

    target_n = int(${params.test_set_target_n})
    max_n = int(${params.test_set_max_n})
    n = 0
    per = []
    for p in sorted(Path('.').glob('*.pvar')):
        c = 0
        with p.open() as fh:
            for line in fh:
                if line and not line.startswith('#'):
                    c += 1
        per.append((p.name, c))
        n += c
    if n <= 0:
        raise SystemExit('compute_abd_global_thin_rate: no variants in A/B/D pvar files')
    rate = 1.0 if n <= max_n else min(1.0, target_n / n)
    with open('global_thin_rate.tsv', 'w') as out:
        out.write('file\\tsites\\n')
        for name, c in per:
            out.write(f'{name}\\t{c}\\n')
        out.write(f'TOTAL_ABD\\t{n}\\n')
        out.write(f'target_n\\t{target_n}\\n')
        out.write(f'max_n\\t{max_n}\\n')
        out.write(f'thin_rate\\t{rate:.8f}\\n')
    Path('thin_rate.txt').write_text(f'{rate:.8f}\\n')
    print(f'ABD_sites={n} target_n={target_n} max_n={max_n} thin_rate={rate:.8f}', file=sys.stderr)
    PY
    """
}

process thin_merged_subgenome_pfile {
    tag "thin merged pfile with global keep rate: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    val thin_rate

    output:
    tuple val(id), val(chr), val("${id}_test.plink2"), path("${id}_test.plink2.pgen"), path("${id}_test.plink2.psam"), path("${id}_test.plink2.pvar"), emit: pfile
    tuple val(id), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > thin_merged_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --thin ${thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}_test.plink2
    """
}

process subsampling_common_variant_pfile_for_test {
    tag "common-thin pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.common.thin"), path("${id}.common.thin.pgen"), path("${id}.common.thin.psam"), path("${id}.common.thin.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > common_thin_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --geno ${params.hf_geno} \\
        --maf ${params.hf_maf} \\
        --thin ${params.hf_thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.common.thin
    """
}

process subsampling_maf_only_pfile_for_test {
    tag "maf-only common pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.common.only"), path("${id}.common.only.pgen"), path("${id}.common.only.psam"), path("${id}.common.only.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > common_only_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --maf ${params.hf_maf} \\
        --thin ${params.hf_thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.common.only
    """
}

process build_chr_rare_inter_extract {
    tag "rare-inter extract: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/extract", mode: 'copy', pattern: "*.txt"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/extract/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/extract/afreq/${id}", mode: 'copy', pattern: "${id}.*.afreq"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'rare_inter_extract'

    input:
    path group_file
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), path("${id}.rare_inter.txt"), emit: pfile_with_extract

    script:
    def published_afreq_dir = "${params.output_dir}/${params.job}/process/${params.mod}/extract/afreq"
    def published_extract = "${params.output_dir}/${params.job}/process/${params.mod}/extract/${id}.rare_inter.txt"
    def published_summary = "${params.output_dir}/${params.job}/process/${params.mod}/extract/info/${id}.rare_inter.summary.tsv"
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("build_rare_inter_${id}.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.group_maf_intersection import write_rare_inter_extract_for_chromosome

    print("Building rare-inter extract for ${id} (all sample_group MAF < ${params.hf_maf}) ...")
    summary = write_rare_inter_extract_for_chromosome(
        "${prefix}",
        "${group_file}",
        "${id}.rare_inter.txt",
        chr_id="${id}",
        maf_threshold=${params.hf_maf},
        summary_path="${id}.rare_inter.summary.tsv",
        published_afreq_dir="${published_afreq_dir}",
        published_extract_path="${published_extract}",
        published_summary_path="${published_summary}",
        afreq_emit_dir=".",
        threads=${task.cpus},
        group_workers=${params.rare_inter_group_workers},
        plink_threads=${params.rare_inter_plink_threads},
    )
    print(summary)
    """
}

process subsampling_common_inter_pfile_for_test {
    tag "common-inter thin pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), path(extract)

    output:
    tuple val(id), val(chr), val("${id}.common.inter"), path("${id}.common.inter.pgen"), path("${id}.common.inter.psam"), path("${id}.common.inter.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > common_inter_pfile_${id}.log 2>&1

    # Keep sites that are NOT all-group-rare (MAF < hf_maf in every sample_group),
    # then randomly thin. Empty rare-inter list => exclude nothing.
    if [ -s "${extract}" ]; then
        echo "Excluding all-group-rare sites listed in ${extract}"
        plink2 --pfile ${prefix} \\
            --exclude ${extract} \\
            --thin ${params.hf_common_inter_thin_rate} \\
            --seed ${params.ld_decay_random_seed} \\
            --allow-extra-chr \\
            --make-pgen \\
            --threads ${task.cpus} \\
            --out ${id}.common.inter
    else
        echo "Empty rare-inter extract for ${id}; thinning all sites (nothing to exclude)."
        plink2 --pfile ${prefix} \\
            --thin ${params.hf_common_inter_thin_rate} \\
            --seed ${params.ld_decay_random_seed} \\
            --allow-extra-chr \\
            --make-pgen \\
            --threads ${task.cpus} \\
            --out ${id}.common.inter
    fi
    """
}

process extract_collinear_pfile_for_test {
    tag "collinear extract pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_collinear'}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_collinear'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), path(extract_list)

    output:
    tuple val(id), val(chr), val("${id}.collinear"), path("${id}.collinear.pgen"), path("${id}.collinear.psam"), path("${id}.collinear.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > collinear_extract_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --extract ${extract_list} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.collinear
    """
}

process publish_collinear_root {
    tag "publish collinear root"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_collinear'}", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    path collinear

    output:
    path collinear, emit: collinear_root

    script:
    """
    set -euo pipefail
    test -d "${collinear}/extract"
    n=\$(ls -1 "${collinear}/extract"/*.txt 2>/dev/null | wc -l)
    echo "collinear extract lists: \${n}"
    test "\${n}" -gt 0
    """
}

process collinear_sites_prep {
    tag "collinear sites prep"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_collinear'}/collinear", mode: 'copy', pattern: "**/*"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_8'

    output:
    path "collinear", emit: collinear_root

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("collinear_sites_prep.log", "w", buffering=1)
    sys.stderr = sys.stdout

    from genetics.genomics.alignment.collinear_sites import prepare_collinear_test_dataset

    print("Preparing collinear test dataset (cached site lists + join + sample)...")
    summary = prepare_collinear_test_dataset(
        syntenic_matrix_dir="${params.syntenic_matrix_dir}",
        main_raw_dir="${params.main_raw_dir ?: params.popdep_dir}",
        full_lineage_summary_tsv="${params.full_lineage_summary_tsv}",
        output_dir="collinear",
        target_n=${params.collinear_target_n},
        maf_threshold=${params.collinear_maf_threshold},
        seed=${params.ld_decay_random_seed},
        workers=${params.collinear_workers},
        force=${params.collinear_force ? 'True' : 'False'},
    )
    print(summary)
    """
}

process prepare_test_site_panel {
    tag "test_site panel prep"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_site'}", mode: 'copy', pattern: "**/*"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_8'

    output:
    path "test_site", emit: site_root

    script:
    def chr_exclude_py = params.test_site_chr_exclude
        ? params.test_site_chr_exclude.toString().split(/[,\\s]+/).collect { it.trim() }.findAll { it }.collect { it.toInteger() }
        : []
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("prepare_test_site_panel.log", "w", buffering=1)
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep_site_panel import prepare_test_site_panel

    chr_exclude = set(${chr_exclude_py})
    print("Preparing test_site random genome panel + popdep extract...")
    summary = prepare_test_site_panel(
        popdep_dir="${params.popdep_dir}",
        main_raw_dir="${params.main_raw_dir ?: params.popdep_dir}",
        output_dir="test_site",
        target_bp=${params.test_site_target_bp},
        seed=${params.test_site_seed ?: params.ld_decay_random_seed},
        chr_exclude=chr_exclude,
        force=${params.test_site_force ? 'True' : 'False'},
    )
    print(summary)
    """
}

process extract_test_site_pfile_for_test {
    tag "test_site extract pfile: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_site'}/variant", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_site'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), path(extract_list)

    output:
    tuple val(id), val(chr), val("${id}.site"), path("${id}.site.pgen"), path("${id}.site.psam"), path("${id}.site.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > test_site_extract_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --extract ${extract_list} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.site
    """
}

process publish_test_site_root {
    tag "publish test_site root"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_site'}", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    path site_root

    output:
    path site_root, emit: site_root

    script:
    """
    set -euo pipefail
    test -f "${site_root}/sites/panel.sites.tsv"
    n=\$(wc -l < "${site_root}/sites/panel.sites.tsv")
    echo "test_site panel rows (incl header): \${n}"
    test "\${n}" -gt 1
    """
}

process prepare_test_wap_keep {
    tag "prepare WAP panel keep list"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.keep"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    path group_file

    output:
    path "test_wap.keep", emit: keep

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("prepare_test_wap_keep.log", "w", buffering=1)
    sys.stderr = sys.stdout

    from genetics.genomics.variant.group_maf_intersection import write_sample_group_keep_file

    include = frozenset(
        g.strip() for g in "${params.test_wap_include_groups}".split(",") if g.strip()
    )
    n = write_sample_group_keep_file("${group_file}", "test_wap.keep", include_groups=include)
    print(f"Wrote test_wap.keep with {n} samples (included groups: {sorted(include)})")
    """
}

process gff_feature_bed {
    tag "GFF ${params.test_wap_exon_gff_type} BED"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/annotate", mode: 'copy', pattern: "*.bed"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    path gff

    output:
    path "${params.mod}.extract.bed", emit: bed
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("gff_feature_bed.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.gene.gff_bed import parse_gff_feature_types, write_merged_gff_bed

    types = parse_gff_feature_types("${params.test_wap_exon_gff_type}")
    n = write_merged_gff_bed("${gff}", "${params.mod}.extract.bed", feature_types=types)
    print(f"wrote {n} merged {','.join(types)} intervals from ${gff}")
    """
}

process extract_wap_exon_pfile {
    tag "WAP exon pfile: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'wap_exon_extract'

    input:
    path keep_file
    path exon_bed
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.wap.exon"), path("${id}.wap.exon.pgen"), path("${id}.wap.exon.psam"), path("${id}.wap.exon.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > wap_exon_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --keep ${keep_file} \\
        --extract bed0 ${exon_bed} \\
        --mac ${params.test_wap_min_mac} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.wap.exon
    """
}

process subsampling_wap_pfile_for_test {
    tag "WAP panel pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    path keep_file
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.wap"), path("${id}.wap.pgen"), path("${id}.wap.psam"), path("${id}.wap.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > wap_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --keep ${keep_file} \\
        --mac ${params.test_wap_min_mac} \\
        --thin ${params.test_wap_thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.wap
    """
}

process subsampling_wap_common_pfile_for_test {
    tag "WAP common pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    path keep_file
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.wapc"), path("${id}.wapc.pgen"), path("${id}.wapc.psam"), path("${id}.wapc.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > wap_common_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --keep ${keep_file} \\
        --maf ${params.test_wap_common_maf} \\
        --thin ${params.test_wap_common_thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.wapc
    """
}

process prepare_test_hex_keep {
    tag "prepare hex panel keep list"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.keep"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    path group_file

    output:
    path "test_hex.keep", emit: keep

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("prepare_test_hex_keep.log", "w", buffering=1)
    sys.stderr = sys.stdout

    from genetics.genomics.variant.group_maf_intersection import write_sample_group_keep_file

    exclude = frozenset(
        g.strip() for g in "${params.test_hex_exclude_groups}".split(",") if g.strip()
    )
    n = write_sample_group_keep_file("${group_file}", "test_hex.keep", exclude_groups=exclude)
    print(f"Wrote test_hex.keep with {n} samples (excluded groups: {sorted(exclude)})")
    """
}

process subsampling_hex_pfile_for_test {
    tag "hex panel pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    path keep_file
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.hex"), path("${id}.hex.pgen"), path("${id}.hex.psam"), path("${id}.hex.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > hex_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --keep ${keep_file} \\
        --mac 2 \\
        --thin ${params.test_hex_thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.hex
    """
}

process subsampling_rare_only_pfile_for_test {
    tag "rare-only pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.rare.only"), path("${id}.rare.only.pgen"), path("${id}.rare.only.psam"), path("${id}.rare.only.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > rare_only_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --max-maf ${params.hf_maf} \\
        --thin ${params.thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.rare.only
    """
}

process merge_subgenome_test_pfile {
    tag "merge plink2 test subgenome pfile: ${subgenome}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{bed,bim,fam,pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    // Receives lists of files/values grouped by subgenome
    tuple val(subgenome), val(ids), val(chrs), val(prefixes), path(pgens), path(psams), path(pvars)

    output:
    tuple val(subgenome), val("sub_${subgenome}"), val("${subgenome}_test.plink"), path("${subgenome}_test.plink.bed"), path("${subgenome}_test.plink.bim"), path("${subgenome}_test.plink.fam"), emit: bfile
    tuple val(subgenome), val("sub_${subgenome}"), val("${subgenome}_test.plink2"), path("${subgenome}_test.plink2.pgen"), path("${subgenome}_test.plink2.psam"), path("${subgenome}_test.plink2.pvar"), emit: pfile
    tuple val(subgenome), path("*.log"), emit: logs

    script:
    def has_multiple_prefixes = prefixes.size() > 1
    def merge_list_block = has_multiple_prefixes
        ? """echo "${prefixes[1..-1].join('\n')}" > ${subgenome}.merge_list.txt

    plink2 --pfile ${prefixes[0]} \\
        --pmerge-list ${subgenome}.merge_list.txt \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${subgenome}_test.plink2"""
        : """plink2 --pfile ${prefixes[0]} \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${subgenome}_test.plink2"""
    """
    set -euo pipefail
    exec > merge_subgenome_test.${subgenome}.log 2>&1
    
    ${merge_list_block}

    plink2 --pfile ${subgenome}_test.plink2 \\
        --make-bed \\
        --threads ${task.cpus} \\
        --out ${subgenome}_test.plink
    """
}

process prepare_test_popdep_annotations {
    tag "test_popdep region-class beds"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}", mode: 'copy', pattern: "annotation/**"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    output:
    path "annotation", emit: annotation_dir
    path "annotation/chr_alloc.tsv", emit: alloc

    script:
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path

    Path("annotation").mkdir(parents=True, exist_ok=True)
    sys.stdout = open("prepare_test_popdep_annotations.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.popdep_region_class import build_region_class_beds
    from genetics.genomics.variant.popdep_site_panel import allocate_test_popdep_sites
    from infra.utils.io import save_df_to_tsv

    print("Building gap/repeat/GC/CNV segment BEDs...")
    counts = build_region_class_beds(
        gap_bed="${params.test_popdep_gap_bed}",
        te_gff="${params.test_popdep_te_gff}",
        fasta_path="${params.test_popdep_fasta}",
        mapq_popdep_dir="${params.test_popdep_mapq_dir}",
        output_dir="annotation",
        window_bp=${params.test_popdep_gc_window_bp ?: 1000},
        gc_high=${params.test_popdep_gc_high ?: 0.55},
        gc_low=${params.test_popdep_gc_low ?: 0.35},
        cnv_frac_min=${params.test_popdep_cnv_frac_min ?: 0.2},
        mapq_prefix="${params.test_popdep_mapq_prefix ?: 'allmapq_1mb'}",
    )
    print(counts)
    alloc = allocate_test_popdep_sites(${params.test_popdep_target_bp ?: 100000000})
    save_df_to_tsv(alloc, "annotation/chr_alloc.tsv")
    print(alloc.to_string(index=False))
    """
}

process prepare_test_popdep_chr_panel {
    tag "test_popdep chr panel: ${chr_id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/sites", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'test_popdep_chr'

    input:
    tuple val(chr_id), val(n_sites)
    path annotation_dir

    output:
    path "${chr_id}.popdep.sites.tsv", emit: sites
    path "*.log", emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("prepare_test_popdep_${chr_id}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.popdep_site_panel import prepare_test_popdep_chr_panel

    print("Sampling ${n_sites} sites on ${chr_id}...")
    prepare_test_popdep_chr_panel(
        "${chr_id}",
        int(${n_sites}),
        popdep_dir="${params.popdep_dir}",
        annotation_dir="${annotation_dir}",
        output_tsv="${chr_id}.popdep.sites.tsv",
        seed=${params.test_popdep_seed ?: 42},
    )
    """
}

process rescale_test_popdep_chr_sites {
    tag "rescale ntaxa500 ${chr_id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/sites", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'test_popdep_chr'

    input:
    val chr_id

    output:
    path "${chr_id}.popdep.sites.tsv", emit: sites
    path "*.log", emit: logs

    script:
    def sites_dir = params.test_popdep_sites_dir ?: "${params.output_dir}/${params.job}/process/${params.mod ?: 'test_popdep'}/sites"
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path

    sys.stdout = open("rescale_test_popdep_${chr_id}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.popdep_site_panel import rescale_test_popdep_sites

    src = Path("${sites_dir}") / "${chr_id}.popdep.sites.tsv"
    print(f"Rescaling {src} to nTaxa=500...")
    rescale_test_popdep_sites(src, "${chr_id}.popdep.sites.tsv")
    """
}

