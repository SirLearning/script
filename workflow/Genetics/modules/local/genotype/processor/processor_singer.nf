nextflow.enable.dsl=2

/*
 * Extract 2 Mb windows from common-inter refchr pfiles, rewrite GT '/' to '|',
 * impute missing genotypes with Beagle, then sample ARGs with SINGER.
 *
 * Does not overwrite test_common_inter or test_common_inter_refchr.
 */


process prepare_wap_china_singer_inputs {
    tag "WAP China keep + ${n_windows} windows"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "WAP_China.keep"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "singer_windows.tsv"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '8.GB'

    input:
    path germplasm_tsv
    path psam
    path pvar_a
    path pvar_b
    path pvar_d
    val n_windows
    val window_bp
    val min_snps
    val country

    output:
    path "WAP_China.keep", emit: keep
    path "singer_windows.tsv", emit: windows
    path "prepare_wap_china.log", emit: log

    script:
    """
    set -euo pipefail
    exec > prepare_wap_china.log 2>&1

    python - <<'PY'
import sys

if "${params.src_dir}" and "${params.src_dir}" != "null":
    sys.path.insert(0, "${params.src_dir}/python")

from genetics.genomics.variant.singer_arg import prepare_wap_china_singer_inputs

n_taxa, n_win = prepare_wap_china_singer_inputs(
    germplasm_tsv="${germplasm_tsv}",
    psam_path="${psam}",
    pvar_by_subgenome={"A": "${pvar_a}", "B": "${pvar_b}", "D": "${pvar_d}"},
    keep_path="WAP_China.keep",
    windows_path="singer_windows.tsv",
    country="${country}",
    n_windows=${n_windows},
    window_bp=${window_bp},
    min_snps=${min_snps},
)
print(f"n_taxa={n_taxa} n_windows={n_win}")
PY
    """
}


process extract_common_inter_window_pfile {
    tag "2Mb ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/windows", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/export", mode: 'copy', pattern: "*.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'singer_window'

    input:
    tuple val(id), val(plink_chr), val(ref_chrom), val(start_bp), val(end_bp), path(pgen), path(psam), path(pvar)
    path keep

    output:
    tuple val(id), val(ref_chrom), val(start_bp), val(end_bp), path("${id}.vcf.gz"), emit: vcf
    tuple val(id), path("${id}.plink2.pgen"), path("${id}.plink2.psam"), path("${id}.plink2.pvar"), emit: pfile
    path "extract_window_${id}.log", emit: logs

    script:
    def prefix = pgen.baseName
    """
    set -euo pipefail
    exec > extract_window_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --chr ${plink_chr} \\
        --from-bp ${start_bp} \\
        --to-bp ${end_bp} \\
        --snps-only just-acgt \\
        --max-alleles 2 \\
        --keep ${keep} \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out window

    printf '%s\\t%s\\n' ${plink_chr} ${ref_chrom} > update.chr
    plink2 --pfile window \\
        --allow-extra-chr \\
        --update-chr update.chr \\
        --sort-vars \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out recoded

    plink2 --pfile recoded \\
        --allow-extra-chr \\
        --export vcf-4.2 bgz \\
        --threads ${task.cpus} \\
        --out recoded_vcf

    mv recoded.pgen ${id}.plink2.pgen
    mv recoded.psam ${id}.plink2.psam
    mv recoded.pvar ${id}.plink2.pvar
    mv recoded_vcf.vcf.gz ${id}.vcf.gz
    """
}


process vcf_selfing_phase_slash_to_pipe {
    tag "selfing GT phase: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/export", mode: 'copy', pattern: "*.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '8.GB'

    input:
    tuple val(id), val(ref_chrom), val(start_bp), val(end_bp), path(vcf)

    output:
    tuple val(id), val(ref_chrom), val(start_bp), val(end_bp), path("${id}.selfing_phased.vcf.gz"), emit: vcf
    path "selfing_phase_${id}.log", emit: log

    script:
    """
    set -euo pipefail
    exec > selfing_phase_${id}.log 2>&1

    python - <<'PY'
import sys

if "${params.src_dir}" and "${params.src_dir}" != "null":
    sys.path.insert(0, "${params.src_dir}/python")

from genetics.genomics.variant.selfing_gt import rewrite_gt_slash_to_pipe

n = rewrite_gt_slash_to_pipe("${vcf}", "${id}.selfing_phased.vcf.gz")
print(f"rewrote_sites={n}")
PY
    """
}


process beagle_impute_vcf {
    tag "beagle impute: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/export", mode: 'copy', pattern: "*.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'singer_window'

    input:
    tuple val(id), val(ref_chrom), val(start_bp), val(end_bp), path(vcf)
    path beagle_jar

    output:
    tuple val(id), val(ref_chrom), val(start_bp), val(end_bp), path("${id}.beagle.vcf.gz"), emit: vcf
    path "beagle_${id}.log", emit: log

    script:
    def jar = beagle_jar
    def xmx = Math.max(2, task.memory.toGiga() - 2)
    """
    set -euo pipefail
    exec > beagle_${id}.log 2>&1

    java -Xmx${xmx}g -jar ${jar} \\
        gt=${vcf} \\
        out=${id}.beagle_raw \\
        nthreads=${task.cpus} \\
        impute=true

    mv ${id}.beagle_raw.vcf.gz ${id}.beagle.vcf.gz
    """
}


process singer_sample_arg {
    tag "SINGER ARG: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/singer", mode: 'copy', pattern: "singer_${id}*"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'singer_arg'

    input:
    tuple val(id), val(ref_chrom), val(start_bp), val(end_bp), path(vcf)
    path singer_dir

    output:
    path "singer_${id}*", emit: artefacts
    path "singer_${id}.run.log", emit: log

    script:
    def bin_dir = singer_dir
    def mut = params.singer_mut_rate
    def n = params.singer_n
    def thin = params.singer_thin
    def polar = params.singer_polar
    def seed_attempts = params.singer_seed_attempts ?: 3
    """
    set -euo pipefail
    exec > singer_${id}.run.log 2>&1

    gzip -dc ${vcf} > ${id}.imputed.vcf
    prefix=${id}.imputed

    python - <<'PY'
import sys
from pathlib import Path

if "${params.src_dir}" and "${params.src_dir}" != "null":
    sys.path.insert(0, "${params.src_dir}/python")

from genetics.genomics.variant.singer_arg import (
    estimate_ne_from_vcf,
    random_singer_seeds,
    write_ne_file,
)

ne = estimate_ne_from_vcf(
    "${id}.imputed.vcf",
    start=${start_bp},
    end=${end_bp},
    mut_rate=${mut},
)
write_ne_file("singer_ne.txt", ne)
print(f"estimated_Ne={ne}")

fixed = "${params.singer_seed}"
if fixed and fixed != "null":
    seeds = [int(fixed)]
else:
    seeds = random_singer_seeds(${seed_attempts})
Path("singer_seeds.txt").write_text("\\n".join(str(s) for s in seeds) + "\\n")
print("singer_seeds=" + ",".join(str(s) for s in seeds))
PY

    Ne=\$(cat singer_ne.txt)
    ok=0
    while read -r seed; do
        echo "trying seed=\${seed} (engine singer; not singer_master -debug)"
        rm -f singer_${id}_nodes_* singer_${id}_branches_* singer_${id}_recombs_* singer_${id}_muts_* singer_${id}_flats_* singer_${id}*.trees
        set +e
        ${bin_dir}/singer \\
            -Ne \${Ne} \\
            -m ${mut} \\
            -r ${mut} \\
            -input \${prefix} \\
            -output singer_${id} \\
            -start ${start_bp} \\
            -end ${end_bp} \\
            -n ${n} \\
            -thin ${thin} \\
            -polar ${polar} \\
            -seed \${seed}
        rc=\$?
        set -e
        if [ "\${rc}" -eq 0 ]; then
            echo "singer succeeded seed=\${seed}"
            ok=1
            break
        fi
        echo "singer failed rc=\${rc} seed=\${seed}"
    done < singer_seeds.txt
    if [ "\${ok}" -ne 1 ]; then
        echo "all seed attempts failed"
        exit 1
    fi

    if [ -x ${bin_dir}/convert_to_tskit ]; then
        ${bin_dir}/convert_to_tskit \\
            -input singer_${id} \\
            -output singer_${id} \\
            -start 0 \\
            -end \$((${n} - 1)) \\
            -step 1 || echo "convert_to_tskit skipped (tskit missing or convert failed)"
    fi
    """
}
