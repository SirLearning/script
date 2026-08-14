nextflow.enable.dsl=2

include { getPopDepTaxaBamFile_v1 } from '../../infra/infra_ref_v1.nf'
include { getPopDepTaxaBamFileAll_v1 } from '../../infra/infra_ref_v1.nf'
include { getPopDepNTaxaForChr_v1 } from '../../infra/infra_ref_v1.nf'
include { getRefV1ChrLength } from '../../infra/infra_ref_v1.nf'

def popdepBenchChrLength(chr) {
    def full = getRefV1ChrLength(chr.toString()).toLong()
    def capParam = params.popdep_bench_max_bp_per_chr ?: params.popdep_loo_bench_max_bp
    if (capParam) {
        def cap = capParam.toLong()
        return Math.min(full, cap)
    }
    return full
}

def popdepVariantPublishDir() {
    def root = params.popdep_publish_dir ?: params.popdep_dir
    return root
        ? "${root}/variant"
        : "${params.output_dir}/${params.job}/process/variant"
}

def popdepLogPublishDir() {
    def root = params.popdep_publish_dir ?: params.popdep_dir
    return root
        ? "${root}/logs"
        : "${params.output_dir}/${params.job}/process"
}

def popdepLooPublishRoot() {
    return params.popdep_loo_publish_root
        ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"
}

def popdepCrossChrTigerExtraArgs() {
    def order = params.popdep_crosschr_taxa_order.toString()
    if (!(order in ['interleave', 'shuffle', 'sorted'])) {
        throw new IllegalArgumentException(
            "popdep_crosschr_taxa_order must be interleave, shuffle, or sorted (got: ${order})"
        )
    }
    def lines = ["-o ${order}"]
    if (params.popdep_crosschr_checkpoint_dir) {
        lines << "-k ${params.popdep_crosschr_checkpoint_dir}"
        lines << "-ci ${params.popdep_crosschr_checkpoint_interval}"
    }
    return lines.collect { "        ${it} \\" }.join('\n')
}

process calc_population_depth {
    tag "prepare popdepth: ${chr}"
    publishDir "${popdepLogPublishDir()}", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'popdep_tiger'

    input:
    tuple val(id), val(chr), path(vcf)
    tuple path(tiger_jar), val(app_name), val(java_version)

    output:
    tuple val(id), val(chr), path("${id}.popdep.txt.gz"), emit: tiger_gz

    script:
    def tb_file = params.popdep_taxa_bam_file ?: getPopDepTaxaBamFile_v1(chr, params.home_dir)
    def chr_length = popdepBenchChrLength(chr)
    """
    set -euo pipefail
    exec > popdep_${chr}.log 2>&1

    echo "Starting population depth analysis for chromosome ${chr}..."
    echo "Using ${java_version} for ${app_name} process"
    echo "The environment java version is:"
    which java
    java -version

    echo "System resources before TIGER execution:"
    echo "Memory: \$(free -h | grep Mem)"
    echo "CPU cores (allocated): ${task.cpus}"
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
    echo "TIGER jar: ${tiger_jar}"

    java -Xmx${task.memory.toGiga()}G -jar ${tiger_jar} \\
        -app ${app_name} \\
        -a ${tb_file} \\
        -b ${chr} \\
        -c ${chr_length} \\
        -d samtools \\
        -e ${task.cpus} \\
        -f ${id}.popdep.txt.gz

    echo "TIGER population depth completed for chromosome ${chr}."
    """
}

process calc_population_depth_crosschr {
    tag "PopDepCrossChr"
    publishDir "${popdepLogPublishDir()}", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'popdep_tiger_crosschr'

    input:
    val(chr_list)
    tuple path(tiger_jar), val(app_name), val(java_version)

    output:
    path("crosschr_out/*.popdep.txt.gz"), emit: tiger_gz

    script:
    def tb_file = params.popdep_taxa_bam_file ?: getPopDepTaxaBamFileAll_v1(params.home_dir)
    def chr_length_lines = chr_list.collect { c ->
        def n_taxa = getPopDepNTaxaForChr_v1(c, params.home_dir)
        "echo -e '${c}\\t${popdepBenchChrLength(c.toString())}\\t${n_taxa}' >> chromosomeLength.txt"
    }.join('\n    ')
    def chr_bash = chr_list.collect { c -> c.toString() }.join(' ')
    """
    set -euo pipefail
    exec > popdep_crosschr.log 2>&1

    echo "Starting PopDepCrossChr (all chromosomes in one pass)..."
    echo "Using ${java_version} for ${app_name} process"
    echo "Chromosomes: ${chr_bash}"
    which java
    java -version

    echo "System resources before TIGER execution:"
    echo "Memory: \$(free -h | grep Mem)"
    echo "CPU cores (allocated): ${task.cpus}"
    echo "Java memory allocation (Xmx): ${task.memory.toGiga()}g"
    echo "TIGER jar: ${tiger_jar}"
    echo "Taxa-BAM map: ${tb_file}"

    echo -e 'Chr\\tLength\\tnTaxa' > chromosomeLength.txt
    ${chr_length_lines}
    echo "Chromosome length file (Chr, Length, nTaxa per subgenome taxaBamMap):"
    cat chromosomeLength.txt
    echo "PopDepCrossChr TIGER extras: taxa_order=${params.popdep_crosschr_taxa_order} checkpoint=${params.popdep_crosschr_checkpoint_dir ?: 'off'} interval=${params.popdep_crosschr_checkpoint_interval}"

    java -Xmx${task.memory.toGiga()}G -jar ${tiger_jar} \\
        -app PopDepCrossChr \\
        -a ${tb_file} \\
        -b chromosomeLength.txt \\
        -d samtools \\
        -e ${task.cpus} \\
${popdepCrossChrTigerExtraArgs()}
        -f crosschr_out/

    for chr in ${chr_bash}; do
        gz="crosschr_out/\${chr}.popdep.txt.gz"
        if [ ! -f "\$gz" ]; then
            echo "Missing TIGER output: \$gz" >&2
            exit 1
        fi
    done

    echo "PopDepCrossChr TIGER completed for chromosomes: ${chr_bash}"
    """
}

process popdep_tiger_gz_to_bgzip_tabix {
    tag "popdep bgz: ${id}"
    publishDir "${popdepVariantPublishDir()}", mode: 'copy', pattern: "*.popdep.txt.bgz*"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'popdep_bgz_tabix'

    input:
    tuple val(id), val(chr), path(gz)

    output:
    tuple val(id), val(chr), path("${id}.popdep.txt.bgz"), path("${id}.popdep.txt.bgz.tbi"), emit: popdep

    script:
    def bgzip = "${params.user_dir}/miniconda3/envs/stats/bin/bgzip"
    def tabix = "${params.user_dir}/miniconda3/envs/stats/bin/tabix"
    """
    set -euo pipefail

    echo "Converting ${id} TIGER gzip to BGZF + tabix (preserve TIGER header; prepend Chrom)..."
    gunzip -c ${gz} | awk -v chr="${chr}" '
        NR == 1 { print "Chrom\\t" \$0; next }
        { print chr "\\t" \$0 }
    ' | ${bgzip} -c -@ ${task.cpus} > ${id}.popdep.txt.bgz
    ${tabix} -S 1 -s 1 -b 2 -e 2 -f ${id}.popdep.txt.bgz
    echo "BGZF + tabix completed for ${id}."
    """
}

process prepare_popdep_loo_taxa_maps {
    tag "popdep loo prepare maps"
    publishDir "${publish_root}/meta", mode: 'copy', pattern: "cases.tsv"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    val publish_root
    path source_map
    val max_loo_cases
    val target_taxa_path

    output:
    path "cases.tsv", emit: manifest
    path "taxa_maps", emit: taxa_maps

    script:
    """
    set -euo pipefail
    mkdir -p taxa_maps
    header=\$(head -1 ${source_map})
    cp ${source_map} taxa_maps/baseline.taxaBamMap.txt
    echo -e "case_id\\tchr\\ttaxa_map" > cases.tsv
    echo -e "baseline\\t${params.popdep_loo_chr}\\ttaxa_maps/baseline.taxaBamMap.txt" >> cases.tsv
    n=0
    max=${max_loo_cases}
    target_file='${target_taxa_path}'

    emit_loo_case() {
      local taxa="\$1"
      [[ -z "\$taxa" ]] && return 0
      local out="taxa_maps/loo_\${taxa}.taxaBamMap.txt"
      {
        echo "\$header"
        awk -F'\\t' -v skip="\$taxa" 'NR>1 && \$1 != skip {print}' ${source_map}
      } > "\$out"
      echo -e "loo_\${taxa}\\t${params.popdep_loo_chr}\\t\$out" >> cases.tsv
      n=\$((n + 1))
    }

    if [[ -n "\$target_file" && -f "\$target_file" ]]; then
      while IFS= read -r taxa || [[ -n "\$taxa" ]]; do
        taxa="\${taxa%%#*}"
        taxa="\$(echo "\$taxa" | xargs)"
        [[ -z "\$taxa" ]] && continue
        if [[ "\$max" -gt 0 && "\$n" -ge "\$max" ]]; then
          break
        fi
        emit_loo_case "\$taxa"
      done < "\$target_file"
    else
      while IFS=\$'\\t' read -r taxa _cov _bam; do
        [[ -z "\$taxa" ]] && continue
        if [[ "\$max" -gt 0 && "\$n" -ge "\$max" ]]; then
          break
        fi
        emit_loo_case "\$taxa"
      done < <(tail -n +2 ${source_map})
    fi
    echo "Prepared \$(awk 'NR>1' cases.tsv | wc -l) PopDep LOO cases (1 baseline + \${n} LOO maps)"
    """
}

process calc_population_depth_loo_case {
    tag "popdep loo: ${case_id} chr${chr}"
    publishDir "${publish_root}/${case_id}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/tiger"
    label 'popdep_loo_tiger'

    input:
    tuple val(publish_root), val(case_id), val(chr), path(taxa_map), path(vcf)
    tuple path(tiger_jar), val(app_name), val(java_version)

    output:
    tuple val(case_id), val(chr), path("${String.format('chr%03d', chr.toInteger())}.popdep.txt.gz"), emit: tiger_gz

    script:
    def grid_id = String.format('chr%03d', chr.toString().toInteger())
    def chr_length = popdepBenchChrLength(chr)
    """
    set -euo pipefail
    exec > popdep_loo_${case_id}_${chr}.log 2>&1

    echo "PopDep LOO case=${case_id} chr=${chr} length=${chr_length} taxa_map=${taxa_map}"
    echo "Using ${java_version} for ${app_name}; jar=${tiger_jar}"
    which java
    java -version
    echo "CPU cores (allocated): ${task.cpus}; Xmx: ${task.memory.toGiga()}g"

    java -Xmx${task.memory.toGiga()}G -jar ${tiger_jar} \\
        -app ${app_name} \\
        -a ${taxa_map} \\
        -b ${chr} \\
        -c ${chr_length} \\
        -d samtools \\
        -e ${task.cpus} \\
        -f ${grid_id}.popdep.txt.gz

    echo "TIGER PopDep LOO completed for ${case_id} chr ${chr}."
    """
}

process popdep_tiger_gz_to_bgzip_tabix_loo {
    tag "popdep loo bgz: ${case_id}"
    publishDir "${publish_root}/${case_id}/variant", mode: 'copy', pattern: "*.popdep.txt.bgz*"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'popdep_bgz_tabix'

    input:
    tuple val(publish_root), val(case_id), val(chr), path(gz)

    output:
    tuple val(case_id), val(chr), path("${String.format('chr%03d', chr.toInteger())}.popdep.txt.bgz"), path("${String.format('chr%03d', chr.toInteger())}.popdep.txt.bgz.tbi"), emit: popdep

    script:
    def grid_id = String.format('chr%03d', chr.toString().toInteger())
    def bgzip = "${params.user_dir}/miniconda3/envs/stats/bin/bgzip"
    def tabix = "${params.user_dir}/miniconda3/envs/stats/bin/tabix"
    """
    set -euo pipefail

    echo "LOO BGZF+tabix case=${case_id} chr=${chr} grid=${grid_id}"
    gunzip -c ${gz} | awk -v chr="${chr}" '
        NR == 1 { print "Chrom\\t" \$0; next }
        { print chr "\\t" \$0 }
    ' | ${bgzip} -c -@ ${task.cpus} > ${grid_id}.popdep.txt.bgz
    ${tabix} -S 1 -s 1 -b 2 -e 2 -f ${grid_id}.popdep.txt.bgz
    echo "BGZF + tabix completed for ${case_id}/${grid_id}."
    """
}

process annotate_subgenome_variant_popdep {
    tag "annotate variant popdep: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.popdep.info.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.popdep.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.popdep.info.tsv"), emit: popdep
    path "${chr}.popdep.log", emit: log

    script:
    def popdep_workers = params.popdep_lookup_workers ?: 0
    def popdep_workers_py = (popdep_workers as int) > 0 ? "${popdep_workers}" : "None"
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.popdep.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.popdep import annotate_variants_popdep_from_pvar

    workers = ${popdep_workers_py}
    print(f"Annotating variant popdep for ${id} (${chr}) from ${params.popdep_dir} (tabix, workers={workers}) ...")
    annotate_variants_popdep_from_pvar(
        "${pvar}",
        "${params.popdep_dir}",
        "${id}.popdep.info.tsv",
        max_workers=workers,
        use_tabix=True,
    )
    """
}
