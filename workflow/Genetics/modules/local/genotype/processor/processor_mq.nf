nextflow.enable.dsl=2

/*
 * Full-chromosome site MQ from bcftools mpileup on a BAM list (no bcftools call).
 *
 * Mean MAPQ per pileup site from INFO/I16:
 *   MQ = (I16[9] + I16[11]) / (I16[1] + I16[2] + I16[3] + I16[4])
 * (1-based I16 indices: ref/non-ref forward/reverse Q13+ counts; ref/non-ref MQ sums)
 *
 * Two published products:
 *   1) *.site_mq.calls.tsv.gz — mpileup row: CHROM POS REF ALT MQ (float)
 *   2) *.site_mq.ref.txt.gz — optional PopDep-like grid (one MQ per reference POS)
 */

process calc_site_mq_bcftools {
    tag "site_mq: ${id}"
    publishDir "${params.mq_dir}/reference", mode: 'copy', pattern: "*.site_mq.*"
    publishDir "${params.mq_dir}/logs", mode: 'copy', pattern: "*.site_mq.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'site_mq_bcftools'

    input:
    tuple val(id), val(chr), val(chr_len), path(reference), path(bam_list)

    output:
    tuple val(id), val(chr), path("${id}.site_mq.ref.*"), emit: site_mq_ref
    tuple val(id), val(chr), path("${id}.site_mq.calls.tsv.gz"), emit: site_mq_calls
    path "${id}.site_mq.ref.info.tsv", emit: info

    script:
    def pad_flag = params.mq_pad_all_positions ? 'true' : 'false'
    def gzip_flag = params.mq_ref_gzip ? 'true' : 'false'
    def mq_dec = params.mq_decimal_places ?: 6
    """
    set -euo pipefail
    exec > ${id}.site_mq.log 2>&1

    echo "Site MQ: chr=${chr} id=${id} len=${chr_len}"
    echo "BAM list: ${bam_list}"
    echo "Reference FASTA: ${reference}"
    echo "mpileup max_depth=${params.mq_mpileup_max_depth}; pad_full_chr=${pad_flag}; gzip_ref=${gzip_flag}; mq_decimals=${mq_dec}"

    bcftools mpileup \\
        -f ${reference} \\
        -b ${bam_list} \\
        -r ${chr}:1-${chr_len} \\
        -d ${params.mq_mpileup_max_depth} \\
        -g 0 \\
        -Ou \\
        --threads ${task.cpus} \\
    | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/I16\\n' \\
    | LC_ALL=C sort -k2,2n -k3,3 -k4,4 \\
    | awk -F'\\t' -v dec=${mq_dec} '
        function mq_from_i16(s,    a, n, cnt, smq) {
            if (s == "." || s == "") return "."
            n = split(s, a, ",")
            if (n < 12) return "."
            cnt = a[1] + a[2] + a[3] + a[4]
            if (cnt <= 0) return "."
            smq = a[9] + a[11]
            return smq / cnt
        }
        {
            mq = mq_from_i16(\$5)
            if (mq == ".") {
                print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t."
            } else {
                printf "%s\\t%s\\t%s\\t%s\\t%.*f\\n", \$1, \$2, \$3, \$4, dec, mq
            }
        }
    ' > ${id}.site_mq.calls.tsv

    gzip -f ${id}.site_mq.calls.tsv

    n_mpileup=\$(zcat ${id}.site_mq.calls.tsv.gz | wc -l)
    echo "Mpileup sites with REF/ALT/MQ: \${n_mpileup}"

    zcat ${id}.site_mq.calls.tsv.gz \\
    | awk -F'\\t' '{print \$1 "\\t" \$2 "\\t" \$5}' \\
    | LC_ALL=C sort -k2,2n \\
    > ${id}.site_mq.perpos.tsv

    n_perpos=\$(wc -l < ${id}.site_mq.perpos.tsv)
    echo "Unique POS with MQ row (for pad): \${n_perpos}"

    if [ "${pad_flag}" = "true" ]; then
        awk -v chr=${chr} -v len=${chr_len} '
            BEGIN { prev = 0 }
            {
                pos = \$2 + 0
                mq = \$3
                while (prev + 1 < pos) {
                    prev++
                    print chr "\\t" prev "\\t."
                }
                prev = pos
                print chr "\\t" pos "\\t" mq
            }
            END {
                while (prev < len) {
                    prev++
                    print chr "\\t" prev "\\t."
                }
            }
        ' ${id}.site_mq.perpos.tsv > ${id}.site_mq.ref.txt
    else
        awk -F'\\t' -v chr=${chr} '{print chr "\\t" \$2 "\\t" \$3}' ${id}.site_mq.perpos.tsv > ${id}.site_mq.ref.txt
    fi

    rm -f ${id}.site_mq.perpos.tsv

    n_ref=\$(wc -l < ${id}.site_mq.ref.txt)
    n_dot=\$(awk -F'\\t' '\$3=="." {c++} END {print c+0}' ${id}.site_mq.ref.txt)
    n_numeric=\$((n_ref - n_dot))

    if [ "${gzip_flag}" = "true" ]; then
        gzip -f ${id}.site_mq.ref.txt
        ref_path="${id}.site_mq.ref.txt.gz"
    else
        ref_path="${id}.site_mq.ref.txt"
    fi

    printf 'Metric\\tValue\\n' > ${id}.site_mq.ref.info.tsv
    printf 'Chrom\\t${chr}\\n' >> ${id}.site_mq.ref.info.tsv
    printf 'Id\\t${id}\\n' >> ${id}.site_mq.ref.info.tsv
    printf 'ChrLength\\t${chr_len}\\n' >> ${id}.site_mq.ref.info.tsv
    printf 'BamList\\t${bam_list}\\n' >> ${id}.site_mq.ref.info.tsv
    printf 'MpileupSites\\t%s\\n' "\${n_mpileup}" >> ${id}.site_mq.ref.info.tsv
    printf 'UniquePosWithMQ\\t%s\\n' "\${n_perpos}" >> ${id}.site_mq.ref.info.tsv
    printf 'ReferenceRows\\t%s\\n' "\${n_ref}" >> ${id}.site_mq.ref.info.tsv
    printf 'NumericMQ\\t%s\\n' "\${n_numeric}" >> ${id}.site_mq.ref.info.tsv
    printf 'DotMQ\\t%s\\n' "\${n_dot}" >> ${id}.site_mq.ref.info.tsv
    printf 'CallsFile\\t${id}.site_mq.calls.tsv.gz\\n' >> ${id}.site_mq.ref.info.tsv
    printf 'ReferenceFile\\t%s\\n' "\${ref_path}" >> ${id}.site_mq.ref.info.tsv
    printf 'CallsColumns\tCHROM;POS;REF;ALT;MQ\\n' >> ${id}.site_mq.ref.info.tsv
    printf 'MqSource\\tI16_mean_MAPQ_no_call\\n' >> ${id}.site_mq.ref.info.tsv
    printf 'MqFormula\\t(I16[9]+I16[11])/(I16[1]+I16[2]+I16[3]+I16[4])\\n' >> ${id}.site_mq.ref.info.tsv
    printf 'MqDecimalPlaces\\t${mq_dec}\\n' >> ${id}.site_mq.ref.info.tsv

    echo "Site MQ completed: calls=${id}.site_mq.calls.tsv.gz ref=\${ref_path}"
    """
}

process annotate_subgenome_variant_mq {
    tag "annotate variant MQ: ${chr}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.mq.info.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.mq.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), path("${id}.mq.info.tsv"), emit: mq
    path "${chr}.mq.log", emit: log

    script:
    def mq_workers = params.mq_tabix_workers ?: 0
    def mq_workers_py = (mq_workers as int) > 0 ? "${mq_workers}" : "None"
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${chr}.mq.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.mq import annotate_variants_mq_from_pvar

    workers = ${mq_workers_py}
    print(f"Annotating variant MQ for ${id} (${chr}) from ${params.mq_dir} (tabix, workers={workers}) ...")
    annotate_variants_mq_from_pvar(
        "${pvar}",
        "${params.mq_dir}",
        "${id}.mq.info.tsv",
        max_workers=workers,
        use_tabix=True,
    )
    """
}
