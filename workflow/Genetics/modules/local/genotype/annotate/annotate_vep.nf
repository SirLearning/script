nextflow.enable.dsl=2

/*
 * Ensembl VEP against a local GFF (no Ensembl cache for wheat).
 * Runs in Conda env `anno` (environment_anno.yml) — never `stats`, whose Hail
 * pin holds JDK 11 and whose compilers break the VEP post-link script.
 *
 * Coordinate freeze: RefSeq v1.0 assembly with IWGSC v1.1 gene models. Use the
 * Lulab GFF (numeric chromosome IDs) against numeric main_raw VCFs; the official
 * HC GFF uses chr1A-style names and needs a chromosome map first.
 *
 * Do not add --offline: it makes VEP look for an Ensembl cache directory, which
 * does not exist for wheat. --gff plus --fasta is already an offline run.
 *
 * Output is a VEP `--tab` table consumed by coding_mask_from_vep.
 */

process vep_prepare_gff {
    tag "vep gff: ${gff.name}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/annotate", mode: 'copy', pattern: "*.gff3.gz*"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/anno"
    label 'cpus_8'

    input:
    path gff

    output:
    tuple path("${params.mod}.vep.gff3.gz"), path("${params.mod}.vep.gff3.gz.tbi"), emit: gff
    path "*.log", emit: logs

    script:
    """
    set -euo pipefail
    exec > ${params.mod}.vep_gff.log 2>&1

    # Segment 0 is renamed because VEP does `unless \$chr`, and Perl reads the
    # string "0" as false, so every unplaced variant aborts the run.
    grep -v '^#' ${gff} \\
        | awk -F'\\t' -v OFS='\\t' -v unplaced='${params.vep_unplaced_name}' \\
            '\$1 == "0" { \$1 = unplaced } { print }' \\
        | sort -t\$'\\t' -k1,1 -k4,4n -k5,5n \\
        | bgzip -c > ${params.mod}.vep.gff3.gz
    tabix -p gff ${params.mod}.vep.gff3.gz
    """
}

/*
 * Sites-only VCF for VEP, built from a pfile's variant IDs.
 *
 * The *_refchr recode rewrote #CHROM/POS onto linear whole chromosomes, but left
 * the variant ID as <segment>-<segment_pos> — the coordinate system the Lulab GFF
 * and abd_iwgscV1.fa.gz use. Annotating the pvar's own #CHROM/POS would silently
 * place every variant against the wrong gene models, so read coordinates back out
 * of the ID instead. See domain/ref_v1_resources.yaml.
 *
 * Segments absent from the GFF are dropped. VEP aborts the whole run on the first
 * variant whose chromosome it cannot find, and the `Others` pfile carries chrUn
 * plus the two organelles (segments 0, 43, 44) while the GFF stops at 42.
 */
process pvar_sites_vcf {
    tag "sites vcf: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/annotate", mode: 'copy', pattern: "*.sites.vcf.gz"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/anno"
    label 'cpus_8'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    tuple path(gff), path(gff_tbi)

    output:
    tuple val(id), path("${id}.sites.vcf.gz"), emit: vcf
    path "*.log", emit: logs

    script:
    """
    set -euo pipefail
    exec > ${id}.sites_vcf.log 2>&1

    tabix --list-chroms ${gff} | sort -u > annotatable_segments.txt
    echo "segments in GFF: \$(wc -l < annotatable_segments.txt)"

    { echo '##fileformat=VCFv4.2'
      printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n'
      grep -v '^#' ${pvar} \\
        | awk -F'\\t' -v OFS='\\t' -v unplaced='${params.vep_unplaced_name}' '
            NR == FNR { annotatable[\$1] = 1; next }
            {
                if (split(\$3, seg, "-") != 2 || seg[1] !~ /^[0-9]+\$/ || seg[2] !~ /^[0-9]+\$/) {
                    bad++
                    next
                }
                chrom = (seg[1] == "0") ? unplaced : seg[1]
                if (!(chrom in annotatable)) {
                    dropped[seg[1]]++
                    next
                }
                print chrom, seg[2], \$3, \$4, \$5, ".", ".", "."
            }
            END {
                if (bad > 0) print "skipped " bad " ids not shaped <segment>-<pos>" > "/dev/stderr"
                for (s in dropped) print "dropped segment " s ": " dropped[s] " variants absent from GFF" > "/dev/stderr"
            }' annotatable_segments.txt - \\
        | sort -t\$'\\t' -k1,1 -k2,2n -S 4G
    } | bgzip -c -@ ${task.cpus} > ${id}.sites.vcf.gz

    echo "sites written: \$(bgzip -dc ${id}.sites.vcf.gz | grep -vc '^#')"
    """
}

process vep_annotate_vcf {
    tag "vep: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.vep.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.{log,html}"
    conda "${params.user_dir}/miniconda3/envs/anno"
    label 'cpus_32'

    input:
    tuple val(id), path(vcf)
    tuple path(gff), path(gff_tbi)
    // The .fai and .gzi must be staged next to the BGZF FASTA for VEP to seek in it.
    tuple path(fasta), path(fasta_fai), path(fasta_gzi)

    output:
    tuple val(id), path("${id}.vep.tsv"), emit: table
    path "*.vep_summary.html", emit: summary, optional: true
    path "*.log", emit: logs

    script:
    def fields = params.vep_fields ? "--fields ${params.vep_fields}" : ''
    def extra = params.vep_extra_args ?: ''
    """
    set -euo pipefail
    exec > ${id}.vep.log 2>&1

    vep --input_file ${vcf} \\
        --output_file ${id}.vep.tsv \\
        --stats_file ${id}.vep_summary.html \\
        --gff ${gff} \\
        --fasta ${fasta} \\
        --species ${params.vep_species} \\
        --tab \\
        --force_overwrite \\
        --buffer_size ${params.vep_buffer_size} \\
        --fork ${task.cpus} \\
        ${fields} ${extra}
    """
}
