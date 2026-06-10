nextflow.enable.dsl=2

/*
 * Upload large raw sequencing files to BIG GSA FTP (submit.big.ac.cn) with resume
 * and retry. One Nextflow task per file; uses curl -C - to continue partial uploads.
 *
 * Source defaults: KSXY-YCS-DO25090403 HiFi BAMs, Illumina PE, and Hi-C FASTQs on s243.
 * Remote target: /GSA/subCRA071073
 *
 * Credentials: export GSA_FTP_PASSWORD before launch (do not pass on the CLI).
 * Config: -c /path/to/workflow/Genetics/nextflow.config (see GENETICS_WORKFLOW.md).
 */

params.ftp_host = 'submit.big.ac.cn'
params.ftp_remote_dir = '/GSA/subCRA071073'
params.ftp_user = 'flu@genetics.ac.cn'
params.ftp_password = null

params.hifi_dir = '/data/home/assembly/01_data/KSXY-YCS-DO25090403-1P'
params.pe_dir = '/data/home/assembly/01_data/KSXY-YCS-DO25090403-2P'
params.hic_dir = '/data/home/assembly/01_data/KSXY-YCS-DO25090403-3P'

params.ftp_max_attempts = 50
params.ftp_retry_sleep_sec = 120
params.ftp_stall_time_sec = 600
params.ftp_stall_limit_bps = 1000
params.upload_only = null
params.test_upload = false
params.verify_md5 = true
params.verify_only = false

def helpMessage() {
    log.info """
    GSA FTP upload (resume + retry)

    Uploads HiFi BAM, Illumina PE, and Hi-C FASTQ files to BIG submit FTP.
    Each file runs in its own task; curl resumes from partial remote copies.
    Optional MD5 verification compares local md5sum with remote content streamed via curl.

    Required:
      params.output_dir          Log / summary publish root
      GSA_FTP_PASSWORD env var   FTP password (preferred; not logged in NF_CMD)

    Optional overrides:
      --ftp_host, --ftp_remote_dir, --ftp_user
      --hifi_dir, --pe_dir, --hic_dir
      --ftp_max_attempts, --ftp_retry_sleep_sec
      --verify_md5               Run MD5 check after upload (default: true)
      --verify_only              Skip upload; only compare local vs remote MD5
      --test_upload              Upload or verify only the smallest HiFi BAM (smoke test)

    Example (s243, run from a dedicated run folder, not the repo):
      export GSA_FTP_PASSWORD='YOUR_PASSWORD'
      screen -dmS gsa_upload bash -c "\\
        cd /data/home/assembly/01_data/gsa_upload/01run && \\
        source ~/.bashrc && conda activate run && \\
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/ops/gsa_ftp_upload.nf \\
          -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \\
          --output_dir /data/home/assembly/01_data/gsa_upload/01run \\
          -resume "

    screen -dmS gsa_upload bash -c "\
        cd /data/home/tusr1/01projects/gsa_upload/02run && \
        source ~/.bashrc && conda activate run && \
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/ops/gsa_ftp_upload.nf \
          -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
          --output_dir /data/home/tusr1/01projects/gsa_upload/02run \
          -resume "

    MD5 verify only (no upload):
      nextflow run .../gsa_ftp_upload.nf ... --verify_only
    """
}

def gsaFileChannel() {
    def pe_dir = params.pe_dir
    def hic_dir = params.hic_dir
    def hifi_dir = params.hifi_dir

    if (params.upload_only) {
        def only_path = file(params.upload_only as String, checkIfExists: true)
        return channel.of(tuple('manual', only_path.toString(), only_path.name))
    }
    if (params.test_upload) {
        def test_bam = file("${hifi_dir}/m84178_250925_102944_s3.hifi_reads.bc2020.bam", checkIfExists: true)
        return channel.of(tuple('hifi', test_bam.toString(), test_bam.name))
    }

    def ch_hifi = channel.fromPath("${hifi_dir}/*bam", checkIfExists: true)
        .map { f -> tuple('hifi', f.toString(), f.name) }

    def ch_pe = channel.of(
        tuple('pe', file("${pe_dir}/JM229_pe_1.fq.gz", checkIfExists: true).toString(), 'JM229_pe_1.fq.gz'),
        tuple('pe', file("${pe_dir}/JM229_pe_2.fq.gz", checkIfExists: true).toString(), 'JM229_pe_2.fq.gz'),
    )

    def ch_hic = channel.of(
        tuple('hic', file("${hic_dir}/JM229-HiC_1.fq.gz", checkIfExists: true).toString(), 'JM229-HiC_1.fq.gz'),
        tuple('hic', file("${hic_dir}/JM229-HiC_2.fq.gz", checkIfExists: true).toString(), 'JM229-HiC_2.fq.gz'),
    )

    return ch_hifi.mix(ch_pe).mix(ch_hic)
}

process gsa_ftp_upload_file {
    tag "upload_${file_name}"
    cpus 1
    memory '4.GB'
    maxForks 1
    errorStrategy 'retry'
    maxRetries 5
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(category), val(local_path), val(file_name)

    output:
    tuple val(category), path("${file_name}.upload.txt"), emit: results

    script:
    def ftp_host = params.ftp_host
    def ftp_user = params.ftp_user
    def remote_dir = params.ftp_remote_dir
    def ftp_password = params.ftp_password ?: System.getenv('GSA_FTP_PASSWORD') ?: ''
    """
    #!/usr/bin/env bash
    set -euo pipefail
    exec > "upload_${file_name}.log" 2>&1

    LOCAL_PATH="${local_path}"
    FILE_NAME="${file_name}"
    FTP_HOST="${ftp_host}"
    FTP_USER="${ftp_user}"
    FTP_REMOTE_DIR="${remote_dir}"
    GSA_FTP_PASSWORD='${ftp_password.replace("'", "'\"'\"'")}'
    MAX_ATTEMPTS=${params.ftp_max_attempts}
    RETRY_SLEEP=${params.ftp_retry_sleep_sec}
    STALL_TIME=${params.ftp_stall_time_sec}
    STALL_LIMIT=${params.ftp_stall_limit_bps}

    if [[ -z "\${GSA_FTP_PASSWORD}" ]]; then
        echo "ERROR: FTP password is not set (params.ftp_password or GSA_FTP_PASSWORD)."
        exit 1
    fi
    if [[ ! -f "\${LOCAL_PATH}" ]]; then
        echo "ERROR: local file not found: \${LOCAL_PATH}"
        exit 1
    fi

    LOCAL_SIZE=\$(stat -c%s "\${LOCAL_PATH}")
    REMOTE_URL="ftp://\${FTP_HOST}\${FTP_REMOTE_DIR}/\${FILE_NAME}"
    LIST_URL="ftp://\${FTP_HOST}\${FTP_REMOTE_DIR}/"

    get_remote_size() {
        curl -s \\
            -u "\${FTP_USER}:\${GSA_FTP_PASSWORD}" \\
            "\${LIST_URL}" 2>/dev/null | awk -v fn="\${FILE_NAME}" '\$NF == fn {print \$5; found=1} END {if (!found) print 0}'
    }

    echo "Category: ${category}"
    echo "Local file: \${LOCAL_PATH} (\${LOCAL_SIZE} bytes)"
    echo "Remote URL: \${REMOTE_URL}"

    REMOTE_SIZE=\$(get_remote_size || echo 0)
    echo "Remote size before upload: \${REMOTE_SIZE} bytes"

    if [[ "\${REMOTE_SIZE}" == "\${LOCAL_SIZE}" ]]; then
        echo "Remote file already complete; skipping upload."
        printf "%s\\t%s\\t%s\\t%s\\tSKIP\\n" "${category}" "\${LOCAL_PATH}" "\${FILE_NAME}" "\${LOCAL_SIZE}" > "${file_name}.upload.txt"
        exit 0
    fi

    attempt=1
    while [[ \${attempt} -le \${MAX_ATTEMPTS} ]]; do
        echo "Upload attempt \${attempt}/\${MAX_ATTEMPTS} (curl resume: -C -)"
        set +e
        curl -f \\
            -u "\${FTP_USER}:\${GSA_FTP_PASSWORD}" \\
            -T "\${LOCAL_PATH}" \\
            -C - \\
            --connect-timeout 60 \\
            --speed-time "\${STALL_TIME}" \\
            --speed-limit "\${STALL_LIMIT}" \\
            "\${REMOTE_URL}"
        curl_rc=\$?
        set -e

        REMOTE_SIZE=\$(get_remote_size || echo 0)
        echo "Remote size after attempt \${attempt}: \${REMOTE_SIZE} bytes"

        if [[ \${curl_rc} -eq 0 && "\${REMOTE_SIZE}" == "\${LOCAL_SIZE}" ]]; then
            echo "Upload complete."
            printf "%s\\t%s\\t%s\\t%s\\tPASS\\n" "${category}" "\${LOCAL_PATH}" "\${FILE_NAME}" "\${LOCAL_SIZE}" > "${file_name}.upload.txt"
            exit 0
        fi

        echo "Attempt \${attempt} incomplete (curl exit \${curl_rc}, remote \${REMOTE_SIZE} / local \${LOCAL_SIZE})."
        attempt=\$((attempt + 1))
        if [[ \${attempt} -le \${MAX_ATTEMPTS} ]]; then
            echo "Sleeping \${RETRY_SLEEP}s before retry..."
            sleep "\${RETRY_SLEEP}"
        fi
    done

    echo "ERROR: upload failed after \${MAX_ATTEMPTS} attempts."
    printf "%s\\t%s\\t%s\\t%s\\tFAIL\\n" "${category}" "\${LOCAL_PATH}" "\${FILE_NAME}" "\${LOCAL_SIZE}" > "${file_name}.upload.txt"
    exit 1
    """
}

process gsa_ftp_verify_md5 {
    tag "md5_${file_name}"
    cpus 1
    memory '2.GB'
    maxForks 1
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.output_dir}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(category), val(local_path), val(file_name)

    output:
    tuple val(category), path("${file_name}.md5verify.txt"), emit: results

    script:
    def ftp_host = params.ftp_host
    def ftp_user = params.ftp_user
    def remote_dir = params.ftp_remote_dir
    def ftp_password = params.ftp_password ?: System.getenv('GSA_FTP_PASSWORD') ?: ''
    """
    #!/usr/bin/env bash
    set -euo pipefail
    exec > "md5verify_${file_name}.log" 2>&1

    LOCAL_PATH="${local_path}"
    FILE_NAME="${file_name}"
    FTP_HOST="${ftp_host}"
    FTP_USER="${ftp_user}"
    FTP_REMOTE_DIR="${remote_dir}"
    GSA_FTP_PASSWORD='${ftp_password.replace("'", "'\"'\"'")}'

    if [[ -z "\${GSA_FTP_PASSWORD}" ]]; then
        echo "ERROR: FTP password is not set (params.ftp_password or GSA_FTP_PASSWORD)."
        exit 1
    fi
    if [[ ! -f "\${LOCAL_PATH}" ]]; then
        echo "ERROR: local file not found: \${LOCAL_PATH}"
        exit 1
    fi

    REMOTE_URL="ftp://\${FTP_HOST}\${FTP_REMOTE_DIR}/\${FILE_NAME}"
    LIST_URL="ftp://\${FTP_HOST}\${FTP_REMOTE_DIR}/"

    LOCAL_SIZE=\$(stat -c%s "\${LOCAL_PATH}")
    REMOTE_SIZE=\$(curl -s \\
        -u "\${FTP_USER}:\${GSA_FTP_PASSWORD}" \\
        "\${LIST_URL}" 2>/dev/null | awk -v fn="\${FILE_NAME}" '\$NF == fn {print \$5; found=1} END {if (!found) print 0}')

    echo "Category: ${category}"
    echo "Local file: \${LOCAL_PATH} (\${LOCAL_SIZE} bytes)"
    echo "Remote URL: \${REMOTE_URL}"
    echo "Remote size: \${REMOTE_SIZE} bytes"

    if [[ "\${REMOTE_SIZE}" != "\${LOCAL_SIZE}" ]]; then
        echo "ERROR: size mismatch (local \${LOCAL_SIZE} vs remote \${REMOTE_SIZE}); skipping MD5."
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\tFAIL\\n" \\
            "${category}" "\${LOCAL_PATH}" "\${FILE_NAME}" "\${LOCAL_SIZE}" "\${REMOTE_SIZE}" "-" "-" > "${file_name}.md5verify.txt"
        exit 1
    fi

    echo "Computing local MD5..."
    LOCAL_MD5=\$(md5sum "\${LOCAL_PATH}" | awk '{print \$1}')
    echo "Local MD5: \${LOCAL_MD5}"

    echo "Streaming remote file and computing MD5 (no local download)..."
    REMOTE_MD5=\$(curl -f -s \\
        -u "\${FTP_USER}:\${GSA_FTP_PASSWORD}" \\
        "\${REMOTE_URL}" | md5sum | awk '{print \$1}')
    echo "Remote MD5: \${REMOTE_MD5}"

    if [[ "\${LOCAL_MD5}" == "\${REMOTE_MD5}" ]]; then
        echo "MD5 match."
        STATUS=PASS
    else
        echo "ERROR: MD5 mismatch."
        STATUS=FAIL
    fi

    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${category}" "\${LOCAL_PATH}" "\${FILE_NAME}" "\${LOCAL_SIZE}" "\${REMOTE_SIZE}" "\${LOCAL_MD5}" "\${REMOTE_MD5}" "\${STATUS}" > "${file_name}.md5verify.txt"

    if [[ "\${STATUS}" != "PASS" ]]; then
        exit 1
    fi
    """
}

process collect_gsa_upload_results {
    tag "collect_uploads"
    publishDir "${params.output_dir}", mode: 'copy', pattern: "gsa_upload_summary.tsv"

    input:
    path result_files

    output:
    path "gsa_upload_summary.tsv"

    script:
    """
    set -euo pipefail
    exec > collect_gsa_upload_results.log 2>&1

    echo -e "category\\tlocal_path\\tfile_name\\tbytes\\tstatus" > gsa_upload_summary.tsv
    cat ${result_files} >> gsa_upload_summary.tsv

    fail_count=\$(awk -F'\\t' '\$NF == "FAIL" {c++} END {print c+0}' gsa_upload_summary.tsv)
    skip_count=\$(awk -F'\\t' '\$NF == "SKIP" {c++} END {print c+0}' gsa_upload_summary.tsv)
    pass_count=\$(awk -F'\\t' '\$NF == "PASS" {c++} END {print c+0}' gsa_upload_summary.tsv)
    echo "Summary: PASS=\${pass_count} SKIP=\${skip_count} FAIL=\${fail_count}"
    if [[ "\${fail_count}" != "0" ]]; then
        echo "One or more uploads failed."
        exit 1
    fi
    """
}

process collect_gsa_md5_results {
    tag "collect_md5"
    publishDir "${params.output_dir}", mode: 'copy', pattern: "gsa_md5_summary.tsv"

    input:
    path result_files

    output:
    path "gsa_md5_summary.tsv"

    script:
    """
    set -euo pipefail
    exec > collect_gsa_md5_results.log 2>&1

    echo -e "category\\tlocal_path\\tfile_name\\tlocal_bytes\\tremote_bytes\\tlocal_md5\\tremote_md5\\tstatus" > gsa_md5_summary.tsv
    cat ${result_files} >> gsa_md5_summary.tsv

    fail_count=\$(awk -F'\\t' '\$NF == "FAIL" {c++} END {print c+0}' gsa_md5_summary.tsv)
    pass_count=\$(awk -F'\\t' '\$NF == "PASS" {c++} END {print c+0}' gsa_md5_summary.tsv)
    echo "MD5 summary: PASS=\${pass_count} FAIL=\${fail_count}"
    if [[ "\${fail_count}" != "0" ]]; then
        echo "One or more MD5 checks failed."
        exit 1
    fi
    """
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    def ftp_pass = params.ftp_password ?: System.getenv('GSA_FTP_PASSWORD')
    if (!params.output_dir) {
        error 'gsa_ftp_upload.nf: params.output_dir is required (run folder for logs and summary).'
    }
    if (!params.ftp_user) {
        error 'gsa_ftp_upload.nf: params.ftp_user is required.'
    }
    if (!ftp_pass) {
        error 'gsa_ftp_upload.nf: set params.ftp_password in the script or export GSA_FTP_PASSWORD.'
    }

    def file_ch = gsaFileChannel()

    if (params.verify_only) {
        gsa_ftp_verify_md5(file_ch)
        collect_gsa_md5_results(gsa_ftp_verify_md5.out.results.map { _cat, res -> res }.collect())
    } else {
        gsa_ftp_upload_file(file_ch)
        collect_gsa_upload_results(gsa_ftp_upload_file.out.results.map { _cat, res -> res }.collect())

        if (params.verify_md5) {
            gsa_ftp_verify_md5(file_ch)
            collect_gsa_md5_results(gsa_ftp_verify_md5.out.results.map { _cat, res -> res }.collect())
        }
    }
}
