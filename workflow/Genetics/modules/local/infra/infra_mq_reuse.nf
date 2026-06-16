nextflow.enable.dsl=2

// Frozen site-MQ reference under test_plink/process/abstract_mq_50_bams (one-off partial build; reuse thereafter).

def hasSiteMqRef(mq_dir, id) {
    file("${mq_dir}/reference/${id}.site_mq.ref.txt.gz").exists() ||
        file("${mq_dir}/reference/${id}.site_mq.ref.txt").exists()
}

def siteMqRefPath(mq_dir, id) {
    def gz = file("${mq_dir}/reference/${id}.site_mq.ref.txt.gz")
    if (gz.exists()) {
        return gz
    }
    return file("${mq_dir}/reference/${id}.site_mq.ref.txt")
}

def siteMqCallsPath(mq_dir, id) {
    return file("${mq_dir}/reference/${id}.site_mq.calls.tsv.gz")
}

def countSiteMqRefs(mq_dir) {
    def ref_dir = new File("${mq_dir}/reference")
    if (!ref_dir.isDirectory()) {
        return 0
    }
    ref_dir.list({ _dir, name -> name.endsWith('.site_mq.ref.txt.gz') } as FilenameFilter)?.length ?: 0
}

// True when process_dir/variant already has per-subgenome MQ annotation (*.mq.info.tsv).
def hasVariantMqForMergedTests(process_dir) {
    if (!process_dir) {
        return false
    }
    def subgenomes = ['A', 'B', 'D', 'Others']
    return subgenomes.every { sg ->
        file("${process_dir}/variant/${sg}.mq.info.tsv").exists()
    }
}

def siteMqRefTabixBgzPath(mq_dir, id) {
    return file("${mq_dir}/reference/${id}.site_mq.ref.txt.bgz")
}

def hasSiteMqRefTabix(mq_dir, id) {
    return siteMqRefTabixBgzPath(mq_dir, id).exists() &&
        file("${mq_dir}/reference/${id}.site_mq.ref.txt.bgz.tbi").exists()
}
