nextflow.enable.dsl=2

/*
 * GWAS Benchmark workflow
 *
 * Runs the Python CLI: python -m gwas.benchmark
 * - Input: a manifest (YAML or CSV) listing GWAS result files from different tools
 * - Optional: ground-truth file and matching window
 * - Output: metrics summary, labeled tables, plots, and best model recommendation
 *
 * Parameters:
 *   --manifest <path>           YAML/CSV manifest; if not provided, build from --methods_dir + --methods_glob
 *   --methods_dir <dir>         Directory with GWAS outputs to auto-build manifest (CSV)
 *   --methods_glob <glob>       Glob pattern to match files in methods_dir (default: *.txt)
 *   --truth <path>              Ground-truth file (overrides manifest truth)
 *   --truth_window <int>        Window flank size around point truths (bp)
 *   --topk <str>                Comma-separated K values for recall/enrichment (default: 10,20,50,100)
 *   --outdir <dir>              Output directory (default: results/gwas_benchmark)
 *   --python <exe>              Python executable (default: python)
 *   --pythonpath <dir>          Python module path to include (default: repo src/python)
 *
 * Usage examples:
 *   nextflow run Association/gwas/run/benchmark.nf --manifest bench.yaml --outdir results/bench
 *   nextflow run Association/gwas/run/benchmark.nf --methods_dir results/gwas --methods_glob "*.assoc.txt" --truth truth.tsv
 */

params.manifest = params.manifest ?: null
params.methods_dir = params.methods_dir ?: null
params.methods_glob = params.methods_glob ?: '*.txt'
params.truth = params.truth ?: null
params.truth_window = params.truth_window instanceof Integer ? params.truth_window : null
params.topk = params.topk ?: '10,20,50,100'
params.outdir = params.outdir ?: 'results/gwas_benchmark'
params.python = params.python ?: 'python'
// Default pythonpath: repo/src/python (relative to this file: ../../../src/python)
params.pythonpath = params.pythonpath ?: null

process BUILD_MANIFEST_FROM_DIR {
	tag "build-manifest"
	publishDir "${params.outdir}", mode: 'copy'

	input:
	tuple val(methods_dir), val(glob)

	output:
	path "manifest.csv", emit: manifest

	script:
	"""
	set -euo pipefail
	shopt -s nullglob
	echo "name,file" > manifest.csv
	for f in ${methods_dir}/${glob}; do
	  base=\$(basename "\$f")
	  name=\${base%.*}
	  echo "\$name,\$f" >> manifest.csv
	done
	echo "[manifest] Wrote \$(wc -l < manifest.csv) lines (incl. header)"
	"""
}

process RUN_GWAS_BENCHMARK {
	tag "gwas-benchmark"
	publishDir "${params.outdir}", mode: 'copy'
	cpus 2
	memory '4 GB'
	errorStrategy 'terminate'

	input:
	path manifest

	output:
	path "metrics_summary.tsv"
	path "best_model.txt"
	path "*.png"
	path "*.tsv"
	path "benchmark.log"

	script:
	// Helper to resolve default PYTHONPATH if not provided
	// TODO: PY is not phenotype
	def defaultPythonPath = new File("${baseDir}/../../../src/python").canonicalPath
	def PY_PATH = params.pythonpath ?: defaultPythonPath
	def truth = params.truth ? "--truth ${params.truth}" : ''
	def tw = params.truth_window ? "--truth-window ${params.truth_window}" : ''
	def topk = params.topk ? "--topk ${params.topk}" : ''
	"""
	set -euo pipefail
	echo "[GWAS] Using PYTHONPATH=${PY_PATH}" | tee benchmark.log
	export PYTHONPATH="${PY_PATH}:\${PYTHONPATH:-}"
	${params.python} -m gwas.benchmark \\
	  --manifest ${manifest} \\
	  --outdir . \\
	  ${truth} ${tw} ${topk} \\
	  >> benchmark.log 2>&1
	"""
}

workflow {
	channel
		.of(1)
		.set { start }

	ch_manifest = channel.empty()

	if (params.manifest) {
		ch_manifest = channel.fromPath(params.manifest)
	} else if (params.methods_dir) {
		def ch_opts = channel.of(tuple(params.methods_dir as String, params.methods_glob as String))
		ch_manifest = BUILD_MANIFEST_FROM_DIR(ch_opts).manifest
	} else {
		log.error "Provide either --manifest or --methods_dir"
		System.exit(1)
	}

	RUN_GWAS_BENCHMARK(ch_manifest)
}

