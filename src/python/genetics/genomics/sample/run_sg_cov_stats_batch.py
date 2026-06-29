"""Batch runner: sample.sg_cov stats — **do not use for publishDir**.

Use ``workflow/Genetics/subworkflows/tmp/sg_cov_replot.nf`` and Nextflow ``publishDir``
instead (workstation-core guardrail 15). This module remains for local /tmp debugging only.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from genetics.genomics.sample.cov import (
    compare_subgenome_ibs_depth_gam,
    sample_coverage_stats_bundle,
)
from genetics.germplasm.sample.mosdepth_sg_depth import sg_depth_taxa_bam_map_path

HOME = "/data/home/tusr1/01projects/vmap4"
PUBLISH_ROOT = Path("/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink")
GROUP_FILE = Path("/data1/dazheng_tusr1/vmap4.VCF.v1/meta_data/sample_group.txt")
JOBS = ("test_thin", "test_common_thin")
MODS = ("A", "B", "D", "Others")
TMP = Path("/tmp/sg_cov_stats_publish")


def _publish(src: Path, stats_dir: Path) -> None:
    name = src.name
    if name.endswith(".info.tsv"):
        dst_dir = stats_dir / "info"
    elif name.endswith(".th.tsv"):
        dst_dir = stats_dir / "thresholds"
    elif name.endswith(".png"):
        dst_dir = stats_dir / "plots"
    elif name.endswith(".log"):
        dst_dir = stats_dir / "logs"
    else:
        return
    dst_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst_dir / name)


def run_job(job: str) -> None:
    proc = PUBLISH_ROOT / "process" / job / "sample"
    stats_dir = PUBLISH_ROOT / "stats" / job
    job_tmp = TMP / job
    job_tmp.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {job} ===")

    ibs_infos: list[str] = []
    for mod in MODS:
        smiss = proc / f"{mod}.info.smiss"
        scount = proc / f"{mod}.info.scount"
        if not smiss.exists() or not scount.exists():
            print(f"[skip] missing smiss/scount for {mod}")
            continue

        prefix = job_tmp / f"{mod}.sample.sg_cov"
        sg_tbm = sg_depth_taxa_bam_map_path(HOME, mod)
        print(f"[run] {mod} sg_tbm={sg_tbm}")
        sample_coverage_stats_bundle(
            str(sg_tbm),
            str(smiss),
            str(scount),
            str(prefix),
            group_file=str(GROUP_FILE),
        )

        for src in job_tmp.glob(f"{mod}.sample.sg_cov*"):
            _publish(src, stats_dir)
        ibs_info = stats_dir / "info" / f"{mod}.sample.sg_cov.ibs_depth_miss.info.tsv"
        if mod in ("A", "B", "D") and ibs_info.exists():
            ibs_infos.append(str(ibs_info))

    abd_prefix = stats_dir / "plots" / "ABD.sample.sg_cov"
    if len(ibs_infos) == 3:
        compare_subgenome_ibs_depth_gam(ibs_infos, output_prefix=str(abd_prefix))
        for stem in (
            "gam_subgenome_summary.info.tsv",
            "gam_te_surface.panels.png",
            "gam_partial_ibs.panels.png",
            "gam_partial_logdepth.panels.png",
        ):
            src = Path(f"{abd_prefix}.{stem}")
            if src.exists():
                _publish(src, stats_dir)
                print(f"[ok] {src}")
    else:
        print(f"[warn] A/B/D sg_cov GAM skipped; ibs infos={ibs_infos}")


def main() -> None:
    for job in JOBS:
        run_job(job)


if __name__ == "__main__":
    main()
