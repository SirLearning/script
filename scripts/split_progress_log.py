#!/usr/bin/env python3
"""Split monolithic doc/TODO_PROGRESS_LOG.md into doc/progress/YYYY-MM-DD.md files."""

from __future__ import annotations

import re
import sys
from collections import OrderedDict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "doc" / "TODO_PROGRESS_LOG.md"
OUT_DIR = REPO_ROOT / "doc" / "progress"

HEADING_RE = re.compile(r"^## (\d{4}-\d{2}-\d{2}) — ", re.MULTILINE)


def split_log(source: Path = SOURCE, out_dir: Path = OUT_DIR) -> dict[str, int]:
    if not source.is_file():
        raise FileNotFoundError(f"Source log not found: {source}")

    text = source.read_text(encoding="utf-8")
    matches = list(HEADING_RE.finditer(text))
    if not matches:
        raise ValueError(f"No dated sections found in {source}")

    by_date: OrderedDict[str, list[str]] = OrderedDict()

    for i, match in enumerate(matches):
        date = match.group(1)
        start = match.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        section = text[start:end].rstrip()
        by_date.setdefault(date, []).append(section)

    out_dir.mkdir(parents=True, exist_ok=True)

    counts: dict[str, int] = {}
    for date, sections in by_date.items():
        body = "\n\n".join(sections)
        out_path = out_dir / f"{date}.md"
        content = f"# Progress — {date}\n\n{body}\n"
        out_path.write_text(content, encoding="utf-8")
        counts[date] = len(sections)

    return counts


def main() -> int:
    counts = split_log()
    total = sum(counts.values())
    print(f"Split {total} sections into {len(counts)} daily files under {OUT_DIR}:")
    for date, n in counts.items():
        print(f"  {date}.md  ({n} sections)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
