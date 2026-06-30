#!/usr/bin/env python3
"""Migrate doc/NF_CMD.md ### blocks into doc/progress/YYYY-MM-DD.md NF replay sections."""

from __future__ import annotations

import re
import sys
from collections import OrderedDict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "doc" / "NF_CMD.md"
OUT_DIR = REPO_ROOT / "doc" / "progress"

HEADING_RE = re.compile(r"^### (\d{4}-\d{2}-\d{2}) — (.+)$", re.MULTILINE)


def _strip_separators(body: str) -> str:
    body = body.strip()
    if body.startswith("---"):
        body = body[3:].lstrip("\n")
    if body.endswith("---"):
        body = body[:-3].rstrip("\n")
    return body.strip()


def parse_nf_cmd(text: str) -> OrderedDict[str, list[tuple[str, str]]]:
    matches = list(HEADING_RE.finditer(text))
    if not matches:
        raise ValueError("No ### YYYY-MM-DD blocks found")

    by_date: OrderedDict[str, list[tuple[str, str]]] = OrderedDict()
    for i, match in enumerate(matches):
        date = match.group(1)
        title = match.group(2).strip()
        start = match.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        body = _strip_separators(text[start:end])
        by_date.setdefault(date, []).append((title, body))
    return by_date


def format_replay_block(title: str, body: str) -> str:
    lines = ["<!-- nf-replay -->", f"#### NF replay — {title}", ""]
    if body:
        lines.append(body)
        lines.append("")
    lines.append("<!-- /nf-replay -->")
    return "\n".join(lines)


def append_nf_replays(
    source: Path = SOURCE,
    out_dir: Path = OUT_DIR,
    *,
    section_suffix: str = "NF replays",
) -> dict[str, int]:
    if not source.is_file():
        raise FileNotFoundError(source)

    by_date = parse_nf_cmd(source.read_text(encoding="utf-8"))
    out_dir.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}

    for date, blocks in by_date.items():
        replay_body = "\n\n".join(format_replay_block(t, b) for t, b in blocks)
        section = f"## {date} — {section_suffix}\n\n{replay_body}\n"
        out_path = out_dir / f"{date}.md"

        if out_path.is_file():
            existing = out_path.read_text(encoding="utf-8").rstrip()
            if f"## {date} — {section_suffix}" in existing:
                print(f"  skip {date}.md — NF replays section already present")
                continue
            content = f"{existing}\n\n---\n\n{section}"
        else:
            content = f"# Progress — {date}\n\n{section}"

        out_path.write_text(content, encoding="utf-8")
        counts[date] = len(blocks)

    return counts


def main() -> int:
    counts = append_nf_replays()
    total = sum(counts.values())
    print(f"Appended {total} NF replay blocks across {len(counts)} daily files:")
    for date, n in sorted(counts.items()):
        print(f"  {date}.md  (+{n})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
