#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


ALGORITHM = "rde26-csop"
TRACK = "CSOP"
RAW_ALGORITHM = "CORDEx_2026_Submission_R39"
DEFAULT_RUNS = 25
DEFAULT_SEED_BASE = 20260611


def parse_id_range(text: str) -> list[int]:
    items: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            left, right = part.split(":", 1)
            start = int(left)
            end = int(right)
            step = 1 if start <= end else -1
            items.extend(range(start, end + step, step))
        else:
            items.append(int(part))
    return items


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    parser = argparse.ArgumentParser(description="Run the RDE26-CSOP algorithm.")
    parser.add_argument("--out-dir", type=Path, default=root / "outputs" / f"{ALGORITHM}_{stamp}")
    parser.add_argument("--functions", default="1:28")
    parser.add_argument("--runs", type=int, default=DEFAULT_RUNS)
    parser.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE)
    parser.add_argument("--paper-id", default="RDE26")
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--skip-build", action="store_true")
    return parser.parse_args()


def build_binary(root: Path, cxx: str) -> Path:
    build_dir = root / "build"
    build_dir.mkdir(parents=True, exist_ok=True)
    binary = build_dir / "rde26_csop"
    command = [
        cxx,
        "-O3",
        "-std=c++11",
        "RDEx.cpp",
        "-o",
        str(binary),
    ]
    subprocess.run(command, cwd=str(root), check=True)
    return binary


def read_matrix(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.strip().split()
            if parts:
                rows.append([float(part) for part in parts])
    return rows


def write_matrix(path: Path, raw_rows: list[list[float]], runs: int) -> None:
    if len(raw_rows) < 2 * runs:
        raise RuntimeError(f"{path.name}: expected at least {2 * runs} rows, got {len(raw_rows)}")
    width = len(raw_rows[0])
    if width == 0:
        raise RuntimeError(f"{path.name}: empty raw row")
    if any(len(row) != width for row in raw_rows[: 2 * runs]):
        raise RuntimeError(f"{path.name}: inconsistent raw row width")
    start_col = 1 if width == 2001 else 0
    with path.open("w", encoding="utf-8") as handle:
        for col in range(start_col, width):
            values: list[str] = []
            for run_id in range(runs):
                values.append(f"{raw_rows[2 * run_id][col]:.17g}")
                values.append(f"{raw_rows[2 * run_id + 1][col]:.17g}")
            handle.write("\t".join(values))
            handle.write("\n")


def main() -> int:
    args = parse_args()
    if args.runs < 1:
        raise SystemExit("--runs must be positive")

    root = Path(__file__).resolve().parent
    binary = root / "build" / "rde26_csop"
    if not args.skip_build or not binary.is_file():
        binary = build_binary(root, args.cxx)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.out_dir / "_raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    input_dst = raw_dir / "inputData"
    if input_dst.exists():
        shutil.rmtree(input_dst)
    shutil.copytree(root / "inputData", input_dst)

    functions = parse_id_range(args.functions)
    if not functions:
        raise SystemExit("no CSOP functions selected")
    start = min(functions)
    end = max(functions)
    subprocess.run(
        [str(binary), str(start), str(end), str(args.runs), str(args.seed_base)],
        cwd=str(raw_dir),
        check=True,
    )

    for function_id in functions:
        if function_id < 1 or function_id > 28:
            raise SystemExit(f"CSOP function id out of range: {function_id}")
        raw_path = raw_dir / f"{RAW_ALGORITHM}_F{function_id}_D30.txt"
        target = args.out_dir / f"{args.paper_id}_CEC26_{TRACK}_F{function_id}.txt"
        write_matrix(target, read_matrix(raw_path), args.runs)
        print(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
