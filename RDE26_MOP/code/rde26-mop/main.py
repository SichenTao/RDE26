#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from datetime import datetime
from pathlib import Path


ALGORITHM = "rde26-mop"
TRACK = "MOP"
DEFAULT_RUNS = 30
DEFAULT_SEED_BASE = 2025170000


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
    track_root = root.parent.parent
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    parser = argparse.ArgumentParser(description="Run the RDE26-MOP algorithm.")
    parser.add_argument("--out-dir", type=Path, default=root / "outputs" / f"{ALGORITHM}_{stamp}")
    parser.add_argument("--problems", default="1:10")
    parser.add_argument("--runs", type=int, default=DEFAULT_RUNS)
    parser.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE)
    parser.add_argument("--N", type=int, default=100)
    parser.add_argument("--M", type=int, default=3)
    parser.add_argument("--D", type=int, default=10)
    parser.add_argument("--maxFE", type=int, default=100000)
    parser.add_argument("--save", type=int, default=501)
    parser.add_argument("--pf-root", type=Path, default=track_root / "benchmark" / "SEC2018_MaOP_M3_D10" / "POF")
    parser.add_argument("--paper-id", default="RDE26")
    parser.add_argument("--skip-build", action="store_true")
    return parser.parse_args()


def build_binary(root: Path) -> Path:
    subprocess.run(["make", "-C", "src"], cwd=str(root), check=True)
    return root / "src" / "rde2026_mop_r05_t0070"


def read_igd(path: Path) -> list[float]:
    values: list[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("index"):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                values.append(float(parts[1]))
    return values


def normalize(values: list[float], rows: int) -> list[float]:
    if not values:
        raise RuntimeError("empty IGD trace")
    if len(values) >= rows:
        return values[-rows:]
    return values + [values[-1]] * (rows - len(values))


def write_matrix(path: Path, columns: list[list[float]]) -> None:
    rows = len(columns[0])
    if any(len(col) != rows for col in columns):
        raise RuntimeError(f"inconsistent trace length for {path.name}")
    with path.open("w", encoding="utf-8") as handle:
        for row in range(rows):
            handle.write("\t".join(f"{col[row]:.17g}" for col in columns))
            handle.write("\n")


def main() -> int:
    args = parse_args()
    if args.runs < 1:
        raise SystemExit("--runs must be positive")

    root = Path(__file__).resolve().parent
    binary = root / "src" / "rde2026_mop_r05_t0070"
    if not args.skip_build or not binary.is_file():
        binary = build_binary(root)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.out_dir / "_raw"
    problems = parse_id_range(args.problems)
    rows = max(1, args.save - 1)

    for problem_id in problems:
        if problem_id < 1 or problem_id > 10:
            raise SystemExit(f"MOP problem id out of range: {problem_id}")
        columns: list[list[float]] = []
        for run_id in range(1, args.runs + 1):
            seed = args.seed_base + problem_id * 1000 + run_id
            subprocess.run(
                [
                    str(binary),
                    "--problem",
                    str(problem_id),
                    "--run",
                    str(run_id),
                    "--N",
                    str(args.N),
                    "--M",
                    str(args.M),
                    "--D",
                    str(args.D),
                    "--maxFE",
                    str(args.maxFE),
                    "--save",
                    str(args.save),
                    "--seed",
                    str(seed),
                    "--pf-root",
                    str(args.pf_root),
                    "--out-dir",
                    str(raw_dir),
                    "--alg-name",
                    ALGORITHM,
                ],
                cwd=str(root),
                check=True,
            )
            trace_path = raw_dir / "Data" / ALGORITHM / (
                f"{ALGORITHM}_MaOP{problem_id}_M{args.M}_D{args.D}_{run_id}_igd.tsv"
            )
            columns.append(normalize(read_igd(trace_path), rows))
        target = args.out_dir / f"{args.paper_id}_CEC26_{TRACK}_F{problem_id}.txt"
        write_matrix(target, columns)
        print(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
