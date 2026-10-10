#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import subprocess
from datetime import datetime
from pathlib import Path


ALGORITHM = "rde26-cmop"
TRACK = "CMOP"
DEFAULT_RUNS = 15
DEFAULT_SEED_BASE = 20260608
METRICS = ("igd", "score_cv", "mean_cv", "feasible_rate")


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
    parser = argparse.ArgumentParser(description="Run the RDE26-CMOP algorithm.")
    parser.add_argument("--out-dir", type=Path, default=root / "outputs" / f"{ALGORITHM}_{stamp}")
    parser.add_argument("--problems", default="1:15")
    parser.add_argument("--runs", type=int, default=DEFAULT_RUNS)
    parser.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE)
    parser.add_argument("--N", type=int, default=100)
    parser.add_argument("--D", type=int, default=30)
    parser.add_argument("--maxFE", type=int, default=200000)
    parser.add_argument("--save", type=int, default=1000)
    parser.add_argument("--pf-root", type=Path)
    parser.add_argument("--paper-id", default="RDE26")
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--skip-build", action="store_true")
    return parser.parse_args()


def build_binary(root: Path, cxx: str) -> Path:
    native = root / "rde_cmop_native"
    build = native / "build"
    build.mkdir(parents=True, exist_ok=True)
    binary = build / "rde_cmop_native"
    subprocess.run(
        [
            cxx,
            "-O3",
            "-std=c++17",
            "rde_cmop_native.cpp",
            "../../../benchmark/cpp/cmop_cec2026/cmop_cec2026_test.cpp",
            "-o",
            str(binary),
        ],
        cwd=str(native),
        check=True,
    )
    return binary


def parse_value(text: str) -> float:
    if text.lower() in {"nan", "-nan"}:
        return math.nan
    return float(text)


def parse_trace(log_text: str, rows: int) -> list[list[float]]:
    series: list[list[float]] = []
    pattern = re.compile(r"([A-Za-z_]+)=([+-]?(?:nan|inf|[0-9.eE+-]+))")
    for line in log_text.splitlines():
        if not line.startswith("TRACE "):
            continue
        fields = {key: value for key, value in pattern.findall(line)}
        series.append([parse_value(fields.get(metric, "nan")) for metric in METRICS])
    if not series:
        raise RuntimeError("empty CMOP trace")
    if len(series) >= rows:
        return series[:rows]
    return series + [series[-1]] * (rows - len(series))


def format_float(value: float) -> str:
    if math.isnan(value):
        return "nan"
    return f"{value:.17g}"


def write_matrix(path: Path, run_series: list[list[list[float]]]) -> None:
    rows = len(run_series[0])
    if any(len(series) != rows for series in run_series):
        raise RuntimeError(f"inconsistent trace length for {path.name}")
    with path.open("w", encoding="utf-8") as handle:
        for row in range(rows):
            values: list[str] = []
            for series in run_series:
                values.extend(format_float(value) for value in series[row])
            handle.write("\t".join(values))
            handle.write("\n")


def main() -> int:
    args = parse_args()
    if args.runs < 1:
        raise SystemExit("--runs must be positive")

    root = Path(__file__).resolve().parent
    binary = root / "rde_cmop_native" / "build" / "rde_cmop_native"
    if not args.skip_build or not binary.is_file():
        binary = build_binary(root, args.cxx)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.out_dir / "_raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    problems = parse_id_range(args.problems)

    for problem_id in problems:
        if problem_id < 1 or problem_id > 15:
            raise SystemExit(f"CMOP problem id out of range: {problem_id}")
        run_series: list[list[list[float]]] = []
        for run_id in range(1, args.runs + 1):
            seed = args.seed_base + problem_id * 1000 + run_id
            command = [
                str(binary),
                "--problem",
                str(problem_id),
                "--run",
                str(run_id),
                "--seed",
                str(seed),
                "--N",
                str(args.N),
                "--D",
                str(args.D),
                "--maxFE",
                str(args.maxFE),
                "--save",
                str(args.save),
            ]
            if args.pf_root is None:
                command.append("--noIGD")
            else:
                command.extend(["--pfRoot", str(args.pf_root)])
            completed = subprocess.run(
                command,
                cwd=str(root / "rde_cmop_native"),
                text=True,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            log_path = raw_dir / f"CMOP{problem_id}_run{run_id}.log"
            log_path.write_text(completed.stdout, encoding="utf-8")
            run_series.append(parse_trace(completed.stdout, args.save))
        target = args.out_dir / f"{args.paper_id}_CEC26_{TRACK}_F{problem_id}.txt"
        write_matrix(target, run_series)
        print(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
