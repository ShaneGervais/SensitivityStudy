#!/usr/bin/env python3
"""Run prepared PPN sweep directories in parallel."""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Execute every ppn.exe below nova_case/runs with bounded parallelism."
    )
    parser.add_argument("nova_case", type=Path, help="Nova case directory containing runs/")
    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=4,
        help="Maximum concurrent PPN jobs, default: 4",
    )
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=Path("runs"),
        help="Runs directory, absolute or relative to nova_case",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print run directories without executing ppn.exe",
    )
    return parser.parse_args()


def resolve_relative_to(base: Path, path: Path) -> Path:
    path = path.expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def ppn_env() -> dict[str, str]:
    env = os.environ.copy()
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        env.setdefault(key, "1")
    return env


def discover_runs(runs_dir: Path) -> list[Path]:
    return sorted(path.parent for path in runs_dir.rglob("ppn.exe") if path.is_file())


def run_one(run_dir: Path, env: dict[str, str]) -> tuple[Path, int, float]:
    start = time.time()
    log_path = run_dir / "run.log"
    with log_path.open("w") as log:
        process = subprocess.run(
            ["./ppn.exe"],
            cwd=run_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=False,
        )
    elapsed = time.time() - start
    return run_dir, process.returncode, elapsed


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")

    nova_case = args.nova_case.expanduser().resolve()
    runs_dir = resolve_relative_to(nova_case, args.runs_dir)
    if not runs_dir.is_dir():
        raise SystemExit(f"runs directory does not exist: {runs_dir}")

    run_dirs = discover_runs(runs_dir)
    if not run_dirs:
        raise SystemExit(f"no ppn.exe files found below {runs_dir}")

    print(f"runs_dir: {runs_dir}")
    print(f"runs: {len(run_dirs)}")
    print(f"jobs: {args.jobs}")

    if args.dry_run:
        for run_dir in run_dirs:
            print(run_dir)
        return 0

    env = ppn_env()
    failed: list[Path] = []
    start = time.time()
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_one, run_dir, env) for run_dir in run_dirs]
        for future in as_completed(futures):
            run_dir, returncode, elapsed = future.result()
            rel = run_dir.relative_to(runs_dir)
            status = "ok" if returncode == 0 else f"failed:{returncode}"
            print(f"{status}\t{elapsed:.1f}s\t{rel}")
            if returncode != 0:
                failed.append(run_dir)

    total = time.time() - start
    print(f"total elapsed: {total:.1f}s")
    if failed:
        print("failed runs:")
        for run_dir in failed:
            print(f"  - {run_dir}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
