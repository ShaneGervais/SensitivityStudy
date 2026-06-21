#!/usr/bin/env python3
"""Build one-rate-at-a-time PPN run directories from a reaction plan."""

from __future__ import annotations

import argparse
import json
import re
import shutil
from decimal import Decimal
from pathlib import Path


REQUIRED_RUNTIME_FILES = (
    "ppn.exe",
    "ppn_frame.input",
    "ppn_physics.input",
    "ppn_solver.input",
    "networksetup.txt",
    "isotopedatabase.txt",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create nova_case/runs from a PPN template and reaction_plan.json. "
            "Each leaf run directory is self-contained except for ../NPDATA."
        )
    )
    parser.add_argument("nova_case", type=Path, help="Nova case directory containing ppn/")
    parser.add_argument(
        "--reaction-plan",
        type=Path,
        default=Path("config/reaction_plan.json"),
        help="Path to reaction_plan.json, absolute or relative to nova_case",
    )
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=Path("runs"),
        help="Output runs directory, absolute or relative to nova_case",
    )
    parser.add_argument(
        "--baseline-only",
        action="store_true",
        help="Only create runs/baseline",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and print what would be built without writing files",
    )
    return parser.parse_args()


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_case_path(path: Path) -> Path:
    return path.expanduser().resolve()


def resolve_relative_to(base: Path, path: Path) -> Path:
    path = path.expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def read_plan(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def parse_ppn_frame_paths(ppn_frame: Path) -> list[Path]:
    text = ppn_frame.read_text()
    paths: list[Path] = []
    for key in ("trajectory_fn", "ini_filename"):
        match = re.search(rf"(?im)^\s*{key}\s*=\s*['\"]([^'\"]+)['\"]", text)
        if match:
            paths.append(Path(match.group(1)))
    return paths


def find_referenced_file(nova_case: Path, ppn_dir: Path, rel_path: Path) -> Path | None:
    if rel_path.is_absolute():
        return rel_path if rel_path.exists() else None

    for base in (ppn_dir, nova_case):
        candidate = base / rel_path
        if candidate.exists():
            return candidate
    return None


def collect_runtime_sources(nova_case: Path) -> tuple[list[tuple[Path, Path]], list[str]]:
    ppn_dir = nova_case / "ppn"
    missing: list[str] = []
    sources: dict[Path, Path] = {}

    for name in REQUIRED_RUNTIME_FILES:
        src = ppn_dir / name
        if src.exists():
            sources[src] = Path(name)
        else:
            missing.append(str(src))

    for src in sorted(ppn_dir.glob("isotopedatabase*.txt")):
        sources[src] = Path(src.name)

    for src in sorted(ppn_dir.glob("*.input")):
        sources[src] = Path(src.name)

    frame = ppn_dir / "ppn_frame.input"
    if frame.exists():
        for rel_path in parse_ppn_frame_paths(frame):
            src = find_referenced_file(nova_case, ppn_dir, rel_path)
            if src is None:
                missing.append(
                    f"{rel_path} referenced by {frame}; searched {ppn_dir} and {nova_case}"
                )
            else:
                sources[src] = rel_path

    return sorted(sources.items(), key=lambda item: str(item[1])), missing


def ensure_npdata_link(link_path: Path, target: Path, dry_run: bool) -> str:
    if link_path.is_symlink():
        current = link_path.resolve(strict=False)
        if current == target:
            return f"ok: {link_path} -> {target}"
        if not dry_run:
            link_path.unlink()
    elif link_path.exists():
        return f"kept existing non-symlink NPDATA: {link_path}"

    if not dry_run:
        link_path.parent.mkdir(parents=True, exist_ok=True)
        link_path.symlink_to(target, target_is_directory=True)
    return f"linked: {link_path} -> {target}"


def strip_rate_factors(text: str) -> str:
    kept = []
    for line in text.splitlines():
        if re.match(r"\s*rate_(index|factor)\s*\(", line):
            continue
        kept.append(line)
    return "\n".join(kept) + "\n"


def insert_rate_factor(text: str, rate_index: int, rate_factor: str) -> str:
    lines = strip_rate_factors(text).splitlines()
    for i, line in enumerate(lines):
        if line.strip() == "/":
            lines[i:i] = [
                f"        rate_index(1) = {rate_index}",
                f"        rate_factor(1) = {rate_factor}",
            ]
            return "\n".join(lines) + "\n"
    raise ValueError("ppn_physics.input does not contain a closing '/' namelist line")


def factor_label(value: int | float | str) -> str:
    decimal = Decimal(str(value)).normalize()
    return format(decimal, "f")


def fortran_factor(value: int | float | str) -> str:
    label = factor_label(value)
    if "." in label:
        return f"{label}d0"
    return f"{label}.d0"


def copy_runtime_files(
    sources: list[tuple[Path, Path]],
    destination: Path,
    dry_run: bool,
) -> None:
    if dry_run:
        return
    destination.mkdir(parents=True, exist_ok=True)
    for src, rel_dst in sources:
        dst = destination / rel_dst
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def write_ppn_physics(
    run_dir: Path,
    source_ppn_physics: Path,
    rate_index: int | None,
    rate_factor: str | None,
    dry_run: bool,
) -> None:
    if dry_run:
        return
    text = source_ppn_physics.read_text()
    if rate_index is None:
        updated = strip_rate_factors(text)
    else:
        if rate_factor is None:
            raise ValueError("rate_factor is required when rate_index is set")
        updated = insert_rate_factor(text, rate_index, rate_factor)
    (run_dir / "ppn_physics.input").write_text(updated)


def reaction_dir_name(reaction_data: dict, fallback_name: str) -> str:
    reaction_id = reaction_data.get("reaction_id")
    if reaction_id:
        return str(reaction_id)
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", fallback_name).strip("_")


def planned_runs(plan: dict, baseline_only: bool) -> list[tuple[str, Path, int | None, str | None]]:
    runs: list[tuple[str, Path, int | None, str | None]] = [
        ("baseline", Path("baseline"), None, None)
    ]
    if baseline_only:
        return runs

    for reaction_name, reaction_data in plan.get("reactions", {}).items():
        selected_index = reaction_data.get("selected_rate_index")
        if selected_index is None:
            continue
        reaction_dir = reaction_dir_name(reaction_data, reaction_name)
        for factor in reaction_data.get("factors", []):
            runs.append(
                (
                    reaction_name,
                    Path(reaction_dir) / f"fact_{factor_label(factor)}",
                    int(selected_index),
                    fortran_factor(factor),
                )
            )
    return runs


def main() -> int:
    args = parse_args()
    nova_case = resolve_case_path(args.nova_case)
    ppn_dir = nova_case / "ppn"
    plan_path = resolve_relative_to(nova_case, args.reaction_plan)
    runs_dir = resolve_relative_to(nova_case, args.runs_dir)
    project_root = project_root_from_script()
    npdata_target = (project_root / "physics" / "NPDATA").resolve()

    if not nova_case.is_dir():
        raise SystemExit(f"nova_case does not exist: {nova_case}")
    if not ppn_dir.is_dir():
        raise SystemExit(f"nova_case is missing ppn/: {ppn_dir}")
    if not plan_path.is_file():
        raise SystemExit(f"reaction plan does not exist: {plan_path}")
    if not npdata_target.exists():
        raise SystemExit(f"NPDATA target does not exist: {npdata_target}")

    plan = read_plan(plan_path)
    sources, missing = collect_runtime_sources(nova_case)
    if missing:
        message = "\n".join(f"  - {item}" for item in missing)
        raise SystemExit(f"Cannot build runs because required runtime files are missing:\n{message}")

    source_ppn_physics = ppn_dir / "ppn_physics.input"
    runs = planned_runs(plan, args.baseline_only)

    skipped = [
        data.get("reaction_id", name)
        for name, data in plan.get("reactions", {}).items()
        if data.get("selected_rate_index") is None
    ]
    if skipped and not args.baseline_only:
        print("Skipping reactions without selected_rate_index:")
        for reaction in skipped:
            print(f"  - {reaction}")

    print(f"nova_case: {nova_case}")
    print(f"reaction_plan: {plan_path}")
    print(f"runs_dir: {runs_dir}")
    print(f"planned leaf runs: {len(runs)}")

    # Preserve the conventional case-level NPDATA link and also create the links
    # needed by PPN's ../NPDATA lookup from each run directory.
    print(ensure_npdata_link(nova_case / "NPDATA", npdata_target, args.dry_run))

    linked_parents: set[Path] = set()
    for _, rel_run_dir, rate_index, rate_factor in runs:
        run_dir = runs_dir / rel_run_dir
        copy_runtime_files(sources, run_dir, args.dry_run)
        write_ppn_physics(run_dir, source_ppn_physics, rate_index, rate_factor, args.dry_run)

        parent = run_dir.parent
        if parent not in linked_parents:
            print(ensure_npdata_link(parent / "NPDATA", npdata_target, args.dry_run))
            linked_parents.add(parent)

    print("dry-run complete" if args.dry_run else "build complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
