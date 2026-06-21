#!/usr/bin/env python3
"""Extract Iliadis et al. (2002) nova sensitivity tables.

This script uses the local PDF in SensitivityStudy/references and writes:
  * normalized .dat files for Tables 1-12,
  * PPN-style abundance files for Tables 2 and 4,
  * nova-case to Iliadis-model match reports.
"""

from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path


STUDY_ROOT = Path(__file__).resolve().parents[1]
PDF_PATH = STUDY_ROOT / "references" / "Iliadis-ApJS-2002.pdf"
OUT_DIR = STUDY_ROOT / "data" / "iliadis2002"
TABLE_DIR = OUT_DIR / "tables"
PPN_INITIAL_DIR = OUT_DIR / "ppn_initial_abundances"
PPN_FINAL_DIR = OUT_DIR / "ppn_final_abundances"
MATCH_DIR = OUT_DIR / "nova_case_matches"

MODELS = ["P1", "P2", "S1", "JCH1", "JCH2", "JH1", "JH2"]
MODEL_LABELS = {
    "P1": "P1",
    "P2": "P2",
    "S1": "S1",
    "JCH1": "JCH 1",
    "JCH2": "JCH 2",
    "JH1": "JH 1",
    "JH2": "JH 2",
}

SENSITIVITY_TABLE_MODELS = {
    5: ("P1", "ONe", 0.290),
    6: ("P2", "ONe", 0.356),
    7: ("S1", "ONe", 0.418),
    8: ("JCH1", "ONe", 0.231),
    9: ("JCH2", "ONe", 0.251),
    10: ("JH1", "CO", 0.145),
    11: ("JH2", "CO", 0.170),
}

ELEMENT_Z = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
}


def run_pdftotext() -> str:
    result = subprocess.run(
        ["pdftotext", "-layout", str(PDF_PATH), "-"],
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout


def normalize_text(text: str) -> str:
    replacements = {
        "\x03": "g",
        "\x04": "a",
        "\x05": "-",
        "þ": "+",
        "¼": "=",
        "": "-",
        "": "g",
        "": "a",
        "": ">=",
        "—": "-",
        "–": "-",
        "−": "-",
        "": "\n",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def table_block(text: str, table_number: int) -> str:
    start = re.search(rf"\bTABLE {table_number}\b", text)
    if not start:
        raise ValueError(f"Could not find TABLE {table_number}")
    if table_number < 12:
        end = re.search(rf"\bTABLE {table_number + 1}\b", text[start.end() :])
        if not end:
            raise ValueError(f"Could not find end marker for TABLE {table_number}")
        return text[start.start() : start.end() + end.start()]
    return text[start.start() :]


def is_numeric_token(token: str) -> bool:
    if token == "...":
        return True
    if re.fullmatch(r"10-\d+", token):
        return True
    try:
        float(token)
    except ValueError:
        return False
    return True


def parse_number(token: str) -> float | None:
    token = token.strip()
    if token == "..." or token == "":
        return None
    if re.fullmatch(r"10-\d+", token):
        return 10 ** (-int(token[3:]))
    return float(token)


def format_float(value: float | None) -> str:
    if value is None:
        return "..."
    return f"{value:.10E}"


def clean_reaction(reaction: str) -> str:
    reaction = normalize_text(reaction)
    reaction = reaction.replace(".", "")
    reaction = re.sub(r"\(([^,]+),\s*([^)]+)\)", r"(\1,\2)", reaction)
    reaction = re.sub(r"\s+", "", reaction)
    return reaction


def reaction_id(reaction: str) -> str:
    match = re.match(r"^(.+)\(([^,]+),([^)]+)\)(.+)$", reaction)
    if not match:
        return reaction
    target, projectile, ejectile, product = match.groups()
    return f"{target}_{projectile}{ejectile}_{product}"


def is_isotope(token: str) -> bool:
    return re.fullmatch(r"\d+[A-Z][a-z]?", token) is not None


def isotope_parts(isotope: str) -> tuple[int, str, int]:
    if isotope == "1H":
        return 1, "H", 1
    match = re.fullmatch(r"(\d+)([A-Z][a-z]?)([gm]?)", isotope)
    if not match:
        raise ValueError(f"Cannot parse isotope {isotope!r}")
    mass = int(match.group(1))
    elem = match.group(2)
    if elem not in ELEMENT_Z:
        raise ValueError(f"No Z mapping for isotope {isotope!r}")
    return ELEMENT_Z[elem], elem, mass


def isotope_label(isotope: str) -> str:
    z, elem, mass = isotope_parts(isotope)
    return f"{elem}-{mass}"


def ppn_species(z: int, elem: str, mass: int) -> str:
    if elem == "H" and mass == 1:
        return "PROT"
    return f"{elem.upper():<2s}{mass:3d}"


def write_dat(path: Path, headers: list[str], rows: list[dict[str, object]], comments: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for comment in comments or []:
            handle.write(f"# {comment}\n")
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(header, "")) for header in headers) + "\n")


def parse_table1(text: str) -> dict[str, dict[str, object]]:
    block = table_block(text, 1)
    rows = {model: {"model": model, "model_label": MODEL_LABELS[model]} for model in MODELS}
    mapping = [
        ("WD mass", "wd_mass_msun", float),
        ("WD composition", "wd_composition", str),
        ("Mixing", "mixing_percent", float),
        ("Tpeak", "tpeak_1e6_k", float),
        ("Lpeak", "lpeak_1e4_lsun", float),
        ("Macc", "macc_1e_minus5_msun", float),
        ("M_ acc", "mdot_acc_1e_minus10_msun_per_yr", float),
        ("Mej", "mej_1e_minus5_msun", float),
    ]
    for marker, key, caster in mapping:
        line = next(line for line in block.splitlines() if marker in line)
        values = line.split()[-7:]
        for model, value in zip(MODELS, values):
            rows[model][key] = caster(value)
    for model in MODELS:
        rows[model]["tpeak_gk"] = float(rows[model]["tpeak_1e6_k"]) / 1000.0
    return rows


def parse_abundance_table(text: str, table_number: int) -> dict[str, dict[str, float]]:
    block = table_block(text, table_number)
    result = {model: {} for model in MODELS}
    value_count = 4 if table_number == 2 else 7
    row_re = re.compile(r"^\s*(\d+[A-Z][a-z]?)\s*\.*\s*(.*)$")
    pending_isotope = ""

    for raw_line in block.splitlines():
        # Table 2 sits in the left column of a two-column page. Cropping avoids
        # interpreting right-column prose as abundance values.
        line = raw_line[:76] if table_number == 2 else raw_line
        match = row_re.match(line)
        isotope = ""
        value_text = ""
        if match:
            isotope = match.group(1)
            value_text = match.group(2)
        elif pending_isotope:
            isotope = pending_isotope
            value_text = line
        else:
            continue

        tokens = value_text.split()
        numeric = [token for token in tokens if is_numeric_token(token) and token != "..."]
        if len(numeric) < value_count:
            if match:
                pending_isotope = isotope
            continue
        pending_isotope = ""
        values = [parse_number(token) for token in numeric[-value_count:]]
        if any(value is None for value in values):
            continue

        if table_number == 2:
            expanded = [
                values[0],
                values[0],
                values[0],
                values[1],
                values[1],
                values[2],
                values[3],
            ]
        else:
            expanded = values
        for model, value in zip(MODELS, expanded):
            result[model][isotope] = float(value)
    return result


def parse_table3(text: str) -> list[dict[str, object]]:
    block = table_block(text, 3)
    rows: list[dict[str, object]] = []
    pattern = re.compile(
        r"(?P<reaction>\d+[A-Z][A-Za-z]*[gm]?\([^)]*\)\d+[A-Z][A-Za-z]*[gm]?)"
        r"(?:\s+(?P<note>[a-e]))?\s*\.*\s+"
        r"(?P<factor>\S+/\S+)"
    )
    for raw_line in block.splitlines():
        line = re.sub(r"\(([^,]+),\s*([^)]+)\)", r"(\1,\2)", raw_line)
        for match in pattern.finditer(line):
            reaction = clean_reaction(match.group("reaction"))
            factor = match.group("factor")
            up, down = factor.split("/", 1)
            rows.append(
                {
                    "reaction_id": reaction_id(reaction),
                    "reaction": reaction,
                    "factor_up": format_float(parse_number(up)),
                    "factor_down": format_float(parse_number(down)),
                    "source_note": match.group("note") or "",
                    "factor_up_source": up,
                    "factor_down_source": down,
                }
            )
    return rows


def parse_sensitivity_table(text: str, table_number: int) -> list[dict[str, object]]:
    block = table_block(text, table_number)
    model, composition, tpeak = SENSITIVITY_TABLE_MODELS[table_number]
    rows: list[dict[str, object]] = []
    current_reaction = ""

    for line in block.splitlines():
        tokens = line.split()
        if len(tokens) < 7:
            continue
        values = tokens[-6:]
        if not all(is_numeric_token(token) for token in values):
            continue
        before = tokens[:-6]
        if not before or not is_isotope(before[-1]):
            continue
        isotope = before[-1]
        reaction_text = " ".join(before[:-1])
        if reaction_text:
            reaction = clean_reaction(reaction_text)
            if "(" not in reaction or ")" not in reaction:
                continue
            current_reaction = reaction
        if not current_reaction:
            continue

        rows.append(
            {
                "model": model,
                "wd_composition": composition,
                "tpeak_gk": f"{tpeak:.3f}",
                "reaction_id": reaction_id(current_reaction),
                "reaction": current_reaction,
                "isotope": isotope,
                "isotope_label": isotope_label(isotope),
                "factor_100": values[0],
                "factor_10": values[1],
                "factor_2": values[2],
                "factor_0.5": values[3],
                "factor_0.1": values[4],
                "factor_0.01": values[5],
            }
        )
    return rows


def parse_table12(text: str) -> list[dict[str, object]]:
    block = table_block(text, 12)
    rows: list[dict[str, object]] = []
    section = ""
    row_re = re.compile(r"^\s*(?P<reaction>\d+[A-Z][A-Za-z]*[gm]?\([^)]*\)\d+[A-Z][A-Za-z]*[gm]?)\.*\s+(?P<isotopes>.+?)\s*$")

    for raw_line in block.splitlines():
        left = raw_line[:85]
        if "CO Nova Models" in left:
            section = "CO"
            continue
        if "ONe Nova Models" in left:
            section = "ONe"
            continue
        if not section:
            continue
        line = re.sub(r"\(([^,]+),\s*([^)]+)\)", r"(\1,\2)", left)
        match = row_re.match(line)
        if not match:
            continue
        reaction = clean_reaction(match.group("reaction"))
        isotopes = [item.strip() for item in match.group("isotopes").split(",") if is_isotope(item.strip())]
        if not isotopes:
            continue
        rows.append(
            {
                "nova_family": section,
                "reaction_id": reaction_id(reaction),
                "reaction": reaction,
                "changed_isotopes": ",".join(isotopes),
                "changed_isotope_labels": ",".join(isotope_label(isotope) for isotope in isotopes),
            }
        )
    return rows


def abundance_rows(data: dict[str, dict[str, float]]) -> list[dict[str, object]]:
    isotopes: list[str] = []
    for model in MODELS:
        for isotope in data[model]:
            if isotope not in isotopes:
                isotopes.append(isotope)

    rows = []
    for isotope in isotopes:
        z, elem, mass = isotope_parts(isotope)
        row: dict[str, object] = {
            "isotope": isotope,
            "isotope_label": isotope_label(isotope),
            "z": z,
            "element": elem,
            "a": mass,
        }
        for model in MODELS:
            row[model] = format_float(data[model].get(isotope))
        rows.append(row)
    return rows


def write_ppn_abundance_files(data: dict[str, dict[str, float]], out_dir: Path, stage: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for model in MODELS:
        rows = list(data[model].items())

        iniab = out_dir / f"iniab_iliadis2002_{stage}_{model}.dat"
        with iniab.open("w") as handle:
            for isotope, value in rows:
                z, elem, mass = isotope_parts(isotope)
                if elem == "H" and mass == 1:
                    handle.write(f"{z:3d} PROT          {value:.10E}\n")
                else:
                    handle.write(f"{z:3d} {elem.lower():<2s} {mass:3d}         {value:.10E}\n")

        iso_massf = out_dir / f"iso_massf_iliadis2002_{stage}_{model}.DAT"
        with iso_massf.open("w") as handle:
            handle.write("H NUM   Z    A    ISOM ABUNDNACE_MF  ISOTP\n")
            handle.write(f" # source Iliadis et al. 2002 Table {'2' if stage == 'initial' else '4'}\n")
            handle.write(f" # model {MODEL_LABELS[model]}\n")
            for idx, (isotope, value) in enumerate(rows, start=1):
                z, elem, mass = isotope_parts(isotope)
                handle.write(
                    f"{idx:5d} {z:5.0f}. {mass:4.0f}. {1:3d}   {value:11.5E}  {ppn_species(z, elem, mass)}\n"
                )


def parse_case_name(case_name: str) -> tuple[str, float | None]:
    lowered = case_name.lower()
    if lowered.startswith("ne_"):
        composition = "ONe"
    elif lowered.startswith("co_"):
        composition = "CO"
    else:
        composition = ""
    mass_match = re.search(r"_nova_(\d+(?:\.\d+)?)_", lowered)
    mass = float(mass_match.group(1)) if mass_match else None
    return composition, mass


def parse_trajectory(path: Path) -> dict[str, object]:
    tunit = "T9K"
    points = 0
    peak_t = -math.inf
    peak_time = None
    peak_rho = None
    for raw_line in path.read_text(errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("TUNIT"):
            tunit = line.split("=", 1)[1].strip()
            continue
        if line.startswith("#") or "=" in line:
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            time = float(parts[0])
            temp = float(parts[1])
            rho = float(parts[2])
        except ValueError:
            continue
        points += 1
        if tunit == "T8K":
            temp_gk = temp / 10.0
        elif tunit == "K":
            temp_gk = temp / 1.0e9
        else:
            temp_gk = temp
        if temp_gk > peak_t:
            peak_t = temp_gk
            peak_time = time
            peak_rho = rho
    return {
        "points": points,
        "tunit": tunit,
        "tpeak_gk": peak_t if peak_t > -math.inf else None,
        "peak_time": peak_time,
        "peak_rho_cgs": peak_rho,
    }


def choose_file(case_dir: Path, patterns: tuple[str, ...]) -> Path | None:
    candidates: list[Path] = []
    skip_parts = {"analysis", "backup", "backup_init_cond", "runs", "logs", ".ipynb_checkpoints"}
    for pattern in patterns:
        for candidate in case_dir.rglob(pattern):
            if any(part in skip_parts for part in candidate.parts):
                continue
            candidates.append(candidate)
    if not candidates:
        return None
    return sorted(candidates, key=lambda path: (len(path.parts), str(path)))[0]


def parse_ppn_abundance_file(path: Path | None) -> dict[str, float]:
    if path is None:
        return {}
    abundances: dict[str, float] = {}
    for raw_line in path.read_text(errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        try:
            if len(parts) == 3 and parts[1].upper() == "PROT":
                abundances["1H"] = float(parts[2])
            elif len(parts) >= 4:
                elem = parts[1].capitalize()
                mass = int(parts[2])
                abundances[f"{mass}{elem}"] = float(parts[3])
        except ValueError:
            continue
    return abundances


def abundance_distance(case_abund: dict[str, float], model_abund: dict[str, float]) -> tuple[float | None, float | None]:
    isotopes = [
        isotope
        for isotope in model_abund
        if isotope in case_abund and max(abs(model_abund[isotope]), abs(case_abund[isotope])) > 1.0e-12
    ]
    if not isotopes:
        return None, None
    floor = 1.0e-30
    log_diffs = [
        math.log10((case_abund[isotope] + floor) / (model_abund[isotope] + floor))
        for isotope in isotopes
    ]
    rms_log = math.sqrt(sum(diff * diff for diff in log_diffs) / len(log_diffs))
    l1 = sum(abs(case_abund[isotope] - model_abund[isotope]) for isotope in isotopes)
    return rms_log, l1


def score_match(
    case_composition: str,
    case_mass: float | None,
    case_tpeak: float | None,
    model_props: dict[str, object],
    rms_log: float | None,
) -> float:
    score = score_property_match(case_composition, case_mass, case_tpeak, model_props)
    if rms_log is None:
        score += 1.5
    else:
        score += min(rms_log, 3.0)
    return score


def score_property_match(
    case_composition: str,
    case_mass: float | None,
    case_tpeak: float | None,
    model_props: dict[str, object],
) -> float:
    score = 0.0
    if case_tpeak is None:
        score += 2.0
    else:
        score += abs(case_tpeak - float(model_props["tpeak_gk"])) / 0.05
    if case_mass is None:
        score += 4.0
    else:
        score += abs(case_mass - float(model_props["wd_mass_msun"])) / 0.10
    if case_composition:
        score += 0.0 if case_composition == model_props["wd_composition"] else 3.0
    else:
        score += 3.0
    return score


def compare_nova_cases(table1: dict[str, dict[str, object]], table2: dict[str, dict[str, float]]) -> None:
    nova_root = STUDY_ROOT / "nova_cases"
    ranking_rows: list[dict[str, object]] = []
    case_best_rows: list[dict[str, object]] = []

    for case_dir in sorted(path for path in nova_root.iterdir() if path.is_dir()):
        trajectory = choose_file(case_dir, ("trajectory*.input",))
        if trajectory is None:
            continue
        abundance_file = choose_file(case_dir, ("iniab*.dat", "initial_abundance.dat"))
        trajectory_info = parse_trajectory(trajectory)
        case_abund = parse_ppn_abundance_file(abundance_file)
        case_composition, case_mass = parse_case_name(case_dir.name)
        case_rows = []

        for model in MODELS:
            model_props = table1[model]
            rms_log, l1 = abundance_distance(case_abund, table2[model])
            tpeak = trajectory_info["tpeak_gk"]
            mass_diff = None if case_mass is None else abs(case_mass - float(model_props["wd_mass_msun"]))
            tpeak_diff = None if tpeak is None else abs(float(tpeak) - float(model_props["tpeak_gk"]))
            property_score = score_property_match(case_composition, case_mass, tpeak, model_props)
            score = score_match(case_composition, case_mass, tpeak, model_props, rms_log)
            row = {
                "case": case_dir.name,
                "trajectory": trajectory.relative_to(STUDY_ROOT),
                "initial_abundance": abundance_file.relative_to(STUDY_ROOT) if abundance_file else "",
                "case_wd_composition": case_composition,
                "case_wd_mass_msun": "" if case_mass is None else f"{case_mass:.3f}",
                "case_tpeak_gk": "" if tpeak is None else f"{float(tpeak):.4f}",
                "case_peak_rho_cgs": "" if trajectory_info["peak_rho_cgs"] is None else f"{float(trajectory_info['peak_rho_cgs']):.6E}",
                "trajectory_points": trajectory_info["points"],
                "iliadis_model": model,
                "iliadis_model_label": MODEL_LABELS[model],
                "iliadis_wd_composition": model_props["wd_composition"],
                "iliadis_wd_mass_msun": f"{float(model_props['wd_mass_msun']):.3f}",
                "iliadis_tpeak_gk": f"{float(model_props['tpeak_gk']):.3f}",
                "abs_tpeak_diff_gk": "" if tpeak_diff is None else f"{tpeak_diff:.4f}",
                "abs_mass_diff_msun": "" if mass_diff is None else f"{mass_diff:.3f}",
                "composition_match": "" if not case_composition else str(case_composition == model_props["wd_composition"]),
                "initial_abundance_rms_log10": "" if rms_log is None else f"{rms_log:.6f}",
                "initial_abundance_l1": "" if l1 is None else f"{l1:.6E}",
                "property_match_score": f"{property_score:.6f}",
                "match_score": f"{score:.6f}",
            }
            ranking_rows.append(row)
            case_rows.append(row)

        best = min(case_rows, key=lambda item: float(item["match_score"]))
        case_best_rows.append(best.copy())

    ranking_rows.sort(key=lambda item: (item["case"], float(item["match_score"])))
    case_best_rows.sort(key=lambda item: item["case"])

    headers = [
        "case",
        "trajectory",
        "initial_abundance",
        "case_wd_composition",
        "case_wd_mass_msun",
        "case_tpeak_gk",
        "case_peak_rho_cgs",
        "trajectory_points",
        "iliadis_model",
        "iliadis_model_label",
        "iliadis_wd_composition",
        "iliadis_wd_mass_msun",
        "iliadis_tpeak_gk",
        "abs_tpeak_diff_gk",
        "abs_mass_diff_msun",
        "composition_match",
        "initial_abundance_rms_log10",
        "initial_abundance_l1",
        "property_match_score",
        "match_score",
    ]
    comments = [
        "Match score is a heuristic combining Tpeak, WD mass, WD composition, and initial abundance distance.",
        "Property match score uses only Iliadis Table 1 properties: Tpeak, WD mass, and WD composition.",
        "Lower score is better. Composition/mass are inferred from case names when available.",
    ]
    write_dat(MATCH_DIR / "nova_case_model_match_rankings.dat", headers, ranking_rows, comments)
    write_dat(MATCH_DIR / "nova_case_best_matches.dat", headers, case_best_rows, comments)
    property_case_best_rows = []
    for case_name in sorted({row["case"] for row in ranking_rows}):
        candidates = [row for row in ranking_rows if row["case"] == case_name]
        property_case_best_rows.append(min(candidates, key=lambda item: float(item["property_match_score"])).copy())
    write_dat(
        MATCH_DIR / "nova_case_best_matches_by_table1_properties.dat",
        headers,
        property_case_best_rows,
        comments,
    )

    model_best_rows = []
    property_model_best_rows = []
    for model in MODELS:
        candidates = [row for row in ranking_rows if row["iliadis_model"] == model]
        if candidates:
            model_best_rows.append(min(candidates, key=lambda item: float(item["match_score"])).copy())
            property_model_best_rows.append(
                min(candidates, key=lambda item: float(item["property_match_score"])).copy()
            )
    model_best_rows.sort(key=lambda item: MODELS.index(str(item["iliadis_model"])))
    property_model_best_rows.sort(key=lambda item: MODELS.index(str(item["iliadis_model"])))
    write_dat(MATCH_DIR / "iliadis_model_best_nova_cases.dat", headers, model_best_rows, comments)
    write_dat(
        MATCH_DIR / "iliadis_model_best_nova_cases_by_table1_properties.dat",
        headers,
        property_model_best_rows,
        comments,
    )


def write_readme(table_counts: dict[str, int]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Iliadis 2002 Extracted Tables",
        "",
        "Source PDF: `SensitivityStudy/references/Iliadis-ApJS-2002.pdf`.",
        "",
        "Generated by `SensitivityStudy/tools/extract_iliadis2002_tables.py` using `pdftotext -layout`.",
        "",
        "Contents:",
        "",
        "- `tables/table_01_model_properties.dat`: normalized Table 1 model metadata.",
        "- `tables/table_02_initial_envelope_composition.dat`: Table 2 initial mass fractions expanded to all seven models.",
        "- `tables/table_03_reaction_rate_uncertainties.dat`: Table 3 selected reaction-rate uncertainties.",
        "- `tables/table_04_final_abundances.dat`: Table 4 final one-zone mass fractions.",
        "- `tables/table_05_sensitivity_P1.dat` through `tables/table_11_sensitivity_JH2.dat`: Tables 5-11 abundance-change ratios.",
        "- `tables/table_12_significant_reactions_summary.dat`: Table 12 qualitative summary.",
        "- `ppn_initial_abundances/`: PPN-style `iniab` and `iso_massf` files from Table 2.",
        "- `ppn_final_abundances/`: PPN-style `iniab` and `iso_massf` files from Table 4.",
        "- `nova_case_matches/`: local nova-case comparisons against Iliadis Table 1 and Table 2.",
        "",
        "Table row counts:",
        "",
    ]
    for name, count in sorted(table_counts.items()):
        lines.append(f"- `{name}`: {count} data rows")
    lines.append("")
    lines.append("Notes:")
    lines.append("")
    lines.append("- Tables 5-11 preserve `...` where the paper reports no value.")
    lines.append("- Tables 2 and 4 are rounded as published; the PPN-style files should be treated as reference comparison inputs, not original NuGrid model output.")
    lines.append("- Nova-case match reports include both a Table-1-only property score and a combined score that also considers Table 2 initial abundances.")
    (OUT_DIR / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    if not PDF_PATH.exists():
        raise SystemExit(f"Missing PDF: {PDF_PATH}")

    text = normalize_text(run_pdftotext())
    table_counts: dict[str, int] = {}

    table1 = parse_table1(text)
    table1_rows = [table1[model] for model in MODELS]
    write_dat(
        TABLE_DIR / "table_01_model_properties.dat",
        [
            "model",
            "model_label",
            "wd_mass_msun",
            "wd_composition",
            "mixing_percent",
            "tpeak_1e6_k",
            "tpeak_gk",
            "lpeak_1e4_lsun",
            "macc_1e_minus5_msun",
            "mdot_acc_1e_minus10_msun_per_yr",
            "mej_1e_minus5_msun",
        ],
        table1_rows,
        ["Iliadis et al. 2002 Table 1."],
    )
    table_counts["table_01_model_properties.dat"] = len(table1_rows)

    table2 = parse_abundance_table(text, 2)
    table2_rows = abundance_rows(table2)
    write_dat(
        TABLE_DIR / "table_02_initial_envelope_composition.dat",
        ["isotope", "isotope_label", "z", "element", "a", *MODELS],
        table2_rows,
        ["Iliadis et al. 2002 Table 2. Grouped source columns are expanded to all seven Table 1 models."],
    )
    table_counts["table_02_initial_envelope_composition.dat"] = len(table2_rows)
    write_ppn_abundance_files(table2, PPN_INITIAL_DIR, "initial")

    table3_rows = parse_table3(text)
    write_dat(
        TABLE_DIR / "table_03_reaction_rate_uncertainties.dat",
        [
            "reaction_id",
            "reaction",
            "factor_up",
            "factor_down",
            "source_note",
            "factor_up_source",
            "factor_down_source",
        ],
        table3_rows,
        ["Iliadis et al. 2002 Table 3."],
    )
    table_counts["table_03_reaction_rate_uncertainties.dat"] = len(table3_rows)

    table4 = parse_abundance_table(text, 4)
    table4_rows = abundance_rows(table4)
    write_dat(
        TABLE_DIR / "table_04_final_abundances.dat",
        ["isotope", "isotope_label", "z", "element", "a", *MODELS],
        table4_rows,
        ["Iliadis et al. 2002 Table 4."],
    )
    table_counts["table_04_final_abundances.dat"] = len(table4_rows)
    write_ppn_abundance_files(table4, PPN_FINAL_DIR, "final")

    all_sensitivity_rows: list[dict[str, object]] = []
    sensitivity_headers = [
        "model",
        "wd_composition",
        "tpeak_gk",
        "reaction_id",
        "reaction",
        "isotope",
        "isotope_label",
        "factor_100",
        "factor_10",
        "factor_2",
        "factor_0.5",
        "factor_0.1",
        "factor_0.01",
    ]
    for table_number, (model, _, _) in SENSITIVITY_TABLE_MODELS.items():
        rows = parse_sensitivity_table(text, table_number)
        all_sensitivity_rows.extend(rows)
        filename = f"table_{table_number:02d}_sensitivity_{model}.dat"
        write_dat(TABLE_DIR / filename, sensitivity_headers, rows, [f"Iliadis et al. 2002 Table {table_number}."])
        table_counts[filename] = len(rows)

    write_dat(
        TABLE_DIR / "tables_05_to_11_sensitivity_all_models.dat",
        sensitivity_headers,
        all_sensitivity_rows,
        ["Combined Iliadis et al. 2002 Tables 5-11."],
    )
    table_counts["tables_05_to_11_sensitivity_all_models.dat"] = len(all_sensitivity_rows)

    table12_rows = parse_table12(text)
    write_dat(
        TABLE_DIR / "table_12_significant_reactions_summary.dat",
        ["nova_family", "reaction_id", "reaction", "changed_isotopes", "changed_isotope_labels"],
        table12_rows,
        ["Iliadis et al. 2002 Table 12."],
    )
    table_counts["table_12_significant_reactions_summary.dat"] = len(table12_rows)

    compare_nova_cases(table1, table2)
    write_readme(table_counts)


if __name__ == "__main__":
    main()
