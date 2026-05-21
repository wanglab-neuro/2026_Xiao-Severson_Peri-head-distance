"""Generate CSV helper tables for manuscript figure reproduction.

Outputs:
    figure_unit_manifest.csv
        Explicit mapping from manuscript figure units to NWB files/rows.
    unit_tuning_summary.csv
        Dataset-level one-row-per-unit summary regenerated from NWB files.

Example:
    python generate_deposition_tables.py E:/XiaoSeversonEtAl2026/NWB/local_nwb_files
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO


FIGURE2_UNITS = [
    {
        "figure": "Figure 2",
        "panel": "A-C",
        "display_unit": "Unit 1",
        "matlab_global_row_1based": 142,
        "matlab_row_label": "20230208_sp_session_1_unit_5",
        "source_table": "PSTH_naive.mat:psth.dataTable; unitTable_naive.mat:unitTable",
    },
    {
        "figure": "Figure 2",
        "panel": "A-C",
        "display_unit": "Unit 2",
        "matlab_global_row_1based": 993,
        "matlab_row_label": "20230626_sp_session_1_unit_8",
        "source_table": "PSTH_naive.mat:psth.dataTable; unitTable_naive.mat:unitTable",
    },
    {
        "figure": "Figure 2",
        "panel": "A-C",
        "display_unit": "Unit 3",
        "matlab_global_row_1based": 1223,
        "matlab_row_label": "20240905_sp_session_1_unit_4",
        "source_table": "PSTH_naive.mat:psth.dataTable; unitTable_naive.mat:unitTable",
    },
    {
        "figure": "Figure 2",
        "panel": "A-C",
        "display_unit": "Unit 4",
        "matlab_global_row_1based": 1128,
        "matlab_row_label": "20240703_sp_session_1_unit_10",
        "source_table": "PSTH_naive.mat:psth.dataTable; unitTable_naive.mat:unitTable",
    },
    {
        "figure": "Figure 2",
        "panel": "A-C",
        "display_unit": "Unit 5",
        "matlab_global_row_1based": 18,
        "matlab_row_label": "20221212_sp_session_2_unit_18",
        "source_table": "PSTH_naive.mat:psth.dataTable; unitTable_naive.mat:unitTable",
    },
    {
        "figure": "Figure 2",
        "panel": "A-C",
        "display_unit": "Unit 6",
        "matlab_global_row_1based": 1276,
        "matlab_row_label": "20240924_sp_session_1_unit_7",
        "source_table": "PSTH_naive.mat:psth.dataTable; unitTable_naive.mat:unitTable",
    },
]

UNIT_COLUMNS = [
    "unit_label",
    "wall_responsive",
    "wall_tuned",
    "activated",
    "suppressed",
    "tuning_type",
    "pref_dist_index",
    "pref_dist_mm",
    "low_firing",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nwb_dir", type=Path, help="Directory containing one *.nwb file per session.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for generated CSV files. Default: this script's folder.",
    )
    parser.add_argument("--float-digits", type=int, default=6, help="Decimal places for floating-point CSV values.")
    return parser.parse_args()


def nwb_files(nwb_dir: Path) -> list[Path]:
    files = sorted(nwb_dir.glob("*.nwb"))
    if not files:
        raise FileNotFoundError(f"No .nwb files found in {nwb_dir}")
    return files


def as_list(value) -> list:
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def scalar_or_empty(value, digits: int):
    if value is None:
        return ""
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return "" if np.isnan(float(value)) else round(float(value), digits)
    return str(value)


def format_float_list(values: Iterable[float], digits: int) -> str:
    return ";".join("" if np.isnan(float(value)) else str(round(float(value), digits)) for value in values)


def session_distance_axis(nwb, digits: int) -> str:
    if "wall_tuning" not in nwb.processing:
        return ""
    module = nwb.processing["wall_tuning"]
    if "distance_axis" not in module.data_interfaces:
        return ""
    distances = module["distance_axis"].to_dataframe()["distance_mm"].to_numpy(float)
    return format_float_list(distances, digits)


def trial_include_counts(unit: pd.Series) -> tuple[int, int]:
    flags = [bool(value) for value in as_list(unit.get("trial_include", []))]
    return len(flags), int(sum(flags))


def dandi_session_id(session_id: str) -> str:
    return session_id.replace("_", "-")


def dandi_asset_path(subject_id: str, session_id: str) -> str:
    if not subject_id or not session_id:
        return ""
    return f"sub-{subject_id}/sub-{subject_id}_ses-{dandi_session_id(session_id)}.nwb"


def unit_summary_rows(nwb_path: Path, digits: int) -> list[dict]:
    rows = []
    with NWBHDF5IO(str(nwb_path), "r", load_namespaces=True) as io:
        nwb = io.read()
        units = nwb.units.to_dataframe()
        subject = nwb.subject
        subject_id = subject.subject_id if subject else ""
        session_id = nwb.session_id or ""
        distance_axis = session_distance_axis(nwb, digits)

        for unit_index, unit in units.iterrows():
            trial_count, trial_include_count = trial_include_counts(unit)
            row = {
                "nwb_file": nwb_path.name,
                "dandi_asset_path": dandi_asset_path(subject_id, session_id),
                "session_id": session_id,
                "session_start_time": str(nwb.session_start_time) if nwb.session_start_time else "",
                "subject_id": subject_id,
                "subject_sex": subject.sex if subject else "",
                "subject_age": subject.age if subject else "",
                "subject_species": subject.species if subject else "",
                "nwb_unit_index_0based": int(unit_index),
                "spike_count": len(as_list(unit["spike_times"])) if "spike_times" in units.columns else 0,
                "trial_count": trial_count,
                "trial_include_count": trial_include_count,
                "tuning_distance_axis_mm": distance_axis,
            }
            for column in UNIT_COLUMNS:
                row[column] = scalar_or_empty(unit[column], digits) if column in units.columns else ""
            rows.append(row)
    return rows


def index_units(files: list[Path]) -> dict[str, dict]:
    index = {}
    for nwb_path in files:
        with NWBHDF5IO(str(nwb_path), "r", load_namespaces=True) as io:
            nwb = io.read()
            units = nwb.units.to_dataframe()
            subject = nwb.subject
            subject_id = subject.subject_id if subject else ""
            session_id = nwb.session_id or ""
            for unit_index, unit in units.iterrows():
                label = str(unit.get("unit_label", "")).strip()
                if label:
                    index[label] = {
                        "session_id": session_id,
                        "subject_id": subject_id,
                        "nwb_file": nwb_path.name,
                        "dandi_asset_path": dandi_asset_path(subject_id, session_id),
                        "nwb_unit_index_0based": int(unit_index),
                        "unit_label": label,
                        "tuning_type": str(unit.get("tuning_type", "")),
                    }
    return index


def figure_manifest(files: list[Path]) -> pd.DataFrame:
    unit_index = index_units(files)
    rows = []
    missing = []
    for entry in FIGURE2_UNITS:
        label = entry["matlab_row_label"]
        if label not in unit_index:
            missing.append(label)
            continue
        rows.append({**entry, **unit_index[label]})
    if missing:
        raise ValueError(f"Figure unit labels not found in NWB files: {', '.join(missing)}")
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    files = nwb_files(args.nwb_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    unit_rows = []
    for nwb_path in files:
        unit_rows.extend(unit_summary_rows(nwb_path, args.float_digits))

    unit_table = pd.DataFrame(unit_rows)
    unit_path = args.output_dir / "unit_tuning_summary.csv"
    unit_table.to_csv(unit_path, index=False, float_format=f"%.{args.float_digits}f")

    manifest = figure_manifest(files)
    manifest_path = args.output_dir / "figure_unit_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    print(f"Wrote {len(unit_table)} unit rows to {unit_path}")
    print(f"Wrote {len(manifest)} figure-unit rows to {manifest_path}")


if __name__ == "__main__":
    main()
