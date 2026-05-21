"""Plot Figure 2-style rasters, PETHs, and tuning curves from NWB files.

The script reads processed PrV peri-head distance NWB files and produces a
multi-row figure with one row per selected unit:

    A. spike rasters aligned to wall-pass onset, grouped by wall distance
    B. PETHs computed from spike times with one trace per wall distance
    C. stored tuning curves from processing/wall_tuning

Example:
    python plot_figure2_from_nwb.py E:/XiaoSeversonEtAl2026/NWB/local_nwb_files \
        --output E:/XiaoSeversonEtAl2026/NWB/figures/figure2_style.png

By default, when given an NWB directory, the script plots the six manuscript
Figure 2 units by their 1-based global row IDs in PSTH_naive.mat dataTable:
142, 993, 1223, 1128, 18, and 1276. Use --global-rows to plot other global
MATLAB table row IDs. Use --units with a single NWB file for 0-based NWB unit
table row indices.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from pynwb import NWBHDF5IO


FIGURE2_GLOBAL_ROW_IDS = [142, 993, 1223, 1128, 18, 1276]
FIGURE2_ROW_LABELS = {
    142: "20230208_sp_session_1_unit_5",
    993: "20230626_sp_session_1_unit_8",
    1223: "20240905_sp_session_1_unit_4",
    1128: "20240703_sp_session_1_unit_10",
    18: "20221212_sp_session_2_unit_18",
    1276: "20240924_sp_session_1_unit_7",
}
PETH_TICK_STEPS = [10, 10, 20, 10, 2, 10]


@dataclass
class PlotData:
    session_id: str
    nwb_path: Path
    trials: pd.DataFrame
    units: pd.DataFrame
    tuning: pd.DataFrame | None
    tuning_distances: np.ndarray | None


@dataclass
class UnitSelection:
    data: PlotData
    unit_idx: int
    global_row_id: int | None = None
    row_label: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make Figure 2-style raster/PETH/tuning plots from processed NWB files."
    )
    parser.add_argument(
        "nwb_input",
        type=Path,
        help="Input NWB file, or a directory containing one NWB file per session.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output image path. Defaults to <nwb stem>_figure2_style.png",
    )
    parser.add_argument(
        "--units",
        type=int,
        nargs="+",
        default=None,
        help=(
            "0-based NWB unit table row indices to plot from a single NWB file. "
            "Overrides --global-rows."
        ),
    )
    parser.add_argument(
        "--global-rows",
        type=int,
        nargs="+",
        default=None,
        help=(
            "1-based global row IDs in PSTH_naive.mat dataTable/unitTable. "
            "Default for NWB directory input: manuscript Figure 2 rows."
        ),
    )
    parser.add_argument(
        "--unit-labels",
        type=str,
        nargs="+",
        default=None,
        help=(
            "Full MATLAB/NWB unit labels to plot, for example "
            "20230208_sp_session_1_unit_5."
        ),
    )
    parser.add_argument(
        "--max-units",
        type=int,
        default=6,
        help="Maximum number of auto-selected units to plot when --auto-select is used.",
    )
    parser.add_argument(
        "--auto-select",
        action="store_true",
        help="Auto-select representative units instead of using the manuscript Figure 2 IDs.",
    )
    parser.add_argument(
        "--window",
        type=float,
        nargs=2,
        default=(0.0, 12.0),
        metavar=("START", "STOP"),
        help="Time window in seconds relative to wall-pass onset.",
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.2,
        help="PETH bin width in seconds.",
    )
    parser.add_argument(
        "--smooth-bins",
        type=int,
        default=3,
        help="Boxcar smoothing width in bins for PETH traces. Use 1 for no smoothing.",
    )
    parser.add_argument(
        "--tuning-scale",
        type=float,
        default=2.0,
        help="Scale factor applied to stored z-scored tuning curves. Default matches manuscript Figure 2.",
    )
    parser.add_argument(
        "--sort-distances",
        choices=("near-to-far", "far-to-near"),
        default="far-to-near",
        help="Distance order within raster rows.",
    )
    parser.add_argument(
        "--include-all-trial-distances",
        action="store_true",
        help=(
            "Use every wall distance present in the trials table. By default, when "
            "available, the script uses processing/wall_tuning/distance_axis to match "
            "the figure tuning convention."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output figure DPI.",
    )
    return parser.parse_args()


def load_plot_data(nwb_path: Path) -> PlotData:
    with NWBHDF5IO(str(nwb_path), "r", load_namespaces=True) as io:
        nwb = io.read()
        trials = nwb.trials.to_dataframe()
        units = nwb.units.to_dataframe()

        tuning = None
        tuning_distances = None
        if "wall_tuning" in nwb.processing:
            module = nwb.processing["wall_tuning"]
            if "tuning_curves" in module.data_interfaces:
                tuning = module["tuning_curves"].to_dataframe()
            if "distance_axis" in module.data_interfaces:
                tuning_distances = module["distance_axis"].to_dataframe()["distance_mm"].to_numpy(float)

        return PlotData(
            session_id=nwb.session_id or nwb.identifier,
            nwb_path=nwb_path,
            trials=trials,
            units=units,
            tuning=tuning,
            tuning_distances=tuning_distances,
        )


def normalize_unit_label(value: object) -> str:
    return str(value).strip()


def unit_indices_from_labels(units: pd.DataFrame, requested_labels: list[str]) -> list[int]:
    if "unit_label" not in units:
        raise ValueError("Cannot select by source unit ID: NWB units table has no unit_label column")

    label_to_index: dict[str, int] = {}
    for row_idx, label in enumerate(units["unit_label"]):
        normalized = normalize_unit_label(label)
        if normalized not in label_to_index:
            label_to_index[normalized] = row_idx

    selected = []
    missing = []
    for label in requested_labels:
        normalized = normalize_unit_label(label)
        if normalized in label_to_index:
            selected.append(label_to_index[normalized])
        else:
            missing.append(label)

    if missing:
        available = ", ".join(str(label) for label in units["unit_label"].head(10))
        raise ValueError(
            f"Source unit ID(s) not found in unit_label: {', '.join(missing)}. "
            f"First available labels: {available}"
        )
    return selected


def choose_units(
    units: pd.DataFrame,
    requested_indices: list[int] | None,
    requested_labels: list[str] | None,
    max_units: int,
    auto_select: bool,
) -> list[int]:
    if requested_indices is not None:
        bad = [idx for idx in requested_indices if idx < 0 or idx >= len(units)]
        if bad:
            raise ValueError(f"Unit index out of range: {bad}; NWB has {len(units)} units")
        return requested_indices

    if requested_labels is not None:
        return unit_indices_from_labels(units, requested_labels)

    if "tuning_type" not in units:
        return list(range(min(max_units, len(units))))

    preferred_types = ["proximity", "map", "ambiguous", "suppressed", "untuned"]
    selected: list[int] = []
    for tuning_type in preferred_types:
        matches = units.index[units["tuning_type"].astype(str).str.lower() == tuning_type].tolist()
        if tuning_type == "map":
            matches = sorted(matches, key=lambda idx: -float(units.iloc[idx].get("wall_responsive", False)))
        selected.extend(idx for idx in matches if idx not in selected)
        if len(selected) >= max_units:
            break

    if len(selected) < max_units:
        selected.extend(idx for idx in range(len(units)) if idx not in selected)

    return selected[: min(max_units, len(units))]


def iter_nwb_files(nwb_dir: Path) -> list[Path]:
    files = sorted(nwb_dir.glob("*.nwb"))
    if not files:
        raise FileNotFoundError(f"No .nwb files found in {nwb_dir}")
    return files


def index_nwb_directory(nwb_dir: Path) -> tuple[list[tuple[int, str, Path, int]], dict[str, tuple[Path, int]]]:
    global_rows = []
    label_to_location = {}
    row_id = 1
    for nwb_path in iter_nwb_files(nwb_dir):
        data = load_plot_data(nwb_path)
        if "unit_label" not in data.units:
            raise ValueError(f"{nwb_path} has no unit_label column")
        for unit_idx, unit_label in enumerate(data.units["unit_label"].astype(str)):
            label = normalize_unit_label(unit_label)
            global_rows.append((row_id, label, nwb_path, unit_idx))
            label_to_location[label] = (nwb_path, unit_idx)
            row_id += 1
    return global_rows, label_to_location


def selections_from_nwb_directory(
    nwb_dir: Path,
    global_row_ids: list[int] | None,
    unit_labels: list[str] | None,
    auto_select: bool,
    max_units: int,
) -> list[UnitSelection]:
    if auto_select:
        first_file = iter_nwb_files(nwb_dir)[0]
        data = load_plot_data(first_file)
        return [UnitSelection(data=data, unit_idx=idx) for idx in choose_units(data.units, None, None, max_units, True)]

    global_rows, label_to_location = index_nwb_directory(nwb_dir)

    selected_rows: list[tuple[int | None, str, Path, int]] = []
    if unit_labels is not None:
        missing = [label for label in unit_labels if label not in label_to_location]
        if missing:
            raise ValueError(f"Unit label(s) not found in NWB directory: {', '.join(missing)}")
        for label in unit_labels:
            nwb_path, unit_idx = label_to_location[label]
            selected_rows.append((None, label, nwb_path, unit_idx))
    else:
        requested_rows = global_row_ids or FIGURE2_GLOBAL_ROW_IDS
        max_row = len(global_rows)
        bad = [row_id for row_id in requested_rows if row_id < 1 or row_id > max_row]
        if bad:
            raise ValueError(f"Global row ID out of range: {bad}; indexed {max_row} NWB units")
        for row_id in requested_rows:
            if row_id in FIGURE2_ROW_LABELS:
                label = FIGURE2_ROW_LABELS[row_id]
                if label not in label_to_location:
                    raise ValueError(f"Figure 2 row {row_id} label not found in NWB directory: {label}")
                nwb_path, unit_idx = label_to_location[label]
            else:
                _, label, nwb_path, unit_idx = global_rows[row_id - 1]
            selected_rows.append((row_id, label, nwb_path, unit_idx))

    data_by_path: dict[Path, PlotData] = {}
    selections = []
    for row_id, label, nwb_path, unit_idx in selected_rows:
        if nwb_path not in data_by_path:
            data_by_path[nwb_path] = load_plot_data(nwb_path)
        selections.append(
            UnitSelection(
                data=data_by_path[nwb_path],
                unit_idx=unit_idx,
                global_row_id=row_id,
                row_label=label,
            )
        )
    return selections


def selections_from_nwb_file(
    nwb_path: Path,
    requested_indices: list[int] | None,
    unit_labels: list[str] | None,
    auto_select: bool,
    max_units: int,
) -> list[UnitSelection]:
    if requested_indices is None and unit_labels is None and not auto_select:
        raise ValueError(
            "Default Figure 2 row selection needs an NWB directory. For a single NWB file, "
            "use --units, --unit-labels, or --auto-select."
        )
    data = load_plot_data(nwb_path)
    unit_indices = choose_units(
        data.units,
        requested_indices=requested_indices,
        requested_labels=unit_labels,
        max_units=max_units,
        auto_select=auto_select,
    )
    return [
        UnitSelection(
            data=data,
            unit_idx=idx,
            row_label=str(data.units.iloc[idx].get("unit_label", idx)),
        )
        for idx in unit_indices
    ]


def distance_values(data: PlotData, order: str, include_all_trial_distances: bool) -> np.ndarray:
    if data.tuning_distances is not None and not include_all_trial_distances:
        distances = np.array(sorted(data.tuning_distances), dtype=float)
    else:
        distances = np.array(sorted(data.trials["wall_distance_mm"].dropna().unique()), dtype=float)
    if order == "far-to-near":
        distances = distances[::-1]
    return distances


def make_distance_colors(distances: np.ndarray) -> dict[float, tuple[float, float, float, float]]:
    cmap = plt.get_cmap("turbo_r")
    norm = Normalize(vmin=float(np.nanmin(distances)), vmax=float(np.nanmax(distances)))
    return {float(distance): cmap(norm(float(distance))) for distance in distances}


def smooth_trace(values: np.ndarray, width: int) -> np.ndarray:
    if width <= 1:
        return values
    kernel = np.ones(width, dtype=float) / width
    return np.convolve(values, kernel, mode="same")


def spikes_for_trial(spike_times: np.ndarray, start_time: float, window: tuple[float, float]) -> np.ndarray:
    t0 = start_time + window[0]
    t1 = start_time + window[1]
    left = np.searchsorted(spike_times, t0, side="left")
    right = np.searchsorted(spike_times, t1, side="right")
    return spike_times[left:right] - start_time


def plot_raster(
    ax: plt.Axes,
    spike_times: np.ndarray,
    trials: pd.DataFrame,
    distances: np.ndarray,
    colors: dict[float, tuple[float, float, float, float]],
    window: tuple[float, float],
) -> None:
    y = 0
    yticks = []
    yticklabels = []
    for distance in distances:
        trial_rows = trials[np.isclose(trials["wall_distance_mm"].to_numpy(float), distance)]
        y_start = y
        for _, trial in trial_rows.iterrows():
            rel_spikes = spikes_for_trial(spike_times, float(trial["start_time"]), window)
            if rel_spikes.size:
                ax.vlines(
                    rel_spikes,
                    y + 0.08,
                    y + 0.92,
                    color=colors[float(distance)],
                    linewidth=0.25,
                    alpha=0.85,
                )
            y += 1
        if len(trial_rows):
            yticks.append((y_start + y - 1) / 2)
            yticklabels.append(f"{distance:g}")

    ax.axvspan(3.0, 11.0, color="0.88", zorder=-10)
    ax.set_xlim(window)
    ax.set_ylim(-1, max(y, 1))
    ax.set_yticks([])
    ax.tick_params(axis="both", labelsize=8, length=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_peth(
    ax: plt.Axes,
    spike_times: np.ndarray,
    trials: pd.DataFrame,
    distances: np.ndarray,
    colors: dict[float, tuple[float, float, float, float]],
    window: tuple[float, float],
    bin_width: float,
    smooth_bins: int,
) -> None:
    edges = np.arange(window[0], window[1] + bin_width, bin_width)
    centers = edges[:-1] + bin_width / 2

    for distance in distances:
        trial_rows = trials[np.isclose(trials["wall_distance_mm"].to_numpy(float), distance)]
        if trial_rows.empty:
            continue

        counts = []
        for _, trial in trial_rows.iterrows():
            rel_spikes = spikes_for_trial(spike_times, float(trial["start_time"]), window)
            hist, _ = np.histogram(rel_spikes, bins=edges)
            counts.append(hist / bin_width)
        rates = np.asarray(counts, dtype=float)
        mean = smooth_trace(np.nanmean(rates, axis=0), smooth_bins)
        sem = smooth_trace(np.nanstd(rates, axis=0, ddof=1) / np.sqrt(max(len(rates), 1)), smooth_bins)

        color = colors[float(distance)]
        ax.plot(centers, mean, color=color, linewidth=1.1)
        ax.fill_between(centers, mean - sem, mean + sem, color=color, alpha=0.15, linewidth=0)

    ax.axvspan(3.0, 11.0, color="0.88", zorder=-10)
    ax.set_xlim(window)
    ax.set_ylabel("sp/s", fontsize=8)
    ax.tick_params(axis="both", labelsize=8, length=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def tuning_row_for_unit(data: PlotData, unit_idx: int) -> pd.Series | None:
    if data.tuning is None:
        return None
    unit_label = str(data.units.iloc[unit_idx].get("unit_label", ""))
    if "unit_label" in data.tuning and unit_label:
        match = data.tuning[data.tuning["unit_label"].astype(str) == unit_label]
        if not match.empty:
            return match.iloc[0]
    if unit_idx < len(data.tuning):
        return data.tuning.iloc[unit_idx]
    return None


def plot_tuning(
    ax: plt.Axes,
    data: PlotData,
    unit_idx: int,
    colors: dict[float, tuple[float, float, float, float]],
    tuning_scale: float,
) -> None:
    row = tuning_row_for_unit(data, unit_idx)
    if row is None or data.tuning_distances is None:
        ax.text(0.5, 0.5, "no tuning", ha="center", va="center", transform=ax.transAxes, fontsize=8)
        ax.set_axis_off()
        return

    y = np.asarray(row["zfr_mean"], dtype=float) * tuning_scale
    yerr = (
        np.asarray(row["zfr_sem"], dtype=float) * tuning_scale
        if "zfr_sem" in row
        else np.full_like(y, np.nan)
    )
    x = data.tuning_distances[: len(y)]

    ax.plot(x, y, color="black", linewidth=0.8, zorder=1)
    for xi, yi, ei in zip(x, y, yerr):
        color = colors.get(float(xi), "black")
        ax.errorbar(
            xi,
            yi,
            yerr=None if np.isnan(ei) else ei,
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=0.7,
            capsize=0,
            markersize=3.5,
            zorder=2,
        )

    ax.set_xlabel("Wall distance (mm)", fontsize=8)
    ax.set_ylabel("Firing rate (z-score)", fontsize=8)
    tick_idx = np.linspace(0, len(x) - 1, min(5, len(x)), dtype=int)
    ticks = x[tick_idx]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{v:.0f}" for v in ticks], rotation=0)
    ax.tick_params(axis="both", labelsize=8, length=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def add_distance_legend(
    fig: plt.Figure,
    distances: np.ndarray,
    colors: dict[float, tuple[float, float, float, float]],
) -> None:
    legend_distances = np.array(sorted(distances), dtype=float)
    left = 0.12
    bottom = 0.965
    width = 0.22
    height = 0.018
    ax = fig.add_axes([left, bottom, width, height])
    for i, distance in enumerate(legend_distances):
        ax.bar(i, 1, color=colors[float(distance)], width=1.0, align="edge")
    ax.set_xlim(0, len(legend_distances))
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    fig.text(left, bottom + 0.023, "Wall distance (mm)", fontsize=8, ha="left")
    fig.text(left, bottom - 0.012, "7", fontsize=8, ha="left")
    fig.text(left + width, bottom - 0.012, "23", fontsize=8, ha="right")


def apply_peth_ticks(ax: plt.Axes, row_idx: int) -> None:
    step = PETH_TICK_STEPS[row_idx] if row_idx < len(PETH_TICK_STEPS) else 10
    ymin, ymax = ax.get_ylim()
    upper = max(step, float(np.ceil(ymax / step) * step))
    ax.set_ylim(bottom=0, top=upper)
    ax.set_yticks(np.arange(0, upper + step * 0.5, step))


def make_figure(
    selections: list[UnitSelection],
    window: tuple[float, float],
    bin_width: float,
    smooth_bins: int,
    tuning_scale: float,
    distance_order: str,
    include_all_trial_distances: bool,
) -> plt.Figure:
    distance_sets = [
        distance_values(selection.data, distance_order, include_all_trial_distances)
        for selection in selections
    ]
    distances = np.array(sorted({float(distance) for row_distances in distance_sets for distance in row_distances}))
    if distance_order == "far-to-near":
        distances = distances[::-1]
    colors = make_distance_colors(distances)

    n_rows = len(selections)
    fig_height = max(2.2, 1.45 * n_rows + 0.8)
    fig, axes = plt.subplots(
        n_rows,
        3,
        figsize=(7.0, fig_height),
        sharex="col",
        gridspec_kw={"width_ratios": [1.0, 1.0, 0.9], "wspace": 0.35, "hspace": 0.18},
    )
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    for row_idx, selection in enumerate(selections):
        data = selection.data
        unit_idx = selection.unit_idx
        row_distances = distance_sets[row_idx]
        unit = data.units.iloc[unit_idx]
        spike_times = np.asarray(unit["spike_times"], dtype=float)
        plot_raster(axes[row_idx, 0], spike_times, data.trials, row_distances, colors, window)
        plot_peth(
            axes[row_idx, 1],
            spike_times,
            data.trials,
            row_distances,
            colors,
            window,
            bin_width,
            smooth_bins,
        )
        apply_peth_ticks(axes[row_idx, 1], row_idx)
        plot_tuning(axes[row_idx, 2], data, unit_idx, colors, tuning_scale)

        label = str(unit.get("unit_label", f"unit_{unit_idx}"))
        short_label = label.split("_")[-1] if "_" in label else label
        tuning_type = str(unit.get("tuning_type", ""))
        axes[row_idx, 0].set_ylabel(
            f"Unit {row_idx + 1}\n{tuning_type}",
            fontsize=8,
            rotation=90,
            labelpad=12,
        )
        if selection.global_row_id is None:
            row_prefix = f"{data.session_id}"
        else:
            row_prefix = f"row {selection.global_row_id}"
        axes[row_idx, 0].text(
            0.02,
            0.94,
            f"{row_prefix} / idx {unit_idx} / {short_label}",
            transform=axes[row_idx, 0].transAxes,
            fontsize=7,
            va="top",
            ha="left",
        )

    axes[0, 0].set_title("Rasters", fontsize=10)
    axes[0, 1].set_title("PETHs", fontsize=10)
    axes[0, 2].set_title("Tuning Curves", fontsize=10)
    axes[-1, 0].set_xlabel("Time (s)", fontsize=8)
    axes[-1, 1].set_xlabel("Time (s)", fontsize=8)

    session_ids = [selection.data.session_id for selection in selections]
    if len(set(session_ids)) == 1:
        title = session_ids[0]
    else:
        title = f"{len(selections)} selected units from {len(set(session_ids))} sessions"
    fig.suptitle(title, fontsize=11, y=0.995)
    add_distance_legend(fig, distances, colors)
    fig.subplots_adjust(top=0.93, left=0.10, right=0.98, bottom=0.06)
    return fig


def main() -> None:
    args = parse_args()
    nwb_input = args.nwb_input
    if not nwb_input.exists():
        raise FileNotFoundError(nwb_input)

    if nwb_input.is_dir():
        if args.units is not None:
            raise ValueError("--units can only be used with a single NWB file; use --global-rows for a directory")
        selections = selections_from_nwb_directory(
            nwb_dir=nwb_input,
            global_row_ids=args.global_rows,
            unit_labels=args.unit_labels,
            auto_select=bool(args.auto_select),
            max_units=args.max_units,
        )
        output = args.output or nwb_input.with_name(f"{nwb_input.name}_figure2_style.png")
    else:
        selections = selections_from_nwb_file(
            nwb_path=nwb_input,
            requested_indices=args.units,
            unit_labels=args.unit_labels,
            auto_select=bool(args.auto_select),
            max_units=args.max_units,
        )
        output = args.output or nwb_input.with_name(f"{nwb_input.stem}_figure2_style.png")

    output.parent.mkdir(parents=True, exist_ok=True)

    fig = make_figure(
        selections=selections,
        window=(float(args.window[0]), float(args.window[1])),
        bin_width=float(args.bin_width),
        smooth_bins=int(args.smooth_bins),
        tuning_scale=float(args.tuning_scale),
        distance_order=args.sort_distances,
        include_all_trial_distances=bool(args.include_all_trial_distances),
    )
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)
    print(f"Wrote {output}")
    plotted = []
    for selection in selections:
        label = selection.row_label or str(selection.data.units.iloc[selection.unit_idx].get("unit_label", ""))
        if selection.global_row_id is None:
            plotted.append(f"{selection.data.session_id}:idx{selection.unit_idx}:{label}")
        else:
            plotted.append(f"row{selection.global_row_id}:idx{selection.unit_idx}:{label}")
    print(f"Plotted units: {', '.join(plotted)}")


if __name__ == "__main__":
    main()
