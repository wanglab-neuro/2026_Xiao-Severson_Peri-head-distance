import numpy as np
import matplotlib.pyplot as plt
import math
import time
import scipy.io as sio
from scipy.signal import find_peaks
import matplotlib as mpl
from pathlib import Path

def plot(
    selected_curves, 
    counts, 
    selected_curves_fine,
    path, 
    suffix="", 
    return_suppressed=False, 
    plot_examples=False, 
    save_svg=False
    ):
    '''
    Plot summary figures.

    suffix: string appended to filenames and example directory,
            e.g. "" for control, "_telc" for TeLC condition.
    '''
    mpl.rcParams.update({
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "font.family": "Arial",
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
    })

    path = Path(path)

    if return_suppressed:
        labels = ['Proximity', 'Map', 'Suppressed', 'Untuned']
        colors = ['lightcoral', 'khaki', 'lightblue', 'lightgray']
    else:
        labels = ['Proximity', 'Map', 'Untuned']
        colors = ['lightcoral', 'khaki', 'lightgray']

    totalN = sum(counts)

    # ---- Pie chart ---- #
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    wedges, texts, autotexts = ax.pie(
        counts,
        labels=labels,
        autopct='%1.1f%%',
        colors=colors,
        startangle=90,
        counterclock=False,
        wedgeprops={'edgecolor': 'black', 'linewidth': 0.6}
    )
    plt.setp(texts, size=7, fontname="Arial")
    plt.setp(autotexts, size=7, fontname="Arial", color="black")
    ax.set_title(f"N={totalN}", fontsize=8, fontname="Arial")
    fig.tight_layout()

    # PNG + SVG
    fig.savefig(path / f"ratio{suffix}.png")
    if save_svg:
        fig.savefig(path / f"ratio{suffix}.svg")
    plt.close(fig)



    # ---- Proximity vs Map bar chart (percentage of total) ---- #
    prox_count = counts[0]
    map_count = counts[1]
    tuned_count = sum(counts[:-1])
    prox_pct = 100.0 * prox_count / tuned_count if tuned_count > 0 else 0.0
    map_pct = 100.0 * map_count / tuned_count if tuned_count > 0 else 0.0

    fig, ax = plt.subplots(figsize=(2.0, 2.2))
    x_pos = np.arange(2)
    bar_heights = [prox_pct, map_pct]
    bar_colors = [colors[0], colors[1]]  # coral, yellow

    bars = ax.bar(x_pos, bar_heights, color=bar_colors, edgecolor="black", width=0.7)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(['Proximity', 'Map'], fontsize=8)
    ax.set_ylabel('Percent distance-tuned neurons', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Fixed 0–100% y-axis for all conditions
    ax.set_ylim(0, 100)

    # Text labels with percentages above bars
    for x, h in zip(x_pos, bar_heights):
        ax.text(x, h + 0.5, f"{h:.1f}%", ha='center', va='bottom', fontsize=7)

    fig.tight_layout()
    fig.savefig(path / f"prox_map_bar{suffix}.png")
    if save_svg:
        fig.savefig(path / f"prox_map_bar{suffix}.svg")
    plt.close(fig)



    # B = number of bins for stimulus
    B = selected_curves.shape[1]
    cmap = plt.get_cmap("jet_r")
    bin_colors = [cmap(i / (B - 1)) for i in range(B)]

    # ---- Map-center histogram ---- #
    max_indices = np.argmax(selected_curves, axis=1)
    counts_hist, edges = np.histogram(max_indices, bins=np.arange(B + 1) - 0.5)
    fig, ax = plt.subplots(figsize=(2.5, 2))
    ax.bar(np.arange(B), counts_hist, color=bin_colors, edgecolor="black", width=0.8)
    ax.set_xlabel("Map Center", fontsize=8)
    ax.set_ylabel("Number of Neurons", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_xticks(np.arange(B))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()

    # PNG + SVG
    fig.savefig(path / f"map_center{suffix}.png")
    if save_svg:
        fig.savefig(path / f"map_center{suffix}.svg")
    plt.close(fig)

    # ---- Directory for saving tuning-curve examples ---- #
    if plot_examples:
        outdir = path / f"examples_centers{suffix}"
        outdir.mkdir(exist_ok=True)

        for f in outdir.glob("*.png"):
            f.unlink()
        for f in outdir.glob("*.svg"):
            f.unlink()

        k = 5

        for c in np.unique(max_indices):
            idx = np.where(max_indices == c)[0]
            if idx.size == 0:
                continue

            sel = np.random.choice(idx, size=min(k, idx.size), replace=False)

            for uid in sel:
                x = np.arange(B)
                y = selected_curves[uid]
                fig, ax = plt.subplots(figsize=(2.2, 1.8))
                ax.plot(x, y, color='k', lw=1.2)
                ax.scatter(x, y, c=bin_colors, s=12, linewidths=0.25, zorder=3)

                ax.set_xlabel("Distance Bin", fontsize=8)
                ax.set_ylabel("Net Activity", fontsize=8)
                ax.set_xticks(np.arange(B))
                ax.tick_params(labelsize=7)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                fig.tight_layout()

                fig.savefig(outdir / f"example_bin{c:02d}_unit{int(uid):05d}.png")
                if save_svg:
                    fig.savefig(outdir / f"example_bin{c:02d}_unit{int(uid):05d}.svg")
                plt.close(fig)


    # ---- Heatmap ---- #
    plot_unit_heatmap(selected_curves_fine, np.arange(B), '', path, suffix=suffix, save_svg=save_svg)

    return selected_curves


def plot_unit_heatmap(resp, x_values, title, path, suffix="", save_svg=False):
    """
    resp: ndarray (n_units, n_input) each row for a PrV neuron
    x_values: input (n_input,)
    """
    path = Path(path)
    eps = 1e-12
    resp_norm = (resp - resp.min(axis=1, keepdims=True)) / \
        (resp.max(axis=1, keepdims=True) - resp.min(axis=1, keepdims=True) + eps)
    # resp_norm = resp / (resp.max(axis=1, keepdims=True) + eps)

    peak_idx = np.argmax(resp_norm, axis=1)
    sort_order = np.argsort(peak_idx)
    resp_sorted = resp_norm[sort_order]

    cmap = plt.get_cmap("viridis", 256)

    fig, ax = plt.subplots(figsize=(3, 6))
    im = ax.imshow(
        resp_sorted,
        aspect='auto',
        origin='lower',
        extent=[x_values[0], x_values[-1], 0, resp.shape[0]],
        cmap=cmap,
        vmin=0,
        vmax=1,
        interpolation='bilinear',
    )

    ax.set_xlabel("Input", fontsize=10)
    ax.set_ylabel("Units (sorted)", fontsize=10)
    ax.set_title(title, fontsize=12)
    ax.set_xticks(np.arange(len(x_values)))
    ax.invert_yaxis()

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Norm. Response", fontsize=10)

    plt.tight_layout()
    fig.savefig(path / f"heatmap{suffix}.png")
    if save_svg:
        fig.savefig(path / f"heatmap{suffix}.svg")
    plt.close(fig)

    # Dump source data into CSV for plotting into matlab
    np.savetxt(path / f"heatmap{suffix}.csv", resp_sorted, delimiter=',')
