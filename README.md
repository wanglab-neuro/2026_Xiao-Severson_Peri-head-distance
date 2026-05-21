# Xiao, Severson, et al. Peri-head distance coding in the mouse brainstem - Code Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20330279.svg)](https://doi.org/10.5281/zenodo.20330279)

Original code for Xiao & Severson et al., 2026: **"Peri-head distance coding in the mouse brainstem"**.

This repository contains the analysis, processing, figure-generation, and simulation code used for the Peri-head Distance project. The code is organized around three main workflows:

- MATLAB scripts for generating manuscript figure panels.
- A MATLAB electrophysiology processing pipeline for high-channel-count silicon probe recordings.
- A Python simulation of PrV circuit mechanisms receiving excitatory and inhibitory trigeminal inputs.

The repository is intended as a code release for transparency and reuse. Some analysis and figure scripts reference local data paths from the original analysis environment and may require path updates, project-specific data files, or external toolboxes before they can be run on another machine.

## Repository Layout

| Path | Description |
|---|---|
| `Manuscript figures/` | Figure reproduction assets, including original MATLAB scripts and Python helpers for reproducing Figure 2-style panels from deposited NWB files. |
| `Processing/` | MATLAB high-channel-count electrophysiology processing pipeline, including export, spike sorting, curation helpers, waveform analysis, and visualization utilities. |
| `PrV_simulation/` | Python simulation, analysis, and plotting code for PrV neurons receiving trigeminal excitatory and inhibitory inputs. |
| `materials/` | Supporting project materials, including vector artwork. |
| `LICENSE` | MIT license for this repository. |

## Manuscript Figure Assets

The `Manuscript figures/` directory contains the original analysis scripts, and helper assets intended to make selected figures easier to reproduce from deposited NWB files.

| Path | Description |
|---|---|
| `Manuscript figures/Original Matlab code/` | Original MATLAB scripts used to generate manuscript and supplement figure panels. |
| `Manuscript figures/plot_figure2_from_nwb.py` | Recreates Figure 2-style raster, PETH, and tuning-curve panels from local NWB files downloaded from DANDI. |
| `Manuscript figures/generate_deposition_tables.py` | Regenerates helper CSV tables from a directory of downloaded NWB files. |
| `Manuscript figures/figure_unit_manifest.csv` | Mapping from manuscript Figure 2 display units to original MATLAB table rows, DANDI asset paths, and NWB unit rows. |
| `Manuscript figures/unit_tuning_summary.csv` | One-row-per-unit summary regenerated from deposited NWB files. |

The NWB-based helper scripts are designed for the DANDI-deposited dataset:

```text
https://dandiarchive.org/dandiset/001687
```

Install the Python dependencies for these helpers with:

```bash
pip install pynwb pandas numpy matplotlib
```

Then run the Figure 2 reproduction helper from the manuscript-figures directory:

```bash
cd "Manuscript figures"
python plot_figure2_from_nwb.py "path/to/local_nwb_files" --output "path/to/output/figure2_from_nwb.png"
```

See [`Manuscript figures/readme.md`](Manuscript%20figures/readme.md) for full instructions, including how to regenerate `figure_unit_manifest.csv` and `unit_tuning_summary.csv`.

The original MATLAB scripts remain available under `Manuscript figures/Original Matlab code/`. They were written for the original analysis environment and may require path updates, project-specific processed data structures, and external MATLAB toolboxes before they can be run on another machine. Representative dependencies include project-specific session loading utilities, the Neural Decoding Toolbox, and Wang lab raster/PSTH helper functions.

## Electrophysiology Processing Pipeline

The `Processing/` directory contains the high-channel-count electrophysiology pipeline used for silicon probe recordings. It is based on the Wang lab High Channel Count Ephys pipeline and covers:

1. Raw Blackrock/Open Ephys data ingestion.
2. Spike sorting with Kilosort.
3. Manual curation in Phy.
4. Export of curated spike times.
5. Duplicate-unit detection across shanks.
6. Waveform statistics and visualization.

The main batch-script template is:

```matlab
Processing/HCCE_loop.m
```

The main pipeline entry point is:

```matlab
Processing/HCCE_pipeline.m
```

See [`Processing/README.md`](Processing/README.md) for detailed setup, dependency, and usage notes for this pipeline.

## PrV Circuit Simulation

The `PrV_simulation/` directory contains a Python simulation of PrV neurons receiving excitatory and inhibitory trigeminal inputs, together with analysis and plotting utilities.

Install the Python dependencies:

```bash
pip install numpy scipy matplotlib
```

Run the simulation:

```bash
cd PrV_simulation
python simulate_prv.py
```

Run the TeLC/disinhibition condition:

```bash
python simulate_prv.py --telc
```

Common options include:

```bash
python simulate_prv.py --M1 20 --seed 42 --fast
```

Simulation outputs are written to a timestamped results directory. See [`PrV_simulation/Readme.md`](PrV_simulation/Readme.md) for the full option list and parameter descriptions.

## Getting Started

Clone the repository:

```bash
git clone https://github.com/wanglab-neuro/2026_Xiao-Severson_Peri-head-distance.git
cd 2026_Xiao-Severson_Peri-head-distance
```

If working with the processing pipeline and any submodules are present, initialize them with:

```bash
git submodule update --init --recursive
```

For the Python simulation, create an environment using your preferred environment manager and install:

```bash
pip install numpy scipy matplotlib
```

For MATLAB processing and figure scripts, add the relevant repository folders and any external toolbox folders to the MATLAB path before running scripts.

## Data Availability

This repository contains code and supporting materials. Large raw and processed experimental datasets are not stored in the repository. Figure and analysis scripts may need the original processed data files, local session metadata, or analysis intermediates to reproduce the manuscript panels.

## Citation

If you use this code, please cite:

Xiao & Severson et al., 2026. "Peri-head distance coding in the mouse brainstem".

## License

This repository is released under the [MIT License](LICENSE).

Third-party packages used by the processing pipeline retain their original licenses. See [`Processing/README.md`](Processing/README.md) for details.
