# Figure Reproduction Assets

Manuscript figures were generated from MATLAB scripts in the `Manuscript figures\Original Matlab code` directory.  

Some helper files are provided here to facilitate reproducing the manuscript figures from the NWB files deposited on [DANDI](https://dandiarchive.org/dandiset/001687).  

## Files

| File | Purpose |
| --- | --- |
| `plot_figure2_from_nwb.py` | Recreates the Figure 2-style raster, PETH, and tuning-curve panels from the deposited NWB files. |
| `generate_deposition_tables.py` | Regenerates the two CSV helper tables below from a directory of downloaded NWB files. |
| `figure_unit_manifest.csv` | Explicit mapping from manuscript Figure 2 display units to MATLAB table rows, subject IDs, DANDI asset paths, and NWB unit rows. |
| `unit_tuning_summary.csv` | Dataset-level one-row-per-unit summary table regenerated from the NWB files, including subject IDs and DANDI asset paths. |

## Requirements

Use Python 3.10+ with the project requirements installed. The scripts need at least:

```sh
pip install pynwb pandas numpy matplotlib
```

On Windows from this repository, the existing local environment can be used:

```bat
cd "<path to this repository>\Manuscript figures"
.venv\Scripts\python -m pip install -r requirements.txt
```

## Regenerate Helper Tables

After downloading the DANDI NWB files into one directory, run:

```bat
cd "<path to this repository>\Manuscript figures"
.venv\Scripts\python "generate_deposition_tables.py" ^
    "path\to\local_nwb_files" ^
    --output-dir "for deposition"
```

This writes:

```text
figure_unit_manifest.csv
unit_tuning_summary.csv
```

`figure_unit_manifest.csv` is needed because the manuscript Figure 2 unit numbers refer to 1-based global row indices in the original MATLAB tables (`PSTH_naive.mat:psth.dataTable` and `unitTable_naive.mat:unitTable`). Those global MATLAB row numbers are not native NWB row indices. The manifest gives the reproducible bridge:

```text
matlab_global_row_1based -> matlab_row_label -> nwb_file -> nwb_unit_index_0based
```

For DANDI-organized downloads, use the `subject_id` and `dandi_asset_path` columns. For example:

```text
subject_id = Ephys2
dandi_asset_path = sub-Ephys2/sub-Ephys2_ses-20230208-sp-session-1.nwb
```

The original conversion output names files as `<session_id>.nwb`, for example `20230208_sp_session_1.nwb`. DANDI-organized downloads rename and place assets by subject/session. The manifest includes both forms:

```text
nwb_file         = 20230208_sp_session_1.nwb
dandi_asset_path = sub-Ephys2/sub-Ephys2_ses-20230208-sp-session-1.nwb
```

## Recreate Figure 2 Panels

Run the plotting script on the downloaded NWB directory:

```bat
cd "<path to this repository>\Manuscript figures"
.venv\Scripts\python "plot_figure2_from_nwb.py" ^
    "path\to\local_nwb_files" ^
    --output "path\to\output\figure2_from_nwb.png"
```

By default the script plots the six Figure 2 units listed in `figure_unit_manifest.csv`. Use `--global-rows` to plot other 1-based MATLAB table rows when a row-to-label mapping is known, or `--unit-labels` to plot explicit labels such as `20230208_sp_session_1_unit_5`.