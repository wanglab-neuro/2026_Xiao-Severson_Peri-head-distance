# PrV Circuit Simulation

This repository contains a simple simulation of PrV neurons receiving excitatory and inhibitory TG inputs, plus analysis and plotting code.

## Installation

1. **Clone or download** this repository.
cd /d C:\Code
```
git clone https://github.com/wanglab-neuro/2026_Xiao-Severson_Peri-head-distance
```

3. (Optional but recommended) Create and activate a virtual environment:
```
conda create -n prv
conda activate prv
```

4.  Install dependencies:
```
pip install numpy scipy matplotlib
```

Make sure the following files are in the current directory:
```
cd ./PrV_simulation
```

simulate_prv.py (this script)
plot_figures.py
analysis.py

## Usage
Run the simulation from the command line:
python simulate_prv.py [OPTIONS]

Example
```
python simulate_prv.py --telc
```
This will run the simulation, perform analysis, and save the generated figures/results under ./results_YYYYMMDD_HHMMSS

Options
--M1 <int>
Number of excitatory TG inputs to each PrV neuron
Default: 20

--path <str>
Directory for saving results (figures, etc.)
Default: current directory (.)

--telc
Enable partial ablation of inhibition (reduces inhibitory input).
Default: off (not included)

--fast
For faster testing of parameters, skips output of example tuning curves, vector object files (SVG).
Default: off

--seed <int>
Random seed for reproducibility.
Default: 42

Model parameters

These are the defaults used for the figure (control condition):

--telc_inh_scale <float>
Multiplicative factor applied to long-range inhibitory weights in the TeLC condition (for non-weak neurons).

< 1 reduces long-range inhibition (simulated disinhibition)
Default: 0.5

--p_untuned <float>
Probability that a given neuron is “weakly tuned” (both excitation and long-range inhibition scaled by weak_gain, making the tuning curve low-amplitude and more likely to be classified as untuned).
Default: 0.4

--strong_gain <float>
Global scaling factor for excitatory TG weights of strongly tuned neurons.
Default: 1.0

--weak_gain <float>
Scaling factor applied to both excitation and long-range inhibition for weakly tuned neurons (is_weak = True).
Default: 0.01

--local_exc_gain <float>
Strength of local (TeLC-insensitive) excitatory drive (E_local).
Default: 1.0

--local_inh_gain <float>
Strength of local (TeLC-insensitive) inhibitory drive (I_local).
Default: 1.7

--long_inh_gain <float>
Strength of long-range (TeLC-sensitive) inhibitory input from SpVi-like neurons.
Default: 10
