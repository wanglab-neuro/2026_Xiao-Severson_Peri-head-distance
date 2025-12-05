# PrV Circuit Simulation

This repository contains a simple simulation of PrV neurons receiving excitatory and inhibitory TG inputs, plus analysis and plotting code.

## Installation

1. **Clone or download** this repository.
cd /d C:\Code
git clone https://github.com/wanglab-neuro/2026_Xiao-Severson_Peri-head-distance

2. (Optional but recommended) Create and activate a virtual environment:
conda create -n prv
conda activate prv

Install dependencies:
pip install numpy scipy matplotlib

Make sure the following files are in the current directory:
cd /d C:\Code
simulate_prv.py (this script)
plot_figures.py
analysis.py

## Usage
Run the simulation from the command line:
python simulate_prv.py [OPTIONS]

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

--seed <int>
Random seed for reproducibility.
Default: 42

Example
Run an opto simulation with custom parameters and a fixed seed:

python simulate_prv.py \
  --M1 30 \
  --telc \
  --seed 42
This will run the simulation, perform analysis, and save the generated figures/results under ./results_YYYYMMDD_HHMMSS