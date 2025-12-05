import os
from pathlib import Path
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import math
import time
import scipy.io as sio
from scipy.signal import find_peaks
import matplotlib as mpl

from plot_figures import plot
from analysis import analysis


def sigmoid(x, beta, theta):
    """
    TG neuron's sigmoid tuning curve
    Monotonically decreases from 1 to 0 as distance increases.
    """
    return 1 / (1 + np.exp(beta * (x - theta)))


def exponential_sampling(a, b, N, rate=5.0):
    """
    Exponential sampling for theta parameters, bounded within [a, b].
    """
    u = np.random.rand(N)
    scaled = (np.exp(rate * u) - 1) / (np.exp(rate) - 1)
    x = a + (b - a) * scaled
    return x


def PrV(x1, x, M1, M2, w1, beta1, beta2, theta1, theta2,
        E_local=0.1, I_local=0.0, I_long=None):
    """
    x1       : pre-tiled input for excitation, shape (B, M1)
    x        : original x (B,), only used for tiling inhibition
    M1, M2   : numbers of exc / inh inputs
    w1       : excitatory weights (M1,)
    E_local  : local excitory weights (incr. shifts population toward proximity/map)
    I_local  : local inhibitory weights (incr. shifts population toward suppressed)
    I_long   : long-range inhibitory weights (M2,)
    """
    def relu(x_):
        x_ = np.asarray(x_)
        return np.maximum(x_, 0)

    # Excitatory drive: reuse x1 (no tile per neuron)
    E = (w1 * sigmoid(x1, beta1, theta1)).sum(axis=1)  # (B,)

    # Inhibition: M2 varies, so we still tile x per neuron here
    if I_long is not None and M2 > 0:
        x2 = np.tile(x[:, None], (1, M2))  # (B, M2)
        I_long_drive = (I_long * sigmoid(x2, beta2, theta2)).sum(axis=1)
    else:
        I_long_drive = 0.0

    # Square-root essentially normalizes based on number of inputs
    Z = (E - I_long_drive) / np.sqrt(M1)

    # Linearly sum inputs and outputs
    return Z + E_local - I_local


def sample_discrete_exp_numpy(M, lambd):
    """
    Sample an integer in [0, M] from a discrete exponential distribution
    with parameter 'lambd'.
    """
    values = np.arange(M + 1)
    probs = np.exp(-lambd * values)
    probs /= probs.sum()
    return np.random.choice(values, p=probs)


def main(
    M1=20,
    telc=False,
    path='.',
    seed=42,
    telc_inh_scale=0.5,
    p_untuned=0.4,
    weak_gain=0.01,
    strong_gain=1.0,
    local_exc_gain=1.0,
    local_inh_gain=1.7,
    long_inh_gain=10,
    return_suppressed=True,
    fast=False,
):
    """
    Simulate PrV neurons, analyze tuning curves, and plot results.

    Parameters
    ----------
    M1 : int
        Number of excitatory TG inputs per PrV neuron.
    telc : bool
        If True, also simulate a TeLC manipulation condition (scaled long-range inhibition).
    path : str
        Parent directory for saving results.
    seed : int
        Random seed for reproducibility.
    telc_inh_scale : float
        Multiplicative factor for long-range inhibitory weights in TeLC condition
        (for non-weak neurons). < 1 = reduced inhibition.
    p_untuned : float
        Probability that a neuron is weakly tuned (weights scaled by weak_gain).
    weak_gain : float
        Scaling factor applied to w1 and w2_long for weakly tuned neurons.
    local_inh_gain : float
        Gain factor for local (TeLC-insensitive) inhibition I_local.
    """
    # Fix random seed for reproducibility of the simulation
    np.random.seed(seed)

    # Base output path
    base_path = Path(path)

    ts = time.perf_counter()

    # Timestamped results folder
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = base_path / f"results_{timestamp}"
    outdir.mkdir(parents=True, exist_ok=True)

    # Stimulus space: distance bins
    num_distances = 9
    num_bins = 99
    N = 1225                                   # total number of PrV neurons = 1225
    x_fine = np.linspace(0, 1, num_bins)         # fine grid
    x1_fine = np.tile(x_fine[:, None], (1, M1))  # pre-tiled excitation input

    PrV_record_control = []
    if telc:
        PrV_record_telc = []

    for _ in range(N):
        # Sample number of long-range inhibitory inputs
        M2 = sample_discrete_exp_numpy(int(M1), 0.5) # 0.03

        # Excitatory / long-range inhibitory weights and TG parameters
        w1 = strong_gain * np.random.rand(M1)
        I_long = long_inh_gain * np.random.rand(M2)

        beta = 20  # base slope of TG distance tuning
        k = 5  # range of variability in slopes across inputs
        beta1 = beta + k * np.random.rand(M1)
        beta2 = beta + k * np.random.rand(M2)
        theta1 = exponential_sampling(0.0, 1.0, M1, rate=2.0) # 3
        theta2 = exponential_sampling(0.0, 1.0, M2, rate=4.0) # 6

        # Background excitation and local inhibition
        E_local = local_exc_gain * np.random.rand()
        I_local = local_inh_gain * np.random.rand()

        # Decide if this neuron is weakly tuned / candidate untuned
        is_weak = np.random.rand() < p_untuned
        if is_weak:
            # Shrink both excitation and long-range inhibition so modulation is tiny.
            # Local inhibition remains as-is and is NOT TeLC-scaled.
            w1 *= weak_gain
            I_long *= weak_gain

        # CONTROL: full long-range + local inhibition
        control_resp = PrV(
            x1_fine, x_fine, M1, M2,
            w1,
            beta1, beta2, theta1, theta2,
            E_local=E_local,
            I_local=I_local,
            I_long=I_long,
        )
        PrV_record_control.append(control_resp)

        # TeLC: only scale long-range inhibition for NON-weak neurons
        if telc:
            if is_weak:
                # Weakly tuned neurons are unaffected by TeLC
                I_telc = I_long
            else:
                I_telc = telc_inh_scale * I_long

            telc_resp = PrV(
                x1_fine, x_fine, M1, M2,
                w1,
                beta1, beta2, theta1, theta2,
                E_local=E_local,  # same local excitation
                I_local=I_local,  # same local inhibition
                I_long=I_telc,
            )
            PrV_record_telc.append(telc_resp)

    # Convert to arrays
    print('Analyzing models for wild-type condition ... ')
    PrV_record_control = np.array(PrV_record_control)

    # Downsample fine curves into 9 coarse distance bins
    if num_bins % num_distances != 0:
        raise ValueError("num_distances must be divisible by num_bins for simple reshape binning.")

    bin_size = num_bins // num_distances  # e.g., 99 // 9 = 11

    # reshape: (N, num_distances) -> (N, num_bins, bin_size), then average over the fine samples
    PrV_record_control_binned = PrV_record_control.reshape(
        PrV_record_control.shape[0], num_distances, bin_size
        ).mean(axis=2)  # -> (N, num_bins)

    # CONTROL: get tuning categories, including suppressed count
    selected_curves_ctrl, counts_ctrl = analysis(
        PrV_record_control_binned,
        return_suppressed=return_suppressed
    )

    # CONTROL: get fine tuning curves for heatmap
    selected_curves_ctrl_fine, _ = analysis(
        PrV_record_control,
        return_suppressed=return_suppressed
    )

    if telc:
        PrV_record_telc = np.array(PrV_record_telc)  # (N, num_distances)
        PrV_record_telc_binned = PrV_record_telc.reshape(
            PrV_record_telc.shape[0], num_distances, bin_size
        ).mean(axis=2)

        selected_curves_telc, counts_telc = analysis(
            PrV_record_telc_binned,
            return_suppressed=return_suppressed
        )

        selected_curves_telc_fine, _ = analysis(
            PrV_record_telc,
            return_suppressed=return_suppressed
        )

    print(f"Simulation and analysis finished in {(time.perf_counter() - ts):.2f} seconds.")
    if not fast:
        print(f"Plotting and saving figures. This may take some time...")
    else:
        print(f"Plotting and saving figures.")

    # Plot control condition
    plot(
        selected_curves_ctrl,
        counts_ctrl,
        selected_curves_ctrl_fine,
        outdir,
        suffix="",
        return_suppressed=return_suppressed,
        plot_examples=not fast,
        save_svg=not fast,
    )

    # Plot TeLC condition
    if telc:
        plot(
            selected_curves_telc,
            counts_telc,
            selected_curves_telc_fine,
            outdir,
            suffix="_telc",
            return_suppressed=return_suppressed,
            plot_examples=not fast,
            save_svg=not fast,
        )

    print(f"Done! Figures saved in {outdir}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Simulate PrV neurons and analyze tuning curves.")

    parser.add_argument("--M1", type=int, default=20,
        help="Number of excitatory TG inputs per PrV neuron (default: 20).")

    parser.add_argument("--path", type=str, default=".",
        help="Parent directory for saving results (default: current directory).")

    parser.add_argument("--seed", type=int, default=42,
        help="Random seed (default: 42).")

    parser.add_argument("--telc_inh_scale", type=float, default=0.5,
        help="Multiplicative factor for long-range inhibitory weights in TeLC condition (default: 0.5).")

    parser.add_argument("--p_untuned", type=float, default=0.4,
        help="Probability that a neuron is weakly tuned (default: 0.4).")

    parser.add_argument("--strong_gain", type=float, default=1.0,
        help="Scaling for weights of strongly tuned neurons (default: 1.0).",)

    parser.add_argument("--weak_gain", type=float, default=0.01,
        help="Scaling for weights of weakly tuned neurons (default: 0.01).",)

    parser.add_argument("--local_exc_gain", type=float, default=1.0,
        help="Strength of local (TeLC-insensitive) excitation (default: 1.0).",)

    parser.add_argument("--local_inh_gain", type=float, default=1.7,
        help="Strength of local (TeLC-insensitive) inhibition (default: 1.7).",)

    parser.add_argument("--long_inh_gain", type=float, default=10,
        help="Strength of long (TeLC-sensitive) inhibition (default: 10).",)

    parser.add_argument("--return_suppressed", action="store_false",
        help="Also return number of suppressed neurons (otherwise include in untuned; default: True).")

    parser.add_argument("--fast", action="store_true",
        help="Fast mode: no example tuning curves, no SVG export.")

    parser.add_argument("--telc", action="store_true",
        help="Also simulate TeLC condition (scaled long-range inhibition).")


    args = parser.parse_args()
    main(
        M1=args.M1,
        telc=args.telc,
        path=args.path,
        seed=args.seed,
        telc_inh_scale=args.telc_inh_scale,
        p_untuned=args.p_untuned,
        strong_gain=args.strong_gain, 
        weak_gain=args.weak_gain,
        local_exc_gain=args.local_exc_gain,
        local_inh_gain=args.local_inh_gain,
        long_inh_gain=args.long_inh_gain,
        return_suppressed=args.return_suppressed,
        fast=args.fast,
    )
