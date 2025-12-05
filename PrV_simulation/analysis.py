import numpy as np
from scipy.signal import find_peaks

def not_flat(curves, threshold=0.1):
    peak_to_peak = curves.max(axis=1) - curves.min(axis=1)
    return peak_to_peak > threshold

def is_monotonic(curve):
    return np.all(np.diff(curve) >= 0) or np.all(np.diff(curve) <= 0)

def is_left_peaked(curve):
    x = np.asarray(curve)
    return (x.size > 0 and x[0] == x.max())

def has_significant_peak(curve, prominence=0.1):
    peaks, props = find_peaks(curve, prominence=prominence)
    return len(peaks) > 0

def is_suppressed(curve, threshold=0.1):
    return max(curve) <= threshold

def analysis(PrV_ori, return_suppressed=False):
    PrV = PrV_ori

    mask_amplitude = not_flat(PrV)
    mask_monotonic = np.array([is_left_peaked(PrV[i]) for i in range(len(PrV))])
    mask_peak = np.array([has_significant_peak(PrV[i]) for i in range(len(PrV))])
    mask_suppressed = np.array([is_suppressed(PrV[i]) for i in range(len(PrV))])

    suppressed = mask_amplitude & mask_suppressed
    only_monotonic = mask_amplitude & mask_monotonic & ~suppressed
    only_peak = mask_amplitude & mask_peak & ~mask_monotonic & ~suppressed
    
    if return_suppressed:
        untuned = (~only_monotonic & ~only_peak) & ~suppressed
        counts = [
            only_monotonic.sum(),   # Proximity
            only_peak.sum(),        # Map
            suppressed.sum(),   # Suppressed
            untuned.sum(),          # Untuned (non-suppressed)
        ]
        print(f'{str(sum(counts))} PrV neurons: Proximity={counts[0]}, Map={counts[1]}, Suppressed={counts[2]}, Untuned={counts[3]}')
    else:
        untuned = (~only_monotonic & ~only_peak)
        counts = [
            only_monotonic.sum(),   # Proximity
            only_peak.sum(),        # Map
            untuned.sum(),          # Untuned (or suppressed)
        ]
        print(f'{str(sum(counts))} PrV neurons: Proximity={counts[0]}, Map={counts[1]}, Untuned={counts[2]}')

    selected_curves = PrV[only_monotonic | only_peak]

    if return_suppressed:
        suppressed_count = mask_suppressed.sum()
        return selected_curves, counts
    else:
        return selected_curves, counts

