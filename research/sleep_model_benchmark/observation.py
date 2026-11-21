"""Convert native events or explicit rate observations to count segments."""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Counts:
    segments: tuple
    bin_s: float
    metadata: dict

    def population(self):
        return tuple(block.sum(axis=1, dtype=np.int64) for block in self.segments)

    def rebin(self, factor):
        if not isinstance(factor, int) or factor < 1:
            raise ValueError("factor must be a positive integer")
        blocks = []
        for block in self.segments:
            n = len(block) // factor
            if n == 0:
                raise ValueError("Segment cannot supply a complete output bin")
            blocks.append(block[:n * factor].reshape(n, factor, block.shape[1]).sum(1))
        width = self.bin_s * factor
        return Counts(tuple(blocks), width, self.metadata | {
            "bin_s": width, "observed_seconds": [len(b) * width for b in blocks],
        })


@dataclass(frozen=True)
class Continuous:
    segments: tuple
    sample_hz: float
    unit: str


def _width(bin_s):
    if not np.isfinite(bin_s) or bin_s <= 0:
        raise ValueError("bin_s must be positive seconds")


def _whole(duration_s, bin_s):
    return int(np.floor(np.nextafter(duration_s / bin_s, np.inf)))


def _metadata(blocks, bin_s, **details):
    seconds = sum(map(len, blocks)) * bin_s
    per_unit = sum((block.sum(0) for block in blocks)) / seconds
    return dict(details, units=len(per_unit), bin_s=bin_s,
                observed_seconds=[len(block) * bin_s for block in blocks],
                per_unit_observed_hz=per_unit.tolist())


def spike_counts(segments, bin_s=.04):
    """Each segment is (start_s, stop_s, one ordered spike array per fixed unit)."""
    _width(bin_s)
    blocks, intervals = [], []
    for start, stop, by_unit in segments:
        if (not np.isfinite([start, stop]).all() or start >= stop or
                (intervals and start < intervals[-1][1])):
            raise ValueError("Segments must be ordered and disjoint")
        n = _whole(stop - start, bin_s)
        if not n or not len(by_unit) or (blocks and len(by_unit) != blocks[0].shape[1]):
            raise ValueError("Each segment needs complete bins and the same nonempty unit set")
        # Searchsorted gives half-open bins. The final partial bin is discarded.
        edges = start + np.arange(n + 1) * bin_s
        cols = []
        for train in by_unit:
            t = np.asarray(train, dtype=float)
            if (t.ndim != 1 or not np.isfinite(t).all() or np.any(np.diff(t) < 0) or
                    np.any(t < start) or np.any(t >= stop)):
                raise ValueError("Expected ordered spike times in [start, stop) in seconds")
            cols.append(np.diff(np.searchsorted(t, edges, side="left")))
        blocks.append(np.stack(cols, axis=1))
        intervals.append([start, stop])
    if not blocks:
        raise ValueError("No spike segments")
    return Counts(tuple(blocks), bin_s, _metadata(blocks, bin_s, mode="native_spikes",
                  requested_intervals_s=intervals, seed=None, noise="none"))


def rate_counts(segments_hz, sample_s, *, bin_s=.04, units=None, seed, shared_shape=None):
    """Integrate piecewise-constant Hz; 1-D is a per-neuron population mean."""
    _width(sample_s)
    _width(bin_s)
    if shared_shape is not None and (not np.isfinite(shared_shape) or shared_shape <= 0):
        raise ValueError("Gamma shape must be positive")
    if units is not None and (not isinstance(units, int) or units < 1):
        raise ValueError("units must be a positive integer")
    rng = np.random.default_rng(seed)
    blocks, source_means, samples = [], [], []
    dimension = None
    for rate in segments_hz:
        rate = np.asarray(rate, dtype=float)
        broadcast_mean = rate.ndim == 1
        if dimension is not None and dimension != rate.ndim:
            raise ValueError("Mixed rate formats within one recording")
        dimension = rate.ndim
        if broadcast_mean:
            if units is None:
                raise ValueError("A population mean rate requires an explicit observed unit count")
            rate = rate[:, None]
        elif units is not None:
            raise ValueError("Individual rate channels already define the units")
        if rate.ndim != 2 or not rate.size or not np.isfinite(rate).all() or np.any(rate < 0):
            raise ValueError("Rates must be finite nonnegative Hz over one or more channels")
        n = _whole(len(rate) * sample_s, bin_s)
        if not n:
            raise ValueError("No complete count bin in rate segment")
        # Integrate at all edges, including overlaps when resolutions are not integer multiples.
        prefix = np.vstack((np.zeros(rate.shape[1]), np.cumsum(rate, axis=0) * sample_s))
        x = np.minimum(np.arange(n + 1) * (bin_s / sample_s), len(rate))
        index = x.astype(np.int64)
        integral = prefix[index] + (x-index)[:, None] * sample_s * rate[np.minimum(index, len(rate)-1)]
        lam = np.maximum(0, np.diff(integral, axis=0))
        if broadcast_mean:
            lam = np.broadcast_to(lam, (n, units))
        if blocks and lam.shape[1] != blocks[0].shape[1]:
            raise ValueError("Unit count changed between segments")
        source_means.append(lam.sum(0))
        samples.append(len(rate))
        if shared_shape is not None:
            lam = lam * rng.gamma(shared_shape, 1/shared_shape, size=(n, 1))
        blocks.append(rng.poisson(lam))
    if not blocks:
        raise ValueError("No rate segments")
    mode = "population_mean_rate" if dimension == 1 else "individual_rates"
    meta = _metadata(blocks, bin_s, mode=mode, seed=int(seed), sample_s=sample_s,
                     sampling_hz=1/sample_s, source_samples=samples,
                     noise="poisson" if shared_shape is None else "shared_gamma_poisson",
                     shared_shape=shared_shape, conditionally_independent=True,
                     expected_unit_hz=(sum(source_means) / (sum(map(len, blocks)) * bin_s)).tolist())
    return Counts(tuple(blocks), bin_s, meta)


def continuous_signal(segments, sample_hz, unit):
    """Preserve sampled EEG/LFP in its physical units on a separate branch."""
    _width(sample_hz)
    segments = tuple(np.asarray(s, dtype=float) for s in segments)
    if (not unit or not segments or
            any(s.ndim != 1 or not len(s) or not np.isfinite(s).all() for s in segments)):
        raise ValueError("Expected finite 1-D samples and an explicit physical unit")
    return Continuous(segments, sample_hz, unit)
