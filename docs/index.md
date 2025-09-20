# Phased‑Array Radar Simulation — Guide & Figures

This page explains what each figure shows and how to reproduce them.

> All images in this page are expected under `docs/`. Generate outputs by running the simulation, then copy PNGs from `results/` to `docs/` and commit.

## Figures

### 1) Array face (`face.png`)
- **What:** Hexagonal (triangular lattice) element layout. Default is 3,333 elements at ~λ/2 spacing.
- **Why:** Sanity‑check geometry, aperture size, and steering direction.

### 2) u–v array factor (`uv_beam.png`)
- **What:** Normalized array factor across the visible region in direction‑cosine space (u, v).
- **Why:** Shows mainlobe and sidelobe field. Circle indicates `u^2+v^2=1` boundary.

### 3) Azimuth cut (`fft_az_cut.png`)
- **What:** Principal azimuth cut at a fixed elevation plane (default 0°), taken from the **aperture‑FFT beampattern**.
- **Why:** The −3 dB crossings are used to report **HPBW** (beamwidth).

### 4) Aperture‑FFT beampattern (`fft_uv.png`)
- **What:** Far‑field pattern via 2‑D FFT of the aperture field deposited on a grid with zero‑padding.
- **Why:** A near‑“ground truth” beampattern for a sampled planar array.

### 5) 3‑D beampattern (`beampattern3d.png`)
- **What:** Sampled power over (az, el) grid (visualized as a heatmap).
- **Why:** Quick intuition of pattern asymmetries and scan effects.

### 6) FMCW chirp (`chirp.png`) and Spectrogram (`spectrogram.png`)
- **What:** Time‑domain LFM chirp and a DIY STFT spectrogram (no toolboxes).
- **Why:** Visualize instantaneous frequency and bandwidth coverage.

### 7) Range–Doppler map (`range_doppler.png`)
- **What:** Toy single‑target map from range FFT (fast‑time) and Doppler FFT (slow‑time).
- **Why:** End‑to‑end FMCW signal processing sanity check.

### 8) Summary (`summary.png`)
- **What:** Textual metrics rendered into a figure (also saved as `summary.pdf`).
- **Contains:**
  - Range & velocity resolution
  - HPBW at boresight and at scan, scan loss
  - Unambiguous range/velocity
  - Grating‑lobe spacing safety check

---

## Reproducing the figures

```matlab
cd matlab
run_demo              % quick start (defaults)
% or use presets:
run_preset('high_res_range')
run_preset('wide_beam_low_elements')
run_preset('fast_scan_az45')
run_preset('long_range_high_prf')
```

Outputs are saved to `results/` (PNG/PDF/FIG). To publish example figures:

```bash
mkdir -p docs
cp results/{face.png,uv_beam.png,fft_az_cut.png,range_doppler.png,summary.png} docs/
git add docs/*.png
git commit -m "docs: add example figures"
git push
```

---

## Tips
- Prefer running from a writable folder (so `results/` can be created).
- To save into a timestamped subfolder, set `useTimestampSubfolder = true` in options/presets.
