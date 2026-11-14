#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the spacing of the purely imaginary vector modes in five dimensions.

The gaps are read from regression/vectorVI.m, which records the spacing Delta Omega between
consecutive overdamped modes of the vector perturbations (n = 3, D = 5) measured at several
resolutions; as that script does, the measurements are averaged for each coupling and multipole
and a straight line is fitted to the averages by least squares. This is Fig. 13 of the article,
redrawn. From this folder,
    python3 plot_readme_figure.py
writes gaussbonnet-gaps.png (Python 3 with NumPy and Matplotlib). The style is the stamped copy
manuscript_style.py beside it.
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

SCRIPT = HERE.parent.parent / "regression" / "vectorVI.m"

rows = re.findall(r"addSpectralGapData\(dataTable,\s*([0-9.]+),\s*([0-9]+),\s*\[([^\]]*)\]\)",
                  SCRIPT.read_text())
data = {}
for alpha, ell, values in rows:
    gaps = np.array([float(v) for v in values.replace(",", " ").split()])
    data.setdefault(float(alpha), []).append((int(ell), gaps))

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, axis = plt.subplots(figsize=(hs.WIDTH, 3.4))
for k, alpha in enumerate(sorted(data)):
    points = sorted(data[alpha], key=lambda p: p[0])
    ell = np.array([p[0] for p in points], float)
    means = np.array([p[1].mean() for p in points])
    A = np.vstack([np.ones_like(ell), ell]).T
    (b0, b1), *_ = np.linalg.lstsq(A, means, rcond=None)
    fit = A @ np.array([b0, b1])
    r2 = 1 - np.sum((means - fit) ** 2) / np.sum((means - means.mean()) ** 2)
    print(f"alpha = {alpha}: Delta Omega = {b1:.6f} l + {b0:.6f}, R2 = {r2:.6f}")
    colour = hs.SERIES[k]
    for (l, gaps) in points:
        axis.plot(np.full(gaps.size, l) + 0.07 * (k - 1), gaps, linestyle="none", marker=hs.MARKERS[k],
                  markersize=2.6, markerfacecolor="white", markeredgecolor=colour, markeredgewidth=0.6,
                  alpha=0.75)
    line = np.linspace(0.8, 4.2, 50)
    axis.plot(line, b0 + b1 * line, color=colour, linestyle=hs.DASHES[k], lw=1.1,
              label=rf"$\hat\alpha = {alpha:.1f}$:  $\Delta\Omega = {b1:.6f}\,\ell + {b0:.5f}$,  $R^2 = {r2:.6f}$")
    axis.plot(ell + 0.07 * (k - 1), means, linestyle="none", marker=hs.MARKERS[k], markersize=5.2,
              markerfacecolor="white", markeredgecolor=colour, markeredgewidth=1.0)
axis.set_xlim(0.8, 4.2)
axis.set_xticks([1, 2, 3, 4])
axis.set_xlabel(r"multipole $\ell$ (vector perturbations, $n = 3$, $D = 5$)")
axis.set_ylabel(r"gap $\Delta\Omega$")
hs.legend_outside(axis, where="below", ncol=1)

out = HERE / "gaussbonnet-gaps.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
