"""Compare three anisotropic interaction potentials:
1. Hard ellipse overlap (Perram-Wertheim contact function)
2. Gaussian overlap potential (from particle_model.jl)
3. Berne-Pechukas (LJ with orientation-dependent sigma)
"""

import numpy as np
import matplotlib.pyplot as plt

l, d = 1.0, 0.3
n = 400
betas = np.linspace(0, np.pi, n)
ys = np.linspace(-3.0, 3.0, n)
BB, YY = np.meshgrid(betas, ys)

# --- 1. Hard ellipse overlap (Perram-Wertheim) ---
# Ellipses with semi-axes (l, d). Particle 1 at origin with alpha=0, particle 2 at (0, y).
# Contact function F = max_lambda f(lambda), overlap iff F <= 1.
# For r=(0,y): f(lam) = lam*(1-lam) * y^2 * Y00 / det(Y)
# where Y = (1-lam)*S1 + lam*S2, S_i = R(theta_i) diag(l^2,d^2) R(theta_i)^T
a2, b2 = l**2, d**2
F = np.zeros((n, n))
for lam in np.linspace(0.01, 0.99, 200):
    s2b = np.sin(BB) ** 2
    Y00 = a2 - lam * (a2 - b2) * s2b
    Y11 = b2 + lam * (a2 - b2) * s2b
    Y01 = lam * (a2 - b2) * np.cos(BB) * np.sin(BB)
    det_Y = Y00 * Y11 - Y01**2
    F = np.maximum(F, lam * (1 - lam) * YY**2 * Y00 / det_Y)
V_overlap = np.maximum(0, 1 - F)

# --- 2. Gaussian overlap potential (vectorized for alpha=0, dx=0) ---
cb, sb = np.cos(BB), np.sin(BB)
S00 = l**2 + l**2 * cb**2 + d**2 * sb**2
S11 = d**2 + l**2 * sb**2 + d**2 * cb**2
S01 = (l**2 - d**2) * cb * sb
det_S = S00 * S11 - S01**2
V_gauss = 1 / (4 * np.pi) * np.sqrt(det_S) * np.exp(-YY**2 * S00 / det_S)

# --- 3. Berne-Pechukas (LJ with orientation-dependent contact distance) ---
# sigma(r_hat, u1, u2) from BP 1972 formula; here alpha=0, r=(0,y)
kappa = l / d
chi = (kappa**2 - 1) / (kappa**2 + 1)
sigma0 = 2 * d  # side-by-side contact distance

sin2B = np.sin(BB) ** 2
cos2B = np.cos(BB) ** 2
sigma = sigma0 * (1 - chi * sin2B / (1 - chi**2 * cos2B)) ** (-0.5)

r = np.maximum(np.abs(YY), 1e-10)
sr6 = (sigma / r) ** 6
V_bp = 4 * (sr6**2 - sr6)

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
titles = ["Ellipse overlap", "Gaussian potential", "Berne-Pechukas"]

for ax, title in zip(axes, titles):
    ax.set_xlabel(r"$\beta$ (orientation of particle 2)")
    ax.set_ylabel(r"$y$ (relative position)")
    ax.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    ax.set_xticklabels([r"$0$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
    ax.set_title(title)

im0 = axes[0].pcolormesh(BB, YY, V_overlap, shading="auto", cmap="inferno")
fig.colorbar(im0, ax=axes[0])

im1 = axes[1].pcolormesh(BB, YY, V_gauss, shading="auto", cmap="inferno")
fig.colorbar(im1, ax=axes[1])

im2 = axes[2].pcolormesh(BB, YY, V_bp, shading="auto", cmap="RdBu_r", vmin=-1, vmax=2)
fig.colorbar(im2, ax=axes[2])

fig.suptitle(rf"$\alpha=0,\;\Delta x=0,\; l={l},\; d={d}$", fontsize=14)
plt.savefig("python/potential_comparison.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved to python/potential_comparison.png")
