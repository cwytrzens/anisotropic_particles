"""Plot the anisotropic particle interaction potential."""

import numpy as np
import matplotlib.pyplot as plt


def potential(xy, alpha, beta, l, d, xi=0.0):
    R = np.array(xy)
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    gamma_i = (l**2 - d**2) * np.array([[ca**2, ca*sa], [ca*sa, sa**2]]) + d**2 * np.eye(2)
    gamma_j = (l**2 - d**2) * np.array([[cb**2, cb*sb], [cb*sb, sb**2]]) + d**2 * np.eye(2)
    Sigma = gamma_i + gamma_j
    return 1 / (4 * np.pi) * np.linalg.det(Sigma) ** (0.5 - xi) * np.exp(-R @ np.linalg.inv(Sigma) @ R)


l, d = 1.0, 0.3
n = 400
betas = np.linspace(0, np.pi, n)
ys = np.linspace(-3.0, 3.0, n)
BB, YY = np.meshgrid(betas, ys)

V = np.zeros_like(BB)
for i in range(n):
    for j in range(n):
        V[i, j] = potential((0.0, YY[i, j]), 0.0, BB[i, j], l, d)

fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
im = ax.pcolormesh(BB, YY, V, shading="auto", cmap="inferno")
ax.set_xlabel(r"$\beta$ (orientation of particle 2)")
ax.set_ylabel(r"$y$ (relative position)")
ax.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
ax.set_xticklabels([r"$0$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
ax.set_title(rf"Potential ($\alpha=0,\;\Delta x=0,\; l={l},\; d={d}$)")
fig.colorbar(im, ax=ax)
plt.savefig("python/potential_beta_y.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved to python/potential_beta_y.png")
