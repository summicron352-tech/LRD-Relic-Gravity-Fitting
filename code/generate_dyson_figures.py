#!/usr/bin/env python3
"""
Generate ApJ-quality Dyson entrainment model figures.
All English labels. Pure Newtonian gravity. Corrected physics.

Figures:
  Fig 8: Density profile (double power-law)
  Fig 9: Velocity cascade — Dyson entrainment + geometric dilution
  Fig 10: Collision time vs radius
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.size': 12, 'axes.labelsize': 14, 'xtick.labelsize': 12,
    'ytick.labelsize': 12, 'legend.fontsize': 9, 'font.family': 'serif',
    'axes.linewidth': 1.2,
})
import os

# ===========================
# Constants
# ===========================
G_N = 4.30091e-3          # [pc (km/s)^2 / Msun]
MSUN = 1.989e33; PC = 3.086e18; YR = 3.156e7
output_dir = os.path.dirname(os.path.abspath(__file__))
fig_dir = os.path.join(os.path.dirname(output_dir), 'figures')
os.makedirs(fig_dir, exist_ok=True)

# ===========================
# GLIMPSE-17775 parameters
# ===========================
Mstar = 10**9.51; r_eff = 83.3; z_src = 5.66
r_c = 0.003; r_b = 0.3; alpha_in = 2.0; alpha_out = 1.2
m_star_avg = 0.5; R_star_avg = 1.0 * 6.96e10 / PC

r_grid = np.logspace(-4, 2, 2000)

# ===========================
# Profile functions
# ===========================
def density_profile(r, r_c, r_b, ai, ao, rho_c):
    rho = np.zeros_like(r)
    rho[r < r_b] = rho_c * (r[r < r_b] / r_c)**(-ai)
    outer = r >= r_b
    rho_at_b = rho_c * (r_b / r_c)**(-ai)
    rho[outer] = rho_at_b * (r[outer] / r_b)**(-ao)
    rho[r < r_c * 0.1] = rho_c * (0.1)**(-ai)
    return rho

def enclosed_mass(r, rho):
    dr = np.diff(np.concatenate([[0], r]))
    return np.cumsum(rho * 4 * np.pi * r**2 * dr)

# Build profiles
rho_prof = density_profile(r_grid, r_c, r_b, alpha_in, alpha_out, 1e6)
M_enc = enclosed_mass(r_grid, rho_prof)
idx_eff = np.argmin(np.abs(r_grid - r_eff))
scale_factor = Mstar / M_enc[idx_eff]
rho_prof *= scale_factor; M_enc = enclosed_mass(r_grid, rho_prof)
Sigma_r = M_enc / (np.pi * r_grid**2)

# Pure Newtonian v_grav
v_grav = np.sqrt(G_N * M_enc / r_grid)

# ===========================
# Dyson Entrainment (corrected)
# ===========================
def dyson_entrainment(r, rho, v_grav, r_seed, eta=1.75, alpha_dil=1.0, r_break=0.3):
    """Dyson Air Multiplier entrainment + geometric dilution."""
    n = len(r); v_total = np.zeros(n)
    idx_seed = np.argmin(np.abs(r - r_seed))
    v_seed = v_grav[idx_seed]
    v_total[idx_seed] = v_seed; rho_seed = rho[idx_seed]

    # Outward: entrainment + geometric dilution
    for i in range(idx_seed + 1, n):
        dlnr = np.log(r[i] / r[i-1])
        drho = max(rho[i] / rho_seed, 1e-10)
        coupling = eta * drho**0.3
        entrain = coupling * v_total[i-1]**2 * dlnr
        dilution = alpha_dil * (r[i] / r_break) * v_total[i-1]**2 * dlnr
        grav = v_grav[i]**2 - v_grav[i-1]**2
        v_total[i] = np.sqrt(max(1.0, v_total[i-1]**2 + grav + entrain - dilution))

    # Inward: entrainment only
    for i in range(idx_seed - 1, -1, -1):
        dlnr = np.log(r[i+1] / r[i])
        drho = max(rho[i] / rho_seed, 1e-10)
        coupling = eta * drho**0.3
        entrain = coupling * v_total[i+1]**2 * dlnr
        grav = max(0, v_grav[i]**2 - v_grav[i+1]**2)
        v_total[i] = np.sqrt(max(1.0, v_total[i+1]**2 + grav + entrain))

    return v_total

# Best parameters
M_central = 4e4  # Msun (compact cluster / small IMBH)
M_tot = M_enc + M_central
v_grav_mc = np.sqrt(G_N * M_tot / r_grid)
r_seed = 0.0007
eta_best = 1.75
v_ent = dyson_entrainment(r_grid, rho_prof, v_grav_mc, r_seed, eta_best)

idx_blr = np.argmin(np.abs(r_grid - 0.03))
fwhm_blr = 2.355 * v_ent[idx_blr]
v_seed_val = v_grav_mc[np.argmin(np.abs(r_grid - r_seed))]

print(f"Dyson Entrainment Model")
print(f"  Seed: r={r_seed:.4f} pc, v={v_seed_val:.0f} km/s")
print(f"  BLR (0.03 pc): v={v_ent[idx_blr]:.0f} km/s, FWHM={fwhm_blr:.0f} km/s")
print(f"  M_central = {M_central:.1e} Msun")

# ===========================
# Fig 8: Density Profile
# ===========================
print("\n--- Fig 8: Density Profile ---")

fig8, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.loglog(r_grid, rho_prof, 'b-', linewidth=2)
ax1.axvline(r_c, color='gray', linestyle='--', alpha=0.5, label=f'$r_c$={r_c:.4f} pc')
ax1.axvline(r_b, color='gray', linestyle=':', alpha=0.5, label=f'$r_b$={r_b:.2f} pc')
ax1.axvline(r_eff, color='red', linestyle='--', alpha=0.5, label=f'$r_{{\\rm eff}}$={r_eff:.0f} pc')
ax1.set_xlabel('r [pc]', fontsize=13)
ax1.set_ylabel(r'$\rho$ [$M_\odot$ pc$^{-3}$]', fontsize=13)
ax1.set_title('Density Profile $\\rho(r)$', fontsize=14)
ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3, which='both')

ax2.loglog(r_grid, Sigma_r, 'g-', linewidth=2)
ax2.axvline(r_c, color='gray', linestyle='--', alpha=0.5)
ax2.axvline(r_b, color='gray', linestyle=':', alpha=0.5)
ax2.axvline(r_eff, color='red', linestyle='--', alpha=0.5)
ax2.set_xlabel('r [pc]', fontsize=13)
ax2.set_ylabel(r'$\Sigma$ [$M_\odot$ pc$^{-2}$]', fontsize=13)
ax2.set_title('Surface Density Profile $\\Sigma(r)$', fontsize=14)
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'fig8_density_profile.png'), dpi=200)
plt.close()
print("  -> fig8_density_profile.png")

# ===========================
# Fig 9: Velocity Cascade
# ===========================
print("\n--- Fig 9: Velocity Cascade ---")

fig9, ax9 = plt.subplots(figsize=(10, 7))

# Show several eta values
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
for j, eta in enumerate([0.5, 1.0, 1.5, 1.75, 2.5]):
    v_eta = dyson_entrainment(r_grid, rho_prof, v_grav_mc, r_seed, eta)
    ls = '-' if eta == 1.75 else '--'
    lw = 2.5 if eta == 1.75 else 1.2
    alpha = 1.0 if eta == 1.75 else 0.6
    ax9.loglog(r_grid, v_eta, ls, linewidth=lw, color=colors[j], alpha=alpha,
              label=f'$\\eta={{{eta}}}$' + (' (best)' if eta == 1.75 else ''))

# Pure Newtonian gravitational velocity (stellar mass only + M_c)
ax9.loglog(r_grid, v_grav, 'k:', linewidth=1, alpha=0.5,
          label='$v_{\\rm grav}$ (Newtonian, $M_c$=' + f'{M_central/1e4:.0f}' + r'$\times10^4 M_\odot$)')

ax9.axhline(1274, color='black', linestyle='-.', alpha=0.5, linewidth=1.5,
            label='Target: $v=1274$ km/s (FWHM=3000)')
ax9.axvline(r_seed, color='gray', linestyle=':', alpha=0.3, label=f'Seed: $r={r_seed:.4f}$ pc')
ax9.axvline(0.03, color='purple', linestyle=':', alpha=0.3, label='BLR: $r=0.03$ pc')
ax9.axvline(r_b, color='gray', linestyle='--', alpha=0.3, label=f'$r_b={r_b:.1f}$ pc')

ax9.set_xlabel('r [pc]', fontsize=14)
ax9.set_ylabel('v [km s$^{-1}$]', fontsize=14)
ax9.set_title('Dyson Entrainment Velocity Cascade (Newtonian Gravity)', fontsize=15)
ax9.legend(fontsize=8, loc='lower left', ncol=2)
ax9.grid(True, alpha=0.3, which='both')
ax9.set_xlim(3e-5, 150)

# Inset: zoom on BLR region
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
ax_inset = inset_axes(ax9, width="35%", height="35%", loc='upper right',
                       bbox_to_anchor=(0, 0, 1, 1), bbox_transform=ax9.transAxes)
ax_inset.semilogx(r_grid, v_ent, 'r-', linewidth=2)
ax_inset.semilogx(r_grid, v_grav_mc, 'k:', linewidth=1, alpha=0.5)
ax_inset.axhline(1274, color='black', linestyle='-.', alpha=0.5, linewidth=1)
ax_inset.axvline(0.03, color='purple', linestyle=':', alpha=0.5)
ax_inset.set_xlim(0.005, 0.3)
ax_inset.set_ylim(500, 2500)
ax_inset.set_xlabel('r [pc]', fontsize=8)
ax_inset.set_ylabel('v [km/s]', fontsize=8)
ax_inset.tick_params(labelsize=7)
ax_inset.grid(True, alpha=0.3)
ax_inset.set_title('BLR Zoom', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'fig9_velocity_cascade.png'), dpi=200)
plt.close()
print("  -> fig9_velocity_cascade.png")

# ===========================
# Fig 10: Collision Time
# ===========================
print("\n--- Fig 10: Collision Time ---")

# Compute collision time profile with entrainment velocity
t_coll_profile = np.zeros_like(r_grid)
n_star_profile = rho_prof / m_star_avg

for i in range(len(r_grid)):
    if n_star_profile[i] < 1 or r_grid[i] < r_c * 0.1:
        t_coll_profile[i] = 1e15
        continue
    sigma_coll = (np.pi * R_star_avg**2 *
                  (1 + 2*G_N*(m_star_avg + m_star_avg) /
                   (R_star_avg * max(v_ent[i]**2, 0.01))))
    t_coll_yr = 1.0 / (n_star_profile[i] * sigma_coll * v_ent[i] * 1e5 / PC * YR)
    t_coll_profile[i] = t_coll_yr / 1e6  # Myr

fig10, ax10 = plt.subplots(figsize=(10, 7))

ax10.loglog(r_grid, t_coll_profile, 'b-', linewidth=2, label='Collision time $t_{\\rm coll}$')
ax10.axhline(1e3, color='red', linestyle='--', alpha=0.5, linewidth=1.5,
            label='Hubble time @ $z\\sim5.7$ ($\\sim$1 Gyr)')
ax10.axvline(r_c, color='gray', linestyle='--', alpha=0.3, label=f'$r_c$={r_c:.4f} pc')
ax10.axvline(r_eff, color='red', linestyle='--', alpha=0.3, label=f'$r_{{\\rm eff}}$={r_eff:.0f} pc')
ax10.axvline(0.03, color='purple', linestyle=':', alpha=0.3, label='BLR ($r$=0.03 pc)')

# Shade frequent-collision zone
ax10.fill_between([1e-5, 0.005], 1e-2, 1e6, alpha=0.08, color='green')
ax10.text(0.0003, 30, 'Frequent\ncollisions', fontsize=9, color='green', ha='center')

ax10.set_xlabel('r [pc]', fontsize=14)
ax10.set_ylabel('$t_{\\rm coll}$ [Myr]', fontsize=14)
ax10.set_title('Stellar Collision Timescale (Dyson Entrainment Model)', fontsize=15)
ax10.legend(fontsize=10, loc='lower right')
ax10.grid(True, alpha=0.3, which='both')
ax10.set_xlim(8e-5, 120)
ax10.set_ylim(0.5, 1e9)

# Annotate key values
key_rs = [0.001, 0.003, 0.01, 0.03, 0.1]
for kr in key_rs:
    idx = np.argmin(np.abs(r_grid - kr))
    t_val = t_coll_profile[idx]
    if t_val < 1e10:
        ax10.annotate(f'{t_val:.0f} Myr', (kr, t_val),
                     textcoords='offset points', xytext=(0, 10),
                     fontsize=8, ha='center', color='#185FA5')

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'fig10_collision_time.png'), dpi=200)
plt.close()
print("  -> fig10_collision_time.png")

# ===========================
# Summary
# ===========================
print(f"\n{'='*60}")
print("Dyson model figures generated.")
print(f"  rho_c = {rho_prof[0]:.2e} Msun/pc^3")
print(f"  v_seed = {v_seed_val:.0f} km/s @ r={r_seed:.4f} pc")
print(f"  v_BLR = {v_ent[idx_blr]:.0f} km/s, FWHM = {fwhm_blr:.0f} km/s")
print(f"  M_central = {M_central:.1e} Msun (1/{(5e6/M_central):.0f}x Kokorev virial)")
print(f"{'='*60}")
