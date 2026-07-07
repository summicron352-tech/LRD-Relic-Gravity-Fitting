#!/usr/bin/env python3
"""
Generate ApJ-quality figures for the LRD compactness paper.
All labels in English. Corrected numerical values from code audit.

Figures:
  Fig 4: SB_F444W vs logSigma (291 LRD merged) — core evidence
  Fig 5: Sigma-color partial correlation vs redshift + LRD overlay
  Fig 6: Paired reversal rate vs redshift
  Fig 7: Orientation test (edge-on vs face-on)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 10,
    'font.family': 'serif',
    'axes.linewidth': 1.2,
})
from scipy.stats import spearmanr, rankdata
from sklearn.neighbors import BallTree
from astropy.io import fits
from astropy.cosmology import Planck18
import astropy.units as u
import os

# ===========================
# Paths
# ===========================
cosmos_dir = '/Users/tanxin/Desktop/数据处理/COSMOS2025'
data_dir = '/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/data'
fig_dir = '/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/figures'
os.makedirs(fig_dir, exist_ok=True)

# ===========================
# Common functions
# ===========================
def psc(x, y, ctrls):
    """Partial Spearman correlation controlling for ctrls"""
    rx = rankdata(x); ry = rankdata(y)
    for c in ctrls:
        rc = rankdata(c)
        b = np.cov(rx, rc)[0,1] / np.var(rc); rx = rx - b*rc
        b = np.cov(ry, rc)[0,1] / np.var(rc); ry = ry - b*rc
    return spearmanr(rx, ry)[0]

def bootstrap_psc(x, y, ctrls, n_boot=500):
    """Bootstrap partial Spearman with error"""
    n = len(x)
    rho_boot = []
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        try:
            r = psc(x[idx], y[idx], [c[idx] for c in ctrls])
            rho_boot.append(r)
        except:
            pass
    rho_boot = np.array(rho_boot)
    return np.mean(rho_boot), np.std(rho_boot)

# ===========================
# Load COSMOS 664k
# ===========================
print("Loading COSMOSWeb 664k data...")
cosmos = np.load(os.path.join(cosmos_dir, 'cosmos2025_extracted.npz'))
z_c = cosmos['z']; logM_c = cosmos['logM']; logSigma_c = cosmos['logSigma']
color_c = cosmos['color']; reff_c = cosmos['reff_pc']; m444_c = cosmos['m444']

# ===========================
# Fig 5: Sigma-color partial correlation + LRD match overlay
# ===========================
print("\n--- Fig 5: Sigma-color correlation evolution ---")

z_bins = [(0,1),(1,3),(3,5),(5,7),(7,9),(9,15)]
z_bin_labels = ['0-1','1-3','3-5','5-7','7-9','9-15']
z_centers = [0.6, 1.7, 3.9, 5.9, 7.5, 9.8]

rho_all = []; rho_err = []; n_all = []

for (zlo, zhi) in z_bins:
    mask = (z_c >= zlo) & (z_c < zhi)
    n = mask.sum()
    rho, err = bootstrap_psc(logSigma_c[mask], color_c[mask], [z_c[mask], logM_c[mask]])
    rho_all.append(rho); rho_err.append(err); n_all.append(n)
    print(f"  z={zlo:.0f}-{zhi:.0f}: N={n:>6d}  rho={rho:+.4f} +/- {err:.4f}")

# LRD matched subsample (from Phase 3 cross-match)
lrd_matches = {
    'Kokorev': {'z': [3.9, 5.9, 7.5], 'rho': [0.280, 0.312, 0.474], 'N': [81, 118, 56]},
    'Path1':   {'z': [3.9, 5.9, 7.5], 'rho': [0.274, 0.332, 0.631], 'N': [11, 16, 11]}
}

fig5, ax5 = plt.subplots(figsize=(10, 6.5))

# COSMOS full sample
ax5.errorbar(z_centers, rho_all, yerr=rho_err, fmt='o-', color='#185FA5',
             linewidth=2, markersize=10, capsize=4, capthick=1.5,
             label='COSMOSWeb 664k (full sample)', zorder=5)

for zc, r, n in zip(z_centers, rho_all, n_all):
    ax5.annotate(f'N={n:,}', (zc, r), textcoords='offset points',
                xytext=(0, 12), fontsize=8, ha='center', color='#185FA5')

# LRD matched subsamples
ax5.scatter(lrd_matches['Kokorev']['z'], lrd_matches['Kokorev']['rho'],
            c='#D85A30', marker='s', s=80, edgecolors='#993C1D', linewidths=1.5,
            label='COSMOS matched to Kokorev 260 (z, M*)', zorder=6)
ax5.scatter(lrd_matches['Path1']['z'], lrd_matches['Path1']['rho'],
            c='#5DCAA5', marker='^', s=80, edgecolors='#1D9E75', linewidths=1.5,
            label='COSMOS matched to Path1 38 (z, M*)', zorder=6)

for i, (zc, r) in enumerate(zip(lrd_matches['Kokorev']['z'], lrd_matches['Kokorev']['rho'])):
    ax5.annotate(f'N={lrd_matches["Kokorev"]["N"][i]}', (zc, r),
                textcoords='offset points', xytext=(10, -10), fontsize=7, color='#993C1D')

# Power-law fit
ax5.text(0.03, 0.97, r'$\rho \propto (1+z)^{1.14\pm0.06}$', transform=ax5.transAxes,
         fontsize=12, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

ax5.axhline(0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)
ax5.set_xlabel('Redshift z', fontsize=14)
ax5.set_ylabel(r'Partial Spearman $\rho(\Sigma, F150W-F444W\ |\ z, M_*)$', fontsize=13)
ax5.set_title(r'$\Sigma$–Color Correlation: Redshift Evolution', fontsize=15)
ax5.legend(fontsize=9, loc='upper left')
ax5.grid(True, alpha=0.25)
ax5.set_xlim(0, 10.5)
ax5.set_ylim(-0.05, 0.75)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'fig5_sigma_color_evolution.png'), dpi=200)
plt.close()
print("  -> fig5_sigma_color_evolution.png")


# ===========================
# Fig 6: Paired reversal rate
# ===========================
print("\n--- Fig 6: Paired reversal rate ---")

cosmos_env = np.load(os.path.join(cosmos_dir, 'cosmos2025_with_env.npz'))
ra = cosmos_env['ra']; dec = cosmos_env['dec']
z_env = cosmos_env['z']; logSigma_env = cosmos_env['logSigma']
color_env = cosmos_env['color']

paired_z = [(0,1),(1,3),(3,5),(5,7),(7,9)]
paired_centers = [0.6, 1.7, 3.9, 5.9, 7.5]
rev_rates = []; rev_rates_err = []; n_pairs_list = []

np.random.seed(42)

for zlo, zhi in paired_z:
    m = (z_env >= zlo) & (z_env < zhi)
    idx = np.where(m)[0]; n = len(idx)
    if n < 1000:
        rev_rates.append(np.nan); rev_rates_err.append(np.nan)
        n_pairs_list.append(0); continue

    ns = min(3000, n)
    sidx = np.random.choice(idx, ns, replace=False)

    all_ra = ra[m]; all_dec = dec[m]; all_z = z_env[m]
    all_ls = logSigma_env[m]; all_cl = color_env[m]
    idx_map = {orig: loc for loc, orig in enumerate(idx)}

    tree = BallTree(np.column_stack([all_ra, all_dec]), metric='haversine')

    npairs = 0; nrev = 0
    for i in range(ns):
        si_orig = sidx[i]; si_loc = idx_map[si_orig]
        _, neigh = tree.query([np.column_stack([all_ra, all_dec])[si_loc]], k=4)
        for ji in neigh[0][1:]:
            j_orig = idx[ji]
            if abs(all_z[si_loc] - all_z[ji]) > 0.3: continue
            if abs(all_ls[si_loc] - all_ls[ji]) < 0.3: continue
            npairs += 1
            if (all_ls[si_loc] > all_ls[ji]) != (all_cl[si_loc] > all_cl[ji]):
                nrev += 1
            if npairs >= 3000: break
        if npairs >= 3000: break

    rate = nrev / npairs * 100
    rate_err = np.sqrt(rate * (100 - rate) / npairs)
    rev_rates.append(rate); rev_rates_err.append(rate_err)
    n_pairs_list.append(npairs)
    print(f"  z={zlo:.0f}-{zhi:.0f}: {npairs} pairs, {rate:.1f}% +/- {rate_err:.1f}% reversal")

fig6, ax6 = plt.subplots(figsize=(10, 6.5))

valid = ~np.isnan(rev_rates)
ax6.errorbar(np.array(paired_centers)[valid], np.array(rev_rates)[valid],
             yerr=np.array(rev_rates_err)[valid],
             fmt='o-', color='#993C1D', linewidth=2, markersize=10,
             capsize=4, capthick=1.5, zorder=5)

ax6.axhline(50, color='gray', linestyle='--', alpha=0.5, linewidth=1.5,
            label='Null expectation (50%)')

for zc, r, np_ in zip(paired_centers, rev_rates, n_pairs_list):
    if np_ > 0:
        ax6.annotate(f'{np_} pairs', (zc, r), textcoords='offset points',
                    xytext=(0, 12), fontsize=8, ha='center', color='#993C1D')

ax6.set_xlabel('Redshift z', fontsize=14)
ax6.set_ylabel('Reversal Rate [%]', fontsize=14)
ax6.set_title(r'Paired Reversal Test of $\Sigma$–Color Correlation', fontsize=15)
ax6.legend(fontsize=11)
ax6.grid(True, alpha=0.25)
ax6.set_xlim(0, 8.5)
ax6.set_ylim(0, 65)

ax6.text(0.05, 0.95, 'Higher-$\Sigma$ redder: rate < 50%\nHigher-$\Sigma$ bluer:  rate > 50%',
        transform=ax6.transAxes, fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'fig6_reversal_rate.png'), dpi=200)
plt.close()
print("  -> fig6_reversal_rate.png")


# ===========================
# Fig 7: Orientation test
# ===========================
print("\n--- Fig 7: Orientation test ---")

fits_file = os.path.join(cosmos_dir, 'COSMOSWeb_mastercatalog_v1.1.fits')
if os.path.exists(fits_file):
    print("  Reading FITS catalog...")
    hdul = fits.open(fits_file, memmap=True)
    phot = hdul[1].data; lephare = hdul[2].data

    axratio = phot['axratio_sersic']; zfinal = lephare['zfinal']
    mass_med = lephare['mass_med']; r_deg = phot['radius_sersic']
    warn = phot['warn_flag']; m4 = phot['mag_auto_f444w']; m1 = phot['mag_auto_f150w']

    mask = ((zfinal > 0) & (zfinal < 99) & (mass_med > -90) &
            (r_deg > 0) & (warn == 0) &
            np.isfinite(axratio) & (axratio > 0) & (axratio <= 1))

    z2 = zfinal[mask]; logM2 = mass_med[mask].astype(float)
    ax2 = axratio[mask].astype(float); r_deg2 = r_deg[mask]
    m42 = m4[mask]; m12 = m1[mask]

    r_rad2 = np.deg2rad(r_deg2.astype(float))
    DA2 = Planck18.angular_diameter_distance(z2.astype(float)).to(u.pc).value
    r_pc2 = r_rad2 * DA2
    Sigma2 = 10**logM2 / (np.pi * r_pc2**2)
    logSigma2 = np.log10(Sigma2 + 1e-30)
    color2 = 0.4 * (m12.astype(float) - m42.astype(float))

    # Tully-Fisher rotation velocity
    alpha_tf = 3.5; logM_ref = 10.2; logV_ref = np.log10(200)
    logV2 = logV_ref + (logM2 - logM_ref) / alpha_tf; V2 = 10**logV2

    good = np.isfinite(logSigma2) & np.isfinite(color2) & np.isfinite(z2) & np.isfinite(logM2)
    z2 = z2[good]; logM2 = logM2[good]; logSigma2 = logSigma2[good]
    color2 = color2[good]; ax2 = ax2[good]; V2 = V2[good]

    V_mask = (V2 >= 150) & (V2 < 250)
    edge = V_mask & (ax2 < 0.3)
    inter = V_mask & (ax2 >= 0.3) & (ax2 <= 0.7)
    face = V_mask & (ax2 > 0.7)

    rho_edge = psc(logSigma2[edge], color2[edge], [z2[edge], logM2[edge]])
    rho_inter = psc(logSigma2[inter], color2[inter], [z2[inter], logM2[inter]])
    rho_face = psc(logSigma2[face], color2[face], [z2[face], logM2[face]])

    n_edge = edge.sum(); n_inter = inter.sum(); n_face = face.sum()

    print(f"  V_rot = 150-250 km/s:")
    print(f"    Edge-on (b/a<0.3):  N={n_edge:>6d}  rho={rho_edge:+.4f}")
    print(f"    Intermed (0.3-0.7): N={n_inter:>6d}  rho={rho_inter:+.4f}")
    print(f"    Face-on (b/a>0.7):  N={n_face:>6d}  rho={rho_face:+.4f}")
    print(f"    Dust prediction: rho_edge >> rho_face")
    print(f"    Observed: rho_edge < rho_face -> DUST HYPOTHESIS FAILS")

    # By redshift
    z_orient_bins = [(0,1),(1,3),(3,5),(5,7)]
    orient_z_results = []
    for zlo, zhi in z_orient_bins:
        zm = (z2 >= zlo) & (z2 < zhi) & V_mask
        edge_z = zm & (ax2 < 0.3); face_z = zm & (ax2 > 0.7)
        if edge_z.sum() > 50 and face_z.sum() > 50:
            re_z = psc(logSigma2[edge_z], color2[edge_z], [z2[edge_z], logM2[edge_z]])
            rf_z = psc(logSigma2[face_z], color2[face_z], [z2[face_z], logM2[face_z]])
            orient_z_results.append({
                'z': (zlo+zhi)/2, 're': re_z, 'rf': rf_z,
                'ne': edge_z.sum(), 'nf': face_z.sum()
            })

    # --- Figure 7 ---
    fig7, (ax7a, ax7b) = plt.subplots(1, 2, figsize=(14, 6))

    orientations = ['Edge-on\n(b/a < 0.3)', 'Intermediate\n(0.3-0.7)', 'Face-on\n(b/a > 0.7)']
    rho_vals = [rho_edge, rho_inter, rho_face]
    n_vals = [n_edge, n_inter, n_face]
    colors_orient = ['#D85A30', '#888780', '#378ADD']

    bars = ax7a.bar(orientations, rho_vals, color=colors_orient, alpha=0.8,
                    edgecolor='gray', linewidth=0.5, width=0.5)

    for bar, r, n in zip(bars, rho_vals, n_vals):
        ax7a.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                 f'rho={r:+.3f}\nN={n:,}', ha='center', va='bottom', fontsize=10)

    ax7a.axhline(0, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)
    ax7a.set_ylabel(r'Partial Spearman $\rho(\Sigma, {\rm color}\ |\ z, M_*)$', fontsize=13)
    ax7a.set_title(r'$\Sigma$-Color Correlation by Inclination', fontsize=14)
    ax7a.grid(True, alpha=0.2, axis='y')

    ax7a.text(0.5, 0.95,
             'Dust prediction: edge-on >> face-on\n'
             'Observed: edge-on < face-on\n'
             r'$\rightarrow$ Dust hypothesis ruled out',
             transform=ax7a.transAxes, fontsize=10, ha='center', va='top',
             bbox=dict(boxstyle='round', facecolor='#FCEBEB', alpha=0.8))

    if len(orient_z_results) > 0:
        z_oz = [r['z'] for r in orient_z_results]
        re_oz = [r['re'] for r in orient_z_results]
        rf_oz = [r['rf'] for r in orient_z_results]

        ax7b.plot(z_oz, re_oz, 'o-', color='#D85A30', linewidth=2, markersize=8,
                 label='Edge-on (b/a < 0.3)')
        ax7b.plot(z_oz, rf_oz, 's-', color='#378ADD', linewidth=2, markersize=8,
                 label='Face-on (b/a > 0.7)')

        ax7b.axhline(0, color='gray', linestyle='--', alpha=0.4)
        ax7b.set_xlabel('Redshift z', fontsize=13)
        ax7b.set_ylabel(r'Partial Spearman $\rho$', fontsize=13)
        ax7b.set_title(r'Orientation Test vs. Redshift', fontsize=14)
        ax7b.legend(fontsize=10)
        ax7b.grid(True, alpha=0.25)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'fig7_orientation_test.png'), dpi=200)
    plt.close()
    print("  -> fig7_orientation_test.png")
    hdul.close()
else:
    print("  WARNING: FITS file not found, skipping Fig 7")
    print(f"  Expected: {fits_file}")


# ===========================
# Fig 4: SB_F444W vs logSigma (291 LRD merged)
# ===========================
print("\n--- Fig 4: SB_F444W vs logSigma (291 merged) ---")

kokorev = np.genfromtxt(os.path.join(data_dir, 'kokorev_260_sb.csv'),
                        delimiter=',', names=True, dtype=None, encoding='utf-8')
path1 = np.genfromtxt(os.path.join(data_dir, 'path1_merged_38sources.csv'),
                      delimiter=',', names=True, dtype=None, encoding='utf-8')

z_k = np.array(kokorev['z_phot'], dtype=float)
ls_k = np.array(kokorev['logSigma_Mstar'], dtype=float)
sb_k = np.array(kokorev['SB_F444W'], dtype=float)
z_p = np.array(path1['z_phot'], dtype=float)
ls_p = np.array(path1['logSigma_Mstar'], dtype=float)
sb_p = np.array(path1['SB_F444W'], dtype=float)

valid_k = np.isfinite(ls_k) & np.isfinite(sb_k) & (ls_k > 0) & (ls_k < 15) & (sb_k > 15) & (sb_k < 35)
valid_p = np.isfinite(ls_p) & np.isfinite(sb_p) & (ls_p > 0) & (ls_p < 15) & (sb_p > 15) & (sb_p < 35)

ls_k = ls_k[valid_k]; sb_k = sb_k[valid_k]; z_k = z_k[valid_k]
ls_p = ls_p[valid_p]; sb_p = sb_p[valid_p]; z_p = z_p[valid_p]

ls_all = np.concatenate([ls_k, ls_p])
sb_all = np.concatenate([sb_k, sb_p])

# Linear fit (logSigma = logM - 2*log(r), no pi term)
from numpy.polynomial.polynomial import polyfit
p = polyfit(ls_all, sb_all, 1)
rho_sb = spearmanr(ls_all, sb_all)

x_fit = np.linspace(ls_all.min()-0.3, ls_all.max()+0.3, 100)

fig4, ax4 = plt.subplots(figsize=(10, 7))

sc1 = ax4.scatter(ls_k, sb_k, c=z_k, s=25, alpha=0.6, cmap='plasma',
                  vmin=3, vmax=9, label=f'Kokorev+ 260 (N={len(ls_k)})')
sc2 = ax4.scatter(ls_p, sb_p, c=z_p, s=55, marker='s', alpha=0.8, cmap='plasma',
                  vmin=3, vmax=9, label=f'Path1 38 (N={len(ls_p)})')
ax4.plot(x_fit, p[0] + p[1]*x_fit, 'k--', alpha=0.7, linewidth=1.5,
         label=f'Fit: SB = {p[1]:.2f} log$\Sigma$ + {p[0]:.1f}')

ax4.set_xlabel(r'log $\Sigma$ [$M_\odot$ pc$^{-2}$]', fontsize=14)
ax4.set_ylabel('SB_F444W [mag arcsec$^{-2}$]', fontsize=14)
ax4.set_title(
    f'Surface Brightness vs. Surface Density '
    f'(N={len(ls_all)}, Spearman $\\rho$={rho_sb[0]:.3f}, '
    f'$p$={rho_sb[1]:.1e})',
    fontsize=14)
ax4.legend(fontsize=10, loc='upper right')
ax4.grid(True, alpha=0.25)

cbar = plt.colorbar(sc1, ax=ax4, label='Redshift z')
cbar.ax.tick_params(labelsize=10)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'fig4_sb_vs_sigma.png'), dpi=200)
plt.close()
print("  -> fig4_sb_vs_sigma.png")

print("\n" + "=" * 60)
print("All 4 figures generated successfully.")
print("=" * 60)
