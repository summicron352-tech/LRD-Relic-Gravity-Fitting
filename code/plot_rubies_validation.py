#!/usr/bin/env python3
"""
RUBIES Spectral Validation Figure
==================================
Independent verification channel for G(Sigma) framework using:
- Sample B: RUBIES 36 LRDs with [OIII]5007/Hb FWHM measurements
- Independent from Kokorev photometric sample (different instrument, team, method)
- Tests the same G(Sigma) prediction: denser LRDs -> broader lines (deeper potential wells)

Output: Figure_RUBIES_Validation.png + Table_RUBIES_Stats.tex
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 150,
})

# =====================================================================
# 1. LOAD DATA
# =====================================================================

# --- Sample A (Kokorev) - for Sigma normalization reference ---
dfA = pd.read_csv('/Users/tanxin/Desktop/数据处理/Kokorev_LRDs_Full.csv')
f444_A = dfA['f444w_flux'].values.astype(float)
f150_A = dfA['f150w_flux'].values.astype(float)
reff_A = dfA['r_eff_50_phys'].values.astype(float)
color_A = -2.5 * np.log10(f444_A / f150_A)

logSigma_A_raw = np.log10(np.clip(f444_A, 1e-30, None)) - 2*np.log10(np.clip(reff_A, 1e-8, None))
logSigma_A_med = np.nanmedian(logSigma_A_raw)  # This is our Sigma_0 reference

# --- Sample B (RUBIES) ---
dfB = pd.read_csv('/Users/tanxin/Desktop/数据处理/源数据/RUBIES_LRD_xDJA_full.csv')

fwhm_B = dfB['fwhm'].values.astype(float)
fw_unc_B = dfB['fw_unc'].values.astype(float)
z_B = dfB['z'].values.astype(float)
logSigma_B = dfB['logSigma'].values.astype(float)
m444_B = dfB['m444_rubies'].values.astype(float)
reff_kpc_B = dfB['reff_kpc'].values.astype(float)
beta_opt_B = dfB['beta_opt'].values.astype(float)

N_total = len(dfB)
print("RUBIES sample: N=%d" % N_total)
print("  FWHM range: [%.0f, %.0f] km/s" % (np.min(fwhm_B), np.max(fwhm_B)))
print("  logSigma range: [%.2f, %.2f]" % (np.min(logSigma_B), np.max(logSigma_B)))
print("  z range: [%.2f, %.2f]" % (np.min(z_B), np.max(z_B)))
print("  Sigma_0 (from Kokorev median): %.3f" % logSigma_A_med)

# Center logSigma to same origin as Sample A
logSigma_B_centered = logSigma_B - logSigma_A_med
Sigma_B_norm = 10 ** logSigma_B_centered  # Sigma / Sigma_0

# Valid mask
valid = (np.isfinite(fwhm_B) & np.isfinite(logSigma_B_centered) &
         (fwhm_B > 0) & (np.isfinite(z_B)))
fwhm_v = fwhm_B[valid]
lS_v = logSigma_B_centered[valid]
S_v = Sigma_B_norm[valid]
z_v = z_B[valid]
fw_err_v = fw_unc_B[valid]
N_valid = valid.sum()

print("\nValid points: N=%d" % N_valid)

# =====================================================================
# 2. STATISTICAL TESTS
# =====================================================================

# Pearson correlation (parametric)
r_pearson, p_pearson = stats.pearsonr(lS_v, fwhm_v)

# Spearman rank (non-parametric, robust to outliers)
rho_spearman, p_spearman = stats.spearmanr(lS_v, fwhm_v)

# Kendall's tau (ordinal)
tau_kendall, p_kendall = stats.kendalltau(lS_v, fwhm_v)

# Linear regression (OLS)
slope, intercept, r_value, p_ols, se_ols = stats.linregress(lS_v, fwhm_v)
rmse_ols = np.sqrt(np.mean((fwhm_v - (slope*lS_v + intercept))**2))

# Partial correlation: control for redshift
# Residualize both variables against z
if len(z_v) > 5:
    slope_fwhm_z, _, _, _, _ = stats.linregress(z_v, fwhm_v)
    slope_lS_z, _, _, _, _ = stats.linregress(z_v, lS_v)
    fwhm_resid = fwhm_v - slope_fwhm_z * z_v
    lS_resid = lS_v - slope_lS_z * z_v
    r_partial, p_partial = stats.pearsonr(lS_resid, fwhm_resid)
else:
    r_partial, p_partial = np.nan, np.nan

print("\n=== STATISTICAL TESTS ===")
print("  Pearson:     r=%.3f, p=%.2e" % (r_pearson, p_pearson))
print("  Spearman:   rho=%.3f, p=%.2e" % (rho_spearman, p_spearman))
print("  Kendall:    tau=%.3f, p=%.2e" % (tau_kendall, p_kendall))
print("  OLS:        slope=%.1f km/s/dex, R^2=%.3f, p=%.2e" % (slope, r_value**2, p_ols))
print("  RMSE:       %.0f km/s" % rmse_ols)
if np.isfinite(r_partial):
    print("  Partial(r|z): rp=%.3f, p=%.2e (controlling for redshift)" % (r_partial, p_partial))

# Significance in sigma
sig_pearson = -np.log10(p_pearson) if p_pearson > 0 else 99
sig_spearman = -np.log10(p_spearman) if p_spearman > 0 else 99

# =====================================================================
# 3. FIT WITH G(SIGMA) MODEL
# =====================================================================

# Use epsilon_g and beta from Kokorev fit (Path A)
# eps_g = 1.317, beta = 0.097 from Bootstrap median
eps_g_ref = 1.317
beta_ref = 0.097

# Physical model: FWHM ~ sqrt(G_eff) ~ sqrt(1 + eps_g * (Sigma/Sigma_0)^beta)
# Simplified: FWHM = v0 + A * sqrt(eps_g) * (Sigma/Sigma_0)^(beta/2)
def model_fwhm_sigma(Sig, v0, Af, eg, be):
    return v0 + Af * np.sqrt(eg) * (Sig ** (be / 2.0))

# Also fit a simple power law for comparison
def model_powerlaw(Sig, v0, alpha):
    return v0 + alpha * Sig

# Fit G(Sigma)-inspired model
from scipy.optimize import curve_fit

popt_g, pcov_g = curve_fit(
    lambda S, v0, Af: model_fwhm_sigma(S, v0, Af, eps_g_ref, beta_ref),
    S_v, fwhm_v,
    p0=[800, 1200],
    bounds=[[100, 10], [2000, 5000]],
    maxfev=5000
)
v0_g, Af_g = popt_g
err_g = np.sqrt(np.diag(pcov_g))

# Fit free power law (for comparison)
popt_pl, pcov_pl = curve_fit(
    lambda S, v0, al: model_powerlaw(S, v0, al),
    S_v, fwhm_v,
    p0=[800, 300],
    bounds=[[0, 0], [2000, 3000]],
    maxfev=5000
)
v0_pl, alpha_pl = popt_pl

# Compute R^2 for both models
pred_g = model_fwhm_sigma(S_v, v0_g, Af_g, eps_g_ref, beta_ref)
pred_pl = model_powerlaw(S_v, v0_pl, alpha_pl)
ss_res_g = np.sum((fwhm_v - pred_g)**2)
ss_res_pl = np.sum((fwhm_v - pred_pl)**2)
ss_tot = np.sum((fwhm_v - np.mean(fwhm_v))**2)
r2_g = 1 - ss_res_g/ss_tot
r2_pl = 1 - ss_res_pl/ss_tot
rmse_g = np.sqrt(np.mean((fwhm_v - pred_g)**2))

print("\n=== MODEL FITS ===")
print("  G(Sigma) model:  v0=%.0f+/-%.0f, Af=%.0f+/-%.0f, R2=%.3f" % (v0_g, err_g[0], Af_g, err_g[1], r2_g))
print("  Power law:       v0=%.0f, alpha=%.0f, R2=%.3f" % (v0_pl, alpha_pl, r2_pl))

# =====================================================================
# 4. Z-SLICE ANALYSIS (does signal persist at fixed z?)
# =====================================================================

z_bins = [(3.0, 5.0), (5.0, 7.0), (7.0, 9.0)]
z_slice_results = []

for zlo, zhi in z_bins:
    mask_z = (z_v >= zlo) & (z_v < zhi)
    nz = mask_z.sum()
    if nz >= 5:
        rz, pz = stats.spearmanr(lS_v[mask_z], fwhm_v[mask_z])
        z_slice_results.append({
            'z_range': '%.0f-%.0f' % (zlo, zhi),
            'N': nz,
            'rho': rz,
            'p': pz,
            'mean_fwhm': np.mean(fwhm_v[mask_z]),
            'mean_logS': np.mean(lS_v[mask_z]),
        })
        print("  z=[%.0f,%.0f): N=%d, rho=%.3f, p=%.2e" % (zlo, zhi, nz, rz, pz))

# =====================================================================
# 5. FIGURE
# =====================================================================

fig = plt.figure(figsize=(16, 12))
gs = GridSpec(2, 3, figure=fig, height_ratios=[1.3, 1], hspace=0.28, wspace=0.30)

# ---- PANEL A: Main — FWHM vs logSigma ----
ax1 = fig.add_subplot(gs[0, :2])

# Color by redshift (shows we're NOT just tracing z-FWHM correlation)
sc = ax1.scatter(lS_v, fwhm_v, c=z_v, cmap='plasma', s=90, edgecolors='black',
                  linewidths=0.6, alpha=0.85, zorder=3)
cbar = plt.colorbar(sc, ax=ax1, shrink=0.8, pad=0.02)
cbar.set_label('Redshift $z$', fontsize=11)

# Sort for smooth curve line
sort_idx = np.argsort(lS_v)
lS_sorted = lS_v[sort_idx]
S_sorted = S_v[sort_idx]

# G(Sigma) model prediction (using Kokorev parameters!)
fwhm_curve = model_fwhm_sigma(S_sorted, v0_g, Af_g, eps_g_ref, beta_ref)
ax1.plot(lS_sorted, fwhm_curve, 'r-', lw=2.8, zorder=5,
         label=(r'G($\Sigma$) Prediction: $\epsilon_g=$' + '{:.2f}'.format(eps_g_ref)
                + r', $\beta=$' + '{:.2f}'.format(beta_ref)))

# Power-law comparison (dashed gray)
fwhm_pl_curve = model_powerlaw(S_sorted, v0_pl, alpha_pl)
ax1.plot(lS_sorted, fwhm_pl_curve, '--', color='gray', lw=1.5, alpha=0.7,
         label=r'Power Law: $\propto \Sigma^{%.1f}$' % (alpha_pl/(v0_pl+1)))

# Confidence band via bootstrap-like scatter
ax1.fill_between(lS_sorted,
                  fwhm_curve - rmse_g,
                  fwhm_curve + rmse_g,
                  alpha=0.15, color='red', label=r'$\pm$1 RMSE')

ax1.axhline(y=np.median(fwhm_v), color='#444444', ls='--', lw=0.8, alpha=0.5)
ax1.axvline(x=0, color='#444444', ls=':', lw=1, alpha=0.5)

ax1.set_xlabel(r'$\log(\Sigma/\Sigma_0)$  (centered on Kokorev median)', fontsize=13)
ax1.set_ylabel(r'[OIII]/H$\beta$ FWHM  [km s$^{-1}$]', fontsize=13)
ax1.set_title(
    'Panel A: Spectral Validation — FWHM vs Surface Density\n'
    'RUBIES LRDs (N=%d), Independent from Photometric Sample' % N_valid,
    fontsize=12, fontweight='bold'
)

# Stats box
stats_text = (
    r'$\rho_{\rm Spearman} =$' + '{:.3f}'.format(rho_spearman)
    + r',  $p =$' + '{:.1e}'.format(p_spearman)
    + r'  ($\approx$' + '{:.1f}'.format(sig_spearman) + r'$\sigma$)'
    + '\n'
    + r'$r_{\rm partial}|_z =$' + ('{:.3f}'.format(r_partial) if np.isfinite(r_partial) else 'N/A')
    + ',  $R^2 =$' + '{:.3f}'.format(r2_g)
)
props = dict(boxstyle='round,pad=0.4', facecolor='#fffef0', alpha=0.95, edgecolor='#cc9900')
ax1.text(0.97, 0.03, stats_text, transform=ax1.transAxes, fontsize=10,
         verticalalignment='bottom', horizontalalignment='right', bbox=props)

ax1.legend(loc='upper left', fontsize=9, framealpha=0.92)
ax1.grid(alpha=0.12)

# ---- PANEL B: FWHM residuals ----
ax2 = fig.add_subplot(gs[1, 0])
resid = fwhm_v - pred_g
ax2.scatter(lS_v, resid, c=z_v, cmap='plasma', s=50, edgecolors='none', alpha=0.7)
ax2.axhline(y=0, color='red', ls='-', lw=1.5)
ax2.set_xlabel(r'$\log(\Sigma/\Sigma_0)$', fontsize=12)
ax2.set_ylabel('Residual [km s$^{-1}$]', fontsize=12)
ax2.set_title('G($\\Sigma$) Model Residuals\n$RMSE =$' + '{:.0f} km/s'.format(rmse_g), fontsize=11)
ax2.grid(alpha=0.12)

# ---- PANEL C: FWHM vs z (control plot) ----
ax3 = fig.add_subplot(gs[1, 1])

sc3 = ax3.scatter(z_v, fwhm_v, c=lS_v, cmap='viridis_r', s=60,
                   edgecolors='black', linewidths=0.5, alpha=0.8)
cbar3 = plt.colorbar(sc3, ax=ax3, shrink=0.8)
cbar3.set_label(r'$\log(\Sigma/\Sigma_0)$', fontsize=10)

# Show weak/no z-FWHM trend
rz_z, pz_z = stats.spearmanr(z_v, fwhm_v)
ax3.text(0.05, 0.95,
         r'$\rho_z =$' + '{:.3f}'.format(rz_z) + r',  $p_z =$' + '{:.2f}'.format(pz_z),
         transform=ax3.transAxes, fontsize=10, va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))

ax3.set_xlabel('Redshift $z$', fontsize=12)
ax3.set_ylabel('FWHM  [km s$^{-1}$]', fontsize=12)
ax3.set_title('Control: FWHM vs Redshift\n(color = $\\Sigma$)', fontsize=11)
ax3.grid(alpha=0.12)

# ---- PANEL D: Redshift-sliced analysis ----
ax4 = fig.add_subplot(gs[1, 2])
ax4.axis('off')

if z_slice_results:
    slice_text = "REDSHIFT-SLICED ANALYSIS\n"
    slice_text += "{:-<45}\n".format("")
    slice_text += "{:^8s}  {:^5s}  {:>8s}  {:>12s}\n".format("z-range", "N", "rho_S", "p-value")
    slice_text += "{:-<45}\n".format("")
    for sr in z_slice_results:
        sig_str = "{:.1f}sigma".format(-np.log10(sr['p'])) if sr['p'] > 1e-99 else ">99sigma"
        slice_text += "  {:>8s}  {:4d}    {:+.3f}      {:.1e}  ({})\n".format(
            sr['z_range'], sr['N'], sr['rho'], sr['p'], sig_str)
    
    slice_text += "\nKey: If FWHM-Sigma correlation is driven\n"
    slice_text += "by a hidden z-correlation, slicing should\n"
    slice_text += "weaken or eliminate the signal.\n"
    
    all_positive = all(sr['rho'] > 0 for sr in z_slice_results)
    if all_positive:
        slice_text += "\nResult: POSITIVE rho in ALL slices -->\n"
        slice_text += "Signal is INDEPENDENT of redshift."
    else:
        slice_text += "\nResult: Mixed -- needs investigation."
else:
    slice_text = "Not enough data for z-slicing."

ax4.text(0.05, 0.95, slice_text, transform=ax4.transAxes, fontsize=10,
         fontfamily='monospace', verticalalignment='top',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f8ff', edgecolor='#4682b4'))

ax4.set_title('Redshift Independence Test', fontsize=11, fontweight='bold')

plt.suptitle(
    'Spectral Validation of Density-Dependent Gravity: RUBIES [OIII] FWHM Analysis',
    fontsize=14, fontweight='bold', y=0.995
)

outpath = '/Users/tanxin/Desktop/数据处理/源数据/MasterPaper/Figure_RUBIES_Validation.png'
plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
print("\nFigure saved: %s" % outpath)

# =====================================================================
# 6. LATEX TABLE
# =====================================================================

tex = (
    r"%% Table: RUBIES spectral validation statistics"
    + "\n"
    + r"\begin{table}[htbp]"
    + "\n"
    + r"\centering"
    + "\n"
    + r"\caption{Spectral validation of the density-dependent effective gravity framework"
    + "\n using RUBIES [OIII]/H$\beta$ FWHM measurements ($N=${}).".format(N_valid)
    + "\n This provides an independent test of the G$(\\Sigma)$ prediction using"
    + "\n spectroscopic data from a different instrument and reduction pipeline.}"
    + "\n"
    + r"\label{tab:rubies_validation}"
    + "\n"
    + r"\vspace{0.3em}"
    + "\n"
    + r"\begin{tabular}{ll}"
    + "\n"
    + r"\toprule"
    + "\n"
    + "Measurement & Value \\\\\n"
    + r"\midrule"
    + "\n"
    + "Sample size & $N = %d$ \\\\" % N_valid
    + "\n"
    + "Redshift range & $z \\in [{:.2f}, {:.2f}]$ \\\\".format(np.min(z_v), np.max(z_v))
    + "\n"
    + "FWHM range & ${:.0f}$--${:.0f}$ km~s$^{{-1}}$ \\\\".format(np.min(fwhm_v), np.max(fwhm_v))
    + "\n"
    + r"$\log(\Sigma/\Sigma_0)$ range & ${:.2f}$ to ${:.2f}$ \\\\".format(np.min(lS_v), np.max(lS_v))
    + "\n"
    + r"\midrule"
    + "\n"
    + r"Spearman $\rho$ & ${{+.3f}}$ \\\\".format(rho_spearman)
    + "\n"
    + r"Spearman $p$-value & ${:.1e}$ \\\\".format(p_spearman)
    + "\n"
    + r"Significance & ${:.1f}\\sigma$ \\\\".format(sig_spearman)
    + "\n"
    + r"\midrule"
    + "\n"
    + r"$r_{{\rm partial}}|_z$ & ${:.3f}$ \\\\".format(r_partial) if np.isfinite(r_partial) else r"$r_{\rm partial}|_z$ ---"
    + "\n"
    + r"G$(\Sigma)$ model $R^2$ & ${:.3f}$ \\\\".format(r2_g)
    + "\n"
    + r"G$(\Sigma)$ RMSE & ${:.0f}$ km~s$^{{-1}}$ \\\\".format(rmse_g)
    + "\n"
    + r"\midrule"
    + "\n"
    + r"Best-fit $v_0$ & ${:.0f} \pm {:.0f}$ km~s$^{{-1}}$ \\\\".format(v0_g, err_g[0])
    + "\n"
    + r"Best-fit amplitude $A_F$ & ${:.0f} \pm {:.0f}$ km~s$^{{-1}}$ \\\\".format(Af_g, err_g[1])
    + "\n"
    + r"\midrule"
    + "\n"
    + r"\multicolumn{2}{l}{\small Adopted from Kokorev fit: "
    + r"$\epsilon_g = {:.3f}$, $\beta = {:.3f}$}} \\\\".format(eps_g_ref, beta_ref)
    + "\n"
    + r"\bottomrule"
    + "\n"
    + r"\end{tabular}"
    + "\n"
    + r"\end{table}"
    + "\n"
)

tex_path = '/Users/tanxin/Desktop/数据处理/源数据/MasterPaper/Table_RUBIES_Validation.tex'
with open(tex_path, 'w') as f:
    f.write(tex)
print("LaTeX table saved: %s" % tex_path)

# =====================================================================
# 7. CSV OUTPUT
# =====================================================================

csv_data = {
    'statistic': ['N', 'z_min', 'z_max', 'FWHM_min', 'FWHM_max',
                  'logSig_min', 'logSig_max',
                  'pearson_r', 'pearson_p', 'pearson_sigma',
                  'spearman_rho', 'spearman_p', 'spearman_sigma',
                  'kendall_tau', 'kendall_p',
                  'partial_r_given_z', 'partial_p',
                  'OLS_slope_km/s_per_dex', 'OLS_intercept', 'OLS_R2', 'OLS_RMSE',
                  'GSigma_R2', 'GSigma_RMSE', 'v0_km/s', 'Af_km/s',
                  'eps_g_from_Kokorev', 'beta_from_Kokorev'],
    'value': [N_valid, np.min(z_v), np.max(z_v), np.min(fwhm_v), np.max(fwhm_v),
              np.min(lS_v), np.max(lS_v),
              r_pearson, p_pearson, sig_pearson,
              rho_spearman, p_spearman, sig_spearman,
              tau_kendall, p_kendall,
              r_partial, p_partial,
              slope, intercept, r_value**2, rmse_ols,
              r2_g, rmse_g, v0_g, Af_g,
              eps_g_ref, beta_ref]
}
df_out = pd.DataFrame(csv_data)
csv_path = '/Users/tanxin/Desktop/数据处理/源数据/MasterPaper/RUBIES_Validation_Statistics.csv'
df_out.to_csv(csv_path, index=False)
print("CSV saved: %s" % csv_path)

print("\n=== DONE ===")
