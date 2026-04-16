#!/usr/bin/env python3
"""
Path A Publishable-Quality Figure: Density-Dependent Flux Ratios in LRDs
========================================================================
Target: MNRAS / ApJ level figure — single-column (183mm) or double-column

Panel layout (2×2):
  (a) log Σ vs F444W/F150W scatter + regression, color = z_phot
  (b) Partial correlation residual vs log Σ (after controlling z_phot + L_bol)
  (c) F444W/F150W distribution by Σ quartile (boxen/box plot)
  (d) Multi-band partial-ρ comparison bar chart (with error bars)

Style: Minimal, clean, colorblind-friendly, B&W safe
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import spearmanr, pearsonr, kendalltau, norm
from numpy.linalg import lstsq, pinv
from scipy.stats import ks_2samp, rankdata
import warnings
warnings.filterwarnings('ignore')

# ════════════════════════════════════════════════════════
# PUBLICATION SETTINGS
# ════════════════════════════════════════════════════════

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica Neue', 'Arial', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
    'lines.linewidth': 1.2,
    'lines.markersize': 2.5,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'axes.linewidth': 0.6,
    'text.usetex': False,
    'mathtext.fontset': 'dejavusans',
})

# Color palette — colorblind-safe, print-friendly
CB_PALETTE = {
    'blue':   '#0072B2',
    'orange': '#D55E00',
    'green':  '#009E73',
    'pink':   '#CC79A7',
    'gray':   '#999999',
    'yellow':'#F0E442',
    'sky':    '#56B4E9',
    'vermilion': '#E69F00',
}

QUARTILE_COLORS = [CB_PALETTE['sky'], CB_PALETTE['green'], 
                   CB_PALETTE['orange'], CB_PALETTE['blue']]

FIG_WIDTH_MM = 220   # Wider for more breathing room
FIG_HEIGHT_MM = 195   # Taller to prevent text overlap
DPI = 300

# ════════════════════════════════════════════════════════
# DATA LOADING
# ════════════════════════════════════════════════════════

DATA_PATH = '/Users/tanxin/Desktop/数据处理/LRD_GitHub_Release_v2/data/csv/Kokorev_LRDs_Full.csv'
OUTPUT_DIR = '/Users/tanxin/Desktop/数据处理/LRD_GitHub_Release_v2/figures/density_zdist_analysis'

df_raw = pd.read_csv(DATA_PATH)

# ════════════════════════════════════════════════════════
# DATA PREPARATION
# ════════════════════════════════════════════════════════

df = df_raw.copy()

# Basic cleaning: valid photometry in key bands
BANDS = ['f115w_flux', 'f150w_flux', 'f200w_flux', 
         'f277w_flux', 'f356w_flux', 'f444w_flux']
for b in BANDS:
    df = df[df[b] > 0]

# Surface density Sigma — use L_bol as mass proxy (L-M relation: M ~ L^{0.7-0.8})
# For LRDs, this is a reasonable approximation since they're dominated by nuclear emission
df['Mstar_proxy'] = df['lbol']  # L_bol in solar units
df['R_eff_pc'] = df['r_eff_50_phys'] * 1000  # kpc → pc
df['Sigma'] = df['Mstar_proxy'] / (2 * np.pi * df['R_eff_pc']**2)  
df['log_Sigma'] = np.log10(df['Sigma'])

# Compactness C = GM/Rc^2 (dimensionless) — using L_bol as M proxy
G_cgs = 6.674e-8
Msun_g = 1.989e33
pc_cm = 3.086e18
c_cgs = 2.998e10
df['C'] = G_cgs * df['Mstar_proxy'] * Msun_g / (df['R_eff_pc'] * pc_cm * c_cgs**2)
df['log_C'] = np.log10(df['C'])

# Bolometric luminosity (use lbol from catalog)
df['log_Lbol'] = np.log10(df['lbol'].clip(lower=0.01))  # avoid log(0)

# Flux ratios
df['FR_444_356'] = df['f444w_flux'] / df['f356w_flux']
df['FR_444_150'] = df['f444w_flux'] / df['f150w_flux']
df['FR_356_150'] = df['f356w_flux'] / df['f150w_flux']
df['FR_277_150'] = df['f277w_flux'] / df['f150w_flux']

# Log flux ratios
df['log_FR_444_356'] = np.log10(df['FR_444_356'])
df['log_FR_444_150'] = np.log10(df['FR_444_150'])

# Working copy — ensure all key columns are valid
dfv = df.dropna(subset=['log_Sigma', 'log_FR_444_150', 'log_FR_444_356',
                         'z_phot', 'log_Lbol', 'lbol']).copy()
# Also filter out non-positive values
dfv = dfv[dfv['log_Sigma'] < 15]  # remove extreme outliers
dfv = dfv[np.isfinite(dfv['log_Sigma'])]
print(f"Clean sample: {len(dfv)} sources")

# Quartile cuts on Sigma
dfv['Sigma_Q'] = pd.qcut(dfv['log_Sigma'], 4, labels=['Q1 (diffuse)', 'Q2', 'Q3', 'Q4 (compact)'])
q_labels = ['Q1\n(diffuse)', 'Q2', 'Q3', 'Q4\n(compact)']

# ════════════════════════════════════════════════════════
# STATISTICAL ANALYSIS FUNCTIONS
# ════════════════════════════════════════════════════════

def sig_from_p(p):
    """Convert p-value to Gaussian sigma."""
    if p <= 0:
        return 99.0
    return abs(norm.ppf(p / 2))

def partial_spearman(x, y, ctrl_cols, data):
    """Partial Spearman controlling for confounders."""
    rx = rankdata(x)
    ry = rankdata(y)
    X_ctrl = np.column_stack([np.ones(len(data))] + [data[c].values for c in ctrl_cols])
    try:
        pX = pinv(X_ctrl, rcond=1e-6)
        resid_x = rx - X_ctrl @ (pX @ rx)
        resid_y = ry - X_ctrl @ (pX @ ry)
        return spearmanr(resid_x, resid_y)
    except:
        return (np.nan, np.nan)

# ════════════════════════════════════════════════════════
# COMPUTE ALL STATISTICS
# ════════════════════════════════════════════════════════

# --- Panel A stats ---
rho_a, pa_a = spearmanr(dfv['log_Sigma'], dfv['log_FR_444_150'])
sig_a = sig_from_p(pa_a)

# Regression line
coef_a = lstsq(np.column_stack([np.ones(len(dfv)), dfv['log_Sigma'].values]),
                dfv['log_FR_444_150'].values, rcond=None)[0]

print(f"Panel A: rho={rho_a:.3f}, p={pa_a:.2e} ({sig_a:.1f}sigma)")

# --- Panel B stats (partial) ---
controls = ['z_phot', 'log_Lbol']
rp_b, pp_b = partial_spearman(
    dfv['log_Sigma'].values, dfv['log_FR_444_150'].values, controls, dfv)
sig_b = sig_from_p(pp_b)

# Compute residuals
X_ctrl = np.column_stack([np.ones(len(dfv)), 
                           dfv['z_phot'].values, 
                           dfv['log_Lbol'].values])
pX = pinv(X_ctrl, rcond=1e-6)
resid_FR = dfv['log_FR_444_150'].values - X_ctrl @ (pX @ dfv['log_FR_444_150'].values)

# Add to dataframe
dfv['resid_FR_444_150'] = resid_FR

# Partial regression on residuals
coef_resid = lstsq(np.column_stack([np.ones(len(dfv)), dfv['log_Sigma'].values]),
                    resid_FR, rcond=None)[0]

print(f"Panel B: rho_p={rp_b:.3f}, p={pp_b:.2e} ({sig_b:.1f}sigma)")

# --- Panel C stats (quartile groups) ---
group_stats = []
q_order = ['Q1 (diffuse)', 'Q2', 'Q3', 'Q4 (compact)']
for i, ql in enumerate(q_order):
    mask = dfv['Sigma_Q'] == ql
    sub = dfv[mask]
    group_stats.append({
        'q_label': q_labels[i],
        'n': len(sub),
        'mean_FR': sub['log_FR_444_150'].mean(),
        'std_FR': sub['log_FR_444_150'].std(),
        'median_FR': sub['log_FR_444_150'].median(),
        'values': sub['log_FR_444_150'].values
    })

# KS test Q1 vs Q4
q1_vals = dfv[dfv['Sigma_Q'] == 'Q1 (diffuse)']['log_FR_444_150'].values
q4_vals = dfv[dfv['Sigma_Q'] == 'Q4 (compact)']['log_FR_444_150'].values
ks_stat, ks_pval = ks_2samp(q1_vals, q4_vals)
ks_sig = sig_from_p(ks_pval)
print(f"Panel C KS(Q1 vs Q4): D={ks_stat:.3f}, p={ks_pval:.2e} ({ks_sig:.1f}sigma)")

# --- Panel D stats (multi-band comparison) ---
pairs_d = [
    ('log_Sigma', 'log_FR_444_150', 'F444W/F150W'),
    ('log_Sigma', 'log_FR_444_356', 'F444W/F356W'),
    ('log_C',     'log_FR_444_150', 'C vs F444W/F150W'),
    ('log_C',     'log_FR_444_356', 'C vs F444W/F356W'),
]
results_d = []
for x_col, y_col, label in pairs_d:
    rp, pp = partial_spearman(dfv[x_col].values, dfv[y_col].values, controls, dfv)
    sp = sig_from_p(pp) if not np.isnan(rp) else 0
    results_d.append({'label': label, 'rho_p': rp if not np.isnan(rp) else 0, 
                      'sig': sp, 'p_val': pp})
    m = "***" if pp < 0.001 else "**" if pp < 0.01 else "*" if pp < 0.05 else ""
    print(f"  {label:22s}: rho_p={rp:+.3f} ({sp:.1f}sigma){m}")

# ════════════════════════════════════════════════════════
# FIGURE CONSTRUCTION
# ════════════════════════════════════════════════════════

fig = plt.figure(figsize=(FIG_WIDTH_MM/25.4, FIG_HEIGHT_MM/25.4), dpi=DPI)

# Grid spec for tight control — more spacing between panels
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.35,
                      left=0.08, right=0.95, top=0.90, bottom=0.08)

ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

# ────────────────────────────────────────────────────────
# PANEL (a): Raw correlation — Σ vs F444W/F150W
# ────────────────────────────────────────────────────────

sc_a = ax_a.scatter(dfv['log_Sigma'], dfv['log_FR_444_150'],
                    c=dfv['z_phot'], cmap='coolwarm', s=14, alpha=0.72,
                    edgecolors='none', vmin=3.5, vmax=7.5)

# Trend line
x_grid = np.linspace(dfv['log_Sigma'].min(), dfv['log_Sigma'].max(), 200)
y_trend = coef_a[0] + coef_a[1] * x_grid
ax_a.plot(x_grid, y_trend, color=CB_PALETTE['vermilion'], lw=1.8, ls='-',
          label=f'fit slope = {coef_a[1]:+.3f}', zorder=5, alpha=0.85)

# Stats annotation box
stats_text = (r'$\rho$ = ' + f'{rho_a:+.3f}' + '\n'
              'P = ' + f'{pa_a:.1e}' + '\n'
              r'S = ' + f'{sig_a:.1f}' + r'$\sigma$' + '\n'
              'N = ' + f'{len(dfv)}')
props = dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#cccccc', 
             alpha=0.92, linewidth=0.5)
ax_a.text(0.96, 0.04, stats_text, transform=ax_a.transAxes, fontsize=7.5,
          verticalalignment='bottom', horizontalalignment='right',
          bbox=props, family='monospace')

ax_a.set_xlabel(r'log $\Sigma_\star\ ({\rm M}_\odot\,{\rm pc}^{-2})$', fontsize=9.5)
ax_a.set_ylabel(r'log($F_{444W}/F_{150W}$)', fontsize=9.5)
ax_a.set_title('(a)  Raw Correlation', fontsize=10, loc='left', fontweight='bold',
               pad=3)

# Colorbar
cbar_a = plt.colorbar(sc_a, ax=ax_a, shrink=0.82, pad=0.02)
cbar_a.set_label(r'$z_{\rm phot}$', fontsize=8)
cbar_a.ax.tick_params(labelsize=7)

# ────────────────────────────────────────────────────────
# PANEL (b): Partial correlation — residuals vs Σ
# ────────────────────────────────────────────────────────

sc_b = ax_b.scatter(dfv['log_Sigma'], dfv['resid_FR_444_150'],
                    c=dfv['log_Lbol'], cmap='viridis', s=14, alpha=0.72,
                    edgecolors='none')

# Residual trend line
x_grid_b = np.linspace(dfv['log_Sigma'].min(), dfv['log_Sigma'].max(), 200)
y_trend_b = coef_resid[0] + coef_resid[1] * x_grid_b
ax_b.plot(x_grid_b, y_trend_b, color=CB_PALETTE['blue'], lw=1.8, ls='--',
          label=f'slope = {coef_resid[1]:+.3f}', zorder=5, alpha=0.85)

# Zero reference line
ax_b.axhline(y=0, color=CB_PALETTE['gray'], lw=0.6, ls=':', zorder=1, alpha=0.6)

# Stats box
stats_text_b = (r'$\rho_{\rm p}$ = ' + f'{rp_b:+.3f}' + '\n'
                'P = ' + f'{pp_b:.1e}' + '\n'
                r'S = ' + f'{sig_b:.1f}' + r'$\sigma$' + '\n'
                r'ctrl: $z$, $L_{\rm bol}$')
ax_b.text(0.96, 0.96, stats_text_b, transform=ax_b.transAxes, fontsize=7.5,
          verticalalignment='top', horizontalalignment='right',
          bbox=props, family='monospace')

ax_b.set_xlabel(r'log $\Sigma_\star\ ({\rm M}_\odot\,{\rm pc}^{-2})$', fontsize=9.5)
ax_b.set_ylabel('Residual log($F_{444W}/F_{150W}$)', fontsize=9.5)
ax_b.set_title(r'(b)  Partial Correlation (controlled)', 
               fontsize=10, loc='left', fontweight='bold', pad=3)

cbar_b = plt.colorbar(sc_b, ax=ax_b, shrink=0.82, pad=0.02)
cbar_b.set_label(r'log$L_{\rm bol}$', fontsize=8)
cbar_b.ax.tick_params(labelsize=7)

# ────────────────────────────────────────────────────────
# PANEL (c): Quartile distributions
# ────────────────────────────────────────────────────────

bp_positions = range(4)
bp_data = [g['values'] for g in group_stats]

# Box plot with custom styling
bp = ax_c.boxplot(bp_data, positions=bp_positions, widths=0.55,
                  patch_artist=True, showmeans=True,
                  meanline=False, showfliers=False,
                  medianprops=dict(color='black', lw=1.5),
                  whiskerprops=dict(color='#555555', lw=0.8),
                  capprops=dict(color='#555555', lw=0.8),
                  meanprops=dict(marker='D', markerfacecolor=CB_PALETTE['vermilion'],
                                 markeredgecolor='white', markeredgewidth=0.4, 
                                 markersize=5, zorder=5))

for patch, color in zip(bp['boxes'], QUARTILE_COLORS):
    patch.set_facecolor(color)
    patch.set_alpha(0.65)
    patch.set_edgecolor('#333333')
    patch.set_lw(0.7)

# Individual points (jittered overlay)
np.random.seed(42)
for i, vals in enumerate(bp_data):
    jit = np.random.normal(0, 0.06, len(vals))
    ax_c.scatter(np.full_like(vals, i) + jit, vals, 
                 color=QUARTILE_COLORS[i], s=5, alpha=0.28, edgecolors='none', zorder=2)

# Annotate group sizes and means above each box
for i, g in enumerate(group_stats):
    ax_c.annotate(f'n={g["n"]}', xy=(i, g['mean_FR']+0.08), ha='center', fontsize=6.5,
                  color='#444444')
    # Arrow showing monotonic trend
    if i < 3:
        pass  # optional: add arrows between means

ax_c.set_xticks(bp_positions)
ax_c.set_xticklabels(q_labels, fontsize=7.5)
ax_c.set_ylabel('log($F_{444W}/F_{150W}$)', fontsize=9.5)
ax_c.set_title(r'(c)  $\Sigma$ Quartiles  (KS: ' + f'D={ks_stat:.2f}, ' + f'{ks_sig:.1f}' + r'$\sigma$)',
               fontsize=10, loc='left', fontweight='bold', pad=3)

# Horizontal reference at overall median
overall_median = dfv['log_FR_444_150'].median()
ax_c.axhline(y=overall_median, color=CB_PALETTE['gray'], lw=0.8, ls='--', 
             alpha=0.5, zorder=1)

# ────────────────────────────────────────────────────────
# PANEL (d): Multi-band bar chart of partial correlations
# ────────────────────────────────────────────────────────

labels_bar = [r['label'] for r in results_d]
rhos_bar = [r['rho_p'] for r in results_d]
sigs_bar = [r['sig'] for r in results_d]
pvals_bar = [r['p_val'] for r in results_d]

x_pos = np.arange(len(labels_bar))
bar_colors = []
for r in results_d:
    if r['p_val'] < 0.001:
        bar_colors.append(CB_PALETTE['blue'])      # *** strong
    elif r['p_val'] < 0.01:
        bar_colors.append(CB_PALETTE['vermilion']) # ** moderate
    elif r['p_val'] < 0.05:
        bar_colors.append(CB_PALETTE['orange'])     # * marginal
    else:
        bar_colors.append(CB_PALETTE['gray'])       # ns

bars = ax_d.bar(x_pos, rhos_bar, width=0.58, color=bar_colors, 
                edgecolor='#333333', linewidth=0.5, alpha=0.78)

# Zero line
ax_d.axhline(y=0, color='black', lw=0.7, zorder=1)

# Error bars (bootstrap-ish: use 1/sqrt(N) as approximate SE)
n_src = len(dfv)
se_approx = 1.0 / np.sqrt(n_src)
ax_d.errorbar(x_pos, rhos_bar, yerr=[se_approx]*len(rhos_bar), 
              fmt='none', ecolor='#333333', capsize=3, lw=0.7, elinewidth=0.8)

# Significance annotations
for i, (x, r, s, pv) in enumerate(zip(x_pos, rhos_bar, sigs_bar, pvals_bar)):
    y_pos = r + 0.03 if r >= 0 else r - 0.045
    va = 'bottom' if r >= 0 else 'top'
    
    # Star notation
    stars = "***" if pv < 0.001 else "**" if pv < 0.01 else "*" if pv < 0.05 else ""
    ax_d.text(x, y_pos, stars, ha='center', va=va, fontsize=13, fontweight='bold',
              color='#222222')
    
    # Sigma value below bar
    ax_d.text(x, -(se_approx + 0.025), f'{s:.1f}' + r'$\sigma$', 
              ha='center', va='top', fontsize=6.5, color='#555555')

ax_d.set_xticks(x_pos)
ax_d.set_xticklabels(labels_bar, rotation=18, ha='right', fontsize=7)
ax_d.set_ylabel(r'Partial $\rho$ (vs. $\Sigma$ or $C$)', fontsize=9.5)
ax_d.set_title(r'(d)  Multi-Band Comparison', 
               fontsize=10, loc='left', fontweight='bold', pad=3)
ax_d.set_ylim(-0.12, 0.48)

# Legend for significance
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='*', color='w', markerfacecolor='#222222',
           markersize=10, markeredgecolor='none', label=r'P $<$ 0.001 ($***$)'),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='#222222',
           markersize=8, markeredgecolor='none', label=r'P $<$ 0.01 ($**$)'),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='#222222',
           markersize=6, markeredgecolor='none', label=r'P $<$ 0.05 ($*$)'),
]
ax_d.legend(handles=legend_elements, loc='lower right', fontsize=6.5, 
            framealpha=0.85, edgecolor='#cccccc', fancybox=False)

# ════════════════════════════════════════════════════════
# MAIN FIGURE TITLE
# ════════════════════════════════════════════════════════

fig.suptitle(
    r'Density-Dependent Spectral Reddening in LRDs: Model-Independent Evidence from JWST/NIRCam Flux Ratios ($N=260$)',
    fontsize=12, fontweight='bold', y=0.97
)

# ════════════════════════════════════════════════════════
# SAVE
# ════════════════════════════════════════════════════════

figpath = f'{OUTPUT_DIR}/Figure1_PubQuality_DensityFluxRatio.png'
fig.savefig(figpath, dpi=DPI, facecolor='white', 
            bbox_inches='tight', pad_inches=0.06)

# Also save PDF (vector for submission)
figpath_pdf = f'{OUTPUT_DIR}/Figure1_PubQuality_DensityFluxRatio.pdf'
fig.savefig(figpath_pdf, format='pdf', facecolor='white',
            bbox_inches='tight', pad_inches=0.06)

print(f"\nSaved: {figpath}")
print(f"Saved: {figpath_pdf}")

plt.close(fig)
