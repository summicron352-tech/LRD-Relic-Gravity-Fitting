#!/usr/bin/env python3
"""
Figure 1 — ApJL Letter Version
=================================
Optimized for ApJL single-column width (~152mm) or double-column (~252mm).
More compact than MNRAS version; fonts and spacing tuned for Letter format.

Panel layout (2×2):
  (a) log Σ vs F444W/F150W scatter + regression
  (b) Partial correlation residual vs log Σ  
  (c) F444W/F150W distribution by Σ quartile
  (d) Multi-band partial-ρ comparison bar chart
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import spearmanr, norm
from numpy.linalg import lstsq, pinv
from scipy.stats import ks_2samp, rankdata
import warnings
warnings.filterwarnings('ignore')

# ════════════════════════════════════════════════════════
# ApJL-SPECIFIC SETTINGS
# ════════════════════════════════════════════════════════
# ApJ guidelines: single column = 5.0 inches (127mm), double = 7.3 inches (185mm)
# We use double-column for readability of 2x2 layout
FIG_WIDTH_INCH = 7.3    # ApJ double-column
FIG_HEIGHT_INCH = 6.0   # Slightly taller for breathing room
DPI = 300

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 8.5,
    'axes.labelsize': 9,
    'axes.titlesize': 9.5,
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 7.5,
    'legend.fontsize': 7,
    'lines.linewidth': 1.1,
    'lines.markersize': 2.2,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 3.5,
    'ytick.major.size': 3.5,
    'axes.linewidth': 0.5,
    'text.usetex': False,
    'mathtext.fontset': 'dejavusans',
})

CB_PALETTE = {
    'blue': '#0072B2', 'orange': '#D55E00', 'green': '#009E73',
    'pink': '#CC79A7', 'gray': '#777777', 'sky': '#56B4E9',
    'vermilion': '#E69F00', 'dark': '#222222',
}
QUARTILE_COLORS = [CB_PALETTE['sky'], CB_PALETTE['green'],
                   CB_PALETTE['orange'], CB_PALETTE['blue']]

# ════════════════════════════════════════════════════════
# DATA
# ════════════════════════════════════════════════════════
DATA_PATH = '/Users/tanxin/Desktop/数据处理/LRD_GitHub_Release_v2/data/csv/Kokorev_LRDs_Full.csv'
OUTPUT_DIR = '/Users/tanxin/Desktop/数据处理/LRD_GitHub_Release_v2/APJL_Letter_DensityDependentRedshift'

df_raw = pd.read_csv(DATA_PATH)

df = df_raw.copy()
BANDS = ['f115w_flux','f150w_flux','f200w_flux','f277w_flux','f356w_flux','f444w_flux']
for b in BANDS:
    df = df[df[b] > 0]

df['Mstar_proxy'] = df['lbol']
df['R_eff_pc'] = df['r_eff_50_phys'] * 1000
df['Sigma'] = df['Mstar_proxy'] / (2 * np.pi * df['R_eff_pc']**2)
df['log_Sigma'] = np.log10(df['Sigma'])
df['log_Lbol'] = np.log10(df['lbol'].clip(lower=0.01))
df['FR_444_356'] = df['f444w_flux'] / df['f356w_flux']
df['FR_444_150'] = df['f444w_flux'] / df['f150w_flux']
df['log_FR_444_356'] = np.log10(df['FR_444_356'])
df['log_FR_444_150'] = np.log10(df['FR_444_150'])

# Compactness
G_cgs, Msun_g, pc_cm, c_cgs = 6.674e-8, 1.989e33, 3.086e18, 2.998e10
df['C'] = G_cgs*df['Mstar_proxy']*Msun_g/(df['R_eff_pc']*pc_cm*c_cgs**2)
df['log_C'] = np.log10(df['C'])

# Clean sample
dfv = df.dropna(subset=['log_Sigma','log_FR_444_150','log_FR_444_356','z_phot','log_Lbol','lbol']).copy()
dfv = dfv[(dfv['log_Sigma'] < 15) & np.isfinite(dfv['log_Sigma'])]
print(f"Clean sample: {len(dfv)} sources")

dfv['Sigma_Q'] = pd.qcut(dfv['log_Sigma'], 4, labels=['Q1(diffuse)','Q2','Q3','Q4(compact)'])
q_labels = ['Q1\n(diffuse)', 'Q2', 'Q3', 'Q4\n(compact)']

# ════════════════════════════════════════════════════════
# STATISTICS
# ════════════════════════════════════════════════════════
def sig_from_p(p):
    return abs(norm.ppf(p/2)) if p > 0 else 99.0

def partial_spearman(x, y, ctrl_cols, data):
    rx, ry = rankdata(x), rankdata(y)
    X_ctrl = np.column_stack([np.ones(len(data))] + [data[c].values for c in ctrl_cols])
    try:
        pX = pinv(X_ctrl, rcond=1e-6)
        return spearmanr(rx - X_ctrl @ (pX @ rx), ry - X_ctrl @ (pX @ ry))
    except:
        return np.nan, np.nan

# Panel A
rho_a, pa_a = spearmanr(dfv['log_Sigma'], dfv['log_FR_444_150'])
sig_a = sig_from_p(pa_a)
coef_a = lstsq(np.c_[np.ones(len(dfv)), dfv['log_Sigma']], dfv['log_FR_444_150'].values, rcond=None)[0]
print(f"Panel A: rho={rho_a:+.3f} ({sig_a:.1f}sigma)")

# Panel B
controls = ['z_phot','log_Lbol']
rp_b, pp_b = partial_spearman(dfv['log_Sigma'].values, dfv['log_FR_444_150'].values, controls, dfv)
sig_b = sig_from_p(pp_b)
X_ctrl = np.c_[np.ones(len(dfv)), dfv['z_phot'].values, dfv['log_Lbol'].values]
pX = pinv(X_ctrl, rcond=1e-6)
resid_FR = dfv['log_FR_444_150'].values - X_ctrl @ (pX @ dfv['log_FR_444_150'].values)
dfv['resid_FR'] = resid_FR
coef_resid = lstsq(np.c_[np.ones(len(dfv)), dfv['log_Sigma']], resid_FR, rcond=None)[0]
print(f"Panel B: rho_p={rp_b:+.3f} ({sig_b:.1f}sigma)")

# Panel C
q_order = ['Q1(diffuse)','Q2','Q3','Q4(compact)']
group_stats = []
for i, ql in enumerate(q_order):
    sub = dfv[dfv['Sigma_Q']==ql]
    group_stats.append({'q':q_labels[i], 'n':len(sub),
                        'mean':sub['log_FR_444_150'].mean(),
                        'vals':sub['log_FR_444_150'].values})

q1_vals = dfv[dfv['Sigma_Q']=='Q1(diffuse)']['log_FR_444_150'].values
q4_vals = dfv[dfv['Sigma_Q']=='Q4(compact)']['log_FR_444_150'].values
ks_stat, ks_pval = ks_2samp(q1_vals, q4_vals)
ks_sig = sig_from_p(ks_pval)
print(f"Panel C KS: D={ks_stat:.3f} ({ks_sig:.1f}sigma)")

# Panel D
pairs_d = [('log_Sigma','log_FR_444_150','F$_{\\rm 444W}$ / F$_{\\rm 150W}$'),
           ('log_Sigma','log_FR_444_356','F$_{\\rm 444W}$ / F$_{\\rm 356W}$')]
results_d = []
for xc, yc, lab in pairs_d:
    rp, pp = partial_spearman(dfv[xc].values, dfv[yc].values, controls, dfv)
    sp = sig_from_p(pp) if not np.isnan(rp) else 0
    results_d.append({'label':lab, 'rp':rp if not np.isnan(rp) else 0, 'sig':sp, 'p':pp})
    print(f"  {lab}: rho_p={rp:+.3f} ({sp:.1f}sigma)")

# ════════════════════════════════════════════════════════
# FIGURE — ApJL compact 2x2
# ════════════════════════════════════════════════════════
fig = plt.figure(figsize=(FIG_WIDTH_INCH, FIG_HEIGHT_INCH), dpi=DPI)
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.32,
                      left=0.07, right=0.96, top=0.91, bottom=0.08)

ax = [fig.add_subplot(gs[i//2, i%2]) for i in range(4)]
props_box = dict(boxstyle='round,pad=0.30', facecolor='white', edgecolor='#cccccc',
                 alpha=0.92, linewidth=0.4)

# ──── Panel (a): Raw ────
sc = ax[0].scatter(dfv['log_Sigma'], dfv['log_FR_444_150'],
                    c=dfv['z_phot'], cmap='coolwarm', s=11, alpha=0.70,
                    edgecolors='none', vmin=3.5, vmax=7.5)
xg = np.linspace(dfv['log_Sigma'].min(), dfv['log_Sigma'].max(), 200)
ax[0].plot(xg, coef_a[0]+coef_a[1]*xg, color=CB_PALETTE['vermilion'], lw=1.5, alpha=0.85, zorder=5)
ax[0].text(0.96, 0.04, r'$\rho$='+f'{rho_a:+.3f}'+'\nS='+f'{sig_a:.1f}$\sigma$\nN='+f'{len(dfv)}',
            transform=ax[0].transAxes, fontsize=6.8, va='bottom', ha='right', bbox=props_box, family='monospace')
ax[0].set_xlabel(r'log $\Sigma_\star$ (M$_\odot$ pc$^{-2}$)', fontsize=9)
ax[0].set_ylabel(r'log($F_{444W}/F_{150W}$)', fontsize=9)
ax[0].set_title('(a) Raw Correlation', loc='left', fontweight='bold', fontsize=9.5, pad=2.5)
cb0 = plt.colorbar(sc, ax=ax[0], shrink=0.80, pad=0.02)
cb0.set_label(r'$z_{\rm phot}$', fontsize=7.5)
cb0.ax.tick_params(labelsize=6.5)

# ──── Panel (b): Partial ────
sc_b = ax[1].scatter(dfv['log_Sigma'], dfv['resid_FR'],
                      c=dfv['log_Lbol'], cmap='viridis', s=11, alpha=0.70, edgecolors='none')
ax[1].plot(xg, coef_resid[0]+coef_resid[1]*xg, color=CB_PALETTE['blue'],
           lw=1.5, ls='--', alpha=0.85, zorder=5)
ax[1].axhline(0, color=CB_PALETTE['gray'], lw=0.5, ls=':', alpha=0.5, zorder=1)
ax[1].text(0.96, 0.96, r'$\rho_{\rm p}$='+f'{rp_b:+.3f}'+'\nS='+f'{sig_b:.1f}'+r'$\sigma$'+'\nctrl: ' + r'$z$, $L_{\rm bol}$',
            transform=ax[1].transAxes, fontsize=6.8, va='top', ha='right', bbox=props_box, family='monospace')
ax[1].set_xlabel(r'log $\Sigma_\star$ (M$_\odot$ pc$^{-2}$)', fontsize=9)
ax[1].set_ylabel('Residual log($F_{444W}/F_{150W}$)', fontsize=9)
ax[1].set_title(r'(b) Partial Correlation', loc='left', fontweight='bold', fontsize=9.5, pad=2.5)
cb1 = plt.colorbar(sc_b, ax=ax[1], shrink=0.80, pad=0.02)
cb1.set_label(r'log $L_{\rm bol}$', fontsize=7.5)
cb1.ax.tick_params(labelsize=6.5)

# ──── Panel (c): Quartiles ────
bp_data = [g['vals'] for g in group_stats]
bp = ax[2].boxplot(bp_data, positions=range(4), widths=0.52, patch_artist=True,
                     showmeans=True, showfliers=False,
                     medianprops=dict(color='black', lw=1.3),
                     whiskerprops=dict(color='#555', lw=0.7),
                     capprops=dict(color='#555', lw=0.7),
                     meanprops=dict(marker='D', mfc=CB_PALETTE['vermilion'], mec='white',
                                    mew=0.35, ms=4.5, zorder=5))
for p, col in zip(bp['boxes'], QUARTILE_COLORS):
    p.set_facecolor(col); p.set_alpha(0.60); p.set_edgecolor('#333'); p.set_lw(0.6)
np.random.seed(42)
for i, v in enumerate(bp_data):
    j = np.random.normal(0, 0.055, len(v))
    ax[2].scatter(np.full_like(v,i)+j, v, c=QUARTILE_COLORS[i], s=4, alpha=0.26, edgecolors='none', zorder=2)
for i, g in enumerate(group_stats):
    ax[2].annotate(f'n={g["n"]}', xy=(i,g['mean']+0.07), ha='center', fontsize=6, color='#444')
med_all = dfv['log_FR_444_150'].median()
ax[2].axhline(med_all, color=CB_PALETTE['gray'], lw=0.7, ls='--', alpha=0.45, zorder=1)
ax[2].set_xticks(range(4))
ax[2].set_xticklabels(q_labels, fontsize=7)
ax[2].set_ylabel('log($F_{444W}/F_{150W}$)', fontsize=9)
ax[2].set_title(r'(c) $\Sigma$ Quartiles  (KS D=' + f'{ks_stat:.2f}, ' + f'{ks_sig:.1f}' + r'$\sigma$)',
                loc='left', fontweight='bold', fontsize=9.5, pad=2.5)

# ──── Panel (d): Bar chart ────
labs_bar = [r['label'] for r in results_d]
rhos_bar = [r['rp'] for r in results_d]
sigs_bar = [r['sig'] for r in results_d]
pvs_bar = [r['p'] for r in results_d]
xp = np.arange(len(labs_bar))
cols_bar = [CB_PALETTE['blue'] if p<0.001 else CB_PALETTE['vermilion'] if p<0.01
            else CB_PALETTE['orange'] if p<0.05 else CB_PALETTE['gray'] for p in pvs_bar]

bars = ax[3].bar(xp, rhos_bar, width=0.50, color=cols_bar, edgecolor='#333', linewidth=0.45, alpha=0.78)
ax[3].axhline(0, color='#222', lw=0.6, zorder=1)
se = 1.0/np.sqrt(len(dfv))
ax[3].errorbar(xp, rhos_bar, yerr=[se]*len(rhos_bar), fmt='none', ecolor='#333', capsize=2.5, lw=0.65, elinewidth=0.75)
for i, (x, r, s, p) in enumerate(zip(xp, rhos_bar, sigs_bar, pvs_bar)):
    yp = r+0.03 if r>=0 else r-0.04
    st = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
    ax[3].text(x, yp, st, ha='center', va='bottom' if r>=0 else 'top', fontsize=12, fontweight='bold', color='#222')
    ax[3].text(x, -(se+0.022), f'{s:.1f}' + r'$\sigma$', ha='center', va='top', fontsize=6, color='#555')
ax[3].set_xticks(xp)
ax[3].set_xticklabels(labs_bar, fontsize=8)
ax[3].set_ylabel(r'Partial $\rho$', fontsize=9)
ax[3].set_title('(d) Multi-Band Comparison', loc='left', fontweight='bold', fontsize=9.5, pad=2.5)
ax[3].set_ylim(-0.12, 0.48)

from matplotlib.lines import Line2D
leg_elems = [Line2D([0],[0], marker='*', color='w', mfc='#222', ms=9, mec='none', label=r'$P<0.001$ ($***$)'),
             Line2D([0],[0], marker='*', color='w', mfc='#222', ms=7, mec='none', label=r'$P<0.01$ ($**$)'),
             Line2D([0],[0], marker='*', color='w', mfc='#222', ms=5, mec='none', label=r'$P<0.05$ ($*$)')]
ax[3].legend(handles=leg_elems, loc='lower right', fontsize=6, framealpha=0.85, edgecolor='#ccc')

# Title
fig.suptitle(
    r'Density-Dependent Spectral Reddening in LRDs: JWST/NIRCam Flux-Ratio Analysis ($N$=' + f'{len(dfv)}' + ')',
    fontsize=11, fontweight='bold', y=0.97)

# Save
for fmt in ['png', 'pdf']:
    fpath = f'{OUTPUT_DIR}/Figure1_ApJL_DensityFluxRatio.{fmt}'
    fig.savefig(fpath, dpi=DPI, facecolor='white', bbox_inches='tight', pad_inches=0.04)
    print(f'Saved: {fpath}')

plt.close(fig)
