#!/usr/bin/env python3
"""
Generate JADES DR5 Cross-Validation figures and tables for paper appendix.
Uses the existing JADES LRD candidate catalog from jades_cross_validation_v2.py output.

Output:
  - figB1_jades_sigma_color.png       : Σ-color scatter + partial correlation
  - figB2_jades_selection_cascade.png : Selection waterfall
  - tableB1_selection_cascade.tex     : Selection step table
  - tableB2_jades_cosmos_comparison.tex: Statistical comparison table
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, rankdata, norm
from numpy.linalg import lstsq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import warnings
warnings.filterwarnings('ignore')
import os

# =============================================================================
# Configuration
# =============================================================================
DATA_CSV = '/Users/tanxin/Desktop/数据处理/重新整理投稿APJ/源数据/JADES_CrossValidation/JADES_LRD_candidates_v2.csv'
OUT_DIR = '/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/paper/figures'
os.makedirs(OUT_DIR, exist_ok=True)

# Matplotlib style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# COSMOS baseline
COSMOS_BASELINE = {
    'N': 260,
    'spearman_rho': 0.341,
    'spearman_sigma': 5.6,
    'ks_D': 0.574,
    'ks_sigma': 6.2,
    'partial_rho': 0.341,
    'partial_sigma': 5.6,
}

# Selection cascade data (from jades_cross_validation_v2.py output)
SELECTION_CASCADE = {
    'GOODS-S': {
        'Total sources': 304366,
        'Valid photo-z': 304366,
        'z > 4': 77121,
        r'$R_{\rm eff} \leq 0.06^{\prime\prime}$': 9170,
        r'${\rm F277W}-{\rm F444W} > 0.5$': 2478,
        r'${\rm F115W}-{\rm F200W} > -0.5$ (BD cut)': 1241,
        'Detection (key bands)': 1219,
        'Valid F150W+F356W': 1074,
    },
    'GOODS-N': {
        'Total sources': 181144,
        'Valid photo-z': 181144,
        'z > 4': 38216,
        r'$R_{\rm eff} \leq 0.06^{\prime\prime}$': 5960,
        r'${\rm F277W}-{\rm F444W} > 0.5$': 966,
        r'${\rm F115W}-{\rm F200W} > -0.5$ (BD cut)': 499,
        'Detection (key bands)': 484,
        'Valid F150W+F356W': 393,
    }
}


def partial_spearman(x, y, ctrl_cols, data):
    """Spearman partial correlation: regress controls from ranks, correlate residuals."""
    rx = rankdata(x)
    ry = rankdata(y)
    X_ctrl = np.column_stack([np.ones(len(data))] + [data[c].values for c in ctrl_cols])
    resid_x = rx - X_ctrl @ lstsq(X_ctrl, rx, rcond=None)[0]
    resid_y = ry - X_ctrl @ lstsq(X_ctrl, ry, rcond=None)[0]
    rho, p = spearmanr(resid_x, resid_y)
    sigma = abs(norm.ppf(p/2)) if p > 0 else np.inf
    return rho, p, sigma


def sigma_from_norm_ppf(p, rho):
    """Convert p-value to Gaussian sigma."""
    if p <= 0 or p >= 1:
        return np.inf
    return abs(norm.ppf(p / 2))


# =============================================================================
# Figure B1: Σ-color scatter + partial correlation (JADES combined)
# =============================================================================
def make_fig_b1(df):
    """Σ-color scatter plot for JADES GOODS-N+S combined."""
    print("Generating Figure B1: JADES Σ-color correlation...")

    # Filter valid data
    mask = (np.isfinite(df['logSigma'])) & (np.isfinite(df['color_444_150']))
    d = df[mask].copy()
    N = len(d)
    print(f"  N_valid = {N}")

    # Raw Spearman
    rho_raw, p_raw = spearmanr(d['logSigma'], d['color_444_150'])
    sig_raw = sigma_from_norm_ppf(p_raw, rho_raw)

    # Partial Spearman (ctrl z + L_bol)
    rho_p, p_p, sig_p = partial_spearman(
        d['logSigma'].values, d['color_444_150'].values,
        ['z_phot', 'log_Lbol_proxy'], d
    )

    # KS test
    med_s = np.median(d['logSigma'])
    high = d[d['logSigma'] >= med_s]['color_444_150'].values
    low = d[d['logSigma'] < med_s]['color_444_150'].values
    D_ks, p_ks = stats.ks_2samp(high, low)
    sig_ks = sigma_from_norm_ppf(p_ks, D_ks)

    # Per-field stats
    stats_per_field = {}
    for field in ['GOODS-S', 'GOODS-N']:
        df_f = d[d['field'] == field]
        if len(df_f) > 10:
            rho_f, p_f = spearmanr(df_f['logSigma'], df_f['color_444_150'])
            sig_f = sigma_from_norm_ppf(p_f, rho_f)
            rho_fp, p_fp, sig_fp = partial_spearman(
                df_f['logSigma'].values, df_f['color_444_150'].values,
                ['z_phot', 'log_Lbol_proxy'], df_f
            )
            stats_per_field[field] = {
                'N': len(df_f), 'rho_raw': rho_f, 'sig_raw': sig_f,
                'rho_p': rho_fp, 'sig_p': sig_fp
            }

    # ── Create figure ──
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: Σ vs color scatter, color-coded by z
    ax = axes[0]
    z = d['z_phot'].values
    sc = ax.scatter(d['logSigma'], d['color_444_150'], c=z, s=8,
                    cmap='plasma', alpha=0.6, edgecolors='none', rasterized=True)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('$z_{\\rm phot}$', fontsize=10)

    # Fit line
    coeffs = np.polyfit(d['logSigma'], d['color_444_150'], 1)
    xs = np.linspace(d['logSigma'].min(), d['logSigma'].max(), 100)
    ax.plot(xs, np.polyval(coeffs, xs), 'k--', lw=1.2, alpha=0.7,
            label=f'Slope={coeffs[0]:+.3f}')

    ax.set_xlabel(r'$\log \Sigma_{\rm phys}$ [$M_\odot$\,pc$^{-2}$]')
    ax.set_ylabel(r'$\log(F_{\rm F444W}/F_{\rm F150W})$')
    ax.set_title(f'JADES GOODS-N+S\n$N={N}$, $\\rho={{{rho_raw:+.3f}}}$ (${{{sig_raw:.1f}}}\\sigma$)')
    ax.legend(fontsize=8)
    ax.text(0.95, 0.05, '(a)', transform=ax.transAxes, fontsize=11, fontweight='bold',
            va='bottom', ha='right')

    # Panel 2: By-field comparison (GOODS-S vs GOODS-N density contours)
    ax = axes[1]
    colors_field = {'GOODS-S': '#2196F3', 'GOODS-N': '#FF9800'}
    for field in ['GOODS-S', 'GOODS-N']:
        df_f = d[d['field'] == field]
        if len(df_f) > 10:
            ax.scatter(df_f['logSigma'], df_f['color_444_150'],
                      c=colors_field[field], s=4, alpha=0.4,
                      label=f"{field} ($N={{{len(df_f)}}}$, $\\rho={{{stats_per_field[field]['rho_raw']:+.3f}}}$)",
                      edgecolors='none', rasterized=True)

    ax.set_xlabel(r'$\log \Sigma_{\rm phys}$ [$M_\odot$\,pc$^{-2}$]')
    ax.set_ylabel(r'$\log(F_{\rm F444W}/F_{\rm F150W})$')
    ax.set_title('By Field')
    ax.legend(fontsize=8, markerscale=3)
    ax.text(0.95, 0.05, '(b)', transform=ax.transAxes, fontsize=11, fontweight='bold',
            va='bottom', ha='right')

    # Panel 3: Δcolor (high-Σ - low-Σ) bar + KS
    ax = axes[2]
    delta_color = np.mean(high) - np.mean(low)
    sigma_delta = delta_color / np.sqrt(np.var(high)/len(high) + np.var(low)/len(low))

    # Bar chart: mean color for low vs high Σ
    labels = [f'Low Σ\n(< median)\nN={len(low)}', f'High Σ\n(> median)\nN={len(high)}']
    means = [np.mean(low), np.mean(high)]
    errors = [np.std(low)/np.sqrt(len(low)), np.std(high)/np.sqrt(len(high))]
    bars = ax.bar(labels, means, yerr=errors, color=['#64B5F6', '#E57373'],
                  edgecolor='k', linewidth=0.8, capsize=5)
    ax.set_ylabel(r'$\langle \log(F_{\rm F444W}/F_{\rm F150W}) \rangle$')
    ax.set_title(f'KS $D={{{D_ks:.3f}}}$ (${{{sig_ks:.1f}}}\\sigma$)\n'
                 f'$\\Delta={{{delta_color:+.3f}}}$ dex (${{{sigma_delta:.1f}}}\\sigma$)')
    ax.text(0.95, 0.05, '(c)', transform=ax.transAxes, fontsize=11, fontweight='bold',
            va='bottom', ha='right')

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'figB1_jades_sigma_color.png')
    fig.savefig(outpath)
    print(f"  → Saved: {outpath}")
    plt.close(fig)

    return {
        'N': N, 'N_GS': (d['field']=='GOODS-S').sum(), 'N_GN': (d['field']=='GOODS-N').sum(),
        'rho_raw': rho_raw, 'sig_raw': sig_raw, 'p_raw': p_raw,
        'rho_partial': rho_p, 'sig_partial': sig_p, 'p_partial': p_p,
        'KS_D': D_ks, 'KS_sigma': sig_ks, 'KS_p': p_ks,
        'delta_color': delta_color, 'delta_sigma': sigma_delta,
        'slope': coeffs[0],
        'GOODS_S': stats_per_field.get('GOODS-S', {}),
        'GOODS_N': stats_per_field.get('GOODS-N', {}),
    }


# =============================================================================
# Figure B2: Selection cascade waterfall
# =============================================================================
def make_fig_b2():
    """Waterfall chart showing selection cascade for both GOODS fields."""
    print("Generating Figure B2: Selection cascade...")

    fig, ax = plt.subplots(1, 1, figsize=(10, 5.5))

    steps_short = [
        'Total',
        r'$z_{\rm phot}>4$',
        r'$R_{\rm eff}\leq0.06^{\prime\prime}$',
        r'${\rm F277W-F444W}>0.5$',
        'BD cut',
        'Detection',
        'Valid\ncolors',
    ]

    gs_vals = [304366, 77121, 9170, 2478, 1241, 1219, 1074]
    gn_vals = [181144, 38216, 5960, 966, 499, 484, 393]

    x = np.arange(len(steps_short))
    width = 0.35

    bars_gs = ax.bar(x - width/2, gs_vals, width, label='GOODS-S',
                     color='#2196F3', edgecolor='k', linewidth=0.5)
    bars_gn = ax.bar(x + width/2, gn_vals, width, label='GOODS-N',
                     color='#FF9800', edgecolor='k', linewidth=0.5)

    ax.set_yscale('log')
    ax.set_ylabel('Number of sources')
    ax.set_xticks(x)
    ax.set_xticklabels(steps_short, fontsize=8)
    ax.legend(fontsize=10)
    ax.set_title('JADES DR5 LRD Selection Cascade (Rinaldi+ 2026 Criteria)')
    ax.grid(axis='y', alpha=0.3, lw=0.5)

    # Annotate final counts
    ax.annotate(f'{gs_vals[-1]:,}', xy=(x[-1] - width/2, gs_vals[-1]),
                ha='center', va='bottom', fontsize=8, fontweight='bold', color='#1565C0')
    ax.annotate(f'{gn_vals[-1]:,}', xy=(x[-1] + width/2, gn_vals[-1]),
                ha='center', va='bottom', fontsize=8, fontweight='bold', color='#E65100')

    # Annotate combined
    ax.annotate(f'Combined: {gs_vals[-1] + gn_vals[-1]:,}',
                xy=(x[-1], (gs_vals[-1] + gn_vals[-1]) * 0.5),
                ha='center', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'figB2_jades_selection_cascade.png')
    fig.savefig(outpath)
    print(f"  → Saved: {outpath}")
    plt.close(fig)


# =============================================================================
# Table B1: Selection cascade (LaTeX)
# =============================================================================
def make_table_b1():
    """Generate LaTeX table for selection cascade."""
    print("Generating Table B1: Selection cascade...")

    lines = []
    lines.append(r'\begin{deluxetable}{lrrr}')
    lines.append(r'\tablecaption{JADES DR5 LRD Selection Cascade \label{tab:jades_cascade}}')
    lines.append(r'\tablehead{')
    lines.append(r'\colhead{Selection Step} & \colhead{GOODS-S} & \colhead{GOODS-N} & \colhead{Combined}')
    lines.append(r'}')
    lines.append(r'\startdata')

    cascade = SELECTION_CASCADE
    steps = list(cascade['GOODS-S'].keys())
    for step in steps:
        gs = cascade['GOODS-S'][step]
        gn = cascade['GOODS-N'][step]
        combined = gs + gn
        step_label = step.replace('_', r'\_')
        lines.append(f'{step_label} & {gs:,} & {gn:,} & {combined:,} \\\\')

    lines.append(r'\enddata')
    lines.append(r'\tablecomments{Selection criteria follow \citet{Rinaldi2026}: '
                 r'photometric redshift $z>4$, compactness $R_{\rm eff} \leq 0.06^{\prime\prime}$, '
                 r'red color ${\rm F277W}-{\rm F444W} > 0.5$, BD rejection cut '
                 r'${\rm F115W}-{\rm F200W} > -0.5$, detection in F277W, F444W, F150W, F356W, '
                 r'and valid color measurements in F150W and F356W.}')
    lines.append(r'\end{deluxetable}')

    tex = '\n'.join(lines)
    outpath = os.path.join(OUT_DIR, 'tableB1_selection_cascade.tex')
    with open(outpath, 'w') as f:
        f.write(tex)
    print(f"  → Saved: {outpath}")

    return tex


# =============================================================================
# Table B2: Statistical comparison COSMOS vs JADES (LaTeX)
# =============================================================================
def make_table_b2(stats):
    """Generate LaTeX comparison table."""
    print("Generating Table B2: Statistical comparison...")

    lines = []
    lines.append(r'\begin{deluxetable}{lccc}')
    lines.append(r'\tablecaption{$\Sigma$--Color Correlation: COSMOSWeb vs.\ JADES DR5 '
                 r'\label{tab:jades_comparison}}')
    lines.append(r'\tablehead{')
    lines.append(r'\colhead{Statistic} & \colhead{COSMOSWeb ($z$=7--9, $N$=5,601)} & '
                 r'\colhead{JADES Combined ($N$=' + f'{stats["N"]:,}' + r')} & '
                 r'\colhead{JADES GOODS-S ($N$=' + f'{stats["N_GS"]:,}' + r')}')
    lines.append(r'}')
    lines.append(r'\startdata')

    def fmt_rho_sig(val, sig, ndigits=3):
        """Format rho ± sigma string, handling sign."""
        s = f'{val:+.{ndigits}f}'
        return f'${s}$ (${{{sig:.1f}}}\\sigma$)'

    # Raw Spearman
    lines.append(
        r'Raw Spearman $\rho$ & '
        r'$+0.618$ (${\gg}20\sigma$) & '
        + fmt_rho_sig(stats["rho_raw"], stats["sig_raw"]) + ' & '
        + fmt_rho_sig(stats["GOODS_S"].get("rho_raw", 0), stats["GOODS_S"].get("sig_raw", 0)) + ' \\\\'
    )

    # Partial Spearman
    partial_label = r'Partial Spearman $\rho$ (ctrl $z+L_{\rm bol}$ proxy)'
    lines.append(
        partial_label + ' & '
        r'$+0.516$ (${\gg}20\sigma$) & '
        + fmt_rho_sig(stats["rho_partial"], stats["sig_partial"]) + ' & '
        + fmt_rho_sig(stats["GOODS_S"].get("rho_p", 0), stats["GOODS_S"].get("sig_p", 0)) + r' \\'
    )

    # KS test
    lines.append(
        r'KS $D$ & '
        r'$0.473$ ($36.0\sigma$) & '
        f'${stats["KS_D"]:.3f}$ (${{{stats["KS_sigma"]:.1f}}}\\sigma$) & '
        f'--- \\\\'
    )

    # Δcolor
    lines.append(
        f'$\\Delta\\langle\\log(F_{{444}}/F_{{150}})\\rangle$ '
        f'(high-$\\Sigma$ $-$ low-$\\Sigma$) & '
        r'$+0.449$ ($38.2\sigma$) & '
        + fmt_rho_sig(stats["delta_color"], stats["delta_sigma"]) + ' & '
        f'--- \\\\'
    )

    # Slope
    slope_str = f'{stats["slope"]:+.3f}'
    lines.append(
        r'Linear slope $d({\rm color})/d(\log\Sigma)$ & '
        r'$+0.222$ & '
        f'${slope_str}$ & '
        r'--- \\'
    )

    lines.append(r'\enddata')
    lines.append(r'\tablecomments{Partial Spearman controls for photometric redshift '
                 r'and $L_{\rm bol}$ proxy (estimated from multi-band NIRCam fluxes). '
                 r'COSMOSWeb statistics from the main text (\S\ref{sec:ev4}). '
                 r'KS $D$ compares color distributions above/below median $\Sigma$. '
                 r'GOODS-N alone ($N$=' + f'{stats["N_GN"]:,}' + r') yields partial $\\rho=' +
                 f'{stats["GOODS_N"].get("rho_p", 0):+.3f}$ '
                 f'(${{{stats["GOODS_N"].get("sig_p", 0):.1f}}}\\sigma$); '
                 r'its smaller sample limits independent statistical power.}')
    lines.append(r'\end{deluxetable}')

    tex = '\n'.join(lines)
    outpath = os.path.join(OUT_DIR, 'tableB2_jades_cosmos_comparison.tex')
    with open(outpath, 'w') as f:
        f.write(tex)
    print(f"  → Saved: {outpath}")

    return tex


# =============================================================================
# Main
# =============================================================================
def main():
    print("=" * 65)
    print("  JADES DR5 Cross-Validation — Appendix Figures & Tables")
    print("=" * 65)

    # Load data
    df = pd.read_csv(DATA_CSV)
    print(f"\nLoaded {len(df):,} JADES LRD candidates")
    print(f"  GOODS-S: {(df['field']=='GOODS-S').sum():,}")
    print(f"  GOODS-N: {(df['field']=='GOODS-N').sum():,}")

    # ── Redshift distribution ──
    z = df['z_phot']
    print(f"\nRedshift distribution:")
    for zlo, zhi in [(4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 12)]:
        n = np.sum((z >= zlo) & (z < zhi))
        print(f"  z ∈ [{zlo},{zhi}): {n:>5d}  ({100*n/len(df):.1f}%)")

    # ── Generate figures ──
    print()
    stats = make_fig_b1(df)
    make_fig_b2()

    # ── Generate tables ──
    print()
    make_table_b1()
    print()
    make_table_b2(stats)

    # ── Print stats summary ──
    print(f"\n{'='*65}")
    print(f"  Key Statistics for LaTeX")
    print(f"{'='*65}")
    print(f"  N_combined      = {stats['N']:,}")
    print(f"  N_GOODS-S       = {stats['N_GS']:,}")
    print(f"  N_GOODS-N       = {stats['N_GN']:,}")
    print(f"  Raw ρ           = {stats['rho_raw']:+.4f} ({stats['sig_raw']:.1f}σ)")
    print(f"  Partial ρ       = {stats['rho_partial']:+.4f} ({stats['sig_partial']:.1f}σ)")
    print(f"  KS D            = {stats['KS_D']:.4f} ({stats['KS_sigma']:.1f}σ)")
    print(f"  Δcolor          = {stats['delta_color']:+.4f} dex ({stats['delta_sigma']:.1f}σ)")
    print(f"  Slope           = {stats['slope']:+.4f}")

    print(f"\n  DONE! 🐢")


if __name__ == '__main__':
    main()
