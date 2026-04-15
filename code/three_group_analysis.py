"""
LRD Three-Group Parameter Distribution Analysis
===============================================
Groups: Support(inc.Weak) / Neutral / Null wins
Analysis: Histograms + KS tests + Best discriminator finder
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 0. Load Data & Merge
# ============================================================
DATA_DIR = '/Users/tanxin/Desktop/数据处理'
results = pd.read_csv(f'{DATA_DIR}/Bulletproof_Results.csv')
original = pd.read_csv(f'{DATA_DIR}/Kokorev_LRDs_Full.csv')

# Merge verdict into original data
df = results[['id','verdict','delta_chi2nu','delta_bic',
              'Av_null','Av_geff','C_AGN','alpha_AGN',
              'chi2nu_null','chi2nu_geff']].merge(
    original, on='id', how='left'
)

# ============================================================
# 1. Compute Derived Parameters FIRST
# ============================================================

def compute_uv_beta(row):
    waves = np.array([1.15, 1.50, 2.00])
    cols = ['f115w_flux', 'f150w_flux', 'f200w_flux']
    fluxes = []
    for c in cols:
        v = row.get(c, np.nan)
        if pd.isna(v) or v < -90:
            return np.nan
        fluxes.append(v)
    fluxes = np.array(fluxes)
    mask = fluxes > 0
    if mask.sum() < 2:
        return np.nan
    log_w = np.log(waves[mask])
    log_f = np.log(fluxes[mask])
    slope, _ = np.polyfit(log_w, log_f, 1)
    beta = slope - 1
    return beta

def compute_ir_excess(row):
    f444 = row.get('f444w_flux', np.nan)
    f150 = row.get('f150w_flux', np.nan)
    if pd.isna(f444) or pd.isna(f150) or f444 < -90 or f150 <= 0:
        return np.nan
    return f444 / f150

def compute_color_red(row):
    f444 = row.get('f444w_flux', np.nan)
    f115 = row.get('f115w_flux', np.nan)
    if pd.isna(f444) or pd.isna(f115) or f444 < -90 or f115 < -90:
        return np.nan
    if f444 > 0 and f115 > 0:
        return -2.5 * np.log10(f444 / f115)
    return np.nan

def compute_mstar_proxy(row):
    muv = row.get('muv', np.nan)
    if pd.isna(muv) or muv < -99:
        return np.nan
    log_mstar = 10.5 - 0.4 * (muv + 20)
    return log_mstar

df['uv_beta'] = df.apply(compute_uv_beta, axis=1)
df['ir_excess'] = df.apply(compute_ir_excess, axis=1)
df['color_f444_f115'] = df.apply(compute_color_red, axis=1)
df['log_mstar_proxy'] = df.apply(compute_mstar_proxy, axis=1)

# ============================================================
# 3. Define All Parameters to Analyze
# ============================================================
CD_COL = 'Cohen_d_SN'   # avoid apostrophe in column name

param_config = {
    r'$z_{\rm phot}$': {
        'col': 'z_phot', 'label': 'Photometric Redshift', 'unit': '',
        'bins': np.linspace(4, 9, 15), 'range': (4, 9),
    },
    r'$F_{\rm F444W}$ ($\mu$Jy)': {
        'col': 'f444w_flux', 'label': 'F444W Flux (MIR proxy)', 'unit': r'$\mu$Jy',
        'bins': np.linspace(-0.05, 0.8, 20), 'range': (-0.02, 0.7),
    },
    r'$L_{\rm bol}$ (erg/s)': {
        'col': 'lbol', 'label': 'Bolometric Luminosity (AGN strength proxy)', 'unit': '',
        'bins': np.linspace(42, 47, 20), 'range': (42, 47),
    },
    r'$\beta_{\rm UV}$': {
        'col': 'uv_beta', 'label': 'UV Slope (from F115W/F150W/F200W)', 'unit': '',
        'bins': np.linspace(-3, 3, 24), 'range': (-3, 3),
    },
    r'$A_V$ (SED fit)': {
        'col': 'av', 'label': 'Dust Extinction from SED fitting', 'unit': 'mag',
        'bins': np.linspace(0, 4, 16), 'range': (0, 4),
    },
    r'$r_{\rm eff,50}$ (pc)': {
        'col': 'r_eff_50_phys', 'label': 'Half-light Radius', 'unit': 'pc',
        'bins': np.linspace(40, 200, 16), 'range': (40, 200),
    },
    r'IR$_{\rm excess}$ (F444/F150)': {
        'col': 'ir_excess', 'label': 'IR Excess (MIR brightness)', 'unit': '',
        'bins': np.linspace(0, 30, 20), 'range': (0, 25),
    },
    r'$m_{\rm UV}$ (mag)': {
        'col': 'muv', 'label': 'Absolute UV Magnitude', 'unit': 'mag',
        'bins': np.linspace(-22, -16, 20), 'range': (-22, -16),
    },
    r'$C_{\rm AGN}$': {
        'col': 'C_AGN', 'label': 'AGN Component Amplitude', 'unit': '',
        'bins': np.linspace(0, 1.0, 20), 'range': (0, 1),
    },
}

# ============================================================
# 2. Define Three Groups (AFTER computing derived parameters)
# ============================================================
df['group'] = df['verdict'].map({
    'Support': 'Support', 'Weak': 'Support',
    'Neutral': 'Neutral', 'Null wins': 'Null wins', 'FAILED': 'FAILED'
})

grp_sup = df[df['group'] == 'Support']
grp_neu = df[df['group'] == 'Neutral']
grp_null = df[df['group'] == 'Null wins']

print("="*70)
print("THREE-GROUP DISTRIBUTION ANALYSIS")
print("="*70)
print(f"Support (incl.Weak): {len(grp_sup):3d} sources")
print(f"Neutral:             {len(grp_neu):3d} sources")
print(f"Null wins:           {len(grp_null):3d} sources")

# ============================================================
# 4. Run KS Tests
# ============================================================
ks_results = []
for name, cfg in param_config.items():
    col = cfg['col']
    s = grp_sup[col].dropna()
    n = grp_neu[col].dropna()
    nw = grp_null[col].dropna()

    result = {'Parameter': name, 'Label': cfg['label']}

    # KS: Support vs Null wins
    if len(s) > 5 and len(nw) > 5:
        ks_sn, p_sn = stats.ks_2samp(s, nw)
        result['KS_SN_stat'] = ks_sn
        result['KS_SN_pval'] = p_sn
        result['KS_SN_sig'] = '***' if p_sn < 0.001 else ('**' if p_sn < 0.01 else ('*' if p_sn < 0.05 else ''))
        pooled_std = np.sqrt((s.std()**2 + nw.std()**2) / 2)
        if pooled_std > 0:
            result[CD_COL] = abs(s.mean() - nw.mean()) / pooled_std

    # KS: Support vs Neutral
    if len(s) > 5 and len(n) > 5:
        ks_sne, p_sne = stats.ks_2samp(s, n)
        result['KS_SNe_stat'] = ks_sne
        result['KS_SNe_pval'] = p_sne

    # KS: Neutral vs Null wins
    if len(n) > 5 and len(nw) > 5:
        ks_nn, p_nn = stats.ks_2samp(n, nw)
        result['KS_NtNw_stat'] = ks_nn
        result['KS_NtNw_pval'] = p_nn

    result['Sup_mean'] = round(s.mean(), 3) if len(s) > 0 else np.nan
    result['Neu_mean'] = round(n.mean(), 3) if len(n) > 0 else np.nan
    result['Nul_mean'] = round(nw.mean(), 3) if len(nw) > 0 else np.nan

    ks_results.append(result)

ks_df = pd.DataFrame(ks_results)

print("\n" + "="*70)
print("KOLMOGOROV-SMIRNOV TEST RESULTS")
print("="*70)

if CD_COL in ks_df.columns and 'KS_SN_pval' in ks_df.columns:
    ks_sorted = ks_df.sort_values('KS_SN_pval')

    hdr = '{:<25s}| {:>7s}| {:>12s}| {:>4s}| {:>10s}| {:>8s}| {:>8s}'.format(
        'Parameter','KS stat','p-value','Sig','Cohen d','Sup mean','Nul mean')
    print(hdr)
    print("-"*95)

    for _, row in ks_sorted.iterrows():
        ks_val = row.get('KS_SN_stat', np.nan)
        pv = row.get('KS_SN_pval', np.nan)
        sig = str(row.get('KS_SN_sig', ''))
        cd = row.get(CD_COL, np.nan)
        sm = row.get('Sup_mean', np.nan)
        nm = row.get('Nul_mean', np.nan)

        ks_s = "{:.3f}".format(ks_val) if not pd.isna(ks_val) else "N/A"
        p_s = "{:.2e}".format(pv) if not pd.isna(pv) else "N/A"
        cd_s = "{:.3f}".format(cd) if not pd.isna(cd) else "N/A"
        sm_s = "{:.3f}".format(sm) if not pd.isna(sm) else "N/A"
        nm_s = "{:.3f}".format(nm) if not pd.isna(nm) else "N/A"

        ln = '{:<23s}| {:>7s}| {:>12s}| {:>4s}| {:>10s}| {:>8s}| {:>8s}'.format(
            str(row['Parameter']), ks_s, p_s, sig, cd_s, sm_s, nm_s)
        print(ln)

# ============================================================
# 5. Visualization
# ============================================================
COLORS = {'Support': '#27ae60', 'Neutral': '#3498db', 'Null wins': '#e74c3c'}

avail_params = param_config
n_params = len(avail_params)
n_cols = 3
n_rows = int(np.ceil(n_params / n_cols))

pdf_path = f'{DATA_DIR}/Three_Group_Distribution_Analysis.pdf'

with PdfPages(pdf_path) as pdf:

    # ---- Page 1: Summary Table ----
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_axes([0.05, 0.08, 0.9, 0.88])
    ax.axis('off')
    
    title_str = (
        "LRD Three-Group Parameter Distribution Analysis\n"
        "Support (incl. Weak) / Neutral / Null Wins\n\n"
        "Total={} | Support={} ({:.1f}%) | Neutral={} ({:.1f}%) | Null wins={} ({:.1f}%)\n".format(
            len(df), len(grp_sup), len(grp_sup)/len(df)*100,
            len(grp_neu), len(grp_neu)/len(df)*100,
            len(grp_null), len(grp_null)/len(df)*100
        )
    )
    ax.text(0.5, 0.97, title_str, transform=ax.transAxes,
            ha='center', va='top', fontsize=14, fontweight='bold', family='monospace')

    lines = ["=" * 95]
    header = '{:<22s}| {:>7s}| {:>12s}| {:>4s}| {:>10s}| {:>8s}| {:>8s}| {:>8s}'.format(
        'Parameter','KS(S-N)','p-val','Sig',"Cohen's d",'Sup_mean','Neu_mean','Nul_mean')
    lines.append(header)
    lines.append("-" * 95)

    if 'KS_SN_pval' in ks_df.columns:
        for _, rr in ks_df.sort_values('KS_SN_pval').iterrows():
            ks_v = "{:.3f}".format(rr.get('KS_SN_stat', 0))
            pv_v = "{:.2e}".format(rr.get('KS_SN_pval', 1))
            sg_v = str(rr.get('KS_SN_sig', ''))
            cd_v = "{:.3f}".format(rr.get(CD_COL, 0))
            sv = "{:.3f}".format(rr.get('Sup_mean', 0))
            nv = "{:.3f}".format(rr.get('Neu_mean', 0))
            nu = "{:.3f}".format(rr.get('Nul_mean', 0))
            ln = '{:<21s}| {:>7s}| {:>12s}| {:>4s}| {:>10s}| {:>8s}| {:>8s}| {:>8s}'.format(
                str(rr['Parameter']), ks_v, pv_v, sg_v, cd_v, sv, nv, nu)
            lines.append(ln)

    lines.append("=" * 95)
    body_text = "\n".join(lines)
    ax.text(0.02, 0.87, body_text, transform=ax.transAxes,
            ha='left', va='top', fontsize=8, family='monospace')

    notes = (
        "\nNotes:\n"
        "- X-ray upper limits: NOT available in Kokorev catalog (no Chandra/XMM cross-match)\n"
        "- H-alpha line width: NOT available (requires spectroscopic data)\n"
        "- MIRI F770W flux: NOT available; using F444W (reddest NIRCam filter) as MIR proxy\n"
        "- Stellar mass: PROXY from M_UV relation; direct SED masses not provided in catalog\n"
        "- UV beta slope: COMPUTED from F115W/F150W/F200W power-law fit per source\n"
        "\nStatistical notes:\n"
        "- Significance: *** p<0.001, ** p<0.01, * p<0.05, ns=not significant\n"
        "- Cohen effect size: |d|<0.2 negligible, 0.2-0.5 small, 0.5-0.8 medium, >0.8 large\n"
        "- KS test: non-parametric two-sample test; sensitive to both location and shape differences"
    )
    ax.text(0.02, 0.26, notes, transform=ax.transAxes,
            ha='left', va='top', fontsize=7.5, family='monospace', color='#555555')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ---- Page 2: Main histogram panel ----
    fig2, axes2 = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    axes2_flat = axes2.flatten() if n_params > 1 else [axes2]

    for idx2, (name2, cfg2) in enumerate(avail_params.items()):
        ax2 = axes2_flat[idx2]
        col2 = cfg2['col']

        sd = grp_sup[col2].dropna()
        nd = grp_neu[col2].dropna()
        nwd = grp_null[col2].dropna()

        alld = pd.concat([sd, nd, nwd]).dropna()
        lo, hi = alld.quantile(0.005), alld.quantile(0.995)
        bp = np.linspace(lo, hi, 25)

        for gn, gd, gc, ga in [
            ('Support', sd, COLORS['Support'], 0.65),
            ('Neutral', nd, COLORS['Neutral'], 0.55),
            ('Null wins', nwd, COLORS['Null wins'], 0.7),
        ]:
            ax2.hist(gd, bins=bp, density=True, alpha=ga, color=gc,
                    edgecolor='white', linewidth=1,
                    label='{g} (N={n})'.format(g=gn, n=len(gd)))

        kr = ks_df[ks_df['Parameter'] == name2]
        if len(kr) > 0:
            pv = kr.iloc[0].get('KS_SN_pval', np.nan)
            kst = kr.iloc[0].get('KS_SN_stat', np.nan)
            if not pd.isna(pv):
                sg = '***' if pv < 0.001 else ('**' if pv < 0.01 else ('*' if pv < 0.05 else 'ns'))
                tstr = '{}\nKS={:.2f}, p={:.1e} {}'.format(name2, kst, pv, sg)
                ax2.set_title(tstr, fontsize=10, fontweight='bold')

        ax2.legend(fontsize=6.5, loc='upper right')
        ax2.set_xlabel(cfg2['label'], fontsize=8)
        ax2.set_ylabel('Density', fontsize=8)

    for jdx2 in range(idx2 + 1, len(axes2_flat)):
        axes2_flat[jdx2].set_visible(False)

    plt.suptitle('Distribution Comparison: Support vs Neutral vs Null Wins',
                 fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    pdf.savefig(fig2, bbox_inches='tight')
    plt.close(fig2)

    # ---- Pages 3+: Individual parameter deep dives ----
    for name3, cfg3 in list(avail_params.items()):
        fig3, axes3 = plt.subplots(1, 3, figsize=(15, 4.5))
        col3 = cfg3['col']

        sd3 = grp_sup[col3].dropna()
        nd3 = grp_neu[col3].dropna()
        nwd3 = grp_null[col3].dropna()

        alld3 = pd.concat([sd3, nd3, nwd3]).dropna()
        lo3, hi3 = alld3.quantile(0.01), alld3.quantile(0.99)
        bp3 = np.linspace(lo3, hi3, 25)

        for ai, (gd, gn, gc) in enumerate([
            (sd3, 'Support', COLORS['Support']),
            (nd3, 'Neutral', COLORS['Neutral']),
            (nwd3, 'Null wins', COLORS['Null wins']),
        ]):
            ax3 = axes3[ai]
            ax3.hist(gd, bins=bp3, density=True, alpha=0.7, color=gc,
                    edgecolor='white', linewidth=1)

            mval = gd.mean()
            medval = np.median(gd)
            sval = gd.std()

            ymax3 = ax3.get_ylim()[1] * 0.9
            ax3.axvline(mval, color='black', linestyle='--', linewidth=2,
                       label='Mean={v:.3f}'.format(v=mval))
            ax3.axvline(medval, color='gray', linestyle=':', linewidth=1.5,
                       label='Median={v:.3f}'.format(v=medval))
            ax3.axvspan(mval - 2*sval, mval + 2*sval, alpha=0.15, color=gc,
                       label='+/-2 sigma')

            t3 = '{}: {}\nN={}, Mean={:.3f}, Med={:.3f}'.format(gn, name3, len(gd), mval, medval)
            ax3.set_title(t3, fontsize=10, fontweight='bold', color=gc)
            ax3.set_xlabel(cfg3['label'])
            ax3.set_ylabel('Density')
            ax3.legend(fontsize=6.5, loc='upper right')

            # Mini boxplot
            bx = ax3.inset_axes([0.72, 0.05, 0.26, 0.85])
            bpp = bx.boxplot(gd, vert=True, widths=0.6, patch_artist=True)
            bpp['boxes'][0].set_facecolor(gc)
            bpp['boxes'][0].set_alpha(0.5)
            bx.set_yticks([])
            bx.set_xticks([])
            for spine in bx.spines.values():
                spine.set_visible(True)

        plt.suptitle('Deep Dive: {}'.format(name3), fontsize=13, fontweight='bold')
        plt.tight_layout()
        pdf.savefig(fig3, bbox_inches='tight')
        plt.close(fig3)

    # ---- Last page: Discrimination Power Summary Bar Chart ----
    fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(14, 6))

    if 'KS_SN_pval' in ks_df.columns and CD_COL in ks_df.columns:
        kv = ks_df.dropna(subset=['KS_SN_pval', CD_COL]).copy()
        kv = kv.sort_values('KS_SN_pval')

        ps_list = [p.replace('$','').replace('\\rm ','').replace('\\mathrm{','').replace('}','')[:22]
                   for p in kv['Parameter']]

        # Left panel: KS statistics
        clist_a = ['#27ae60' if p < 0.05 else ('#f39c12' if p < 0.1 else '#e74c3c') for p in kv['KS_SN_pval']]
        bars_a = ax4a.barh(range(len(kv)), kv['KS_SN_stat'], color=clist_a, edgecolor='white')
        ax4a.set_yticks(range(len(kv)))
        ax4a.set_yticklabels(ps_list, fontsize=9)
        ax4a.set_xlabel('KS Statistic (Support vs Null)', fontsize=11)
        ax4a.set_title('Discrimination Power: Higher = More Different', fontsize=12, fontweight='bold')
        ax4a.invert_yaxis()
        for bi, (bar, rv) in enumerate(zip(bars_a, kv.iterrows())):
            _, r = rv
            pv4 = r['KS_SN_pval']
            ax4a.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                     'p={:.1e}'.format(pv4), va='center', fontsize=7)

        # Right panel: Cohen's d
        cd_valid = kv.dropna(subset=[CD_COL]).copy()
        if len(cd_valid) > 0:
            clist_b = ['#27ae60' if d > 0.8 else ('#f39c12' if d > 0.5 else ('#f1c40f' if d > 0.2 else '#bdc3c7'))
                      for d in cd_valid[CD_COL]]
            ps_b = [p.replace('$','').replace('\\rm ','').replace('\\mathrm{','').replace('}','')[:22]
                    for p in cd_valid['Parameter']]
            bars_b = ax4b.barh(range(len(cd_valid)), cd_valid[CD_COL], color=clist_b, edgecolor='white')
            ax4b.set_yticks(range(len(cd_valid)))
            ax4b.set_yticklabels(ps_b, fontsize=9)
            ax4b.set_xlabel("Cohen's d Effect Size", fontsize=11)
            ax4b.set_title("|d|>0.8 large, 0.5-0.8 medium, <0.2 negligible",
                          fontsize=11, fontweight='bold')
            ax4b.invert_yaxis()
            ax4b.axvline(0.8, color='#e74c3c', linestyle='--', alpha=0.5, label='Large threshold')
            ax4b.axvline(0.5, color='#f39c12', linestyle=':', alpha=0.5, label='Medium threshold')
            ax4b.legend(fontsize=8, loc='lower right')

            # Annotate best discriminator
            best_idx_local = cd_valid[CD_COL].idxmax()
            best_param_local = cd_valid.loc[best_idx_local, 'Parameter']
            best_d_local = cd_valid.loc[best_idx_local, CD_COL]
            best_p_local = cd_valid.loc[best_idx_local, 'KS_SN_pval']
            ypos = list(cd_valid.index).index(best_idx_local)
            ax4b.annotate(
                'BEST:\n{}\nd={:.2f}\np={:.1e}'.format(best_param_local, best_d_local, best_p_local),
                xy=(best_d_local, ypos),
                xytext=(best_d_local + 0.3, ypos - 1),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=9, color='red', fontweight='bold'
            )

    plt.suptitle("Which Parameter Best Separates Support from Null Wins?",
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    pdf.savefig(fig4, bbox_inches='tight')
    plt.close(fig4)

print("\nFull report saved to:", pdf_path)

# ============================================================
# 6. Print Final Summary
# ============================================================
print("\n" + "="*70)
print("BEST DISCRIMINATORS (Support vs Null wins)")
print("="*70)

if 'KS_SN_pval' in ks_df.columns and CD_COL in ks_df.columns:
    ranked2 = ks_df.dropna(subset=['KS_SN_pval', CD_COL]).sort_values('KS_SN_pval')

    print('\n{:<4s}| {:<22s}| {:>7s}| {:>12s}| {:>10s}| {}'.format(
        'Rank','Parameter','KS stat','p-value',"Cohen's d",'Interpretation'))
    print("-"*95)

    for rank, (_, rw) in enumerate(ranked2.iterrows(), 1):
        ksv = rw['KS_SN_stat']
        pv = rw['KS_SN_pval']
        cdv = rw[CD_COL]
        if cdv > 0.8:
            interp = "LARGE difference"
        elif cdv > 0.5:
            interp = "Medium diff"
        elif cdv > 0.2:
            interp = "Small but real"
        else:
            interp = "Negligible"
        if pv < 0.001:
            interp += ", HIGHLY SIGNIFICANT"
        elif pv < 0.01:
            interp += ", significant"
        elif pv < 0.05:
            interp += ", marginal"

        print('{:>4d}| {:<21s}| {:>7.3f}| {:>12.1e}| {:>10.3f}| {}'.format(
            rank, str(rw['Parameter']), ksv, pv, cdv, interp))

    # Winner
    best_row = ranked2.iloc[0]
    print('\n' + '*'*55)
    print('  * BEST DISCRIMINATOR: {}'.format(best_row['Parameter']))
    print('  * KS statistic = {:.3f}, p = {:.1e}'.format(best_row['KS_SN_stat'], best_row['KS_SN_pval']))
    print("  * Cohen's d = {:.3f}".format(best_row[CD_COL]))
    print('*'*55)

# Save CSV
csv_out = f'{DATA_DIR}/KS_Test_Results_ThreeGroups.csv'
ks_df.to_csv(csv_out, index=False)
print('\nKS test results saved to:', csv_out)
