"""
G_eff Window Model — Parameter Space Scan for LRD SED Fitting
=============================================================
Scans (A_G, gamma, beta, z_on) to find optimal window parameters
that maximize support for G_eff model across CEERS LRDs.

Outputs:
  1. ParamScan_Heatmap_AGamma.png   — Δχ²_ν heatmap (A_G × γ)
  2. ParamScan_Slices.png           — 1D slices per parameter
  3. ParamScan_TopParams.csv        — Top-N parameter combinations ranked by score
  4. ParamScan_BestFit.png          — SED fits with best-scoring parameters
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import product
import warnings
warnings.filterwarnings('ignore')

# ==========================================
# 0. G_eff Window Model + SED Models (same as figure 3 updat.py)
# ==========================================

def H_z(z, Omega_m=0.3, Omega_L=0.7, H0=70.0):
    return H0 * np.sqrt(Omega_m * (1+z)**3 + Omega_L)


def G_eff_window(z, A_G=0.15, beta=0.5, z_p=10.0, z_on=8.0, z_off=30.0, dz=2.0):
    W = 0.25 * (1 + np.tanh((z - z_on) / dz)) * (1 - np.tanh((z - z_off) / dz))
    Hz_ratio = H_z(z) / H_z(z_p)
    return A_G * (Hz_ratio ** beta) * W


def z_dist_from_phot(z_phot, gamma=1.0, **wp):
    delta_G = G_eff_window(z_phot, **wp)
    return gamma * np.clip(delta_G, 0, None)


def f_uv(wave_rest):
    return (wave_rest / 0.3)**(-0.5)


def f_int(wave_rest):
    break_wave = 0.364
    flux = np.zeros_like(wave_rest)
    mask = wave_rest > break_wave
    flux[~mask] = 0.1 * (wave_rest[~mask] / break_wave)**2
    flux[mask] = ((wave_rest[mask] / break_wave)**0.5 * np.exp(-(wave_rest[mask] - 1.0)**2 / 2.0) + 
                  0.5 * (1.0 / wave_rest[mask]))
    return flux


def sed_model(wave_rest, C_blue, C_red, z_dist):
    distorted_wave = wave_rest / (1 + z_dist)
    return C_blue * f_uv(wave_rest) + C_red * f_int(distorted_wave) / (1 + z_dist)


def fit_two_param(rest_wave, nuFnu_obs, nuFnu_err, z_dist_fixed):
    """Fixed z_dist, fit C_blue + C_red only."""
    def model(wave, C_blue, C_red):
        return sed_model(wave, C_blue, C_red, z_dist_fixed)
    
    try:
        popt, pcov = curve_fit(model, rest_wave, nuFnu_obs,
                               sigma=nuFnu_err,
                               bounds=([0, 0], [10, 10]),
                               maxfev=10000)
        
        n_data = len(rest_wave)
        k_params = 2
        model_flux = model(rest_wave, popt[0], popt[1])
        residuals = (nuFnu_obs - model_flux) / nuFnu_err
        
        chi2 = np.sum(residuals**2)
        dof = n_data - k_params
        red_chi2 = chi2 / dof if dof > 0 else chi2
        bic = k_params * np.log(n_data) + chi2
        
        return {'chi2': chi2, 'red_chi2': red_chi2, 'bic': bic,
                'C_blue': popt[0], 'C_red': popt[1],
                'model_flux': model_flux, 'residuals': residuals,
                'success': True}
    except Exception as e:
        return {'success': False, 'error': str(e)}


# ==========================================
# 1. 数据加载（与主脚本一致）
# ==========================================
print("Loading data...")
df = pd.read_csv('Kokorev_LRDs_Full.csv')
ceers_df = df[df['field'] == 'ceers-full']
target_ids = ceers_df['id'].head(6).tolist()

filters = {
    'f090w': 0.90, 'f115w': 1.15, 'f150w': 1.50, 'f200w': 2.00,
    'f277w': 2.77, 'f356w': 3.56, 'f444w': 4.44, 'f770w': 7.70
}

sources_data = []
for target_id in target_ids:
    source = ceers_df[ceers_df['id'] == target_id].iloc[0]
    z_phot = source['z_phot']
    
    obs_wave, obs_flux, obs_err = [], [], []
    for f, wave in filters.items():
        flux_col = f'{f}_flux'
        err_col = f'{f}_fluxerr'
        if flux_col in source and pd.notnull(source[flux_col]) and source[flux_col] > -90:
            obs_wave.append(wave)
            obs_flux.append(source[flux_col])
            obs_err.append(source[err_col] if source[err_col] > 0 else source[flux_col] * 0.1)
    
    obs_wave = np.array(obs_wave)
    obs_flux = np.array(obs_flux)
    obs_err = np.array(obs_err)
    rest_wave = obs_wave / (1 + z_phot)
    
    nuFnu_obs = obs_flux / obs_wave
    nuFnu_err = obs_err / obs_wave
    norm_factor = np.max(nuFnu_obs) if np.max(nuFnu_obs) > 0 else 1.0
    
    sources_data.append({
        'target_id': target_id, 'z_phot': z_phot,
        'rest_wave': rest_wave,
        'nuFnu_norm': nuFnu_obs / norm_factor,
        'nuFnu_err_norm': nuFnu_err / norm_factor,
    })
    print(f"  ID={target_id} (z={z_phot:.2f}, {len(obs_wave)} bands)")

n_sources = len(sources_data)
print(f"\nLoaded {n_sources} sources.\n")

# ==========================================
# 2. 参数空间定义
# ==========================================
scan_grid = {
    'A_G':   np.linspace(0.03, 0.40, 8),      # 引力增强幅度
    'gamma': np.linspace(0.5,  8.0,  8),       # 比例常数（δG → z_dist 映射）
    'beta':  np.linspace(0.2,  2.5,  6),       # H(z) 幂律指数
    'z_on':  np.array([2.0, 3.0, 4.0, 5.0]),  # 窗口开启红移
}

# 固定参数（不扫描的）
fixed_params = {
    'z_p':   7.0,
    'z_off': 25.0,
    'dz':    3.0,
}

print("=" * 70)
print("PARAMETER SPACE SCAN — G_eff Window Model")
print("=" * 70)
for key, vals in scan_grid.items():
    print(f"  {key:>6s}: {len(vals)} values from {vals.min():.3f} to {vals.max():.3f}")
total_combos = np.prod([len(v) for v in scan_grid.values()])
print(f"  Total combinations: {int(total_combos)}")
print(f"  Total fits needed:  ~{int(total_combos) * n_sources * 2:.0f}")
print()

# ==========================================
# 3. 扫描引擎：先跑 A_G × γ 二维网格（固定 beta, z_on）
# ==========================================

# Phase 1: 粗扫描 A_G × γ（固定 beta=1.0, z_on=3.0)
print("-" * 70)
print("Phase 1: Coarse 2D Scan (A_G × gamma) — fixing beta=1.0, z_on=3.0")
print("-" * 70)

AG_grid = scan_grid['A_G']
gamma_grid = scan_grid['gamma']
beta_fixed = 1.0
z_on_fixed = 3.0

results_2d = []  # (AG_idx, gamma_idx) → stats dict

for i_ag, AG_val in enumerate(AG_grid):
    for i_gam, gam_val in enumerate(gamma_grid):
        
        wp = {
            'A_G': AG_val, 'beta': beta_fixed,
            'z_p': fixed_params['z_p'],
            'z_on': z_on_fixed, 'z_off': fixed_params['z_off'], 'dz': fixed_params['dz']
        }
        
        total_dc2 = 0.0
        total_dbic = 0.0
        n_support = 0
        n_nullwin = 0
        dc2_per_source = []
        valid_count = 0
        
        for src in sources_data:
            z_phot = src['z_phot']
            
            # Null model
            fit_n = fit_two_param(src['rest_wave'], src['nuFnu_norm'], src['nuFnu_err_norm'], 0.0)
            
            # G_eff model
            zd = z_dist_from_phot(z_phot, gamma=gam_val, **wp)
            fit_g = fit_two_param(src['rest_wave'], src['nuFnu_norm'], src['nuFnu_err_norm'], zd)
            
            if fit_n['success'] and fit_g['success']:
                dc2 = fit_n['red_chi2'] - fit_g['red_chi2']
                dbic = fit_n['bic'] - fit_g['bic']
                
                dc2_per_source.append(dc2)
                total_dc2 += dc2
                total_dbic += dbic
                valid_count += 1
                
                if dc2 > 1 and dbic < -2:
                    n_support += 1
                elif dc2 <= -1 and dbic > 2:
                    n_nullwin += 1
        
        mean_dc2 = total_dc2 / valid_count if valid_count > 0 else -999
        med_dbic = np.median([fit_n['bic'] - fit_g['bic'] for _ in [1]])  # placeholder
        min_dc2 = min(dc2_per_source) if dc2_per_source else -999
        
        # 综合得分：mean Δχ²_ν 加权 + 支持/反对惩罚
        # 正分：平均改善；负分：最差源拖后腿
        score = mean_dc2 + min_dc2 * 0.3 + n_support * 5 - n_nullwin * 10
        
        results_2d.append({
            'A_G': AG_val, 'gamma': gam_val,
            'beta': beta_fixed, 'z_on': z_on_fixed,
            'mean_delta_chi2': mean_dc2,
            'min_delta_chi2': min_dc2,
            'n_support': n_support,
            'n_nullwin': n_nullwin,
            'score': score,
            'dc2_per_source': dc2_per_source.copy(),
        })

# 转换为矩阵用于热图绘图
heatmap_mean = np.full((len(gamma_grid), len(AG_grid)), np.nan)
heatmap_min  = np.full((len(gamma_grid), len(AG_grid), n_sources), np.nan)
heatmap_score = np.full((len(gamma_grid), len(AG_grid)), np.nan)

for r in results_2d:
    i_gam = list(gamma_grid).index(r['gamma'])
    i_ag = list(AG_grid).index(r['A_G'])
    heatmap_mean[i_gam, i_ag] = r['mean_delta_chi2']
    heatmap_score[i_gam, i_ag] = r['score']
    for isrc, dc2 in enumerate(r['dc2_per_source']):
        heatmap_min[i_gam, i_ag, isrc] = dc2

print(f"Phase 1 done: {len(results_2d)} parameter points evaluated.")

# 找最佳点
best_idx = np.argmax([r['score'] for r in results_2d])
best_2d = results_2d[best_idx]
print(f"  BEST (Phase 1): A_G={best_2d['A_G']:.3f}, gamma={best_2d['gamma']:.2f}, "
      f"score={best_2d['score']:.1f}, ⟨Δχ²ν⟩={best_2d['mean_delta_chi2']:+.1f}, "
      f"min Δχ²ν={best_2d['min_delta_chi2']:+.1f}")

# ==========================================
# 4. Phase 2: 围绕最佳点的精细扫描（加 beta 和 z_on）
# ==========================================
print("\n" + "-" * 70)
print("Phase 2: Fine Scan around best (A_G, gamma) — varying beta & z_on")
print("-" * 70)

AG_best = best_2d['A_G']
gam_best = best_2d['gamma']

beta_fine = scan_grid['beta']
z_on_fine = scan_grid['z_on']

results_fine = []

for beta_val in beta_fine:
    for zo_val in z_on_fine:
        wp = {
            'A_G': AG_best, 'beta': beta_val,
            'z_p': fixed_params['z_p'],
            'z_on': zo_val, 'z_off': fixed_params['z_off'], 'dz': fixed_params['dz']
        }
        
        total_dc2 = 0.0
        n_support = 0
        n_nullwin = 0
        dc2_ps = []
        
        for src in sources_data:
            z_phot = src['z_phot']
            fit_n = fit_two_param(src['rest_wave'], src['nuFnu_norm'], src['nuFnu_err_norm'], 0.0)
            zd = z_dist_from_phot(z_phot, gamma=gam_best, **wp)
            fit_g = fit_two_param(src['rest_wave'], src['nuFnu_norm'], src['nuFnu_err_norm'], zd)
            
            if fit_n['success'] and fit_g['success']:
                dc2 = fit_n['red_chi2'] - fit_g['red_chi2']
                dbic = fit_n['bic'] - fit_g['bic']
                dc2_ps.append(dc2)
                total_dc2 += dc2
                if dc2 > 1 and dbic < -2:
                    n_support += 1
                elif dc2 <= -1 and dbic > 2:
                    n_nullwin += 1
        
        mean_d = total_dc2 / len(dc2_ps) if dc2_ps else -999
        min_d = min(dc2_ps) if dc2_ps else -999
        score = mean_d + min_d * 0.3 + n_support * 5 - n_nullwin * 10
        
        results_fine.append({
            'A_G': AG_best, 'gamma': gam_best,
            'beta': beta_val, 'z_on': zo_val,
            'mean_delta_chi2': mean_d,
            'min_delta_chi2': min_d,
            'n_support': n_support,
            'n_nullwin': n_nullwin,
            'score': score,
            'dc2_per_source': dc2_ps.copy(),
        })

best_fidx = np.argmax([r['score'] for r in results_fine])
best_overall = results_fine[best_fidx]
print(f"  BEST OVERALL: A_G={best_overall['A_G']:.3f}, gamma={best_overall['gamma']:.2f}, "
      f"beta={best_overall['beta']:.2f}, z_on={best_overall['z_on']:.1f}, "
      f"score={best_overall['score']:.1f}, ⟨Δχ²ν⟩={best_overall['mean_delta_chi2']:+.1f}")

# ==========================================
# 5. 图 1: A_G × γ 热图（核心输出）
# ==========================================
fig_hm, axes_hm = plt.subplots(1, 3, figsize=(20, 7))

# --- Panel A: Mean Δχ²_ν heatmap ---
ax_a = axes_hm[0]
im_a = ax_a.imshow(heatmap_mean, aspect='auto', origin='lower',
                     extent=[AG_grid[0], AG_grid[-1], gamma_grid[0], gamma_grid[-1]],
                     cmap='RdYlGn', vmin=-100, vmax=50)
ax_a.set_xlabel(r'$A_G$ (Gravitational Enhancement Amplitude)', fontsize=12)
ax_a.set_ylabel(r'$\gamma$ (Scaling: $\delta G \to z_{\rm dist}$)', fontsize=12)
ax_a.set_title(r'(A) Mean $\Delta\chi^2_\nu$ Across All Sources', fontsize=12, fontweight='bold')

# 标记 fiducial 值
ax_a.axvline(x=0.15, color='blue', linestyle='--', linewidth=1.5, alpha=0.7, label='Fiducial $A_G=0.15$')
ax_a.axhline(y=3.0, color='purple', linestyle=':', linewidth=1.5, alpha=0.7, label=r'Fiducial $\gamma=3.0$')
ax_a.plot(best_2d['A_G'], best_2d['gamma'], '*', color='yellow', markersize=20,
          markeredgecolor='black', markeredgewidth=1.5, zorder=10, label=f'Best ({best_2d["A_G"]:.3f},{best_2d["gamma"]:.1f})')
ax_a.legend(fontsize=9, loc='upper right')

cbar_a = plt.colorbar(im_a, ax=ax_a, shrink=0.85)
cbar_a.set_label(r'Mean $\Delta\chi^2_\nu$', fontsize=11)

# --- Panel B: Min Δχ²_ν heatmap (worst source) ---
ax_b = axes_hm[1]
im_b = ax_b.imshow(np.min(heatmap_min, axis=2), aspect='auto', origin='lower',
                     extent=[AG_grid[0], AG_grid[-1], gamma_grid[0], gamma_grid[-1]],
                     cmap='RdYlGn_r', vmin=-600, vmax=20)  # reversed: green=good (least negative)
ax_b.set_xlabel(r'$A_G$', fontsize=12)
ax_b.set_ylabel(r'$\gamma$', fontsize=12)
ax_b.set_title(r'(B) Minimum $\Delta\chi^2_\nu$ (Worst Source)', fontsize=12, fontweight='bold')
ax_b.axvline(x=0.15, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
ax_b.axhline(y=3.0, color='purple', linestyle=':', linewidth=1.5, alpha=0.7)
ax_b.plot(best_2d['A_G'], best_2d['gamma'], '*', color='yellow', markersize=20,
          markeredgecolor='black', markeredgewidth=1.5, zorder=10)

cbar_b = plt.colorbar(im_b, ax=ax_b, shrink=0.85)
cbar_b.set_label(r'Min $\Delta\chi^2_\nu$ (higher = better)', fontsize=11)

# --- Panel C: Composite Score ---
ax_c = axes_hm[2]
im_c = ax_c.imshow(heatmap_score, aspect='auto', origin='lower',
                     extent=[AG_grid[0], AG_grid[-1], gamma_grid[0], gamma_grid[-1]],
                     cmap='viridis')
ax_c.set_xlabel(r'$A_G$', fontsize=12)
ax_c.set_ylabel(r'$\gamma$', fontsize=12)
ax_c.set_title('(C) Composite Score\n(mean + 0.3×min + 5×supp − 10×nullwin)', fontsize=12, fontweight='bold')
ax_c.axvline(x=0.15, color='white', linestyle='--', linewidth=1.5, alpha=0.7)
ax_c.axhline(y=3.0, color='cyan', linestyle=':', linewidth=1.5, alpha=0.7)
ax_c.plot(best_2d['A_G'], best_2d['gamma'], '*', color='red', markersize=20,
          markeredgecolor='white', markeredgewidth=1.5, zorder=10, label='Global Best')
ax_c.legend(fontsize=9, loc='upper right')

cbar_c = plt.colorbar(im_c, ax=ax_c, shrink=0.85)
cbar_c.set_label('Score (higher → better)', fontsize=11)

plt.suptitle(
    r'Parameter Scan: $G_{\rm eff}(z)$ Window Model — LRD SED Statistical Test'
    + '\nN=' + str(n_sources) + r' CEERS LRDs  |  Fixing $\beta$=' + str(beta_fixed)
    + r', $z_{\rm on}$=' + str(z_on_fixed),
    fontsize=14, y=1.02
)
plt.tight_layout()
fig_hm.savefig('ParamScan_Heatmap_AGamma.png', dpi=250, bbox_inches='tight')
print("✅ Saved: ParamScan_Heatmap_AGamma.png")

# ==========================================
# 6. 图 2: Beta × Z_on 细扫热图
# ==========================================
hm_beta = np.full((len(z_on_fine), len(beta_fine)), np.nan)
hm_beta_score = np.full((len(z_on_fine), len(beta_fine)), np.nan)

for r in results_fine:
    iz = list(z_on_fine).index(r['z_on'])
    ib = list(beta_fine).index(r['beta'])
    hm_beta[iz, ib] = r['mean_delta_chi2']
    hm_beta_score[iz, ib] = r['score']

fig_fine, (ax_fb, ax_fs) = plt.subplots(1, 2, figsize=(16, 6))

im_fb = ax_fb.imshow(hm_beta, aspect='auto', origin='lower',
                       extent=[beta_fine[0], beta_fine[-1], z_on_fine[0], z_on_fine[-1]],
                       cmap='RdYlGn', vmin=-100, vmax=50)
ax_fb.set_xlabel(r'$\beta$ ($H(z)$ Power Index)', fontsize=12)
ax_fb.set_ylabel(r'$z_{\rm on}$ (Window Turn-on Redshift)', fontsize=12)
ax_fb.set_title(r'(D) Mean $\Delta\chi^2_\nu$' + '\n' + r'($A_G$=' + f'{AG_best:.3f}' + r', $\gamma$=' + f'{gam_best:.1f}' + ')', fontsize=12, fontweight='bold')
ax_fb.plot(best_overall['beta'], best_overall['z_on'], '*r', markersize=18,
           markeredgecolor='yellow', markeredgewidth=1.5, zorder=10, label='Best')
ax_fb.legend(fontsize=9)
plt.colorbar(im_fb, ax=ax_fb, shrink=0.85, label=r'Mean $\Delta\chi^2_\nu$')

im_fs = ax_fs.imshow(hm_beta_score, aspect='auto', origin='lower',
                       extent=[beta_fine[0], beta_fine[-1], z_on_fine[0], z_on_fine[-1]],
                       cmap='viridis')
ax_fs.set_xlabel(r'$\beta$', fontsize=12)
ax_fs.set_ylabel(r'$z_{\rm on}$', fontsize=12)
ax_fs.set_title(r'(E) Composite Score' + '\n' + r'($A_G$=' + f'{AG_best:.3f}' + r', $\gamma$=' + f'{gam_best:.1f}' + ')', fontsize=12, fontweight='bold')
ax_fs.plot(best_overall['beta'], best_overall['z_on'], '*r', markersize=18,
           markeredgecolor='yellow', markeredgewidth=1.5, zorder=10, label='Best')
ax_fs.legend(fontsize=9)
plt.colorbar(im_fs, ax=ax_fs, shrink=0.85, label='Score')

plt.tight_layout()
fig_fine.savefig('ParamScan_Slices.png', dpi=250, bbox_inches='tight')
print("✅ Saved: ParamScan_Slices.png")

# ==========================================
# 7. 图 3: 最佳参数的完整拟合图
# ==========================================
wp_best = {
    'A_G': best_overall['A_G'], 'beta': best_overall['beta'],
    'z_p': fixed_params['z_p'],
    'z_on': best_overall['z_on'], 'z_off': fixed_params['z_off'], 'dz': fixed_params['dz']
}
gamma_best = best_overall['gamma']

print(f"\n--- Running final fits with best params: A_G={best_overall['A_G']:.3f}, "
      f"γ={gamma_best:.2f}, β={best_overall['beta']:.2f}, z_on={best_overall['z_on']} ---")

best_results = []
for src in sources_data:
    z_phot = src['z_phot']
    tid = src['target_id']
    
    fit_n = fit_two_param(src['rest_wave'], src['nuFnu_norm'], src['nuFnu_err_norm'], 0.0)
    zd = z_dist_from_phot(z_phot, gamma=gamma_best, **wp_best)
    fit_g = fit_two_param(src['rest_wave'], src['nuFnu_norm'], src['nuFnu_err_norm'], zd)
    
    if fit_n['success'] and fit_g['success']:
        dc2 = fit_n['red_chi2'] - fit_g['red_chi2']
        dbic = fit_n['bic'] - fit_g['bic']
        verdict = "Support" if (dc2 > 1 and dbic < -2) else ("Weak" if dc2 > 0 else ("Neutral" if abs(dc2)<=1 else "Null wins"))
        print(f"  ID={tid} (z={z_phot:.2f}): z_dist={zd:.3f}, Δχ²ν={dc2:+.1f}, ΔBIC={dbic:+.1f} → {verdict}")
        best_results.append({**src, **{
            'zd': zd, 'fit_n': fit_n, 'fit_g': fit_g,
            'dc2': dc2, 'dbic': dbic, 'verdict': verdict
        }})

# 绘制最佳参数拟合图
if len(best_results) > 0:
    fig_best = plt.figure(figsize=(20, 16))
    gs = fig_best.add_gridspec(4, 3, height_ratios=[3, 1, 3, 1], hspace=0.28, wspace=0.22)
    
    vcolors = {'Support': '#2ca02c', 'Weak': '#98df8a', 'Neutral': '#ff7f0e', 'Null wins': '#d62728'}
    
    for idx, br in enumerate(best_results):
        row_sed = 0 if idx < 3 else 2
        row_res = 1 if idx < 3 else 3
        col = idx % 3
        
        ax_s = fig_best.add_subplot(gs[row_sed, col])
        ax_r = fig_best.add_subplot(gs[row_res, col], sharex=ax_s)
        
        wave_sm = np.logspace(-1, 1.2, 300)
        m_null = sed_model(wave_sm, br['fit_n']['C_blue'], br['fit_n']['C_red'], 0.0)
        m_geff = sed_model(wave_sm, br['fit_g']['C_blue'], br['fit_g']['C_red'], br['zd'])
        
        vc = vcolors.get(br['verdict'], 'gray')
        
        ax_s.plot(wave_sm, m_geff, color='#2ca02c', linewidth=2.5, label=r'$G_{\rm eff}$ Model (Best Params)')
        ax_s.plot(wave_sm, m_null, color='#1f77b4', linewidth=1.8, linestyle='--', label=r'Null ($z_{\rm dist}=0$)')
        ax_s.errorbar(br['rest_wave'], br['nuFnu_norm'], yerr=br['nuFnu_err_norm'],
                      fmt='o', color='black', markersize=6, capsize=2, zorder=10)
        
        ax_s.set_xscale('log'); ax_s.set_xlim(0.1, 15)
        ymax = max(max(m_geff), max(m_null)) * 1.25
        ax_s.set_ylim(0, ymax); plt.setp(ax_s.get_xticklabels(), visible=False)
        if col == 0: ax_s.set_ylabel(r'Norm. $\nu F_\nu$', fontsize=13)
        
        info_t = (
            r"ID: {:d} ($z$={:.2f})".format(br['target_id'], br['z_phot']) + "\n"
            + r"$z_{{dist}}$ = {:.3f}".format(br['zd']) + "\n"
            + r"$\chi^2_\nu$(N) = {:.1f}".format(br['fit_n']['red_chi2']) + "\n"
            + r"$\chi^2_\nu$(G) = {:.1f}".format(br['fit_g']['red_chi2']) + "\n"
            + r"$\Delta\chi^2_\nu$ = {:+.1f}".format(br['dc2']) + "\n"
            + r"$\Delta$BIC = {:+.1f}".format(br['dbic'])
        )
        ax_s.text(0.06, ymax*0.78, info_t, fontsize=9,
                  bbox=dict(facecolor='white', alpha=0.88, edgecolor=vc, linewidth=1.5),
                  family='monospace')
        if idx == 0: ax_s.legend(loc='upper right', fontsize=9, framealpha=0.85)
        
        # Residual
        ax_r.axhline(0, color='black', lw=1, alpha=0.5)
        res_g = br['fit_g']['residuals']
        res_n = (br['nuFnu_norm'] - br['fit_n']['model_flux']) / br['nuFnu_err_norm']
        ax_r.errorbar(br['rest_wave'], res_g, yerr=1.0, fmt='o', color='#2ca02c', ms=6, label='$G_{eff}$')
        ax_r.errorbar(br['rest_wave'], res_n, yerr=1.0, fmt='s', color='#1f77b4', ms=5, alpha=0.6, label='Null')
        ax_r.set_ylim(-4, 4); ax_r.set_xticks([0.1, 0.5, 1, 5, 10]); ax_r.set_xticklabels(['0.1','0.5','1','5','10'])
        if row_res == 3: ax_r.set_xlabel(r'Rest-frame $\lambda$ ($\mu$m)', fontsize=13)
        if col == 0: ax_r.set_ylabel(r'$\Delta \chi$', fontsize=13)
        if idx == len(best_results)-1 or idx == 2: ax_r.legend(loc='upper right', fontsize=8, ncol=2)
    
    fig_best.suptitle(
        r'Best-Fit Parameters: $A_G=$' + f"{best_overall['A_G']:.3f}" + r', $\gamma$=' + f"{gamma_best:.2f}"
        + r', $\beta$=' + f"{best_overall['beta']:.2f}" + r', $z_{\rm on}=$' + f"{best_overall['z_on']:.1f}"
        + '\nLRD SED — G_eff Window Model vs Null Hypothesis',
        fontsize=13, y=0.98
    )
    fig_best.savefig('ParamScan_BestFit.png', dpi=300, bbox_inches='tight')
    print("✅ Saved: ParamScan_BestFit.png")

# ==========================================
# 8. 输出 Top-10 参数组合 CSV
# ==========================================
all_results = results_fine  # fine scan 已经包含了所有组合

# 按 score 排序取 Top 20
sorted_results = sorted(all_results, key=lambda x: x['score'], reverse=True)[:20]

rows_out = []
for rank, r in enumerate(sorted_results, 1):
    row = {
        'rank': rank,
        'A_G': r['A_G'], 'gamma': r['gamma'], 'beta': r['beta'], 'z_on': r['z_on'],
        'mean_delta_chi2': round(r['mean_delta_chi2'], 2),
        'min_delta_chi2': round(r['min_delta_chi2'], 2),
        'n_support': r['n_support'], 'n_nullwin': r['n_nullwin'],
        'score': round(r['score'], 1),
        'per_source_dc2': ';'.join([f'{x:.1f}' for x in r['dc2_per_source']])
    }
    rows_out.append(row)

df_top = pd.DataFrame(rows_out)
df_top.to_csv('ParamScan_TopParams.csv', index=False)
print("✅ Saved: ParamScan_TopParams.csv")

# ==========================================
# 9. 打印汇总报告
# ==========================================
print("\n" + "=" * 75)
print("PARAMETER SCAN SUMMARY REPORT")
print("=" * 75)
print(f"\n  FIDUCIAL (original):  A_G=0.15, gamma=3.0, beta=1.0, z_on=3.0")
print(f"  BEST     (scanned):  A_G={best_overall['A_G']:.3f}, gamma={best_overall['gamma']:.2f}, "
      f"beta={best_overall['beta']:.2f}, z_on={best_overall['z_on']:.1f}")
print(f"\n  Best Score = {best_overall['score']:.1f}")
print(f"  Mean Δχ²_ν = {best_overall['mean_delta_chi2']:+.1f}")
print(f"  Min Δχ²_ν  = {best_overall['min_delta_chi2']:+.1f}")
print(f"  Support: {best_overall['n_support']}/{n_sources}  |  Null wins: {best_overall['n_nullwin']}/{n_sources}")
print(f"\n--- TOP 5 Parameter Combinations ---")
print(df_top[['rank','A_G','gamma','beta','z_on','mean_delta_chi2','min_delta_chi2','n_support','n_nullwin']].head().to_string(index=False))
print("\n" + "=" * 75)
print("✅ Parameter scan complete!")
