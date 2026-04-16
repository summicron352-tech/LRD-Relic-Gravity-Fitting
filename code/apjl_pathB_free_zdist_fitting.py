#!/usr/bin/env python3
"""
=============================================================================
  LRD 路径B：自由参数 z_dist SED 拟合
  
  核心思想：不再对所有源使用固定 z_dist=0.17，
  而是对每个源让 z_dist 作为自由拟合参数，
  看最优 z_dist 是否与面密度 Σ 相关。
  
  物理预期（UID）：
    高 Σ → 深势阱 → 高最优 z_dist
    低 Σ → 浅势阱 → 低最优 z_dist
    
=============================================================================
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize_scalar
from scipy.stats import spearmanr, norm as sp_norm, ks_2samp as sp_ks2samp, rankdata
from numpy.linalg import lstsq
from math import sqrt
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── 全局设置 ──
DATA_DIR = '/Users/tanxin/Desktop/数据处理/LRD_GitHub_Release_v2/data/csv'
OUTPUT_DIR = '/Users/tanxin/Desktop/数据处理/LRD_GitHub_Release_v2/figures/density_zdist_analysis'
import os
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams.update({
    'font.size': 11,
    'axes.unicode_minus': False,
    'figure.facecolor': 'white',
})

print("=" * 72)
print("  LRD 路径B：自由参数 z_dist SED 拟合")
print("=" * 72)

# ═══════════════════════════════════════════════════════════════════
# 1. SED 物理模型（与 param_scan_geff.py 完全一致）
# ═══════════════════════════════════════════════════════════════════
def f_uv(wave_rest):
    """UV 连续谱: 年轻恒星的幂律谱"""
    return (wave_rest / 0.3)**(-0.5)


def f_int(wave_rest):
    """光学/近红外连续谱: 老年恒星 + Balmer break"""
    break_wave = 0.364
    flux = np.zeros_like(wave_rest)
    mask = wave_rest > break_wave
    flux[~mask] = 0.1 * (wave_rest[~mask] / break_wave)**2
    flux[mask] = ((wave_rest[mask] / break_wave)**0.5 
                 * np.exp(-(wave_rest[mask] - 1.0)**2 / 2.0) 
                 + 0.5 * (1.0 / wave_rest[mask]))
    return flux


def sed_model_3param(wave_rest, C_blue, C_red, z_dist):
    """三参数 SED 模型：C_blue(UV), C_red(光红), z_dist(引力红移)"""
    # 引力红移效应：rest-frame 光学成分被红移
    distorted_wave = wave_rest / (1 + z_dist)
    # 总通量 = UV分量 + 红移后的光红分量（光度守恒除以1+z）
    return C_blue * f_uv(wave_rest) + C_red * f_int(distorted_wave) / (1 + z_dist)


# ═══════════════════════════════════════════════════════════════════
# 2. 数据加载
# ═══════════════════════════════════════════════════════════════════
df = pd.read_csv(f'{DATA_DIR}/Kokorev_LRDs_Full.csv')
print(f"\n[1] 原始数据: {len(df)} 个源")

FILTERS = {
    'f090w': 0.90, 'f115w': 1.15, 'f150w': 1.50, 'f200w': 2.00,
    'f277w': 2.77, 'f356w': 3.56, 'f444w': 4.44, 'f770w': 7.70,
}

# ═══════════════════════════════════════════════════════════════════
# 3. 对每个源做三参数拟合 (C_blue, C_red, z_dist)
# ═══════════════════════════════════════════════════════════════════
print("\n[2] 开始自由参数 SED 拟合 (z_dist ∈ [0, 0.50])...")

results = []

for idx, row in df.iterrows():
    target_id = row['id']
    z_phot = row['z_phot']
    
    # 提取观测数据
    obs_wave, obs_flux, obs_err = [], [], []
    for fname, wave_um in FILTERS.items():
        fcol = f'{fname}_flux'
        ecol = f'{fname}_fluxerr'
        if fcol in df.columns and pd.notna(row[fcol]) and row[fcol] > -90:
            obs_wave.append(wave_um)
            obs_flux.append(row[fcol])
            err_val = row[ecol] if pd.notna(row.get(ecol, np.nan)) and row.get(ecol, 1) > 0 else abs(row[fcol]) * 0.1
            obs_err.append(max(err_val, abs(row[fcol]) * 0.01))  # 最小误差 1%
    
    if len(obs_wave) < 5:  # 至少需要5个波段
        continue
    
    obs_wave = np.array(obs_wave)
    obs_flux = np.array(obs_flux)
    obs_err = np.array(obs_err)
    
    rest_wave = obs_wave / (1 + z_phot)
    nuFnu_obs = obs_flux / obs_wave
    nuFnu_err = obs_err / np.array(obs_wave)
    norm_factor = np.max(nuFnu_obs) if np.max(nuFnu_obs) > 0 else 1.0
    
    nuFnu_norm = nuFnu_obs / norm_factor
    nuFnu_err_norm = nuFnu_err / norm_factor
    
    try:
        # 三参数拟合: C_blue, C_red, z_dist
        popt, pcov = curve_fit(
            sed_model_3param, rest_wave, nuFnu_norm,
            sigma=nuFnu_err_norm,
            p0=[1.0, 1.0, 0.17],  # 初始猜测: z_dist 从论文值开始
            bounds=(
                [0,   0,   0],      # 下界: 无负值
                [10, 10,   0.50]     # 上界: z_dist 不超过 0.5
            ),
            maxfev=20000
        )
        
        C_blue_opt, C_red_opt, zdist_opt = popt
        
        # 计算拟合质量
        n_data = len(rest_wave)
        model_flux = sed_model_3param(rest_wave, C_blue_opt, C_red_opt, zdist_opt)
        residuals = (nuFnu_norm - model_flux) / nuFnu_err_norm
        chi2 = np.sum(residuals**2)
        dof = n_data - 3  # 3个自由参数
        red_chi2 = chi2 / dof if dof > 0 else chi2
        
        # Null model (z_dist=0, 两参数拟合) 用于对比
        from scipy.optimize import curve_fit as cf2
        def sed_null(w, cb, cr):
            return cb * f_uv(w) + cr * f_int(w)
        
        pnull, _ = cf2(sed_null, rest_wave, nuFnu_norm, sigma=nuFnu_err_norm,
                        bounds=([0, 0], [10, 10]), maxfev=10000)
        null_model = sed_null(rest_wave, pnull[0], pnull[1])
        null_residuals = (nuFnu_norm - null_model) / nuFnu_err_norm
        null_chi2 = np.sum(null_residuals**2)
        null_dof = n_data - 2
        null_red_chi2 = null_chi2 / null_dof if null_dof > 0 else null_chi2
        
        delta_chi2 = null_red_chi2 - red_chi2
        
        results.append({
            'id': target_id,
            'field': row.get('field', ''),
            'z_phot': z_phot,
            'zdist_opt': zdist_opt,
            'zdist_err': sqrt(pcov[2, 2]) if pcov[2, 2] > 0 else np.nan,
            'C_blue': C_blue_opt,
            'C_red': C_red_opt,
            'red_chi2_zfree': red_chi2,
            'red_chi2_null': null_red_chi2,
            'delta_chi2_nu': delta_chi2,
            'N_bands': len(obs_wave),
            'lbol': row.get('lbol', np.nan),
            'R_eff_pc': row.get('r_eff_50_phys', np.nan),
        })
        
    except Exception as e:
        # 拟合失败，记录 NaN
        results.append({
            'id': target_id, 'z_phot': z_phot,
            'zdist_opt': np.nan, 'zdist_err': np.nan,
            'C_blue': np.nan, 'C_red': np.nan,
            'red_chi2_zfree': np.nan, 'red_chi2_null': np.nan,
            'delta_chi2_nu': np.nan, 'N_bands': len(obs_wave),
            'lbol': row.get('lbol', np.nan),
            'R_eff_pc': row.get('r_eff_50_phys', np.nan),
        })

df_results = pd.DataFrame(results)
print(f"    成功拟合: {df_results['zdist_opt'].notna().sum()} / {len(df)}")

# ═══════════════════════════════════════════════════════════════════
# 4. 关键检验：z_dist_opt vs 密度
# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("  核心：z_dist(free) vs 面密度 Σ")
print("=" * 72)

dfv = df_results[df_results['zdist_opt'].notna() & 
                  df_results['lbol'].notna() & 
                  df_results['R_eff_pc'].notna() &
                  (df_results['R_eff_pc'] > 0)].copy()

print(f"[3] 有效源数（有 z_dist_opt + Σ）: {len(dfv)}")

# 计算面密度和紧致度
dfv['log_Lbol'] = np.log10(dfv['lbol'])
dfv['log_Reff'] = np.log10(dfv['R_eff_pc'])
dfv['log_Sigma'] = np.log10(dfv['lbol']) - 2 * np.log10(dfv['R_eff_pc'])
dfv['log_C'] = np.log10(dfv['lbol']) - np.log10(dfv['R_eff_pc'])

print(f"\n    z_dist_opt 分布:")
print(f"      均值 = {dfv['zdist_opt'].mean():.4f}")
print(f"      中位数 = {dfv['zdist_opt'].median():.4f}")
print(f"      标准差 = {dfv['zdist_opt'].std():.4f}")
print(f"      范围 = [{dfv['zdist_opt'].min():.4f}, {dfv['zdist_opt'].max():.4f}]")

# ── 原始相关 ──
pairs_b = [
    ('log_Sigma', 'zdist_opt', 'Sigma vs z_dist(opt)'),
    ('log_C',     'zdist_opt', 'Compactness vs z_dist(opt)'),
    ('log_Reff',  'zdist_opt', 'R_eff vs z_dist(opt)'),
]

print("\n  ── 原始 Spearman 相关 ──")
raw_results_b = []
for x_col, y_col, label in pairs_b:
    rho, p_val = spearmanr(dfv[x_col], dfv[y_col])
    sig = abs(sp_norm.ppf(p_val/2)) if p_val > 0 else float('inf')
    raw_results_b.append({'pair': label, 'rho': rho, 'p': p_val, 'sig': sig})
    m = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
    print(f"  {label:30s}: rho = {rho:+.3f} ({sig:+.1f}σ){m}")

# ── 偏相关 ──
print("\n  ── 偏相关 (控制 z_phot + log_Lbol) ──")
partial_b = []

def partial_spearman(x, y, ctrl_cols, data):
    rx = rankdata(x)
    ry = rankdata(y)
    X_ctrl = np.column_stack([np.ones(len(data))] + [data[c].values for c in ctrl_cols])
    # 添加微小扰动避免奇异矩阵
    rng = np.random.default_rng(42)
    X_ctrl = X_ctrl + rng.normal(0, 1e-10, X_ctrl.shape)
    # 用伪逆处理可能的奇异矩阵
    try:
        pinv_X = np.linalg.pinv(X_ctrl, rcond=1e-6)
        resid_x = rx - X_ctrl @ (pinv_X @ rx)
        resid_y = ry - X_ctrl @ (pinv_X @ ry)
        return spearmanr(resid_x, resid_y)
    except np.linalg.LinAlgError:
        # 如果仍然失败，返回 NaN
        return (np.nan, np.nan)

controls_b = ['z_phot', 'log_Lbol']
for x_col, y_col, label in pairs_b[:2]:
    rp, pp = partial_spearman(dfv[x_col].values, dfv[y_col].values, controls_b, dfv)
    if np.isnan(rp):
        sp = 0
        m = "N/A"
        print(f"  {label:30s}: rho_p = N/A (singular matrix)")
    else:
        sp = abs(sp_norm.ppf(pp/2)) if pp > 0 else float('inf')
        m = "***" if pp < 0.001 else "**" if pp < 0.01 else "*" if pp < 0.05 else ""
        print(f"  {label:30s}: rho_p = {rp:+.3f} ({sp:+.1f}σ){m}")
    partial_b.append({'pair': label, 'rho_p': rp if not np.isnan(rp) else 0, 'p_p': pp if not np.isnan(pp) else 1, 'sig_p': sp})

# ── 分组统计 ──
print("\n  ── 按 Σ 四分位数分组的 z_dist(opt) ──")
dfv['Sigma_q'] = pd.qcut(dfv['log_Sigma'], 4, labels=['Q1(Low)', 'Q2', 'Q3', 'Q4(High)'])

group_zd = dfv.groupby('Sigma_q', observed=True).agg(
    zdist_mean=('zdist_opt', 'mean'),
    zdist_std=('zdist_opt', 'std'),
    zdist_median=('zdist_opt', 'median'),
    Sigma_mean=('log_Sigma', 'mean'),
    N=('id', 'count')
).round(4)
print(group_zd.to_string())

# KS 检验
q1_zd = dfv[dfv['Sigma_q'] == 'Q1(Low)']['zdist_opt'].values
q4_zd = dfv[dfv['Sigma_q'] == 'Q4(High)']['zdist_opt'].values
ks_zd, ks_pval_zd = sp_ks2samp(q1_zd, q4_zd)
print(f"\n  KS(Q1 vs Q4): D={ks_zd:.3f}, p={ks_pval_zd:.2e}")

# 与论文固定值 0.17 的比较
above_17 = (dfv['zdist_opt'] > 0.17).sum()
below_17 = (dfv['zdist_opt'] <= 0.17).sum()
print(f"\n  z_dist(opt) vs 固定值 0.17:")
print(f"    > 0.17: {above_17} ({100*above_17/len(dfv):.1f}%)")
print(f"    ≤ 0.17: {below_17} ({100*below_17/len(dfv):.1f}%)")

# ═══════════════════════════════════════════════════════════════════
# 5. 绘图
# ═══════════════════════════════════════════════════════════════════
print("\n[4] 生成图表...")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Panel A: z_dist(opt) vs Σ — 核心散点图！
ax = axes[0, 0]
sc = ax.scatter(dfv['log_Sigma'], dfv['zdist_opt'],
                c=dfv['z_phot'], cmap='plasma', s=28, alpha=0.75, edgecolors='none')

# 趋势线（用稳健拟合避免数值问题）
try:
    zf = np.polyfit(dfv['log_Sigma'], dfv['zdist_opt'], 1)
    slope_label = f'{zf[0]:+.4f}'
except np.linalg.LinAlgError:
    # fallback: 简单线性回归
    x_mean = dfv['log_Sigma'].mean()
    y_mean = dfv['zdist_opt'].mean()
    zf = [((dfv['zdist_opt'] * (dfv['log_Sigma'] - x_mean)).sum() / 
          ((dfv['log_Sigma'] - x_mean)**2).sum()),
         y_mean - ((dfv['zdist_opt'] * (dfv['log_Sigma'] - x_mean)).sum() / 
           ((dfv['log_Sigma'] - x_mean)**2).sum()) * x_mean]
    slope_label = f'{zf[0]:+.4f}'

xt = np.linspace(dfv['log_Sigma'].min(), dfv['log_Sigma'].max(), 100)
ax.plot(xt, np.poly1d(zf)(xt), 'r--', lw=2.5, alpha=0.85,
        label=f'slope = {slope_label}')

# 参考线: 论文固定值 0.17
ax.axhline(y=0.17, color='blue', ls=':', lw=1.8, alpha=0.7, label='Fixed value (paper) = 0.17')

rho_a, pa = spearmanr(dfv['log_Sigma'], dfv['zdist_opt'])
sa = abs(sp_norm.ppf(pa/2)) if pa > 0 else 99

ax.set_xlabel(r'Compactness $\log_{10}(L_{\rm bol}/R^2)$', fontsize=11)
ax.set_ylabel(r'Optimal $z_{\rm dist}$ (free parameter fit)', fontsize=11)
ax.set_title(f'Panel A: Free-$z_{{dist}}$ vs Compactness\n$\\rho$={rho_a:.3f} ({sa:.1f}$\\sigma$) N={len(dfv)}',
              fontsize=12, fontweight='bold')
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('$z_{phot}$')
ax.legend(loc='lower right', fontsize=9)

# 标注 UID 预测方向
if rho_a > 0:
    ax.annotate('UID Prediction:\nHigher $\\Sigma$ → Higher $z_{dist}$',
                xy=(0.97, 0.82), xycoords='axes fraction',
                fontsize=9, ha='right', va='top', color='darkgreen',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.6))

# Panel B: z_dist(opt) 直方图 + 分组分布
ax = axes[0, 1]

# 整体直方图（背景）
bins_global = np.linspace(0, 0.5, 40)
ax.hist(dfv['zdist_opt'], bins=bins_global, alpha=0.25, color='gray', density=True, 
         label='All sources')

# 分组直方图
order_q = ['Q1(Low)', 'Q2', 'Q3', 'Q4(High)']
colors_q = ['#3498db', '#2ecc71', '#f39c12', '#e74c3c']

for q, col in zip(order_q, colors_q):
    sub = dfv[dfv['Sigma_q'] == q]['zdist_opt']
    ax.hist(sub, bins=bins_global, alpha=0.55, color=col, density=True, label=q)

# 论文固定值参考
ax.axvline(x=0.17, color='navy', ls='--', lw=2.5, label='Paper fixed = 0.17')
ax.set_xlabel(r'Optimal $z_{\rm dist}$', fontsize=11)
ax.set_ylabel('Normalized Density', fontsize=11)
ax.set_title('Panel B: $z_{dist}$ Distribution Stratified by Density\n(UID: high-$\\Sigma$ should shift right)',
              fontsize=12, fontweight='bold')
ax.legend(loc='upper right', fontsize=9)
ax.set_xlim(0, 0.5)

# Panel C: Delta-chi2 vs z_dist(opt)
ax = axes[1, 0]
sc_c = ax.scatter(dfv['zdist_opt'], dfv['delta_chi2_nu'],
                  c=dfv['log_Sigma'], cmap='coolwarm_r', s=28, alpha=0.75, edgecolors='none')
ax.axhline(y=0, color='black', ls='-', lw=1, alpha=0.5)
ax.axvline(x=0.17, color='blue', ls=':', lw=1.5, alpha=0.6, label='Fixed = 0.17')

# 标注哪些区域改善最大
best_mask = dfv['delta_chi2_nu'] > dfv['delta_chi2_nu'].quantile(0.9)
if best_mask.sum() > 0:
    ax.scatter(dfv.loc[best_mask, 'zdist_opt'], dfv.loc[best_mask, 'delta_chi2_nu'],
               facecolors='none', edgecolors='gold', s=80, lw=2, label='Top 10% improvers')

ax.set_xlabel(r'Optimal $z_{\rm dist}$', fontsize=11)
ax.set_ylabel(r'$\Delta\chi^2_\nu$ (Null − Free)', fontsize=11)
ax.set_title("Panel C: Fit Improvement vs $z_{dist}$\n(Higher = Relic Model Wins)",
              fontsize=12, fontweight='bold')
cbar_c = plt.colorbar(sc_c, ax=ax)
cbar_c.set_label(r'$\log(\Sigma)$')
ax.legend(loc='upper right', fontsize=9)

# Panel D: 综合——z_dist(opt) 的 Σ 依赖性总结
ax = axes[1, 1]

# 误差棒柱状图: 每组的 mean ± SEM
means_zd = []
sems_zd = []
labels_bar = []

for q in order_q:
    sub = dfv[dfv['Sigma_q'] == q]['zdist_opt']
    means_zd.append(sub.mean())
    sems_zd.append(sub.std() / sqrt(len(sub)))
    labels_bar.append(q.replace('(', '\n('))

x_pos = np.arange(len(order_q))
bars = ax.bar(x_pos, means_zd, yerr=sems_zd, color=colors_q, alpha=0.75,
             edgecolor='black', linewidth=1.5, capsize=5)

# 参考线
ax.axhline(y=dfv['zdist_opt'].mean(), color='purple', ls='--', lw=2, alpha=0.7,
           label=f'Global mean = {dfv["zdist_opt"].mean():.3f}')
ax.axhline(y=0.17, color='blue', ls=':', lw=2, alpha=0.7, label=f'Paper fixed = 0.170')

ax.set_xticks(x_pos)
ax.set_xticklabels(labels_bar, fontsize=9)
ax.set_ylabel(r'Mean Optimal $z_{\rm dist}$', fontsize=11)
ax.set_xlabel('Surface Density Quartile', fontsize=11)
ax.set_title(f'Panel D: Mean $z_{{dist}}$(opt) by Compactness\nKS(Q1vsQ4): D={ks_zd:.3f}, p={ks_pval_zd:.1e}',
              fontsize=12, fontweight='bold')
ax.legend(loc='upper left', fontsize=9)
ax.set_ylim(0, max(means_zd) * 1.4)

# 在柱子上标注数值
for i, (m, s) in enumerate(zip(means_zd, sems_zd)):
    ax.text(i, m + s + 0.005, f'{m:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

figpath_b = f'{OUTPUT_DIR}/Figure_PathB_FreeZdist_vs_Density.png'
plt.savefig(figpath_b, dpi=250, facecolor='white')
print(f"  ✅ 图表保存: {figpath_b}")
plt.close()

# ═══════════════════════════════════════════════════════════════════
# 6. 输出汇总
# ═══════════════════════════════════════════════════════════════════
summary_b = f"""
╔══════════════════════════════════════════════════════════════╗
║     LRD 路径B：自由参数 z_dist SED 拟合结果                ║
╠══════════════════════════════════════════════════════════════╣
║ 有效样本: N = {len(dfv)}                                         ║
╠══════════════════════════════════════════════════════════════╣
║                                                              ║
║ 【z_dist(opt) 统计】                                          ║
║   均值   = {dfv['zdist_opt'].mean():.4f}                                     ║
║   中位数 = {dfv['zdist_opt'].median():.4f}                                    ║
║   标准差 = {dfv['zdist_opt'].std():.4f}                                     ║
║   范围   = [{dfv['zdist_opt'].min():.4f}, {dfv['zdist_opt'].max():.4f}]                          ║
║                                                              ║
║ > 0.17: {above_17} ({100*above_17/len(dfv):.1f}%)  |  <= 0.17: {below_17} ({100*below_17/len(dfv):.1f}%)          ║
║                                                              ║
║ 【原始 Spearman 相关】                                        ║
"""

for r in raw_results_b:
    m = "***" if r['p'] < 0.001 else "**" if r['p'] < 0.01 else "*" if r['p'] < 0.05 else ""
    summary_b += f"║   {r['pair']:28s} ρ = {r['rho']:+.3f} ({r['sig']:5.1f}σ) {m}\n"

summary_b += f"""║                                                              ║
║ 【偏相关 (控制 z_phot + L_bol)】                             ║
"""

for r in partial_b:
    m = "***" if r['p_p'] < 0.001 else "**" if r['p_p'] < 0.01 else "*" if r['p_p'] < 0.05 else ""
    summary_b += f"║   {r['pair']:28s} ρ_p = {r['rho_p']:+.3f} ({r['sig_p']:5.1f}σ) {m}\n"

summary_b += f"""║                                                              ║
║ 【按 Σ 分组均值 z_dist(opt)】                                 ║
"""
for q in order_q:
    row = group_zd.loc[q] if q in group_zd.index else None
    if row is not None:
        summary_b += f"║   {q:12s}: {row['zdist_mean']:.4f} +/- {row['zdist_std']:.4f} (N={int(row['N']):3.0f})\n"

summary_b += f"""║                                                              ║
║ KS(Q1 vs Q4): D = {ks_zd:.3f}, p = {ks_pval_zd:.2e}                      ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
"""

print(summary_b)

with open(f'{OUTPUT_DIR}/PathB_results_summary.txt', 'w') as f:
    f.write(summary_b)

# 保存完整结果 CSV
csv_path = f'{OUTPUT_DIR}/PathB_FreeZdist_AllSources.csv'
dfv.to_csv(csv_path, index=False)
print(f"\n✅ 全部结果已保存到: {OUTPUT_DIR}/")
print(f"  • Figure_PathB_FreeZdist_vs_Density.png")
print(f"  • PathB_FreeZdist_AllSources.csv")
