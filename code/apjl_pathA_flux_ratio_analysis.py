#!/usr/bin/env python3
"""
=============================================================================
  LRD 路径A：通量比代理分析 — 密度依赖引力红移检验
  
  核心物理：高 Σ → 深势阱 → 高 z_dist → SED 整体红移 → 长波通量相对增强
=============================================================================
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, rankdata
from numpy.linalg import lstsq
from math import sqrt
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')  # 非交互后端，避免GUI问题
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
print("  LRD 密度依赖引力红移 — 路径A：通量比代理分析")
print("=" * 72)

# ════════════════════════════════════════════════
# 1. 数据加载与清洗
# ════════════════════════════════════════════════
df = pd.read_csv(f'{DATA_DIR}/Kokorev_LRDs_Full.csv')
print(f"\n[1] 原始数据: {len(df)} 个源")

# 关键列名（已确认）
FLUX_KEYS = ['f150w_flux', 'f277w_flux', 'f356w_flux', 'f444w_flux']

mask = pd.Series(True, index=df.index)
for col in FLUX_KEYS + ['z_phot', 'lbol', 'r_eff_50_phys']:
    mask &= df[col].notna() & (df[col] > 0)

dfv = df[mask].copy()
print(f"[2] 有效源数: {len(dfv)}")

# ════════════════════════════════════════════════
# 2. 物理量计算
# ════════════════════════════════════════════════
dfv['log_Lbol'] = np.log10(dfv['lbol'])
dfv['R_eff_pc'] = dfv['r_eff_50_phys']
dfv['log_Reff'] = np.log10(dfv['R_eff_pc'])

# 面密度代理：Σ ∝ L_bol / R²
dfv['log_Sigma'] = np.log10(dfv['lbol']) - 2 * np.log10(dfv['R_eff_pc'])

# 紧致度代理：C ∝ L_bol / R
dfv['log_C'] = np.log10(dfv['lbol']) - np.log10(dfv['R_eff_pc'])

# 通量比（引力红移的关键观测指标）
dfv['FR_444_356'] = dfv['f444w_flux'] / dfv['f356w_flux']
dfv['FR_444_150'] = dfv['f444w_flux'] / dfv['f150w_flux']
dfv['FR_356_150'] = dfv['f356w_flux'] / dfv['f150w_flux']

dfv['log_FR_444_356'] = np.log10(dfv['FR_444_356'])
dfv['log_FR_444_150'] = np.log10(dfv['FR_444_150'])
dfv['log_FR_356_150'] = np.log10(dfv['FR_356_150'])

print(f"[3] 物理量计算完成:")
print(f"    log(Σ) 范围: [{dfv['log_Sigma'].min():.2f}, {dfv['log_Sigma'].max():.2f}]")
print(f"    log(F444W/F356W): [{dfv['log_FR_444_356'].min():.2f}, {dfv['log_FR_444_356'].max():.2f}]")
print(f"    z_phot 范围: [{dfv['z_phot'].min():.2f}, {dfv['z_phot'].max():.2f}]")

# ════════════════════════════════════════════════
# 3. 统计检验
# ════════════════════════════════════════════════
print("\n" + "=" * 72)
print("  原始 Spearman 相关")
print("=" * 72)

pairs_raw = [
    ('log_Sigma', 'log_FR_444_356', 'Σ vs F444W/F356W'),
    ('log_Sigma', 'log_FR_444_150', 'Σ vs F444W/F150W'),
    ('log_C',     'log_FR_444_356', 'C vs F444W/F356W'),
    ('log_C',     'log_FR_444_150', 'C vs F444W/F150W'),
]

results_raw = []
for x_col, y_col, label in pairs_raw:
    rho, p_val = spearmanr(dfv[x_col], dfv[y_col])
    sig = abs(stats.norm.ppf(p_val / 2)) if p_val > 0 else float('inf')
    results_raw.append({'pair': label, 'rho': rho, 'p': p_val, 'sig': sig})
    m = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
    print(f"  {label:25s}: ρ = {rho:+.3f} ({sig:+.1f}σ){m}")

# ── 偏相关：控制 z_phot + L_bol ──
print("\n" + "=" * 72)
print("  偏相关 (控制 z_phot + L_bol)")
print("=" * 72)

def partial_spearman(x, y, ctrl_cols, data):
    """Spearman偏相关：对控制变量回归后取残差再算相关"""
    rx = rankdata(x)
    ry = rankdata(y)
    X_ctrl = np.column_stack([np.ones(len(data))] + [data[c].values for c in ctrl_cols])
    
    resid_x = rx - X_ctrl @ lstsq(X_ctrl, rx, rcond=None)[0]
    resid_y = ry - X_ctrl @ lstsq(X_ctrl, ry, rcond=None)[0]
    
    return spearmanr(resid_x, resid_y)

controls = ['z_phot', 'log_Lbol']
results_partial = []

for x_col, y_col, label in pairs_raw[:4]:
    rho_p, p_p = partial_spearman(
        dfv[x_col].values, dfv[y_col].values, controls, dfv
    )
    sig_p = abs(stats.norm.ppf(p_p / 2)) if p_p > 0 else float('inf')
    results_partial.append({'pair': label, 'rho_p': rho_p, 'p_p': p_p, 'sig_p': sig_p})
    m = "***" if p_p < 0.001 else "**" if p_p < 0.01 else "*" if p_p < 0.05 else ""
    print(f"  {label:25s}: ρ_p = {rho_p:+.3f} ({sig_p:+.1f}σ){m}")

# ── 分组统计 ──
print("\n" + "=" * 72)
print("  按 Σ 四分位数分组")
print("=" * 72)

dfv['Sigma_q'] = pd.qcut(dfv['log_Sigma'], 4, labels=['Q1(最稀疏)', 'Q2', 'Q3', 'Q4(最致密)'])

group_means = dfv.groupby('Sigma_q', observed=True).agg(
    FR_444_356_mean=('log_FR_444_356', 'mean'),
    FR_444_150_mean=('log_FR_444_150', 'mean'),
    z_mean=('z_phot', 'mean'),
    Sigma_mean=('log_Sigma', 'mean'),
    N=('log_FR_444_356', 'count')
).round(3)
print(group_means.to_string())

# KS 检验 Q1 vs Q4
q1_data = dfv[dfv['Sigma_q'] == 'Q1(最稀疏)']['log_FR_444_356'].values
q4_data = dfv[dfv['Sigma_q'] == 'Q4(最致密)']['log_FR_444_356'].values
ks_d, ks_pval = stats.ks_2samp(q1_data, q4_data)
print(f"\n  KS(Q1 vs Q4): D={ks_d:.3f}, p={ks_pval:.2e}")
if ks_pval < 0.05:
    print(f"  ⭐ Q1 和 Q4 的 F444W/F356W 分布显著不同！")
else:
    print(f"  分布差异不显著——但注意当前所有源共用固定 z_dist=0.17")

# IR Excess 分层
dfv['IRExcess'] = dfv['f444w_flux'] / dfv['f150w_flux']
ire_low = dfv[dfv['IRExcess'] <= 5]['log_Sigma']
ire_high = dfv[dfv['IRExcess'] > 6.67]['log_Sigma']
ks_ire, p_ire = stats.ks_2samp(ire_low, ire_high)
print(f"\n  IR Excess KS (Low≤5 vs High>6.67): D={ks_ire:.3f}, p={p_ire:.3f}")

# ════════════════════════════════════════════════
# 4. 绘图（4面板）
# ════════════════════════════════════════════════
print("\n[4] 生成图表...")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Panel A: F444W/F356W vs Σ, color=z_phot
ax = axes[0, 0]
sc = ax.scatter(dfv['log_Sigma'], dfv['log_FR_444_356'], 
                c=dfv['z_phot'], cmap='plasma', s=25, alpha=0.75, edgecolors='none')

# 趋势线
z_fit = np.polyfit(dfv['log_Sigma'], dfv['log_FR_444_356'], 1)
x_trend = np.linspace(dfv['log_Sigma'].min(), dfv['log_Sigma'].max(), 100)
ax.plot(x_trend, np.poly1d(z_fit)(x_trend), 'r--', lw=2, alpha=0.8,
        label=f'slope={z_fit[0]:+.3f}')

rho_a, pa = spearmanr(dfv['log_Sigma'], dfv['log_FR_444_356'])
sa = abs(stats.norm.ppf(pa/2)) if pa > 0 else 99

ax.set_xlabel(r'Compactness Proxy $\log_{10}(L_{\rm bol}/R^2)$', fontsize=11)
ax.set_ylabel(r'$\log_{10}(F_{\rm 444W}/F_{\rm 356W})$', fontsize=11)
ax.set_title(f'Panel A: IR Excess vs Compactness\n$\\rho$={rho_a:.3f} ({sa:.1f}$\\sigma$) N={len(dfv)}',
              fontsize=12, fontweight='bold')
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('$z_{phot}$')
ax.legend(loc='lower right', fontsize=9)

# UID 预测箭头
if rho_a > 0:
    ax.annotate('UID Prediction:\nHigher $\\Sigma$ → Higher $z_{dist}$ → Higher FR',
                xy=(0.95, 0.85), xycoords='axes fraction',
                fontsize=9, ha='right', color='darkgreen',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.5))

# Panel B: 控制 z 后的残差 vs Σ
ax = axes[0, 1]
X_z = np.column_stack([np.ones(len(dfv)), dfv['z_phot'].values])
coef_z = lstsq(X_z, dfv['log_FR_444_356'].values, rcond=None)[0]
resid_FR = dfv['log_FR_444_356'].values - X_z @ coef_z

sc2 = ax.scatter(dfv['log_Sigma'], resid_FR, 
           c=dfv['log_Lbol'], cmap='coolwarm', s=25, alpha=0.75, edgecolors='none')

z2 = np.polyfit(dfv['log_Sigma'], resid_FR, 1)
ax.plot(x_trend, np.poly1d(z2)(x_trend), 'r--', lw=2, alpha=0.8)

rho_b, pb = spearmanr(dfv['log_Sigma'], resid_FR)
sb = abs(stats.norm.ppf(pb/2)) if pb > 0 else 99

ax.set_xlabel(r'Compactness Proxy $\log_{10}(L_{\rm bol}/R^2)$', fontsize=11)
ax.set_ylabel('Residual $\\log(F_{444W}/F_{356W})$ after $z_{phot}$ regression', fontsize=11)
ax.set_title(f'Panel B: Density Signal After Redshift Control\n$\\rho$={rho_b:.3f} ({sb:.1f}$\\sigma$)',
              fontsize=12, fontweight='bold')
cbar2 = plt.colorbar(sc2, ax=ax)
cbar2.set_label(r'$\log(L_{\rm bol})$')
ax.axhline(y=0, color='gray', ls=':', alpha=0.5)

# Panel C: 分组箱线图
ax = axes[1, 0]
order_labels = ['Q1(最稀疏)', 'Q2', 'Q3', 'Q4(最致密)']
bp_data = [dfv[dfv['Sigma_q'] == q]['log_FR_444_356'].values for q in order_labels]
colors_box = ['#3498db', '#2ecc71', '#f39c12', '#e74c3c']

bp = ax.boxplot(bp_data, patch_artist=True, notch=True)
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

means = [np.mean(d) for d in bp_data]
ax.scatter(range(1, 5), means, marker='D', c='black', s=50, zorder=5, label='Mean')
ax.plot(range(1, 5), means, 'k--', lw=1.5, alpha=0.7)

ax.set_xticklabels(['Q1\n(Low $\\Sigma$)', 'Q2', 'Q3', 'Q4\n(High $\\Sigma$)'], fontsize=10)
ax.set_ylabel(r'$\log_{10}(F_{444W}/F_{356W})$', fontsize=11)
ax.set_title(f'Panel C: Compactness Stratification\nKS(Q1vsQ4): D={ks_d:.3f}, p={ks_pval:.1e}',
              fontsize=12, fontweight='bold')
ax.legend(loc='upper left')

# Panel D: 多通量比综合柱状图
ax = axes[1, 1]
fr_pairs = [
    ('log_FR_444_356', r'$F_{444W}/F_{356W}$', '#e74c3c'),
    ('log_FR_444_150', r'$F_{444W}/F_{150W}$', '#3498db'),
    ('log_FR_356_150', r'$F_{356W}/F_{150W}$', '#2ecc71'),
]

for i, (col, name, color) in enumerate(fr_pairs):
    rt, pt = spearmanr(dfv['log_Sigma'], dfv[col])
    st = abs(stats.norm.ppf(pt/2)) if pt > 0 else 99
    
    x_pos = i + 1
    ax.bar(x_pos, rt, color=color, alpha=0.7, edgecolor='black', lw=1)
    
    m = "***" if pt < 0.001 else "**" if pt < 0.01 else "*" if pt < 0.05 else ""
    y_off = 0.02 if rt >= 0 else -0.06
    va_txt = 'bottom' if rt >= 0 else 'top'
    ax.text(x_pos, rt + y_off, m, ha='center', va=va_txt, fontsize=14, fontweight='bold')
    
    # 数值标注
    y_label_pos = rt - 0.04 if rt >= 0 else rt + 0.03
    va_str = 'bottom' if rt >= 0 else 'top'
    ax.text(x_pos, y_label_pos,
            f'{rt:.3f} ({st:.1f}sig)', ha='center', fontsize=8, color='black')

ax.axhline(y=0, color='gray', ls='-', lw=1)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels([r[1] for r in fr_pairs], fontsize=11)
ax.set_ylabel('Spearman rho (vs Sigma)', fontsize=11)
ax.set_title('Panel D: Multi-Band Flux Ratios Correlation with Density',
              fontsize=12, fontweight='bold')
ax.set_ylim(-0.15, 0.20)
ax.grid(axis='y', alpha=0.3)

figpath = f'{OUTPUT_DIR}/Figure_PathA_FluxRatio_vs_Density.png'
plt.savefig(figpath, dpi=250, facecolor='white')
print(f"  ✅ 图表保存: {figpath}")
plt.close()

# ════════════════════════════════════════════════
# 5. 输出汇总文本
# ════════════════════════════════════════════════
summary = f"""
╔══════════════════════════════════════════════════════════════╗
║      LRD 密度依赖引力红移 — 路径A 通量比代理分析结果         ║
╠══════════════════════════════════════════════════════════════╣
║ 样本: N = {len(dfv)} (原始260个源)                              ║
╠══════════════════════════════════════════════════════════════╣
║                                                              ║
║ 【原始 Spearman 相关】                                        ║
"""

for r in results_raw:
    m = "***" if r['p'] < 0.001 else "**" if r['p'] < 0.01 else "*" if r['p'] < 0.05 else ""
    summary += f"║   {r['pair']:23s} ρ = {r['rho']:+.3f} ({r['sig']:5.1f}σ) {m}\n"

summary += f"""║                                                              ║
║ 【偏相关 (控制 z_phot + L_bol)】                             ║
"""

for r in results_partial:
    m = "***" if r['p_p'] < 0.001 else "**" if r['p_p'] < 0.01 else "*" if r['p_p'] < 0.05 else ""
    summary += f"║   {r['pair']:23s} ρ_p = {r['rho_p']:+.3f} ({r['sig_p']:5.1f}σ) {m}\n"

summary += f"""║                                                              ║
║ 【按 Σ 分组均值 (log F444W/F356W)】                          ║
"""
for q in order_labels:
    row = group_means.loc[q] if q in group_means.index else None
    if row is not None:
        summary += f"║   {q:15s}: mean={row['FR_444_356_mean']:+.3f}, N={int(row['N']):3.0f}\n"

summary += f"""║                                                              ║
║ 【KS 检验】                                                   ║
║   Q1 vs Q4 (F444W/F356W): D = {ks_d:.3f}, p = {ks_pval:.2e}            ║
║   Low IRE vs High IRE (Σ): D = {ks_ire:.3f}, p = {p_ire:.3e}             ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
"""

print(summary)

with open(f'{OUTPUT_DIR}/PathA_results_summary.txt', 'w') as f:
    f.write(summary)

print(f"\n✅ 所有结果已保存到: {OUTPUT_DIR}/")
