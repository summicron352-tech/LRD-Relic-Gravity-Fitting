"""
Loiacono 5源统一 Dyson 应力测试
================================
BiRD + ID1017 + ID3646 + ID8511 + ID9008

对每个源: 扫 (M*, r_eff) → 计算 Dyson M_central
比较: M_central / M_BH(virial)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 11
import os

G_N = 4.30091e-3
output_dir = os.path.dirname(os.path.abspath(__file__))
fig_dir = os.path.join(os.path.dirname(output_dir), 'figures')
os.makedirs(fig_dir, exist_ok=True)

# Dyson 参数 (GLIMPSE-17775 校准值)
r_c, r_b = 0.003, 0.3
alpha_in, alpha_out = 2.0, 1.2
Sigma0 = 1e3
beta_val = 0.7
epsilon_g = 0.15
eta_visc = 0.5
r_grid = np.logspace(-4, 3, 2000)

# ===========================
# 5个源 (from Loiacono+ 2026 Table 4 + BiRD paper)
# ===========================
sources = {
    'BiRD':    {'z': 2.33, 'Lbol': 33e44, 'MBH_vir': 24e7,  'FWHM': 2000, 'b_opt': 2.71, 'cls': 'LRD'},
    'ID 1017': {'z': 2.50, 'Lbol': 3e44,  'MBH_vir': 1.7e7, 'FWHM': 1073, 'b_opt': 0.00, 'cls': '中间'},
    'ID 3646': {'z': 2.41, 'Lbol': 7e44,  'MBH_vir': 3.5e7, 'FWHM': 1237, 'b_opt':-1.79, 'cls': 'LBD'},
    'ID 8511': {'z': 4.50, 'Lbol': 2.8e44,'MBH_vir': 1.6e7, 'FWHM': 1668, 'b_opt': 2.46, 'cls': 'LRD'},
    'ID 9008': {'z': 4.50, 'Lbol': 4.7e44,'MBH_vir': 0.83e7,'FWHM': 1067, 'b_opt': 1.88, 'cls': 'LRD'},
}

def run_dyson(Mstar, r_eff, z, FWHM_target, eps_g=epsilon_g, eta=eta_visc):
    v_target = FWHM_target / 2.355

    rho_c_g = 1e6
    rho = np.zeros_like(r_grid)
    inner = r_grid < r_b
    rho[inner] = rho_c_g * (r_grid[inner] / r_c)**(-alpha_in)
    outer = r_grid >= r_b
    rho_at_b = rho_c_g * (r_b / r_c)**(-alpha_in)
    rho[outer] = rho_at_b * (r_grid[outer] / r_b)**(-alpha_out)
    rho[r_grid < r_c*0.1] = rho_c_g * (0.1)**(-alpha_in)

    dr = np.diff(np.concatenate([[0], r_grid]))
    M_enc_init = np.cumsum(rho * dr * 4 * np.pi * r_grid**2)
    idx_eff = np.argmin(np.abs(r_grid - r_eff))
    rho *= Mstar / M_enc_init[idx_eff]
    M_enc = np.cumsum(rho * dr * 4 * np.pi * r_grid**2)

    Sigma_r = M_enc / (np.pi * r_grid**2)
    G_ratio = 1.0 + (1+z)**3 * eps_g * (Sigma_r / Sigma0)**beta_val

    v_grav = np.sqrt(G_ratio * G_N * M_enc / r_grid)

    v_dyson = np.zeros_like(r_grid)
    v_dyson[0] = v_grav[0]
    for i in range(1, len(r_grid)):
        dr_i = r_grid[i] - r_grid[i-1]
        r_mid = (r_grid[i] + r_grid[i-1]) / 2
        coupling = eta * max(rho[i]/rho[0], 1e-8)**0.3
        visc = coupling * v_dyson[i-1]**2 / r_mid * dr_i
        grav_term = v_grav[i]**2 - v_grav[i-1]**2
        v_dyson[i] = np.sqrt(max(0, v_dyson[i-1]**2 + grav_term + visc))

    idx_blr = np.argmin(np.abs(r_grid - 0.03))
    Menc_blr = M_enc[idx_blr]
    G_blr = G_ratio[idx_blr]
    vg_blr = v_grav[idx_blr]
    vd_blr = v_dyson[idx_blr]

    if vg_blr >= v_target:
        M_central = 0.0
        mechanism = '纯星团核'
    else:
        M_central = max(0, v_target**2 * 0.03 / (G_blr * G_N) - Menc_blr)
        mechanism = '需IMBH' if M_central > 0 else '边缘'

    return {'M_central': M_central, 'v_grav_blr': vg_blr, 'v_dyson_blr': vd_blr,
            'G_ratio_blr': G_blr, 'mechanism': mechanism, 'M_enc_blr': Menc_blr,
            'Sigma_blr': Sigma_r[idx_blr], 'v_target': v_target}

# ===========================
# 扫描参数
# ===========================
Mstar_range = np.logspace(8.5, 11.0, 26)
reff_range = np.logspace(1.3, 2.5, 21)

print("=" * 68)
print("Loiacono 5源 Dyson 应力测试 (ε_g=0.15, GLIMPSE校准)")
print("=" * 68)

results = {}
for name, s in sources.items():
    z = s['z']
    FWHM = s['FWHM']
    v_target = FWHM / 2.355
    MBH_vir = s['MBH_vir']
    Lbol = s['Lbol']
    MBH_edd = Lbol / 1.26e38

    grid = np.zeros((len(Mstar_range), len(reff_range)))
    pure = np.zeros((len(Mstar_range), len(reff_range)), dtype=bool)
    G_blr_grid = np.zeros((len(Mstar_range), len(reff_range)))

    for i, Ms in enumerate(Mstar_range):
        for j, re in enumerate(reff_range):
            r = run_dyson(Ms, re, z, FWHM)
            grid[i, j] = r['M_central']
            pure[i, j] = (r['M_central'] == 0)
            G_blr_grid[i, j] = r['G_ratio_blr']

    n_pure = np.sum(pure)
    if n_pure > 0:
        # 纯星团核可行的参数
        pi, pj = np.where(pure)
        min_re_pure = reff_range[pj.min()]
        min_Ms_pure = Mstar_range[pi.min()]
    else:
        min_re_pure = min_Ms_pure = None

    # 找到 M_central < Naidu 上限 (10^5.1) 的区域
    naidu_mask = (grid > 0) & (grid < 10**5.1)
    n_naidu = np.sum(naidu_mask)
    # M_central < MBH_edd
    edd_mask = (grid > 0) & (grid < MBH_edd)
    n_edd = np.sum(edd_mask)

    # 最佳情况
    best_mc = np.min(grid)
    best_i, best_j = np.unravel_index(np.argmin(grid), grid.shape)

    results[name] = {
        'grid': grid, 'pure': pure, 'G_blr_grid': G_blr_grid,
        'n_pure': n_pure, 'n_naidu': n_naidu, 'n_edd': n_edd,
        'min_re_pure': min_re_pure, 'min_Ms_pure': min_Ms_pure,
        'best_mc': best_mc,
        'best_Ms': Mstar_range[best_i], 'best_re': reff_range[best_j],
        'v_target': v_target, 'MBH_vir': MBH_vir, 'MBH_edd': MBH_edd,
    }

    print(f"\n{'─'*68}")
    print(f"  {name}  z={z}  FWHM={FWHM}  MBH_vir=10^{np.log10(MBH_vir):.1f}  {s['cls']}")
    print(f"  Lbol={Lbol/1e44:.1f}e44  MBH_edd=10^{np.log10(MBH_edd):.1f}  β_opt={s['b_opt']}")
    print(f"{'─'*68}")
    print(f"  纯星团核: {n_pure}/{len(Mstar_range)*len(reff_range)} 组合可行")
    if min_Ms_pure:
        print(f"    最小 M*={min_Ms_pure/1e9:.1f}e9, r_eff={min_re_pure:.0f}pc")
    print(f"  M_c < Naidu(10^5.1): {n_naidu} 组合")
    print(f"  M_c < Eddington:     {n_edd} 组合")
    print(f"  最佳 M_c = 10^{np.log10(max(best_mc,1)):.1f} (M*={Mstar_range[best_i]/1e9:.1f}e9, r_eff={reff_range[best_j]:.0f}pc)")

# ===========================
# 图: 5源并列热力图
# ===========================
fig, axes = plt.subplots(2, 3, figsize=(20, 12))
axes = axes.flatten()

X, Y = np.meshgrid(np.log10(reff_range), np.log10(Mstar_range))

for idx, (name, s) in enumerate(sources.items()):
    ax = axes[idx]
    r = results[name]

    # 显示 log(M_central/M_BH_vir)
    ratio = r['grid'] / r['MBH_vir']
    ratio_log = np.where(r['grid'] > 0, np.log10(np.maximum(ratio, 1e-8)), -8)

    levels = [-8, -6, -4, -3, -2, -1, 0]
    cmap = plt.cm.RdYlBu_r
    c = ax.pcolormesh(X, Y, ratio_log, cmap=cmap, shading='auto',
                       vmin=-8, vmax=0)
    ax.contour(X, Y, ratio_log, levels=levels, colors='black', linewidths=0.4, alpha=0.5)

    # 纯星团核区
    if r['n_pure'] > 0:
        pi, pj = np.where(r['pure'])
        ax.scatter(X[pi, pj], Y[pi, pj], marker='.', color='lime', s=2, alpha=0.5)

    # 参考线
    ax.axhline(np.log10(r['MBH_edd']/0.1), color='green', ls=':', lw=1, alpha=0.5)
    ax.axhline(np.log10(r['MBH_edd']/0.01), color='green', ls='--', lw=1, alpha=0.5)

    ax.set_xlabel('log r_eff [pc]', fontsize=9)
    ax.set_ylabel('log M* [Msun]', fontsize=9)
    ax.set_title(f"{name} (z={s['z']}, {s['cls']})\n"
                 f"FWHM={s['FWHM']}, MBH_vir=10^{np.log10(r['MBH_vir']):.1f}",
                 fontsize=10)

plt.colorbar(c, ax=axes[4], label='log(M_central / M_BH_virial)', shrink=0.8)
axes[5].set_visible(False)

plt.suptitle('Dyson Model Stress Test: 5 Loiacono Sources (ε_g=0.15)',
             fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'loiacono_5sources_heatmap.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"\n  ✓ loiacono_5sources_heatmap.png")

# ===========================
# 图: M_central vs MBH_virial 对比
# ===========================
fig, ax = plt.subplots(figsize=(10, 8))

colors = {'BiRD': '#E53935', 'ID 1017': '#FF9800', 'ID 3646': '#4CAF50',
          'ID 8511': '#2196F3', 'ID 9008': '#9C27B0'}

for name, s in sources.items():
    r = results[name]
    MBH_vir = r['MBH_vir']
    z = s['z']
    FWHM = s['FWHM']

    # 扫几个代表性 M* 的 M_central 范围
    for mstar_idx in [5, 10, 15, 20]:
        mc_vals = r['grid'][mstar_idx, :]
        Ms = Mstar_range[mstar_idx]
        # 只标合理值
        valid = mc_vals > 0
        if np.any(valid):
            ax.plot(reff_range[valid], mc_vals[valid], '-', color=colors[name],
                   alpha=0.3, lw=0.5)

    # 标记最佳点
    ax.scatter(r['best_re'], r['best_mc'], color=colors[name], s=100, zorder=5,
              edgecolors='black', linewidths=1.5)
    ax.annotate(f"{name}\nFWHM={FWHM}", (r['best_re'], r['best_mc']),
               xytext=(10, 10), textcoords='offset points', fontsize=8,
               color=colors[name], fontweight='bold')

    # MBH_vir 线
    ax.axhline(MBH_vir, color=colors[name], ls=':', lw=1, alpha=0.4)

ax.axhline(10**7.4, color='gray', ls='--', lw=1, alpha=0.3, label='Eddington limit (BiRD)')
ax.axhline(10**5.1, color='blue', ls='--', lw=1, alpha=0.3, label='Naidu upper 10^5.1')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('r_eff [pc]', fontsize=13)
ax.set_ylabel('M_central [Msun]', fontsize=13)
ax.set_title('Dyson M_central for 5 Loiacono Sources\n(ε_g=0.15, swept over M* ∈ [3e8, 1e11] Msun)',
             fontsize=13)
ax.legend(fontsize=9, loc='lower left')
ax.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'loiacono_5sources_comparison.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  ✓ loiacono_5sources_comparison.png")

# ===========================
# 汇总表
# ===========================
print(f"\n{'='*68}")
print("汇总: 5源 Dyson 应力测试")
print(f"{'='*68}")
print(f"{'源':>10s} {'z':>5s} {'FWHM':>6s} {'cls':>6s} "
      f"{'MBH_vir':>10s} {'MBH_edd':>10s} {'纯核%':>6s} {'最佳Mc':>10s} {'vir/最佳':>8s}")
print(f"{'─'*68}")

for name, s in sources.items():
    r = results[name]
    ratio = r['MBH_vir'] / max(r['best_mc'], 1) if r['best_mc'] > 0 else float('inf')
    print(f"{name:>10s} {s['z']:>5.2f} {s['FWHM']:>6d} {s['cls']:>6s} "
          f"10^{np.log10(r['MBH_vir']):>4.1f}  10^{np.log10(r['MBH_edd']):>4.1f}   "
          f"{r['n_pure']/len(Mstar_range)/len(reff_range)*100:>5.0f}% "
          f"10^{np.log10(max(r['best_mc'],1)):>4.1f}  "
          f"{ratio:>6.0f}x{' (∞)' if ratio == float('inf') else ''}")

print(f"\n  结论:")
print(f"  ──────────────────────────────────────────────")
# 哪些是纯星团核
pure_sources = [name for name in sources if results[name]['n_pure'] > 0]
need_imbh = [name for name in sources if results[name]['n_pure'] == 0]

for name in pure_sources:
    r = results[name]
    pct = r['n_pure']/len(Mstar_range)/len(reff_range)*100
    print(f"  ✅ {name}: {pct:.0f}% 参数空间 → 纯星团核可行")
    if pct > 50:
        print(f"     → 大概率无中心SMBH/IMBH, X射线沉默天然解释")

for name in need_imbh:
    r = results[name]
    best = r['best_mc']
    ratio = r['MBH_vir'] / max(best, 1)
    print(f"  ⚠️ {name}: 需要 IMBH, 最佳 M_c=10^{np.log10(max(best,1)):.1f}, "
          f"维里修正 {ratio:.0f}x")
    if best < r['MBH_edd']:
        print(f"     M_c < Eddington 极限 → 亚爱丁顿, 但仍可能无强X射线")

print(f"\n  🔑 关键模式:")
print(f"  - 高z (ID 8511/9008 @ z=4.5): 需要更大的M*或更紧凑结构")
print(f"    由于(1+z)³=166, G_eff中等强度")
print(f"  - 低z (BiRD/1017/3646 @ z=2.3-2.5): 质量大→Σ大→G_eff强")
print(f"    (1+z)³=37, 但Σ补偿了背景减弱")
print(f"  - LBD (ID 3646) vs LRD: 3646 β_opt=-1.79(蓝) vs BiRD β_opt=+2.71(红)")
print(f"    → LBD的低Σ可能意味着它需要中心引擎, 而LRD的致密性=引擎无关性")
