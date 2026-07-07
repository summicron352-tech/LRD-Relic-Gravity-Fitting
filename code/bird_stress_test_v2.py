"""
BiRD 应力测试 v2: 扫描 ε_g 参数
===============================
目的: 找到使 BiRD 的 Dyson 模型产生合理 M_central 的 ε_g 范围

合理区间: M_central ∈ [10^4.1 (Naidu best), 10^7.4 (Eddington)]
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 12
import os

G_N = 4.30091e-3
output_dir = os.path.dirname(os.path.abspath(__file__))
fig_dir = os.path.join(os.path.dirname(output_dir), 'figures')
os.makedirs(fig_dir, exist_ok=True)

# BiRD 固定参数
z_bird = 2.33
FWHM_target = 2000
v_target = FWHM_target / 2.355

# Dyson 结构参数
r_c, r_b = 0.003, 0.3
alpha_in, alpha_out = 2.0, 1.2
Sigma0 = 1e3
beta_val = 0.7
eta_visc = 0.5
r_grid = np.logspace(-4, 3, 3000)

# 三个 M* 场景
scenarios = [
    (1e10, 150, 'M*=10^10, r_eff=150 pc'),
    (5e9, 100, 'M*=5×10^9, r_eff=100 pc'),
    (2e10, 200, 'M*=2×10^10, r_eff=200 pc'),
]

def run_dyson(Mstar, r_eff, z, epsilon_g):
    rho_c_g = 1e6
    rho = np.zeros_like(r_grid)
    inner = r_grid < r_b
    rho[inner] = rho_c_g * (r_grid[inner] / r_c)**(-alpha_in)
    outer = r_grid >= r_b
    rho_at_b = rho_c_g * (r_b / r_c)**(-alpha_in)
    rho[outer] = rho_at_b * (r_grid[outer] / r_b)**(-alpha_out)
    rho[r_grid < r_c*0.1] = rho_c_g * (0.1)**(-alpha_in)

    dr = np.diff(np.concatenate([[0], r_grid]))
    dV = 4 * np.pi * r_grid**2 * dr
    M_enc_init = np.cumsum(rho * dV)

    idx_eff = np.argmin(np.abs(r_grid - r_eff))
    scale = Mstar / M_enc_init[idx_eff]
    rho *= scale
    M_enc = np.cumsum(rho * dV)

    Sigma_r = M_enc / (np.pi * r_grid**2)
    G_ratio_z = 1.0 + (1+z)**3 * epsilon_g * (Sigma_r / Sigma0)**beta_val

    v_grav = np.sqrt(G_ratio_z * G_N * M_enc / r_grid)

    # Dyson cascade
    v_dyson = np.zeros_like(r_grid)
    v_dyson[0] = v_grav[0]
    for i in range(1, len(r_grid)):
        dr_i = r_grid[i] - r_grid[i-1]
        density_ratio = max(rho[i] / rho[0], 1e-8)
        coupling = eta_visc * density_ratio**0.3
        visc_term = coupling * v_dyson[i-1]**2 / ((r_grid[i]+r_grid[i-1])/2) * dr_i
        grav_term = v_grav[i]**2 - v_grav[i-1]**2
        v_dyson[i] = np.sqrt(max(0, v_dyson[i-1]**2 + grav_term + visc_term))

    idx_blr = np.argmin(np.abs(r_grid - 0.03))
    Menc_blr = M_enc[idx_blr]
    Gratio_blr = G_ratio_z[idx_blr]
    vg_blr = v_grav[idx_blr]
    vd_blr = v_dyson[idx_blr]

    # 如果 v_grav 已达到目标, 不需要中心质量
    if vg_blr >= v_target:
        M_central = 0.0
        mechanism = '纯星团核'
    else:
        # 需要多少中心质量使 v_dyson 达到 v_target
        M_central = max(0, v_target**2 * 0.03 / (Gratio_blr * G_N) - Menc_blr)
        mechanism = '需IMBH' if M_central > 0 else '边缘'

    return {
        'M_central': M_central,
        'v_grav_blr': vg_blr,
        'v_dyson_blr': vd_blr,
        'G_ratio_blr': Gratio_blr,
        'mechanism': mechanism,
        'M_enc_blr': Menc_blr,
    }

# 扫描 ε_g
eps_vals = np.logspace(-4, 0, 41)  # 0.0001 to 1.0

MBH_virial = 2.4e8
MBH_edd = 3.3e45 / (1.26e38)
MBH_naidu_upper = 10**5.1
MBH_naidu_best = 10**4.1

print(f"BiRD 应力测试 v2: ε_g 参数扫描")
print(f"目标: FWHM={FWHM_target}, z={z_bird}")
print(f"{'='*80}")

for Mstar, reff, label in scenarios:
    print(f"\n{'─'*80}")
    print(f"场景: {label}")
    print(f"{'─'*80}")
    print(f"{'ε_g':>10s} {'G_eff/G_N':>10s} {'v_grav':>8s} {'v_dyson':>9s} "
          f"{'M_central':>12s} {'log M_c':>8s} {'机制':>10s}")
    print(f"{'─'*65}")

    for eps in eps_vals:
        r = run_dyson(Mstar, reff, z_bird, eps)
        log_Mc = np.log10(r['M_central']) if r['M_central'] > 1 else 0.0
        print(f"{eps:>10.5f} {r['G_ratio_blr']:>10.1f} {r['v_grav_blr']:>8.0f} "
              f"{r['v_dyson_blr']:>9.0f} {r['M_central']:>12.1e} {log_Mc:>8.1f} {r['mechanism']:>10s}")

# 汇总: 各场景达到合理的 ε_g
print(f"\n{'='*80}")
print(f"汇总: 合理的 ε_g 范围")
print(f"  维里: 10^{np.log10(MBH_virial):.1f}, Eddington: 10^{np.log10(MBH_edd):.1f}")
print(f"  Naidu 上限: 10^5.1, 最佳: 10^4.1")
print(f"{'='*80}")

for Mstar, reff, label in scenarios:
    eps_for_naidu_range = []
    eps_for_edd_range = []
    eps_for_pure = None

    for eps in eps_vals:
        r = run_dyson(Mstar, reff, z_bird, eps)
        mc = r['M_central']
        if mc == 0:
            if eps_for_pure is None:
                eps_for_pure = eps
        elif mc <= MBH_naidu_upper:
            eps_for_naidu_range.append(eps)
        elif mc <= MBH_edd:
            eps_for_edd_range.append(eps)

    print(f"\n  {label}:")
    if eps_for_pure is not None:
        print(f"    ε_g ≲ {eps_for_pure:.4f} → 纯星团核 (v_grav > v_target)")
    if eps_for_naidu_range:
        print(f"    ε_g ∈ [{min(eps_for_naidu_range):.4f}, {max(eps_for_naidu_range):.4f}] → M_c < Naidu 上限 10^5.1")
    if eps_for_edd_range:
        print(f"    ε_g ∈ [{min(eps_for_edd_range):.4f}, {max(eps_for_edd_range):.4f}] → M_c < Eddington 10^{np.log10(MBH_edd):.1f}")

# 推荐参数
print(f"\n{'='*80}")
print(f"推荐: GLIMPSE-17775 (z=5.66) 的 ε_g=0.15 对应 BiRD (z=2.33)")

# 对典型场景,找到给出类似 G_eff 增强的 ε_g
Mstar_typ, reff_typ = 1e10, 150
# GLIMPSE @ z=5.66: G_eff/G_N @ 0.03pc 约多少?
# 从 lrd_dyson_gradient_model.py: M*=10^9.51, r_eff=83.3, z=5.66
# G_ratio_z @ 0.03pc 需要跑一下

# 快速估算 GLIMPSE 的 G_ratio
rho_g = np.zeros_like(r_grid)
inner = r_grid < r_b
rho_g[inner] = 1e6 * (r_grid[inner] / r_c)**(-alpha_in)
outer = r_grid >= r_b
rho_at_b_g = 1e6 * (r_b / r_c)**(-alpha_in)
rho_g[outer] = rho_at_b_g * (r_grid[outer] / r_b)**(-alpha_out)
rho_g[r_grid < r_c*0.1] = 1e6 * (0.1)**(-alpha_in)
dr = np.diff(np.concatenate([[0], r_grid]))
M_enc_g = np.cumsum(rho_g * dr * 4 * np.pi * r_grid**2)
idx_eff_g = np.argmin(np.abs(r_grid - 83.3))
scale_g = 10**9.51 / M_enc_g[idx_eff_g]
M_enc_g *= scale_g
Sigma_g = M_enc_g / (np.pi * r_grid**2)
idx_blr_g = np.argmin(np.abs(r_grid - 0.03))
Gratio_glimpse = 1.0 + (1+5.66)**3 * 0.15 * (Sigma_g[idx_blr_g] / 1000)**0.7

# 找匹配的 ε_g for BiRD
for eps in eps_vals:
    r = run_dyson(Mstar_typ, reff_typ, z_bird, eps)
    if r['G_ratio_blr'] <= Gratio_glimpse:
        eps_matched = eps
        break

print(f"  GLIMPSE G_eff/G_N @ 0.03pc ≈ {Gratio_glimpse:.0f}")
print(f"  BiRD 要达到相同 G_eff, 需 ε_g ≈ {eps_matched:.4f}")
print(f"  → BiRD 的合理 ε_g 约 = 0.15 × [{Gratio_glimpse:.0f}匹配值] ≈ 需要校准")

# 找出"最佳"ε_g: M_central 在 Naidu 和 Eddington 之间
print(f"\n  推荐 ε_g 范围 (M_central 在 Naidu 和 Eddington 之间):")
for Mstar, reff, label in scenarios:
    eps_candidates = []
    for eps in eps_vals:
        r = run_dyson(Mstar, reff, z_bird, eps)
        if MBH_naidu_best < r['M_central'] < MBH_edd and r['M_central'] > 0:
            eps_candidates.append(eps)
    if eps_candidates:
        print(f"    {label}: ε_g ∈ [{min(eps_candidates):.4f}, {max(eps_candidates):.4f}]")
    else:
        print(f"    {label}: 无解 — 纯星团核覆盖太广")

# 图: ε_g vs M_central
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

colors = ['#2196F3', '#4CAF50', '#FF9800']
for (Mstar, reff, label), color in zip(scenarios, colors):
    mc_vals = [run_dyson(Mstar, reff, z_bird, eps)['M_central'] for eps in eps_vals]
    mc_log = [np.log10(max(m, 1)) for m in mc_vals]
    ax1.plot(eps_vals, mc_log, 'o-', color=color, markersize=3, label=label, alpha=0.8)

ax1.axhline(np.log10(MBH_virial), color='red', ls='--', lw=1.5, label=f'Virial (8.4)')
ax1.axhline(np.log10(MBH_edd), color='orange', ls='--', lw=1.5, label=f'Eddington (7.4)')
ax1.axhline(np.log10(MBH_naidu_upper), color='blue', ls=':', lw=1.5, label=f'Naidu upper (5.1)')
ax1.axhline(0, color='gray', ls='-', lw=1, alpha=0.5, label='Pure cluster (0)')
ax1.set_xscale('log')
ax1.set_xlabel('epsilon_g', fontsize=13)
ax1.set_ylabel('log M_central [Msun]', fontsize=13)
ax1.set_title('BiRD: M_central vs epsilon_g', fontsize=14)
ax1.legend(fontsize=9, loc='upper left')
ax1.grid(True, alpha=0.3)

# G_eff vs eps_g
for (Mstar, reff, label), color in zip(scenarios, colors):
    gr_vals = [run_dyson(Mstar, reff, z_bird, eps)['G_ratio_blr'] for eps in eps_vals]
    ax2.plot(eps_vals, gr_vals, 'o-', color=color, markersize=3, label=label, alpha=0.8)

ax2.axhline(1, color='gray', ls='-', lw=1, alpha=0.5, label='G_eff = G_N')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('epsilon_g', fontsize=13)
ax2.set_ylabel('G_eff/G_N @ 0.03 pc', fontsize=13)
ax2.set_title('BiRD: G_eff enhancement vs epsilon_g', fontsize=14)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'bird_eps_scan.png'), dpi=150)
plt.close()
print(f"\n图表已保存: bird_eps_scan.png")
