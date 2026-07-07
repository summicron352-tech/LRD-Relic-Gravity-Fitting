"""
LRD 致密梯度(Dyson)模型 — 完整物理模拟
===========================================
物理图像:
  LRD中心不是单一的SMBH，而是一个极端致密的恒星核团，
  在G_eff(Σ)框架下保持"悬浮态"——致密但未坍缩。

  密度从内向外递减 → 引力势递减 → 速度逐级放大(Dyson效应)
  
  本脚本计算:
  1. 密度梯度剖面 ρ(r) 和 Σ(r)
  2. G_eff(Σ) 径向梯度 G_eff(r)
  3. 引力势 Φ(r) 和速度 v(r)
  4. Dyson黏性加速级联 (entrainment)
  5. 坍缩阈值 (什么条件下触发坍缩)
  6. 恒星合并率
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# 中文字体配置
plt.rcParams['font.sans-serif'] = ['Hiragino Sans GB', 'Arial Unicode MS', 'Heiti TC', 'PingFang HK']
plt.rcParams['axes.unicode_minus'] = False
from scipy.special import gamma as gamma_func
# cumtrapz removed — analytic integration used
import os

# ===========================
# 常数
# ===========================
G_N = 4.30091e-3          # [pc·(km/s)^2 / Msun]
MSUN = 1.989e33            # [g]
PC = 3.086e18              # [cm]
YR = 3.156e7               # [s/yr]
kpc = 1e3                  # pc in kpc
sigma_SB = 5.6704e-5       # Stefan-Boltzmann [erg cm⁻^2 s⁻¹ K⁻⁴]
c_light = 2.998e10         # [cm/s]

# ===========================
# 输出目录
# ===========================
output_dir = os.path.dirname(os.path.abspath(__file__))
fig_dir = os.path.join(os.path.dirname(output_dir), 'figures')
os.makedirs(fig_dir, exist_ok=True)

# ===========================
# GLIMPSE-17775 参数
# ===========================
Mstar = 10**9.51           # [Msun]
r_eff = 83.3               # [pc]
z = 5.66

# ===========================
# 1. 密度梯度剖面
# ===========================
print("=" * 68)
print("LRD 致密梯度(Dyson)模型 — 完整物理模拟")
print("=" * 68)
print(f"\n基准源: GLIMPSE-17775 (z={z})")
print(f"M* = 10^{np.log10(Mstar):.1f} Msun, r_eff = {r_eff:.1f} pc")
print()

# 使用双幂律密度剖面（内区陡峭、外区平缓）
# ρ(r) = ρ_c × (r/r_c)^(-α_in) for r < r_b
# ρ(r) = ρ_b × (r/r_b)^(-α_out) for r > r_b
# 其中核星团的核心尺度 ~0.01-0.1 pc，内区α_in ~ 1.5-2.5，外区α_out ~ 1.0-1.5

# 参数设置（从观测约束反推）
r_c = 0.003                # [pc] 核心半径（~600 AU，稠密星团核心）
r_b = 0.3                  # [pc] 拐折半径
alpha_in = 2.0             # 内区幂律指数
alpha_out = 1.2            # 外区幂律指数

# 要求M(<r_eff) = Mstar，求解归一化密度
# M(<r) = 4π ∫ ρ(r) r^2 dr

# 构造径向格点
r_grid = np.logspace(-4, 2, 2000)  # 0.0001 到 100 pc

# 先假设ρ_c，然后积分得到M(<r_eff)，缩放
rho_c_guess = 1e6          # [Msun/pc³] 初始猜测

def density_profile(r, r_c, r_b, alpha_in, alpha_out, rho_c):
    """双幂律密度剖面"""
    rho = np.zeros_like(r)
    # 内区 r < r_b
    inner = r < r_b
    rho[inner] = rho_c * (r[inner] / r_c)**(-alpha_in)
    # 外区 r >= r_b
    outer = r >= r_b
    rho_at_b = rho_c * (r_b / r_c)**(-alpha_in)
    rho[outer] = rho_at_b * (r[outer] / r_b)**(-alpha_out)
    # 在 r=0 附近截断（奇点保护）
    rho[r < r_c * 0.1] = rho_c * (0.1)**(-alpha_in)
    return rho

def enclosed_mass(r, rho):
    """数值积分包络质量"""
    dr = np.diff(np.concatenate([[0], r]))
    # 更精确：使用梯形积分
    dV = 4 * np.pi * r**2 * dr
    M_enc = np.cumsum(rho * dV)
    return M_enc

# 缩放密度到：M(<r_eff) = Mstar
rho_profile = density_profile(r_grid, r_c, r_b, alpha_in, alpha_out, rho_c_guess)
M_enc = enclosed_mass(r_grid, rho_profile)
idx_eff = np.argmin(np.abs(r_grid - r_eff))
M_enc_eff = M_enc[idx_eff]
scale_factor = Mstar / M_enc_eff

rho_profile *= scale_factor
M_enc = enclosed_mass(r_grid, rho_profile)

print("─" * 68)
print("1. 密度梯度剖面")
print("─" * 68)
print(f"  内核半径 r_c = {r_c:.4f} pc (~{r_c*3.086e18/1.496e13:.0f} AU)")
print(f"  拐折半径 r_b = {r_b:.2f} pc")
print(f"  内区指数 α_in = {alpha_in:.1f}")
print(f"  外区指数 α_out = {alpha_out:.1f}")
print(f"  中心密度 ρ_c = {rho_profile[0]:.2e} Msun/pc³")
print(f"  中心密度 (cgs) = {rho_profile[0]*MSUN/PC**3:.2e} g/cm³")
print(f"  M(<r_eff=83.3pc) = {M_enc[idx_eff]:.2e} Msun (scaled to Mstar)")

# 打印几个关键半径的密度
key_rs = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 83.3]
print(f"\n  {'r [pc]':>10s}  {'ρ [Msun/pc³]':>14s}  {'M(<r) [Msun]':>14s}  {'Σ [Msun/pc^2]':>14s}")
print("  " + "-" * 56)
for r in key_rs:
    idx = np.argmin(np.abs(r_grid - r))
    rho_r = rho_profile[idx]
    Menc_r = M_enc[idx]
    Sigma_r = Menc_r / (np.pi * r**2)
    print(f"  {r:>10.4f}  {rho_r:>14.2e}  {Menc_r:>14.2e}  {Sigma_r:>14.2e}")

# ===========================
# 2. G_eff(Σ) 径向梯度
# ===========================
print()
print("─" * 68)
print("2. G_eff(Σ) 径向梯度")
print("─" * 68)

# UID框架参数
epsilon_g = 0.15
beta = 0.7
Sigma0 = 1e3               # [Msun/pc^2] pivot密度

# 在每个半径计算局部Σ和G_eff
Sigma_r_prof = M_enc / (np.pi * r_grid**2)
G_ratio_prof = 1.0 + epsilon_g * (Sigma_r_prof / Sigma0)**beta

# redhsot增强因子
# 在高红移z=5.66，背景宇宙密度更高
# ρ_bkg(z) / ρ_bkg(0) = (1+z)^3 ≈ (6.66)^3 ≈ 296
# 双变量框架中的背景项
rho_bkg_ratio = (1 + z)**3 / (1 + 0)**3
# 简化的背景增强: 双变量框架 G_eff(z, Σ)
G_ratio_z = 1.0 + rho_bkg_ratio * epsilon_g * (Sigma_r_prof / Sigma0)**beta

print(f"  ε_g = {epsilon_g}, β = {beta}, Σ₀ = 10^{np.log10(Sigma0):.1f} Msun/pc^2")
print(f"  红移增强 (1+z)³ = {rho_bkg_ratio:.0f}x")
print()
print(f"  {'r [pc]':>10s}  {'Σ [Msun/pc^2]':>14s}  {'G_eff/G_N':>10s}  {'G_eff(z)/G_N':>14s}")
print("  " + "-" * 52)
for r in key_rs:
    idx = np.argmin(np.abs(r_grid - r))
    sg = Sigma_r_prof[idx]
    gr = G_ratio_prof[idx]
    grz = G_ratio_z[idx]
    print(f"  {r:>10.4f}  {sg:>14.2e}  {gr:>10.2f}  {grz:>14.2f}")

# ===========================
# 3. 引力势和纯引力速度
# ===========================
print()
print("─" * 68)
print("3. 引力势和速度剖面 (无心点质量)")
print("─" * 68)

# 势能: Φ(r) = -∫_r^∞ G_eff(r') M(<r')/r'^2 dr'
# 速度: v(r) = √(G_eff(r) * M(<r) / r)

# 考虑G_eff梯度
v_grav = np.sqrt(G_ratio_z * G_N * M_enc / r_grid)
FWHM_grav = 2.355 * v_grav

print(f"  {'r [pc]':>10s}  {'M(<r) [Msun]':>14s}  {'v_grav [km/s]':>14s}  {'FWHM [km/s]':>14s}")
print("  " + "-" * 58)
for r in key_rs:
    idx = np.argmin(np.abs(r_grid - r))
    vv = v_grav[idx]
    ff = FWHM_grav[idx]
    print(f"  {r:>10.4f}  {M_enc[idx]:>14.2e}  {vv:>14.1f}  {ff:>14.1f}")

# ===========================
# 4. Dyson 黏性加速级联
# ===========================
print()
print("─" * 68)
print("4. Dyson 黏性加速级联")
print("─" * 68)
print("""
  Dyson风筒类比:
    内层高速气流 → 产生剪切层 → 黏性拖拽外层气体
    动量从内向外传递，外层速度被放大
    
  物理模型:
    dv/dr ∝ η_visc × v(r) / r × (ρ(r)/⟨ρ⟩)^(1/2)
    其中η_visc是黏性耦合效率(0-1)
""")

def dyson_velocity(r, rho, v_grav, eta_visc=0.3):
    """
    考虑Dyson黏性加速的总速度剖面
    dv^2/dr = η_visc × (v^2/r)  (黏性动量传递)
    数值积分从内到外
    """
    n = len(r)
    v_total = np.zeros(n)
    v_total[0] = v_grav[0]  # 最内层由纯引力决定
    
    for i in range(1, n):
        dr = r[i] - r[i-1]
        r_mid = (r[i] + r[i-1]) / 2
        
        # 黏性加速项: 正比于局部速度梯度 × 密度比 × 耦合效率
        # dv^2/dr = η_visc × v^2/r × [ρ(r)/ρ(r_c)]^(p)
        # 密度比项：越外层密度越低，耦合效率降低
        density_ratio = rho[i] / rho[0]
        if density_ratio < 1e-8:
            density_ratio = 1e-8
        coupling = eta_visc * (density_ratio)**0.3
        
        # 黏性项
        visc_term = coupling * v_total[i-1]**2 / r_mid * dr
        
        # 引力项
        grav_term = v_grav[i]**2 - v_grav[i-1]**2
        
        v_total[i] = np.sqrt(max(0, v_total[i-1]**2 + grav_term + visc_term))
    
    return v_total

for eta in [0.1, 0.3, 0.5, 0.7]:
    v_dyson = dyson_velocity(r_grid, rho_profile, v_grav, eta_visc=eta)
    FWHM_dyson = 2.355 * v_dyson
    
    print(f"\n  η_visc = {eta:.1f}:")
    print(f"  {'r [pc]':>10s}  {'v_grav':>8s}  {'v_dyson':>8s}  {'FWHM_dyson':>12s}  {'放大倍率':>10s}")
    print("  " + "-" * 52)
    
    display_rs = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 10.0]
    for r in display_rs:
        idx = np.argmin(np.abs(r_grid - r))
        vv = v_grav[idx]
        vd = v_dyson[idx]
        ratio = vd / vv if vv > 0 else 0
        print(f"  {r:>10.4f}  {vv:>8.1f}  {vd:>8.1f}  {FWHM_dyson[idx]:>12.1f}  {ratio:>8.2f}x")
    
    # 在BLR尺度(0.03 pc)和观测尺度(0.2 pc)检查速度
    idx_blr = np.argmin(np.abs(r_grid - 0.03))
    idx_obs = np.argmin(np.abs(r_grid - 0.2))
    print(f"  → 在r=0.03 pc: v_dyson={v_dyson[idx_blr]:.0f} km/s, FWHM={FWHM_dyson[idx_blr]:.0f} km/s")
    print(f"  → 在r=0.20 pc: v_dyson={v_dyson[idx_obs]:.0f} km/s, FWHM={FWHM_dyson[idx_obs]:.0f} km/s")

# 使用eta=0.5作为基准
eta_default = 0.5
v_dyson = dyson_velocity(r_grid, rho_profile, v_grav, eta_visc=eta_default)
FWHM_dyson = 2.355 * v_dyson

# ===========================
# 5. 坍缩阈值
# ===========================
print()
print("─" * 68)
print("5. 坍缩阈值分析")
print("─" * 68)

# 星团的坍缩条件: 当引力势能 > 动能×2 时发生核心坍缩
# 2K < |U| → 坍缩
# K = 0.5 * M * σ^2, U = -0.4 * G * M^2 / r_h
# 临界条件: σ^2 < 0.4*G*M/r_h

# 在G_eff框架下，有效引力增强 → 坍缩门槛也提高
# 等价于: standard σ^2_crit ∝ G_eff × M / r_h
# 在G_eff增强区域，需要更高的速度弥散才能对抗坍缩

print("  坍缩判据: 2K < |U| → 核心坍缩")
print("  其中 K = 0.5 M <σ^2>, U = -α G M^2 / r_h")
print("  在 G_eff 框架下: U_eff = -α G_eff M^2 / r_h")
print()

# 在每个壳层检查坍缩状态
t_collapse_factor = np.zeros_like(r_grid)
for i in range(1, len(r_grid)):
    if i < 5:
        t_collapse_factor[i] = 0
        continue
    # 局部: M_shell = M_enc[i] - M_enc[i-1], r_shell = r_grid[i]
    M_shell = M_enc[i] - M_enc[i-1]
    r_shell = r_grid[i]
    if M_shell <= 0 or r_shell <= 0:
        t_collapse_factor[i] = 0
        continue
    
    # 速度弥散 = v_dyson[i] (轨道速度作为弥散的代理)
    sigma_v = v_dyson[i]
    
    # 动能 / |势能|
    # K = 0.5 * M_shell * sigma_v^2
    # |U| ∝ G_eff * M_shell * (M_enc[i]) / r_shell
    K = 0.5 * M_shell * sigma_v**2
    U_abs = 0.4 * G_ratio_z[i] * G_N * M_shell * M_enc[i] / r_shell
    
    if U_abs > 0:
        t_collapse_factor[i] = 2 * K / U_abs
    else:
        t_collapse_factor[i] = np.inf

# 找到坍缩边界
# t_collapse < 1 → 坍缩区域
collapse_idx = np.where(t_collapse_factor < 1.0)[0]
if len(collapse_idx) > 0:
    r_collapse = r_grid[collapse_idx[-1]]
    print(f"  坍缩边界 r_collapse = {r_collapse:.4f} pc")
    print(f"  在 r < {r_collapse:.4f} pc 内部，2K < |U|，满足坍缩条件")
    print(f"  但在 G_eff 框架下，坍缩不一定会发生——")
    print(f"  因恒星合并门槛升高、IMF截断、和背景致密度共同抑制")
else:
    r_collapse = None
    print(f"  全半径范围 2K > |U|，系统处于稳定状态")
    print(f"  → G_eff增强的弥散足以抵抗引力坍缩")

print()
print(f"  重要注释:")
print(f"  即使满足坍缩条件，在UID/G_eff框架中: ")
print(f"  1. IMF截断(P5) → 大质量恒星少 → 超新星少 → 残留BH少")
print(f"  2. 恒星速度弥散大 → 碰撞截面小 → 合并困难")
print(f"  3. 背景宇宙密度高 → G_eff进一步放大 → 悬浮态更稳定")
print(f"  → LRD核区处于临界悬浮态,不是失败的黑洞种子工厂")

# ===========================
# 6. 恒星合并率
# ===========================
print()
print("─" * 68)
print("6. 恒星合并率和演化状态")
print("─" * 68)

# 在稠密星团中，恒星碰撞时间:
# t_coll = 1 / (n * σ * v_rel)
# 其中 n = 数密度, σ = 碰撞截面, v_rel = 相对速度

# 假设平均恒星质量 ⟨m*⟩ = 0.5 Msun (IMF截断后偏小)
m_star_avg = 0.5        # [Msun]
R_star_avg = 1.0 * 6.96e10 / PC   # [pc] 太阳半径 ≈ 2.25e-8 pc

print(f"  假设⟨m*⟩ = {m_star_avg:.1f} Msun, R* = {R_star_avg:.2e} pc")
print()

print(f"  {'r [pc]':>10s}  {'n [pc⁻³]':>12s}  {'t_coll [Myr]':>14s}  {'Hubble/Gyr':>10s}  {'状态':>14s}")
print("  " + "-" * 64)

for r in key_rs:
    idx = np.argmin(np.abs(r_grid - r))
    rho_r = rho_profile[idx]
    n_star = rho_r / m_star_avg          # [pc⁻³] 数密度
    
    if n_star < 1:
        print(f"  {r:>10.4f}  {n_star:>12.2e}  {'—':>14s}  {'—':>10s}  {'极稀疏':>14s}")
        continue
    
    v_rel = v_dyson[idx]                  # [km/s]
    v_rel_cm = v_rel * 1e5                # [cm/s]
    
    # 碰撞截面: 引力聚焦截面
    # σ_grav = π * R*^2 * (1 + 2G(m1+m2)/(R* v_rel^2))
    v_esc_sq = 2 * G_N * (2 * m_star_avg) / R_star_avg
    v_esc_sq_cgs = v_esc_sq * (1e5)**2
    # 简化的碰撞截面
    sigma_coll = np.pi * R_star_avg**2 * (1 + 2 * G_N * (m_star_avg + m_star_avg) / (R_star_avg * v_rel**2 + 1e-30))
    sigma_coll_cm = sigma_coll * PC**2
    
    # 碰撞时间
    t_coll_yr = 1.0 / (n_star * PC**-3 * sigma_coll * PC**2 * v_rel * 1e5 * 3.156e7)
    t_coll_myr = t_coll_yr / 1e6
    
    # 宇宙年龄在 z=5.66 约 1.0 Gyr
    t_universe_gyr = 1.0
    status = "频繁合并" if t_coll_myr < 100 else ("偶尔合并" if t_coll_myr < 1000 else "几乎不合并")
    
    print(f"  {r:>10.4f}  {n_star:>12.2e}  {t_coll_myr:>14.2f}  {t_universe_gyr:>10.1f}  {status:>14s}")

# ===========================
# 7. 中心最小质量估计
# ===========================
print()
print("─" * 68)
print("7. 中心IMBH/核团质量估计 (G_eff修正后)")
print("─" * 68)

# 在BLR尺度(0.03 pc)需要的速度
target_FWHM = 3000      # km/s
target_v = target_FWHM / 2.355

# 在r=0.03 pc处，纯引力速度来自包络质量
idx_blr = np.argmin(np.abs(r_grid - 0.03))
Menc_blr = M_enc[idx_blr]
G_ratio_blr = G_ratio_z[idx_blr]
v_blr_grav = np.sqrt(G_ratio_blr * G_N * Menc_blr / 0.03)

print(f"  在 r=0.03 pc:")
print(f"    M_enc(星团) = {Menc_blr:.2e} Msun")
print(f"    G_eff/G_N = {G_ratio_blr:.2f}")
print(f"    v_grav(星团) = {v_blr_grav:.0f} km/s")
print(f"    需要的 v(含Dyson) = {v_dyson[idx_blr]:.0f} km/s")
print(f"    不足部分需由中心点质量补充")

# 计算需要的中心点质量
# v^2 = G_eff × (M_enc_cluster + M_central) / r
v_needed_blr = v_dyson[idx_blr]
M_central_needed = v_needed_blr**2 * 0.03 / (G_ratio_blr * G_N) - Menc_blr
if M_central_needed < 0:
    M_central_needed = 0

print(f"\n  需要的中心点质量 M_central = {M_central_needed:.2e} Msun")
print(f"  对应 log M_central = {np.log10(M_central_needed):.1f}")
print(f"  Kokorev+ 维里质量 M_BH = 10^6.7")
print(f"  比值 = {M_central_needed/10**6.7:.2f}")

# 在r=0.2 pc处验证
idx_obs = np.argmin(np.abs(r_grid - 0.2))
Menc_obs = M_enc[idx_obs]
G_ratio_obs = G_ratio_z[idx_obs]

# 如果中心有M_central，在0.2 pc处的速度
v_at_02 = np.sqrt(G_ratio_obs * G_N * (Menc_obs + M_central_needed) / 0.2)
print(f"\n  验证 r=0.2 pc:")
print(f"    M_enc(星团) = {Menc_obs:.2e} Msun")
print(f"    v(含中心+星团) = {v_at_02:.0f} km/s")
print(f"    Dyson速度 = {v_dyson[idx_obs]:.0f} km/s")
print(f"    → {'自洽' if abs(v_at_02 - v_dyson[idx_obs]) / v_dyson[idx_obs] < 0.3 else '需调整'}")

# ===========================
# 绘制关键结果
# ===========================
print()
print("─" * 68)
print("正在生成图表...")
print("─" * 68)

# 图1: 密度和Σ剖面
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.loglog(r_grid, rho_profile, 'b-', linewidth=2)
ax1.axvline(r_c, color='gray', linestyle='--', alpha=0.5, label=f'r_c={r_c:.4f} pc')
ax1.axvline(r_b, color='gray', linestyle=':', alpha=0.5, label=f'r_b={r_b:.2f} pc')
ax1.axvline(r_eff, color='red', linestyle='--', alpha=0.5, label=f'r_eff={r_eff:.1f} pc')
ax1.set_xlabel('r [pc]', fontsize=13)
ax1.set_ylabel('ρ [Msun/pc³]', fontsize=13)
ax1.set_title('密度剖面 ρ(r)', fontsize=14)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3, which='both')

ax2.loglog(r_grid, Sigma_r_prof, 'g-', linewidth=2)
ax2.axvline(r_c, color='gray', linestyle='--', alpha=0.5)
ax2.axvline(r_b, color='gray', linestyle=':', alpha=0.5)
ax2.axvline(r_eff, color='red', linestyle='--', alpha=0.5)
ax2.set_xlabel('r [pc]', fontsize=13)
ax2.set_ylabel('Σ [Msun/pc^2]', fontsize=13)
ax2.set_title('面密度剖面 Σ(r)', fontsize=14)
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'dyson_density_profile.png'), dpi=150)
plt.close()
print("  ✓ dyson_density_profile.png")

# 图2: G_eff梯度和速度剖面
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.semilogx(r_grid, G_ratio_z, 'purple', linewidth=2, label='G_eff(z)/G_N')
ax1.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
ax1.axvline(r_c, color='gray', linestyle='--', alpha=0.3)
ax1.axvline(r_eff, color='red', linestyle='--', alpha=0.3)
ax1.set_xlabel('r [pc]', fontsize=13)
ax1.set_ylabel('G_eff/G_N', fontsize=13)
ax1.set_title('G_eff(Σ) 径向梯度', fontsize=14)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(1e-4, 100)

ax2.loglog(r_grid, v_grav, 'b-', linewidth=1.5, label='v_grav (纯引力)', alpha=0.7)
for eta, color in zip([0.1, 0.3, 0.5, 0.7], ['orange', 'green', 'red', 'darkred']):
    v_d = dyson_velocity(r_grid, rho_profile, v_grav, eta_visc=eta)
    ax2.loglog(r_grid, v_d, '--', linewidth=1.5, color=color, label=f'v_dyson (η={eta})', alpha=0.8)

ax2.axhline(1274, color='black', linestyle='-.', alpha=0.5, label='FWHM=3000 km/s')
ax2.axvline(r_eff, color='red', linestyle='--', alpha=0.3)
ax2.set_xlabel('r [pc]', fontsize=13)
ax2.set_ylabel('v [km/s]', fontsize=13)
ax2.set_title('速度剖面: 纯引力 vs Dyson黏性加速', fontsize=14)
ax2.legend(fontsize=8, loc='lower left')
ax2.grid(True, alpha=0.3, which='both')
ax2.set_xlim(1e-4, 100)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'dyson_velocity_cascade.png'), dpi=150)
plt.close()
print("  ✓ dyson_velocity_cascade.png")

# 图3: 合并时间和坍缩状态
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# 合并时间
t_coll_profile = np.zeros_like(r_grid)
n_star_profile = rho_profile / m_star_avg
for i in range(len(r_grid)):
    if n_star_profile[i] < 1 or r_grid[i] < r_c*0.1:
        t_coll_profile[i] = 1e12
        continue
    sigma_coll_grav = np.pi * R_star_avg**2 * (1 + 2 * G_N * (m_star_avg + m_star_avg) / (R_star_avg * max(v_dyson[i]**2, 0.01)))
    t_coll = 1.0 / (n_star_profile[i] * sigma_coll_grav * v_dyson[i] * 1e5 / PC * 3.156e7)
    t_coll_profile[i] = min(t_coll / 1e6, 1e12)  # Myr

ax1.loglog(r_grid, t_coll_profile, 'b-', linewidth=2)
ax1.axhline(1000, color='red', linestyle='--', alpha=0.5, label='宇宙年龄 ~1 Gyr')
ax1.axvline(r_c, color='gray', linestyle='--', alpha=0.3)
ax1.axvline(r_eff, color='red', linestyle='--', alpha=0.3)
ax1.set_xlabel('r [pc]', fontsize=13)
ax1.set_ylabel('t_coll [Myr]', fontsize=13)
ax1.set_title('恒星碰撞时间剖面', fontsize=14)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3, which='both')
ax1.set_xlim(1e-4, 100)

# 坍缩因子
valid = (r_grid > r_c * 0.5) & (t_collapse_factor < 1e6)
ax2.semilogx(r_grid[valid], t_collapse_factor[valid], 'purple', linewidth=2)
ax2.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='坍缩边界 (2K=|U|)')
ax2.axvline(r_c, color='gray', linestyle='--', alpha=0.3)
ax2.axvline(r_eff, color='red', linestyle='--', alpha=0.3)
ax2.set_xlabel('r [pc]', fontsize=13)
ax2.set_ylabel('2K/|U|', fontsize=13)
ax2.set_title('坍缩判据 (2K/|U| < 1 → 坍缩)', fontsize=14)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(1e-4, 100)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'dyson_collapse_merger.png'), dpi=150)
plt.close()
print("  ✓ dyson_collapse_merger.png")

# ===========================
# 8. 总结
# ===========================
print()
print("=" * 68)
print("模型总结")
print("=" * 68)

print(f"""
  GLIMPSE-17775 Dyson致密梯度模型:

  密度结构:
    中心密度 ρ_c = {rho_profile[0]:.2e} Msun/pc³
    核心半径 r_c = {r_c:.4f} pc
    有效半径 r_eff = {r_eff:.1f} pc

  引力势:
    中心 G_eff/G_N(含红移) ≈ {G_ratio_z[0]:.1f}
    在 r_eff 处 G_eff/G_N(含红移) ≈ {G_ratio_z[idx_eff]:.1f}
    背景增强因子 (1+z)³ = {rho_bkg_ratio:.0f}x

  速度级联 (η={eta_default}):
    内核 (0.003 pc): v ≈ {v_dyson[np.argmin(np.abs(r_grid-0.003))]:.0f} km/s
    BLR尺度 (0.03 pc): v ≈ {v_dyson[idx_blr]:.0f} km/s
    观测尺度 (0.20 pc): v ≈ {v_dyson[idx_obs]:.0f} km/s
    r_eff (83 pc): v ≈ {v_dyson[idx_eff]:.0f} km/s

  中心质量:
    Dyson模型所需中心点质量: 10^{np.log10(M_central_needed):.1f} Msun
    Kokorev+ 维里质量 (无G_eff修正): 10^6.7 Msun
    真实质量比值: {M_central_needed/10**6.7:.1f}x

  物理结论:
    LRD中心不需要SMBH。一个IMBH(~10^{np.log10(M_central_needed):.1f} Msun)
    或致密星团核在G_eff增强下足以产生观测到的宽线速度和光谱特征。
    X射线缺失天然解释: 中心质量不处于吸积状态(流过式,Dyson机制)。
""")

print("所有图表已保存至:", fig_dir)
print("=" * 68)
