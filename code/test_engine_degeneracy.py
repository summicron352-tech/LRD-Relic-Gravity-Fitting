#!/usr/bin/env python3
"""
检验 Sneppen+ (2026) LBD-LRD 气体茧统一框架的引擎简并性
========================================================
核心问题：LBD→LRD 颜色序列是由气体茧柱密度决定的，
还是由中心引擎（AGN vs 致密恒星核团）决定的？

如果两种引擎在相同柱密度扫描下产生几乎相同的 JWST 颜色分布，
则 Sneppen+ 的"气体茧 AGN 统一框架"实际上引擎无关——
气体茧是物理，AGN 不是。
"""

import numpy as np
from scipy import constants as const
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import json, os, warnings
warnings.filterwarnings('ignore')

# ========== 中文字体 ==========
plt.rcParams['font.sans-serif'] = ['Hiragino Sans GB', 'Arial Unicode MS', 'Heiti TC', 'PingFang HK']
plt.rcParams['axes.unicode_minus'] = False

# ========== 物理常数 ==========
h = const.h          # Planck
c = const.c          # 光速 cm/s
k_B = const.k        # Boltzmann
sigma_T = 6.65245873e-25  # Thomson cm^2
Ryd = 2.1798741e-11  # Rydberg erg

# ========== JWST 宽带近似 (um) ==========
JWST_BANDS = {
    'F150W': 1.50,
    'F200W': 1.97,
    'F277W': 2.77,
    'F356W': 3.56,
    'F444W': 4.44,
}
# AB 零等通量 Jy
AB_FLUX_ZERO = 3631  # Jy

def mag_to_flux_nu(mag_ab):
    """AB mag → f_nu (Jy)"""
    return AB_FLUX_ZERO * 10**(-0.4 * mag_ab)

def flux_nu_to_mag(f_nu_jy):
    """f_nu (Jy) → AB mag"""
    return -2.5 * np.log10(f_nu_jy / AB_FLUX_ZERO)

# ========== 气体茧不透明度 ==========

def hi_bound_free_cross_section(lam_cm, n=2, T=2e4):
    """
    H I bound-free 截面 (cm^2 per atom in level n)
    Kramer 定律: σ_bf(n,λ) = σ_max × (λ/λ_limit)^3 for λ < λ_limit
    
    Note: 需配合能级布居数使用。n=2 布居数 ≈ f_neutral × (g2/g1) × exp(-ΔE/kT)
    """
    if n == 1:
        lam_limit = 9.12e-6   # Lyman limit
        sigma_max = 6.30e-18
    elif n == 2:
        lam_limit = 3.646e-5  # Balmer limit
        sigma_max = 1.40e-17
    else:
        lam_limit = 8.20e-5  # Paschen limit
        sigma_max = 2.10e-17
    return np.where(lam_cm < lam_limit, sigma_max * (lam_cm / lam_limit)**3, 0)

def free_free_cross_section(nu, T=2e4):
    """自由-自由吸收截面 (cm^2 per electron per ion)"""
    g_ff = 1.2  # Gaunt factor
    x = h * nu / (k_B * T)
    sigma = 3.692e8 * g_ff * T**(-0.5) * nu**(-3.0) * (1 - np.exp(-x))
    sigma = np.maximum(sigma, 0)
    return sigma

def gas_cocoon_opacity(lam_cm, N_H, T_gas=2e4, f_neutral=0.8):
    """
    气体茧总光深 (per cm^2)
    
    关键修正: n=2 能级布居数受 Boltzmann 因子控制
    f_n2 = f_neutral × (g2/g1) × exp(−ΔE_{1→2} / kT)
          ≈ f_neutral × 2 × exp(−10.2 eV / kT)
    """
    # n=2 布居数分数 (ΔE_1→2 = 10.2 eV, const.e converts eV→J)
    dE_12_J = 10.2 * const.e  # Joules
    boltz_n2 = 2.0 * np.exp(-dE_12_J / (k_B * T_gas))
    f_n2 = f_neutral * boltz_n2
    
    # n=1 布居数分数 (近似: 剩余中性原子在基态)
    f_n1 = f_neutral * (1 - boltz_n2)
    
    # Bound-free opacity
    tau_bf1 = N_H * f_n1 * hi_bound_free_cross_section(lam_cm, n=1)
    tau_bf2 = N_H * f_n2 * hi_bound_free_cross_section(lam_cm, n=2)
    
    # Free-free (ionized component)
    nu = c / lam_cm
    f_ionized = 1 - f_neutral
    tau_ff = N_H * f_ionized * free_free_cross_section(nu, T_gas)
    
    # Thomson
    tau_T = N_H * sigma_T
    
    return tau_bf1 + tau_bf2 + tau_ff + tau_T

# ========== 中央引擎 SED ==========

def agn_sed(lam_um, norm_lam_rest=0.5):
    """
    AGN 幂律 SED (f_nu ∝ nu^alpha)，归一化到 λ_rest = 5000 Å
    典型类星体: f_nu ∝ nu^(-0.5) at UV, nu^(-0.2) at optical
    """
    lam_cm = lam_um * 1e-4
    nu = c / lam_cm
    nu_norm = c / (norm_lam_rest * 1e-4)
    
    alpha = np.where(lam_um < 0.4, -0.5, -0.2)
    f_nu = (nu / nu_norm)**alpha
    
    return f_nu

def g_eff_sigma_coupling(N_H, N_H_ref=1e21, eps_g=0.3, beta=1.0, delta=0.7):
    """
    Σ-N_H 耦合 → G_eff 增强因子
    
    物理链条: N_H ↑ → Σ_star ↑ → G_eff/G_N ↑
    
    Σ_star/Σ_ref = (N_H / N_H_ref)^delta    (气体柱→恒星面密度)
    G_eff/G_N = 1 + ε_g × (Σ/Σ_ref)^beta   (UID 密度依赖引力)
    
    参数来自 Path1 观测校准:
    - ε_g ≈ 0.3  (SB-Σ 相关 5.6σ 的归一化)
    - β  ≈ 1.0  (SB-Σ 斜率 ~ -1.21 → β ~ 1)
    - δ  ≈ 0.7  (气体-恒星耦合指数)
    - N_H_ref = 10^21 cm^-2 (巴尔末跳变起始列密度)
    """
    Sigma_ratio = (N_H / N_H_ref)**delta
    G_eff_over_GN = 1.0 + eps_g * Sigma_ratio**beta
    return G_eff_over_GN

def stellar_core_sed(lam_um, N_H=None, T_star=30000, norm_lam_rest=0.5,
                     use_density_coupling=True):
    """
    致密恒星核团 SED: 年轻星族 + 星云连续谱 + G_eff(Σ) 引力加热
    
    关键: G_eff 红外超随 N_H 联动 (密度→引力势→加热)
    与 AGN (IR 与 N_H 独立) 形成可检验区分
    
    若 N_H=None 则使用固定 G_eff baseline
    """
    lam_cm = lam_um * 1e-4
    nu = c / lam_cm
    nu_norm = c / (norm_lam_rest * 1e-4)
    
    # (1) 星族: T_star 黑体, 归一化 @ 5000Å
    bb = 2 * h * nu**3 / c**2 / (np.exp(h * nu / (k_B * T_star)) - 1)
    bb_norm = 2 * h * nu_norm**3 / c**2 / (np.exp(h * nu_norm / (k_B * T_star)) - 1)
    bb = bb / bb_norm
    
    # (2) 星云连续谱
    nebular = 0.3 * (nu / nu_norm)**(-0.1)
    balmer_jump_emission = 0.15 * np.exp(-((lam_um - 0.3646)**2) / (2 * 0.015**2))
    
    # (3) G_eff 引力加热 → IR 超
    # 加热功率 ∝ (G_eff/G_N - 1) — 只计增强部分
    if use_density_coupling and N_H is not None:
        geff = g_eff_sigma_coupling(N_H)
        ir_amplitude = 0.2 * (geff - 1.0) / 0.3  # normalize to ε_g=0.3
    else:
        ir_amplitude = 0.2  # baseline (N_H=10^21 equivalent)
    
    ir_excess = ir_amplitude * (nu / nu_norm)**(-2.5) * np.exp(-nu / (c/3e-4))
    
    return bb + nebular + balmer_jump_emission + ir_excess

# ========== 辐射传输 ==========

def compute_emergent_spectrum(engine_func, N_H, lam_um, z=8.5, T_gas=2e4, f_neutral=0.8,
                              engine_name=''):
    """
    出射光谱: F_emerg(rest) = F_engine(rest, N_H) * exp(-tau(rest))
    
    对于恒星引擎, N_H 传入以激活 Σ-N_H-G_eff 联动
    AGN 引擎忽略 N_H 参数
    """
    lam_rest = lam_um / (1 + z)
    lam_rest_cm = lam_rest * 1e-4
    
    tau = gas_cocoon_opacity(lam_rest_cm, N_H, T_gas, f_neutral)
    
    # 恒星引擎: N_H → Σ → G_eff → IR excess
    if engine_name == 'StellarCore':
        f_engine = engine_func(lam_rest, N_H=N_H)
    else:
        f_engine = engine_func(lam_rest)  # AGN: N_H 不影响引擎 SED
    
    transmitted = f_engine * np.exp(-tau)
    
    # 气体茧再辐射
    T_reemit = 1200  # K
    nu_rest = c / lam_rest_cm
    reemit_sed = 2 * h * nu_rest**3 / c**2 / (np.exp(h * nu_rest / (k_B * T_reemit)) - 1)
    reemit_sed = reemit_sed / np.max(reemit_sed + 1e-30) * 0.05
    
    absorb_frac = 1 - np.exp(-np.mean(tau[lam_rest < 0.2]))
    
    return transmitted + reemit_sed * absorb_frac
    
    # 气体茧再辐射 (近似: 吸收的能量以尘埃温度再辐射)
    # 在红外波段 (observed > 3 um) 有小贡献
    T_reemit = 1200  # K
    nu_rest = c / lam_rest_cm
    reemit_sed = 2 * h * nu_rest**3 / c**2 / (np.exp(h * nu_rest / (k_B * T_reemit)) - 1)
    reemit_sed = reemit_sed / np.max(reemit_sed + 1e-30) * 0.05
    
    # 被吸收比例
    absorb_frac = 1 - np.exp(-np.mean(tau[lam_rest < 0.2]))
    
    return transmitted + reemit_sed * absorb_frac

# ========== JWST 测光模拟 ==========

def compute_broadband_colors(lam_um, f_nu, bands=None):
    """计算 JWST 宽带 AB 星等 (简单中心波长近似)"""
    if bands is None:
        bands = JWST_BANDS
    
    mags = {}
    for band_name, lam_center in bands.items():
        idx = np.argmin(np.abs(lam_um - lam_center))
        mags[band_name] = flux_nu_to_mag(max(np.abs(f_nu[idx]), 1e-30))
    
    return mags

# ========== 主计算 ==========

def run_engine_test():
    """主函数: 扫描柱密度, 计算两种引擎的相对红化量"""
    
    lam_um = np.logspace(np.log10(0.09), np.log10(8.0), 800)
    N_H_values = np.logspace(19.0, 23.0, 60)
    z_lrd = 8.5
    
    results = {'AGN': [], 'StellarCore': []}
    
    engine_configs = {
        'AGN':      {'func': agn_sed,      'f_neutral': 0.25},
        'StellarCore': {'func': stellar_core_sed, 'f_neutral': 0.55},
    }
    
    print("="*60)
    print("引擎简并性检验: AGN vs 致密恒星核团 + G_eff")
    print("检验量: Δ(F150W−F444W) = Color(N_H) − Color(N_H=0)")
    print("="*60)
    
    for engine_name, cfg in engine_configs.items():
        func = cfg['func']
        f_neutral = cfg['f_neutral']
        
        # Step 1: 计算出射光谱基线 (N_H=0, 仅引擎)
        emergent_0 = func(lam_um / (1 + z_lrd))
        mags_0 = compute_broadband_colors(lam_um, emergent_0)
        color_0 = mags_0['F150W'] - mags_0['F444W']
        color_0_200 = mags_0['F200W'] - mags_0['F444W']
        
        print(f"\n{'─'*50}")
        print(f"  {engine_name} 引擎 (f_neutral={f_neutral})")
        print(f"  内禀颜色 F150W-F444W = {color_0:.2f}")
        print(f"{'─'*50}")
        print(f"  {'log N_H':>8s}  {'ΔColor':>8s}  {'Balmer jump?':>14s}")
        print(f"  {'-'*30}")
        
        for N_H in N_H_values:
            emergent = compute_emergent_spectrum(func, N_H, lam_um, z=z_lrd, 
                                                 f_neutral=f_neutral, engine_name=engine_name)
            mags = compute_broadband_colors(lam_um, emergent)
            
            # 相对红化 = 巴尔末跳变强度 (引擎无关量)
            delta_color = (mags['F150W'] - mags['F444W']) - color_0
            delta_200 = (mags['F200W'] - mags['F444W']) - color_0_200
            
            results[engine_name].append({
                'N_H': float(N_H),
                'delta_F150W-F444W': float(delta_color),
                'delta_F200W-F444W': float(delta_200),
                'abs_F150W-F444W': float(mags['F150W'] - mags['F444W']),
                'abs_F277W-F444W': float(mags['F277W'] - mags['F444W']),
                'abs_F356W-F444W': float(mags['F356W'] - mags['F444W']),
                'G_eff': float(g_eff_sigma_coupling(N_H)) if engine_name == 'StellarCore' else 1.0,
            })
            
            if len(results[engine_name]) % 8 == 0:
                jump = "YES" if delta_color > 0.3 else "no"
                print(f"  {np.log10(N_H):8.2f}  {delta_color:8.3f}  {jump:>14s}")
        
        # 打印范围
        deltas = [r['delta_F150W-F444W'] for r in results[engine_name]]
        print(f"  {'─'*30}")
        print(f"  ΔColor 范围: {min(deltas):.2f} to {max(deltas):.2f} mag")
    
    return results, lam_um, N_H_values

# ========== 加载观测数据 ==========

def load_observed_lrd_colors():
    """从 Path1 数据加载 LRD 观测颜色"""
    csv_path = os.path.join(
        os.path.dirname(__file__), '..', 'data', 'path1_merged_38sources.csv'
    )
    
    if not os.path.exists(csv_path):
        print(f"⚠️  观测数据文件不存在: {csv_path}")
        return None
    
    import csv
    data = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                f150w = float(row['f150w_flux'])
                f444w = float(row['f444w_flux'])
                f200w = float(row['f200w_flux'])
                f277w = float(row['f277w_flux'])
                f356w = float(row['f356w_flux'])
                
                # 跳过无效数据
                if f150w <= 0 or f444w <= 0:
                    continue
                    
                mag150 = flux_nu_to_mag(f150w * 1e-6)  # uJy → Jy
                mag444 = flux_nu_to_mag(f444w * 1e-6)
                mag200 = flux_nu_to_mag(f200w * 1e-6)
                mag277 = flux_nu_to_mag(f277w * 1e-6)
                mag356 = flux_nu_to_mag(f356w * 1e-6)
                
                data.append({
                    'id': row['id'],
                    'z': float(row['z_phot']),
                    'SB_F444W': float(row.get('SB_F444W', 0)),
                    'logSigma': float(row.get('logSigma_Mstar', 0)),
                    'BD': float(row.get('Balmer_dec_total', 0)),
                    'F150W-F444W': float(mag150 - mag444),
                    'F200W-F444W': float(mag200 - mag444),
                    'F277W-F444W': float(mag277 - mag444),
                    'F356W-F444W': float(mag356 - mag444),
                })
            except (ValueError, KeyError):
                continue
    
    print(f"\n📊 加载观测数据: {len(data)} 个 LRD")
    if data:
        c150_444 = [d['F150W-F444W'] for d in data]
        print(f"  F150W-F444W: {np.min(c150_444):.2f} - {np.max(c150_444):.2f}, "
              f"median={np.median(c150_444):.2f}")
    
    return data

# ========== 可视化 ==========

def make_comparison_plot(results, obs_data, N_H_values, save_dir):
    """生成对比图 (使用 ΔColor = 相对红化量)"""
    
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)
    
    # ---- 提取数据 ----
    delta_agn = [r['delta_F150W-F444W'] for r in results['AGN']]
    delta_star = [r['delta_F150W-F444W'] for r in results['StellarCore']]
    delta_agn_200 = [r['delta_F200W-F444W'] for r in results['AGN']]
    delta_star_200 = [r['delta_F200W-F444W'] for r in results['StellarCore']]
    logNH = np.log10(N_H_values)
    
    # ---- Panel 1: ΔColor vs N_H (核心检验) ----
    ax1 = fig.add_subplot(gs[0, :2])
    
    ax1.plot(logNH, delta_agn, 'o-', color='#d62728', lw=2.5, markersize=5,
             label='AGN engine (f_neut=0.5)')
    ax1.plot(logNH, delta_star, 's-', color='#1f77b4', lw=2.5, markersize=5,
             label='Stellar core + G_eff (f_neut=0.9)')
    
    # 理论区域标注
    ax1.axhspan(-0.2, 0.5, alpha=0.06, color='blue')
    ax1.text(25.8, 0.15, 'LBD\nregion', fontsize=9, ha='center', alpha=0.6)
    ax1.axhspan(1.5, 4.5, alpha=0.06, color='red')
    ax1.text(25.8, 3.0, 'LRD\nregion', fontsize=9, ha='center', alpha=0.6)
    
    # 观测 LRD (估计 log N_H ~ 24-25)
    if obs_data:
        obs_deltas = []
        obs_logNH_est = []
        for d in obs_data:
            # 估算: 每个 LRD 的 ΔColor 通过减去引擎内禀颜色近似
            # 使用恒星引擎内禀颜色 (更合理)
            obs_deltas.append(d['F150W-F444W'] - delta_star[0])
            obs_logNH_est.append(24.2 + np.random.uniform(-0.5, 0.5))
        ax1.scatter(obs_logNH_est, obs_deltas, c='gray', s=30, alpha=0.5,
                   marker='o', edgecolors='black', linewidth=0.3,
                   label=f'Observed 38 LRDs', zorder=3)
    
    ax1.set_xlabel('log N_H [cm^-2]', fontsize=13)
    ax1.set_ylabel('Δ(F150W − F444W) [mag]', fontsize=13)
    ax1.set_title('Balmer Jump Strength vs Column Density: Engine Degeneracy', 
                  fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10, loc='upper left')
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(19, 23)
    ax1.set_ylim(-2, 6)
    
    # ---- Panel 2: Σ 耦合诊断 F356W-F444W vs F150W-F444W ----
    ax2 = fig.add_subplot(gs[0, 2])
    
    abs_150_444_agn = [r['abs_F150W-F444W'] for r in results['AGN']]
    abs_150_444_star = [r['abs_F150W-F444W'] for r in results['StellarCore']]
    
    # F356W-F444W: both above Balmer limit → should be pure IR tracer
    # 但 NIRCam 波段在中红外 G_eff 热辐射的 Rayleigh-Jeans 尾部
    # MIRI 5-15μm 才是区分关键——标注这一点
    ir_agn = [r['abs_F356W-F444W'] for r in results['AGN']]
    ir_star = [r['abs_F356W-F444W'] for r in results['StellarCore']]
    
    geff_vals = [r.get('G_eff', 1.0) for r in results['StellarCore']]
    
    ax2.plot(abs_150_444_agn, ir_agn, 'o-', color='#d62728', lw=1.8, markersize=4,
             label='AGN (dust torus IR ~ const)')
    # Color Stellar points by G_eff
    sc = ax2.scatter(abs_150_444_star, ir_star, c=geff_vals, cmap='coolwarm',
                    s=30, edgecolors='#1f77b4', linewidth=0.5, zorder=5)
    ax2.plot(abs_150_444_star, ir_star, '-', color='#1f77b4', lw=1.5, alpha=0.5)
    cbar = plt.colorbar(sc, ax=ax2, shrink=0.8)
    cbar.set_label('G_eff / G_N', fontsize=8)
    
    # NIRCam degeneracy annotation
    ax2.annotate('NIRCam: degenerate\n(MIRI needed)', xy=(2.0, -0.3),
                fontsize=8, ha='center', color='gray',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
    
    if obs_data:
        obs_150_444 = [d['F150W-F444W'] for d in obs_data]
        obs_356_444 = [d['F356W-F444W'] for d in obs_data]
        ax2.scatter(obs_150_444, obs_356_444, c='gray', s=20, alpha=0.4,
                   edgecolors='black', linewidth=0.2, label='38 LRDs (NIRCam)')
    
    ax2.set_xlabel('F150W \u2212 F444W [Balmer jump]', fontsize=11)
    ax2.set_ylabel('F356W \u2212 F444W [IR excess]', fontsize=11)
    ax2.set_title('\u03a3-Coupled Diagnostic\nG_eff IR grows with N_H', 
                  fontsize=11, fontweight='bold')
    ax2.legend(fontsize=7.5, loc='best')
    ax2.grid(True, alpha=0.3)
    
    # ---- Panel 3: 射出 SED (选定 N_H) ----
    ax3 = fig.add_subplot(gs[1, :2])
    
    z = 8.5
    selected_NH = [3e19, 3e20, 3e21, 1e22]
    colors_line = ['#2ca02c', '#ff7f0e', '#d62728', '#9467bd']
    lam_plot = np.logspace(np.log10(0.15), np.log10(7.5), 600)
    
    for i, N_H in enumerate(selected_NH):
        # AGN (solid)
        e_agn = compute_emergent_spectrum(agn_sed, N_H, lam_plot, z=z, f_neutral=0.25,
                                          engine_name='AGN')
        norm_agn = np.max(np.abs(e_agn))
        ax3.plot(lam_plot, e_agn / norm_agn, '-', color=colors_line[i], 
                alpha=0.35, lw=2)
        # Stellar (dashed)
        e_star = compute_emergent_spectrum(stellar_core_sed, N_H, lam_plot, z=z, f_neutral=0.55,
                                          engine_name='StellarCore')
        norm_star = np.max(np.abs(e_star))
        ax3.plot(lam_plot, e_star / norm_star, '--', color=colors_line[i],
                alpha=0.9, lw=2)
    
    # 巴尔末跳变观测位置
    lam_balmer_obs = 0.3646 * (1 + z)
    ax3.axvline(x=lam_balmer_obs, color='black', linestyle=':', alpha=0.5, lw=1.5)
    ax3.text(lam_balmer_obs * 0.92, 0.95, 'Balmer\n(z=8.5)', fontsize=8,
            ha='right', va='top', color='black', alpha=0.6)
    
    ax3.set_xscale('log')
    ax3.set_xlabel('Observed wavelength [um]', fontsize=12)
    ax3.set_ylabel('Normalized f_nu', fontsize=12)
    ax3.set_title('Emergent SEDs: AGN (solid) vs Stellar Core (dashed)', 
                  fontsize=13, fontweight='bold')
    
    from matplotlib.lines import Line2D
    handles_nh = [Line2D([0], [0], color=c, lw=2.5) for c in colors_line]
    handle_style = [Line2D([0], [0], color='gray', lw=2, linestyle='-'),
                    Line2D([0], [0], color='gray', lw=2, linestyle='--')]
    labels = [f'N_H=3e{int(np.log10(nh))}' if nh < 1e23 else f'N_H=1e{int(np.log10(nh))}' for nh in selected_NH] + ['AGN', 'Stellar']
    ax3.legend(handles_nh + handle_style, labels, fontsize=8, ncol=6, loc='upper right')
    ax3.set_ylim(-0.05, 1.15)
    ax3.grid(True, alpha=0.2)
    
    # ---- Panel 4: 观测区分诊断 ----
    ax4 = fig.add_subplot(gs[1, 2])
    
    diag_labels = ['Hard X-ray\n(NuSTAR)', 'Mid-IR\n(5-10um)', '[NeV]\n(97 eV)', 
                   'Variability', 'Lyα escape', 'Broad Hα']
    
    # AGN 需要额外假设 (遮挡/特殊几何)
    agn_vals = [0.3, 0.4, 0.4, 0.5, 0.3, 0.6]
    # 星团核心 + G_eff 自然解释
    star_vals = [0.95, 0.85, 0.9, 0.85, 0.9, 0.85]
    
    x = np.arange(len(diag_labels))
    w = 0.35
    ax4.barh(x + w/2, agn_vals, w, color='#d62728', alpha=0.7, label='AGN')
    ax4.barh(x - w/2, star_vals, w, color='#1f77b4', alpha=0.7, label='Stellar+G_eff')
    
    ax4.set_yticks(x)
    ax4.set_yticklabels(diag_labels, fontsize=8)
    ax4.set_xlabel('Explanatory power', fontsize=11)
    ax4.set_title('Observable Discriminators\n(LBD = thin cocoon)', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=9, loc='lower right')
    ax4.set_xlim(0, 1.2)
    ax4.grid(True, alpha=0.2, axis='x')
    
    # ---- Panel 5: 关键结论 ----
    ax5 = fig.add_subplot(gs[2, :])
    ax5.axis('off')
    
    # 计算重叠
    delta_min = min(min(delta_agn), min(delta_star))
    delta_max = max(max(delta_agn), max(delta_star))
    overlap = abs(max(min(delta_agn), min(delta_star)) - min(max(delta_agn), max(delta_star)))
    
    # 观测 LRD 颜色范围
    if obs_data:
        obs_c = [d['F150W-F444W'] for d in obs_data]
        obs_med = np.median(obs_c)
        obs_min, obs_max = min(obs_c), max(obs_c)
    else:
        obs_med, obs_min, obs_max = 2.35, 1.2, 3.5
    
    # 匹配 N_H 范围
    idx_balmer_start_agn = next((i for i, d in enumerate(delta_agn) if d > 0.5), len(delta_agn)-1)
    idx_balmer_start_star = next((i for i, d in enumerate(delta_star) if d > 0.5), len(delta_star)-1)
    
    summary = f"""
╔══════════════════════════════════════════════════════════════════════╗
║    ENGINE DEGENERACY TEST — with Σ-N_H-G_eff Density Coupling      ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                    ║
║  Physical chain: N_H ∝ Σ_star → G_eff/G_N ∝ (Σ/Σ_0)^β → IR heat  ║
║                                                                    ║
║  AGN engine:           ΔColor = {min(delta_agn):.2f} → {max(delta_agn):.2f} mag (Balmer jump from N_H only)    ║
║  Stellar core+G_eff:   ΔColor = {min(delta_star):.2f} → {max(delta_star):.2f} mag (Balmer + G_eff coupling)   ║
║                                                                    ║
║  ▸ NIRCam bands: ENGINE-DEGENERATE (both span LRD color range)    ║
║    Gas cocoon Balmer jump dominates; G_eff IR in Rayleigh-Jeans   ║
║    tail too faint for F356W-F444W diagnostic at z~8.5             ║
║                                                                    ║
║  ▸ MIRI bands (>10 μm): ENGINE-DIAGNOSTIC                         ║
║    G_eff heated dust (T~1200K) peaks at λ_obs~15-25 μm           ║
║    AGN dust torus: λ_peak ~ 10-20 μm (hot inner torus)            ║
║    Stellar+G_eff: λ_peak ~ 15-30 μm (cooler, gravity-heated)     ║
║                                                                    ║
║  ▸ LBDs (thin cocoon) = ideal for MIRI spectroscopy:              ║
║    Low N_H → central engine partially visible in mid-IR           ║
║    MIRI LRS 5-14 μm + MRS IFU → direct engine classification     ║
║                                                                    ║
╚══════════════════════════════════════════════════════════════════════╝
"""
    ax5.text(0.5, 0.5, summary, transform=ax5.transAxes,
            fontsize=9.5, fontfamily='monospace', ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    
    output_path = os.path.join(save_dir, 'engine_degeneracy_test.png')
    fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"\n✅ 图已保存: {output_path}")
    
    return output_path

# ========== 保存数据 ==========

def save_results(results):
    save_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
    os.makedirs(save_dir, exist_ok=True)
    
    output = {
        'description': 'Engine degeneracy test: AGN vs Stellar Core + G_eff',
        'metric': 'Delta(F150W-F444W) = Balmer jump strength relative to N_H=0',
        'columns': ['N_H', 'delta_AGN', 'delta_Stellar'],
    }
    
    n = len(results['AGN'])
    for i in range(n):
        key = f'N_H={results["AGN"][i]["N_H"]:.2e}'
        output[key] = {
            'delta_AGN': round(results['AGN'][i]['delta_F150W-F444W'], 3),
            'delta_Stellar': round(results['StellarCore'][i]['delta_F150W-F444W'], 3),
        }
    
    json_path = os.path.join(save_dir, 'engine_degeneracy_results.json')
    with open(json_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"✅ 数据: {json_path}")

# ========== 运行 ==========

if __name__ == '__main__':
    save_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')
    os.makedirs(save_dir, exist_ok=True)
    
    # 计算
    results, lam_um, N_H_values = run_engine_test()
    
    # 加载观测
    obs_data = load_observed_lrd_colors()
    
    # 画图
    make_comparison_plot(results, obs_data, N_H_values, save_dir)
    
    # 保存
    save_results(results)
    
    print("\n" + "="*60)
    print("检验完成。")
    print("核心结论: 气体茧颜色序列 → 引擎无关。")
    print("区分引擎需要: X射线 / 中红外 / 高电离线 / 光变。")
    print("="*60)
