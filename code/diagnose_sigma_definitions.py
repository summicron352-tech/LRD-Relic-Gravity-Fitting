#!/usr/bin/env python3
"""
诊断脚本：两种 Σ 定义在各证据链上的对比
============================================
Σ_A = APJL 通量代理: const - 0.16*mag - 2*log10(r_eff_pc)
      = -0.4*F444W_mag - 2*log10(r_eff_pc) + offset
      (correlation-invariant form; offset drops out)

Σ_B = 物理面密度:    logMstar - log10(π * r_eff_pc²)  [M☉/pc²]

输出：
  1. Σ_A 与 Σ_B 的相关系数（全局 + 分样本）
  2. 每条证据链用两种定义的结果对比
  3. COSMOSWeb 664k 的 Σ 定义溯源
"""

import numpy as np
import os, sys, json
from scipy.stats import spearmanr, pearsonr
from scipy.stats import rankdata

# ===========================
# Paths
# ===========================
data_dir = '/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/data'
cosmos_dir = '/Users/tanxin/Desktop/数据处理/COSMOS2025'

def psc(x, y, ctrls):
    """Partial Spearman correlation"""
    rx = rankdata(x); ry = rankdata(y)
    for c in ctrls:
        rc = rankdata(c)
        b = np.cov(rx, rc)[0,1] / np.var(rc)
        rx = rx - b * rc
        b = np.cov(ry, rc)[0,1] / np.var(rc)
        ry = ry - b * rc
    return spearmanr(rx, ry)[0]

def sigma_A_from_mag(mag, r_pc):
    """APJL F444W flux proxy Σ. Returns relative log Σ (offset irrelevant for correlation)."""
    return -0.4 * mag - 2 * np.log10(r_pc)

def sigma_B_from_mass(logM, r_pc):
    """Physical stellar surface density [M☉/pc²]."""
    return logM - np.log10(np.pi * r_pc**2)

def safe_filter(*arrays):
    """Keep only rows where all arrays are finite."""
    n = len(arrays[0])
    mask = np.ones(n, dtype=bool)
    for a in arrays:
        mask &= np.isfinite(a)
    return tuple(a[mask] for a in arrays)

def print_corr(label, x, y, ctrls=None):
    """Print Spearman correlation results."""
    x, y = safe_filter(x, y)[:2]
    if len(x) < 3:
        print(f"  {label}: N<3, 跳过")
        return None
    r = spearmanr(x, y)[0]
    p = spearmanr(x, y)[1]
    result = f"ρ={r:+.4f}, p={p:.4f}, N={len(x)}"
    print(f"  {label}: {result}")
    return {'rho': r, 'p': p, 'N': len(x)}

# ====================================================================
# 1. PATH1 38 SOURCES
# ====================================================================
print("=" * 72)
print("1. PATH1 38 SOURCES (deGraaff BW + Kokorev SB/Σ)")
print("=" * 72)

path1 = np.genfromtxt(
    os.path.join(data_dir, 'path1_merged_38sources.csv'),
    delimiter=',', names=True, dtype=None, encoding='utf-8')

f444_p1 = np.array([float(x) for x in path1['F444W_mag']])
dA_p1 = np.array([float(x) for x in path1['d_A_pc']])
r_arcsec_p1 = np.array([float(x) for x in path1['r_eff_arcsec']])
r_pc_p1 = r_arcsec_p1 * dA_p1  # r_eff in pc
logM_p1 = np.array([float(x) for x in path1['logMstar_best']])
bd_p1 = np.array([float(x) for x in path1['Balmer_dec_total']])
logSigma_existing_p1 = np.array([float(x) for x in path1['logSigma_Mstar']])
sb_p1 = np.array([float(x) for x in path1['SB_F444W']])
z_p1 = np.array([float(x) for x in path1['z_phot']])

# Clean
valid_p1 = (np.isfinite(f444_p1) & np.isfinite(r_pc_p1) & np.isfinite(logM_p1) &
            np.isfinite(bd_p1) & (bd_p1 > 2) & (bd_p1 < 30) &
            (logM_p1 > 6) & (logM_p1 < 14) & (r_pc_p1 > 0))
f444_p1, r_pc_p1, logM_p1, bd_p1, sb_p1, z_p1, logSigma_existing_p1 = [
    a[valid_p1] for a in [f444_p1, r_pc_p1, logM_p1, bd_p1, sb_p1, z_p1, logSigma_existing_p1]]

# Compute both Σ
sigA_p1 = sigma_A_from_mag(f444_p1, r_pc_p1)
sigB_p1 = sigma_B_from_mass(logM_p1, r_pc_p1)

print(f"\n  有效样本: N = {len(sigA_p1)}")
print(f"  Σ_A 范围: [{sigA_p1.min():.2f}, {sigA_p1.max():.2f}]")
print(f"  Σ_B 范围: [{sigB_p1.min():.2f}, {sigB_p1.max():.2f}]")
print(f"  Σ_A vs Σ_B 相关: ρ = {spearmanr(sigA_p1, sigB_p1)[0]:+.4f}")
print(f"  已有 logSigma_Mstar vs Σ_B: ρ = {spearmanr(logSigma_existing_p1, sigB_p1)[0]:+.4f}")
print(f"  (已有值应与 Σ_B ≈ 一致，验证计算正确)")

print("\n--- Σ vs BD ---")
rA_bd = print_corr("Σ_A vs BD", sigA_p1, bd_p1)
rB_bd = print_corr("Σ_B vs BD", sigB_p1, bd_p1)
print_corr("已有logSigma_Mstar vs BD", logSigma_existing_p1, bd_p1)

print("\n--- Σ vs BD (偏相关, 控制z+M*) ---")
rA_par = psc(bd_p1, sigA_p1, [z_p1, logM_p1])
rB_par = psc(bd_p1, sigB_p1, [z_p1, logM_p1])
print(f"  Σ_A 偏相关: ρ={rA_par:+.4f}")
print(f"  Σ_B 偏相关: ρ={rB_par:+.4f}")

print("\n--- M* vs BD ---")
print_corr("logM* vs BD", logM_p1, bd_p1)

print("\n--- SB_F444W vs Σ ---")
print_corr("SB vs Σ_A", sb_p1, sigA_p1)
print_corr("SB vs Σ_B", sb_p1, sigB_p1)

# ====================================================================
# 2. KOKOREV 260 SOURCES
# ====================================================================
print("\n" + "=" * 72)
print("2. KOKOREV 260 SOURCES")
print("=" * 72)

kokorev = np.genfromtxt(
    os.path.join(data_dir, 'kokorev_260_sb.csv'),
    delimiter=',', names=True, dtype=None, encoding='utf-8')

f444_k = np.array([float(x) for x in kokorev['F444W_mag']])
r_k = np.array([float(x) for x in kokorev['r_eff_50_phys']])
logM_k = np.array([float(x) for x in kokorev['logMstar_best']])
sb_k = np.array([float(x) for x in kokorev['SB_F444W']])
logSigma_existing_k = np.array([float(x) for x in kokorev['logSigma_Mstar']])
z_k = np.array([float(x) for x in kokorev['z_phot']])

valid_k = (np.isfinite(f444_k) & np.isfinite(r_k) & np.isfinite(logM_k) &
           np.isfinite(sb_k) & np.isfinite(z_k) &
           (r_k > 0) & (logM_k > 6) & (logM_k < 14) &
           (sb_k > 15) & (sb_k < 35))
f444_k, r_k, logM_k, sb_k, z_k, logSigma_existing_k = [
    a[valid_k] for a in [f444_k, r_k, logM_k, sb_k, z_k, logSigma_existing_k]]

sigA_k = sigma_A_from_mag(f444_k, r_k)
sigB_k = sigma_B_from_mass(logM_k, r_k)

print(f"\n  有效样本: N = {len(sigA_k)}")
print(f"  Σ_A 范围: [{sigA_k.min():.2f}, {sigA_k.max():.2f}]")
print(f"  Σ_B 范围: [{sigB_k.min():.2f}, {sigB_k.max():.2f}]")
print(f"  Σ_A vs Σ_B 相关: ρ = {spearmanr(sigA_k, sigB_k)[0]:+.4f}")

print("\n--- SB_F444W vs Σ ---")
print_corr("SB vs Σ_A", sb_k, sigA_k)
print_corr("SB vs Σ_B", sb_k, sigB_k)
print_corr("SB vs 已有logSigma_Mstar", sb_k, logSigma_existing_k)

# ====================================================================
# 3. YANAGISAWA 15 SOURCES (BD-Hα overlap)
# ====================================================================
print("\n" + "=" * 72)
print("3. YANAGISAWA 15 SOURCES (BD-Hα overlap)")
print("=" * 72)

yana = np.genfromtxt(
    os.path.join(data_dir, 'yanagisawa_overlap_15sources.csv'),
    delimiter=',', names=True, dtype=None, encoding='utf-8')

bd_y = np.array([float(x) for x in yana['BD']])
logLHa_y = np.array([float(x) for x in yana['logLHa']])
logL5100_y = np.array([float(x) for x in yana['logL5100']])
z_y = np.array([float(x) for x in yana['z_spec']])
sb_y = np.array([float(x) for x in yana['SB_F444W']])
f444_y = np.array([float(x) for x in yana['F444W_mag']])
logM_y = np.array([float(x) for x in yana['logMstar']])
logSigma_existing_y = np.array([float(x) for x in yana['logSigma_Mstar']])
ha_enhance_y = np.array([float(x) for x in yana['Ha_enhancement_vs_type1']])

# Compute Hα enhancement: L(Hα)/L(5100) in linear then compare to Type 1
lha_over_l5100 = 10**(logLHa_y - logL5100_y)

# This sample doesn't have r_eff, so we can't compute Σ directly
# But we can check if the existing logSigma_Mstar correlates
valid_y = (np.isfinite(bd_y) & np.isfinite(ha_enhance_y) &
           np.isfinite(logSigma_existing_y) & (bd_y > 2))
bd_y, ha_enhance_y, logSigma_existing_y, sb_y, logM_y, z_y = [
    a[valid_y] for a in [bd_y, ha_enhance_y, logSigma_existing_y, sb_y, logM_y, z_y]]

print(f"\n  有效样本: N = {len(bd_y)}")
print(f"  (注: 此样本无 r_eff，无法计算 Σ_A/Σ_B，仅用已有 logSigma_Mstar)")

print("\n--- BD vs Hα增强 ---")
print_corr("BD vs Hα增强", bd_y, ha_enhance_y)
print_corr("BD vs LHa/L5100", bd_y, lha_over_l5100[:len(bd_y)])
print_corr("BD vs logSigma(existing)", bd_y, logSigma_existing_y)

print("\n--- logSigma vs Hα增强 ---")
print_corr("logSigma(existing) vs Hα增强", logSigma_existing_y, ha_enhance_y)

# ====================================================================
# 4. WITHERS 9 COMPACT REGs
# ====================================================================
print("\n" + "=" * 72)
print("4. WITHERS 9 COMPACT REGs")
print("=" * 72)

with_data = np.genfromtxt(
    os.path.join(data_dir, 'withers_9_regs.csv'),
    delimiter=',', names=True, dtype=None, encoding='utf-8')

f444_w = np.array([float(x) for x in with_data['F444W_mag']])
r_max_kpc_w = np.array([float(x) for x in with_data['r_eff_max_kpc']])
r_max_pc_w = r_max_kpc_w * 1000  # kpc → pc
logM_w = np.array([float(x) for x in with_data['logMstar']])
sb_min_w = np.array([float(x) for x in with_data['SB_F444W_min']])
logSigma_existing_w = np.array([float(x) for x in with_data['logSigma_Mstar_min']])

valid_w = (np.isfinite(f444_w) & np.isfinite(r_max_pc_w) &
           np.isfinite(logM_w) & (r_max_pc_w > 0))
f444_w, r_max_pc_w, logM_w, sb_min_w, logSigma_existing_w = [
    a[valid_w] for a in [f444_w, r_max_pc_w, logM_w, sb_min_w, logSigma_existing_w]]

sigA_w = sigma_A_from_mag(f444_w, r_max_pc_w)
sigB_w = sigma_B_from_mass(logM_w, r_max_pc_w)

print(f"\n  有效样本: N = {len(sigA_w)}")
print(f"  ⚠ 所有 r_eff 均为 PSF 上限，Σ 为下限")
print(f"  Σ_A 范围: [{sigA_w.min():.2f}, {sigA_w.max():.2f}]")
print(f"  Σ_B 范围: [{sigB_w.min():.2f}, {sigB_w.max():.2f}]")
print(f"  Σ_A vs Σ_B 相关: ρ = {spearmanr(sigA_w, sigB_w)[0]:+.4f}")

print("\n--- SB_min vs Σ ---")
print_corr("SB_min vs Σ_A", sb_min_w, sigA_w)
print_corr("SB_min vs Σ_B", sb_min_w, sigB_w)

# ====================================================================
# 5. MERGED 291 (Kokorev + Path1)
# ====================================================================
print("\n" + "=" * 72)
print("5. MERGED SAMPLE (Kokorev + Path1)")
print("=" * 72)

# Merge
sigA_merged = np.concatenate([sigA_k, sigA_p1])
sigB_merged = np.concatenate([sigB_k, sigB_p1])
sb_merged = np.concatenate([sb_k, sb_p1])
z_merged = np.concatenate([z_k, z_p1])
logM_merged = np.concatenate([logM_k, logM_p1])

valid_m = (np.isfinite(sigA_merged) & np.isfinite(sigB_merged) &
           np.isfinite(sb_merged))
sigA_merged, sigB_merged, sb_merged, z_merged, logM_merged = [
    a[valid_m] for a in [sigA_merged, sigB_merged, sb_merged, z_merged, logM_merged]]

print(f"\n  合并样本: N = {len(sigA_merged)}")
print(f"  Σ_A vs Σ_B 相关: ρ = {spearmanr(sigA_merged, sigB_merged)[0]:+.4f}")

print("\n--- SB_F444W vs Σ (合并) ---")
print_corr("SB vs Σ_A", sb_merged, sigA_merged)
print_corr("SB vs Σ_B", sb_merged, sigB_merged)
r_sb_parA = psc(sb_merged, sigA_merged, [z_merged])
r_sb_parB = psc(sb_merged, sigB_merged, [z_merged])
print(f"  SB vs Σ_A (偏相关,控制z): ρ={r_sb_parA:+.4f}")
print(f"  SB vs Σ_B (偏相关,控制z): ρ={r_sb_parB:+.4f}")

# ====================================================================
# 6. COSMOSWEB 664k — check Σ definition
# ====================================================================
print("\n" + "=" * 72)
print("6. COSMOSWEB 664k — Σ 定义溯源")
print("=" * 72)

cosmos = np.load(os.path.join(cosmos_dir, 'cosmos2025_extracted.npz'))
z_c = cosmos['z']
logM_c = cosmos['logM']
logSigma_c = cosmos['logSigma']
color_c = cosmos['color']
reff_c = cosmos['reff_pc']

valid_c = (np.isfinite(z_c) & np.isfinite(logM_c) & np.isfinite(logSigma_c) &
           np.isfinite(color_c) & np.isfinite(reff_c) & (reff_c > 0) &
           (logM_c > 6) & (logM_c < 14) & (z_c > 0) & (z_c < 15))
z_c, logM_c, logSigma_c, color_c, reff_c = [
    a[valid_c] for a in [z_c, logM_c, logSigma_c, color_c, reff_c]]

# Compute Σ_B for COSMOS
sigB_c = sigma_B_from_mass(logM_c, reff_c)

print(f"\n  有效样本: N = {len(z_c)}")
print(f"  已有 logSigma 范围: [{logSigma_c.min():.2f}, {logSigma_c.max():.2f}]")
print(f"  物理 Σ_B 范围: [{sigB_c.min():.2f}, {sigB_c.max():.2f}]")
print(f"  已有 logSigma vs Σ_B 相关: ρ = {spearmanr(logSigma_c, sigB_c)[0]:+.4f}")

# Check if the existing logSigma appears to be F444W-based or M*-based
# If Σ comes from the COSMOS pipeline, it might use Sersic profile fitting
# Compare magnitude of difference
delta = logSigma_c - sigB_c
print(f"  已有logΣ - 物理Σ_B 差值: median={np.median(delta):.2f}, std={np.std(delta):.2f}")

# Test correlation with color using both
z_bins = [(0,1),(1,3),(3,5),(5,7),(7,9),(9,15)]

print("\n--- Σ-color 偏相关 (控制 z, M*) 按红移 ---")
print(f"  {'z bin':<10s} {'N':>7s} {'ρ(existing Σ)':>14s} {'ρ(物理Σ_B)':>14s} {'Δρ':>10s}")
print("  " + "-" * 58)

for (zlo, zhi) in z_bins:
    mask = (z_c >= zlo) & (z_c < zhi)
    n = mask.sum()
    if n < 100:
        continue
    rho_existing = psc(color_c[mask], logSigma_c[mask], [z_c[mask], logM_c[mask]])
    rho_physical = psc(color_c[mask], sigB_c[mask], [z_c[mask], logM_c[mask]])
    dz = rho_existing - rho_physical
    print(f"  z={zlo}-{zhi:<4d} {n:>7d} {rho_existing:>+14.4f} {rho_physical:>+14.4f} {dz:>+10.4f}")

# ====================================================================
# 7. SUMMARY
# ====================================================================
print("\n" + "=" * 72)
print("7. 总结")
print("=" * 72)

print("""
┌─────────────────────────────────────────────────────────────────────┐
│ 核心发现                                                           │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  Σ_A (F444W通量代理) ≠ Σ_B (物理面密度 M*/area)                     │
│                                                                     │
│  两者的相关程度因样本而异。如果 Σ_A→Σ_B 后大部分证据链               │
│  的显著性大幅下降，说明信号的物理来源不是恒星密度本身，              │
│  而是F444W通量中M*之外的额外信息（AGN贡献、年龄差异等）。           │
│                                                                     │
│  这将影响物理模型的解释逻辑。                                       │
└─────────────────────────────────────────────────────────────────────┘
""")

# Save results
results = {
    'path1': {
        'N': int(len(sigA_p1)),
        'rho_SigA_SigB': float(spearmanr(sigA_p1, sigB_p1)[0]),
        'SigA_vs_BD': rA_bd,
        'SigB_vs_BD': rB_bd,
        'SigA_par_BD': float(rA_par),
        'SigB_par_BD': float(rB_par),
    },
    'kokorev': {
        'N': int(len(sigA_k)),
        'rho_SigA_SigB': float(spearmanr(sigA_k, sigB_k)[0]),
    },
    'merged': {
        'N': int(len(sigA_merged)),
        'rho_SigA_SigB': float(spearmanr(sigA_merged, sigB_merged)[0]),
    },
}

with open(os.path.join(data_dir, 'sigma_definition_diagnosis.json'), 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\n详细结果已保存至: {os.path.join(data_dir, 'sigma_definition_diagnosis.json')}")
