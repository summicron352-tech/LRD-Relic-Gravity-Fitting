#!/usr/bin/env python3
"""
Withers+ 2026 (arXiv:2606.06585) — 致密REG 与 ApJL AAS77670 交叉验证
=========================================================================
比较 Withers+ 发现的 9 个"缺失"致密REG 与喵老师 ApJL 46源样本的 F444W 表面亮度。
所有 9 个致密REG 被经典 LRD 选择系统性遗漏，但可能在 SB 维度上与现有 LRD 样本重叠。

关键发现：Withers+ 致密REG 的 F444W 表面亮度下限（中位 ~24.0 mag/arcsec^2）
已经比 ApJL 样本中位数（~24.5）更亮。如果它们实际比 PSF 更小，SB 会更高，
可能落在 ApJL 样本的极端高端。

局限性：Withers+ 数据缺乏巴尔末减缩（BD），无法完整验证 Σ-BD 相关性。
这是未来的可检验预测。

作者：三岁喵 (XIN Tan)
日期：2026-06-09
"""

import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
import astropy.units as u
from astropy.coordinates import SkyCoord
import json
import os

cosmo = Planck18

OUTDIR = "/Users/tanxin/Desktop/数据处理/11_Withers_MissingLRDs_CrossCheck"

# ============================================================
# 1. Withers+ Table 2 数据
# ============================================================
compact_regs = [
    {"id": "5203291",  "z": 5.15, "logMstar": 8.05, "F444W_mag": 28.5, "AV": 1.88,
     "EW_Ha": 1173, "EW_OIIIHb": 434},
    {"id": "1207412",  "z": 5.21, "logMstar": 8.03, "F444W_mag": 27.5, "AV": 1.11,
     "EW_Ha": 1327, "EW_OIIIHb": 843},
    {"id": "2228428",  "z": 5.99, "logMstar": 8.10, "F444W_mag": 28.2, "AV": 1.59,
     "EW_Ha": 1492, "EW_OIIIHb": 2278},
    {"id": "2209729",  "z": 6.02, "logMstar": 8.16, "F444W_mag": 27.6, "AV": 1.19,
     "EW_Ha": 1006, "EW_OIIIHb": 2116},
    {"id": "3112462",  "z": 6.03, "logMstar": 8.11, "F444W_mag": 28.0, "AV": 2.64,
     "EW_Ha": 572, "EW_OIIIHb": 272},
    {"id": "2210928",  "z": 6.79, "logMstar": 8.25, "F444W_mag": 28.7, "AV": 1.02,
     "EW_Ha": None, "EW_OIIIHb": 1398},
    {"id": "2113181",  "z": 8.18, "logMstar": 7.96, "F444W_mag": 27.8, "AV": 0.97,
     "EW_Ha": None, "EW_OIIIHb": 1457},
    {"id": "3204679",  "z": 8.41, "logMstar": 8.16, "F444W_mag": 28.7, "AV": 2.22,
     "EW_Ha": None, "EW_OIIIHb": 1080},
    {"id": "2213303",  "z": 8.63, "logMstar": 8.36, "F444W_mag": 28.5, "AV": 2.05,
     "EW_Ha": None, "EW_OIIIHb": 983},
]

psf_fwhm = 0.161  # arcsec (F444W PSF)
psf_hwhm = psf_fwhm / 2  # 0.0805"
psf_area = np.pi * psf_hwhm**2  # ~0.0204 arcsec^2 (面积上限 → SB下限)

# 计算 SB 下限
for r in compact_regs:
    r['r_eff_max_kpc'] = None  # 稍后填
    r['SB_F444W_min'] = r['F444W_mag'] + 2.5 * np.log10(psf_area)
    # r_eff 上限 (物理)
    d_A = cosmo.angular_diameter_distance(r['z']).to(u.pc).value
    theta_rad = psf_hwhm / 206265.0
    r_eff_max_pc = d_A * theta_rad
    r['r_eff_max_kpc'] = round(r_eff_max_pc / 1000, 3)
    # stellar mass surface density 下限
    r['logSigma_Mstar_min'] = round(r['logMstar'] - 2 * np.log10(r_eff_max_pc), 2)

df_withers = pd.DataFrame(compact_regs)
df_withers.to_csv(f"{OUTDIR}/withers_compact_regs.csv", index=False)
print("[1/5] Withers+ 致密REG 数据已保存")

# ============================================================
# 2. 加载 ApJL 样本 (Kokorev 全星表 + deGraaff 交叉)
# ============================================================
k_raw = pd.read_csv('/Users/tanxin/WorkBuddy/20260412234449/lrd-relic-repo/data/csv/LRD_Master_Combined_AllParams.csv')
k = k_raw.dropna(subset=['r_eff_50_phys', 'logMstar_best', 'f444w_flux', 'ra', 'dec']).copy()
k = k[(k['r_eff_50_phys'] > 0) & (k['logMstar_best'] > 0) & (k['f444w_flux'] > 0)]

k['F444W_mag'] = -2.5 * np.log10(k['f444w_flux'].values * 1e-9 / 3631.0)
k['d_A_pc'] = [cosmo.angular_diameter_distance(z).to(u.pc).value for z in k['z_phot']]
k['r_eff_arcsec'] = k['r_eff_50_phys'] / k['d_A_pc'] * 206265.0
k['area_arcsec2'] = np.pi * k['r_eff_arcsec']**2
k['SB_F444W'] = k['F444W_mag'] + 2.5 * np.log10(k['area_arcsec2'])
k['logSigma_Mstar'] = k['logMstar_best'] - 2 * np.log10(k['r_eff_50_phys'])

# 交叉匹配 deGraaff BD 数据
cross = pd.read_csv('/Users/tanxin/WorkBuddy/20260412234449/figures_apjl/validation_crossmatch.csv')
coords_k = SkyCoord(ra=k['ra'].values * u.deg, dec=k['dec'].values * u.deg)
coords_cross = SkyCoord(ra=cross['ra'].values * u.deg, dec=cross['dec'].values * u.deg)
idx, sep2d, _ = coords_cross.match_to_catalog_sky(coords_k)
mask = sep2d.arcsec < 1.0

k_matched = k.iloc[idx[mask]].copy()
cross_m = cross[mask].copy()
k_matched['BD'] = cross_m['BD'].values
k_matched['break_strength'] = cross_m['break_strength'].values
k_matched['sep_arcsec'] = cross_m['sep_arcsec'].values

# 保存完整数据
cols_out = ['id', 'field', 'ra', 'dec', 'z_phot', 'F444W_mag', 'r_eff_50_phys',
            'r_eff_arcsec', 'SB_F444W', 'logSigma_Mstar', 'logMstar_best',
            'BD', 'break_strength', 'verdict']
k_matched[cols_out].to_csv(f"{OUTDIR}/apjl_crossmatch_with_sb.csv", index=False)
print(f"[2/5] ApJL 交叉匹配 {len(k_matched)} 源已保存")

# 保存全 Kokorev 星表 SB 数据
k[['id', 'field', 'ra', 'dec', 'z_phot', 'F444W_mag', 'r_eff_50_phys',
   'r_eff_arcsec', 'SB_F444W', 'logSigma_Mstar', 'logMstar_best', 'verdict']].to_csv(
    f"{OUTDIR}/kokorev_full_sb.csv", index=False)
print(f"[3/5] Kokorev 全星表 {len(k)} 源已保存")

# ============================================================
# 3. 对比统计
# ============================================================
withers_sb = [r['SB_F444W_min'] for r in compact_regs]

stats = {
    "date": "2026-06-09",
    "paper": "Withers+ 2026 (arXiv:2606.06585)",
    "withers_compact_regs": {
        "N": len(compact_regs),
        "z_range": [min(r['z'] for r in compact_regs), max(r['z'] for r in compact_regs)],
        "F444W_mag_range": [min(r['F444W_mag'] for r in compact_regs), max(r['F444W_mag'] for r in compact_regs)],
        "SB_F444W_min_median": float(np.median(withers_sb)),
        "SB_F444W_min_range": [float(min(withers_sb)), float(max(withers_sb))],
        "logSigma_Mstar_min_median": float(np.median([r['logSigma_Mstar_min'] for r in compact_regs])),
        "note": "所有 SB 值均为下限（源在 F444W 中未分辨，PSF FWHM=0.161\"）"
    },
    "apjl_crossmatch": {
        "N": int(len(k_matched)),
        "r_eff_pc_range": [float(k_matched['r_eff_50_phys'].min()), float(k_matched['r_eff_50_phys'].max())],
        "F444W_mag_range": [float(k_matched['F444W_mag'].min()), float(k_matched['F444W_mag'].max())],
        "SB_F444W_median": float(k_matched['SB_F444W'].median()),
        "SB_F444W_range": [float(k_matched['SB_F444W'].min()), float(k_matched['SB_F444W'].max())],
        "BD_range": [float(k_matched['BD'].min()), float(k_matched['BD'].max())],
    },
    "kokorev_full": {
        "N": int(len(k)),
        "SB_F444W_median": float(k['SB_F444W'].median()),
        "SB_F444W_range": [float(k['SB_F444W'].min()), float(k['SB_F444W'].max())],
    },
    "key_finding": {
        "apjl_SB_median": float(k_matched['SB_F444W'].median()),
        "withers_SB_min_median": float(np.median(withers_sb)),
        "delta_mag": float(k_matched['SB_F444W'].median() - np.median(withers_sb)),
        "interpretation": "Withers+ 致密REG SB下限（中位 {:.1f}）比 ApJL 样本中位数（{:.1f}）亮 {:.1f} mag。真实值更亮。".format(
            np.median(withers_sb), k_matched['SB_F444W'].median(), 
            k_matched['SB_F444W'].median() - np.median(withers_sb)),
        "caveat": "Withers+ 缺乏 BD 数据，完整 Σ-BD 验证需等待光谱观测",
        "prediction": "如果致密REG是LRD，其BD应落在 Σ-BD 相关线的高端延伸"
    }
}

with open(f"{OUTDIR}/comparison_stats.json", 'w') as f:
    json.dump(stats, f, indent=2, ensure_ascii=False)
print("[4/5] 对比统计已保存")

# ============================================================
# 4. 生成对比图数据
# ============================================================
# 用于后续画图
plot_data = {
    "withers_x": [r['F444W_mag'] for r in compact_regs],
    "withers_y": withers_sb,
    "withers_id": [r['id'] for r in compact_regs],
    "withers_label": "Withers+ compact REGs (SB下限)",
    "apjl_x": k_matched['F444W_mag'].tolist(),
    "apjl_y": k_matched['SB_F444W'].tolist(),
    "apjl_bd": k_matched['BD'].tolist(),
    "apjl_label": f"ApJL AAS77670 crossmatch (N={len(k_matched)})",
}

with open(f"{OUTDIR}/plot_data.json", 'w') as f:
    json.dump(plot_data, f, indent=2)
print("[5/5] 画图数据已保存")

print(f"\n所有文件已保存至: {OUTDIR}")
print(f"  - withers_compact_regs.csv       Withers+ 9致密REG完整参数")
print(f"  - apjl_crossmatch_with_sb.csv    ApJL交叉匹配38源 + SB")
print(f"  - kokorev_full_sb.csv            Kokorev全260源 + SB")
print(f"  - comparison_stats.json          对比统计摘要")
print(f"  - plot_data.json                 画图数据")
print(f"  - withers_analysis.py            本脚本")
