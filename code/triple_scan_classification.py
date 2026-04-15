"""
LRD G_eff 检验 — 三组分类扫描（Gold / Silver / Noise）
=========================================================
核心思想：LRD 的定义性特征是 X-ray 静默 + 极紧凑 + 高红移。
用物理属性将 260 个源分成三个质量等级，分别检验 G_eff 信号。

分组标准（基于 Kokorev v1.1 星表可获取的物理量）：

  🥇 GOLD:   z_phot ≥ 7    （真正的高红移 LRD 候选体）
              r_eff < 中位数（极紧凑，点源特征）
              f444w/f150w < P75（非 AGN 型红外 excess）
              n_bands ≥ 6      （数据质量好）
              → 这些是最"纯"的 LRD 样本，G_eff 信号应该最清晰

  🥈 SILVER: z_phot 4 ~ 7   （中高红移，辅助样本）
              SED 形态合理
              不满足 Gold 但也不像 Noise
              → 辅助证据组

  🥉 NOISE:  z_phot < 4     （低红移，可能是污染星系）
              或 r_eff > P90  （弥散形态，不是点源）
              或 f444w/f150w > P95（极端红外 excess = AGN 幂律？）
              或 Null χ²_ν > 100（两个模型都拟合很差）
              → 控制组：预期不应该看到 G_eff 信号

输出：
  1. Triple_Scan_Summary.png     — 三组对比大图（6~8 面板）
  2. Triple_Scan_Stats.txt       — 数值统计汇总
  3. Triple_Scan_Gold_SEDs.pdf   — Gold 组全部 SED 对比图
  4. Triple_Scan_Results.csv     — 全部结果含分组标签

Dependencies: pandas, numpy, matplotlib, scipy
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats as sp_stats
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')


# ==========================================
# 0. G_eff Window Model — 使用扫描最佳参数
# ==========================================

def H_z(z, Omega_m=0.3, Omega_L=0.7, H0=70.0):
    return H0 * np.sqrt(Omega_m * (1+z)**3 + Omega_L)

def G_eff_window(z, A_G=0.15, beta=0.5, z_p=10.0, z_on=8.0, z_off=30.0, dz=2.0):
    W = 0.25 * (1 + np.tanh((z - z_on) / dz)) * (1 - np.tanh((z - z_off) / dz))
    Hz_ratio = H_z(z) / H_z(z_p)
    return A_G * (Hz_ratio ** beta) * W

def z_dist_from_phot(z_phot, gamma=1.0, **window_params):
    delta_G = G_eff_window(z_phot, **window_params)
    return gamma * np.clip(delta_G, 0, None)


# ==========================================
# 1. 唯象模型定义
# ==========================================

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


# ==========================================
# 2. 双模型拟合引擎
# ==========================================

def fit_two_param(rest_wave, nuFnu_obs, nuFnu_err, z_dist_fixed):
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
        
        return {
            'C_blue': popt[0], 'C_red': popt[1],
            'perr_Cb': np.sqrt(pcov[0, 0]),
            'perr_Cr': np.sqrt(pcov[1, 1]),
            'chi2': chi2, 'red_chi2': red_chi2,
            'bic': bic, 'dof': dof, 'n_data': n_data,
            'model_flux': model_flux,
            'residuals': residuals,
            'success': True
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}


# ==========================================
# 3. 数据加载与预处理
# ==========================================

print("=" * 70)
print("TRIPLE SCAN: Gold / Silver / Noise Classification")
print("=" * 70)

df = pd.read_csv('Kokorev_LRDs_Full.csv')
print(f"Total sources in catalog: {len(df)}")

filters = {
    'f090w': 0.90, 'f115w': 1.15, 'f150w': 1.50, 'f200w': 2.00,
    'f277w': 2.77, 'f356w': 3.56, 'f444w': 4.44, 'f770w': 7.70
}

# 加载全部源 + 计算辅助量
all_sources = []
for idx, row in df.iterrows():
    target_id = row['id']
    z_phot = row['z_phot']
    field = row.get('field', 'unknown')
    r_eff = row['r_eff_50_phys'] if pd.notnull(row.get('r_eff_50_phys')) else np.nan
    lbol = row['lbol'] if pd.notnull(row.get('lbol')) else np.nan
    muv = row['muv'] if pd.notnull(row.get('muv')) else np.nan
    
    obs_wave, obs_flux, obs_err = [], [], []
    for f, wave in filters.items():
        flux_col = f'{f}_flux'
        err_col = f'{f}_fluxerr'
        if flux_col in row and pd.notnull(row[flux_col]) and row[flux_col] > -90:
            obs_wave.append(wave)
            obs_flux.append(row[flux_col])
            err_val = row[err_col] if err_col in row and row[err_col] > 0 else row[flux_col] * 0.1
            obs_err.append(err_val)
    
    if len(obs_wave) >= 4:
        obs_wave = np.array(obs_wave)
        obs_flux = np.array(obs_flux)
        obs_err = np.array(obs_err)
        
        # 计算 IR excess ratio
        f150_idx = None; f444_idx = None
        for i, w in enumerate(obs_wave):
            if abs(w - 1.50) < 0.01: f150_idx = i
            if abs(w - 4.44) < 0.01: f444_idx = i
        
        ir_excess = np.nan
        if f150_idx is not None and f444_idx is not None:
            if obs_flux[f150_idx] > 0:
                ir_excess = obs_flux[f444_idx] / obs_flux[f150_idx]
        
        rest_wave = obs_wave / (1 + z_phot)
        nuFnu_obs = obs_flux / obs_wave
        nuFnu_err = obs_err / obs_wave
        norm_factor = np.max(nuFnu_obs) if np.max(nuFnu_obs) > 0 else 1.0
        
        all_sources.append({
            'target_id': target_id,
            'z_phot': z_phot,
            'field': field,
            'r_eff': r_eff,
            'lbol': lbol,
            'muv': muv,
            'ir_excess': ir_excess,
            'rest_wave': rest_wave,
            'nuFnu_norm': nuFnu_obs / norm_factor,
            'nuFnu_err_norm': nuFnu_err / norm_factor,
            'norm_factor': norm_factor,
            'obs_wave': obs_wave,
            'n_bands': len(obs_wave),
        })

print(f"Sources with ≥4 valid bands: {len(all_sources)}")

# ==========================================
# 4. 计算全局阈值用于分类
# ==========================================

valid_r_eff = [s['r_eff'] for s in all_sources if not np.isnan(s['r_eff'])]
valid_ir = [s['ir_excess'] for s in all_sources if not np.isnan(s['ir_excess'])]

R_EFF_MEDIAN = np.median(valid_r_eff) if valid_r_eff else 100
R_EFF_P90 = np.percentile(valid_r_eff, 90) if valid_r_eff else 500
IR_P25 = np.percentile(valid_ir, 25) if valid_ir else 3
IR_P75 = np.percentile(valid_ir, 75) if valid_ir else 15
IR_P95 = np.percentile(valid_ir, 95) if valid_ir else 30

print(f"\n--- Classification Thresholds ---")
print(f"  r_eff: median={R_EFF_MEDIAN:.1f} pc, P90={R_EFF_P90:.1f} pc")
print(f"  IR excess (f444w/f150w):")
print(f"    P25={IR_P25:.2f}, P75={IR_P75:.2f}, P95={IR_P95:.2f}")


# ==========================================
# 5. 三组分类逻辑
# ==========================================

def classify_source(src):
    """
    分类优先级（从严格到宽松）：
    
    GOLD 条件（必须同时满足所有）：
      ① z_phot ≥ 7（高红移，真正的 EoR LRD）
      ② r_eff < median（极紧凑）
      ③ IR excess < P75（非 AGN 型红外）
      ④ n_bands ≥ 6（数据质量足够）
    
    NOISE 条件（满足任一即排除）：
      ① z_phot < 4（低红移，大概率是污染星系/低z星暴）
      ② r_eff > P90 或 NaN（弥散形态）
      ③ IR excess > P95 或 > 30（极端红外 = 可能是 AGN/PAH 主导）
    
    其余归入 SILVER
    """
    z = src['z_phot']
    reff = src['r_eff']
    ir = src['ir_excess']
    nb = src['n_bands']
    
    # === GOLD 检查 ===
    gold_checks = {
        'z≥7': z >= 7,
        'compact': (not np.isnan(reff)) and reff < R_EFF_MEDIAN,
        'moderate_IR': (not np.isnan(ir)) and ir < IR_P75,
        'good_data': nb >= 6,
    }
    is_gold = all(gold_checks.values())
    
    # === NOISE 检查 ===
    is_noise = False
    noise_reason = []
    
    if z < 4:
        is_noise = True
        noise_reason.append('low-z')
    if (not np.isnan(reff) and reff > R_EFF_P90) or (np.isnan(reff)):
        # NaN 的 r_eff 也归入 Noise（数据缺失不可信）
        if np.isnan(reff):
            is_noise = True
            noise_reason.append('no_size')
        elif reff > R_EFF_P90:
            is_noise = True
            noise_reason.append('extended')
    if (not np.isnan(ir) and ir > IR_P95) or (not np.isnan(ir) and ir > 30):
        is_noise = True
        noise_reason.append('extreme_IR')
    
    if is_gold:
        return 'Gold', gold_checks, ''
    elif is_noise:
        return 'Noise', {}, '|'.join(noise_reason)
    else:
        return 'Silver', {}, ''


# 执行分类
for src in all_sources:
    group, checks, reason = classify_source(src)
    src['group'] = group
    src['gold_checks'] = checks
    src['noise_reason'] = reason


gold_srcs = [s for s in all_sources if s['group'] == 'Gold']
silver_srcs = [s for s in all_sources if s['group'] == 'Silver']
noise_srcs = [s for s in all_sources if s['group'] == 'Noise']

print(f"\n=== CLASSIFICATION RESULTS ===")
print(f"  🥇 GOLD:   {len(gold_srcs):>3d} sources ({100*len(gold_srcs)/len(all_sources):.1f}%)")
print(f"  🥈 SILVER: {len(silver_srcs):>3d} sources ({100*len(silver_srcs)/len(all_sources):.1f}%)")
print(f"  🥉 NOISE:  {len(noise_srcs):>3d} sources ({100*len(noise_srcs)/len(all_sources):.1f}%)")

# 打印 Gold 列表
if gold_srcs:
    print(f"\n  --- Gold Sample Details ---")
    for s in sorted(gold_srcs, key=lambda x: x['z_phot'], reverse=True):
        gc = s['gold_checks']
        print(f"    ID={s['target_id']:>5d}  z={s['z_phot']:>6.2f}  "
              f"r_eff={s['r_eff']:>6.1f}  IR={s['ir_excess']:>5.2f}  "
              f"N={s['n_bands']}  field={s['field']:<20s}")


# ==========================================
# 6. G_eff 窗口参数（扫描最佳值）
# ==========================================

window_optimal = {
    'A_G': 0.189,
    'beta': 1.58,
    'z_p': 7.0,
    'z_on': 2.0,
    'z_off': 25.0,
    'dz': 3.0
}
GAMMA_OPTIMAL = 0.50


# ==========================================
# 7. 单组运行引擎
# ==========================================

def run_group_test(group_name, sources_list):
    """对一组源运行 G_eff vs Null 双模型检验"""
    
    # 预计算 z_dist
    for src in sources_list:
        src['z_dist_pred'] = z_dist_from_phot(
            src['z_phot'], gamma=GAMMA_OPTIMAL, **window_optimal
        )
        src['delta_G_pred'] = G_eff_window(src['z_phot'], **window_optimal)
    
    results = []
    print(f"\n{'─'*75}")
    print(f"GROUP: {group_name.upper()}  (N={len(sources_list)})")
    print(f"{'─'*75}")
    print(f"{'ID':>7s} {'z':>5s} {'N':>2s} | {'χ²ν(N)':>8s} {'χ²ν(G)':>8s} "
          f"{'Δχ²ν':>8s} | {'ΔBIC':>7s} | {'Verdict':>10s}")
    print(f"{'─'*75}")
    
    for src in sources_list:
        z_phot = src['z_phot']
        target_id = src['target_id']
        
        # Null model (z_dist = 0)
        fit_null = fit_two_param(
            src['rest_wave'], src['nuFnu_norm'], src['nuFnu_err_norm'],
            z_dist_fixed=0.0
        )
        
        # G_eff model (z_dist from theory)
        z_dist_val = src['z_dist_pred']
        fit_geff = fit_two_param(
            src['rest_wave'], src['nuFnu_norm'], src['nuFnu_err_norm'],
            z_dist_fixed=z_dist_val
        )
        
        if fit_null['success'] and fit_geff['success']:
            delta_chi2_nu = fit_null['red_chi2'] - fit_geff['red_chi2']
            delta_bic = fit_null['bic'] - fit_geff['bic']
            
            if delta_chi2_nu > 1 and delta_bic < -2:
                verdict = "Support"
            elif delta_chi2_nu > 0 and delta_bic >= -2:
                verdict = "Weak"
            elif abs(delta_chi2_nu) <= 1:
                verdict = "Neutral"
            else:
                verdict = "Null wins"
            
            results.append({
                'id': target_id,
                'field': src['field'],
                'z_phot': z_phot,
                'z_dist': z_dist_val,
                'delta_G': src['delta_G_pred'],
                'n_bands': src['n_bands'],
                'ir_excess': src['ir_excess'],
                'r_eff': src['r_eff'],
                'group': group_name,
                'Cb_null': fit_null['C_blue'],
                'Cr_null': fit_null['C_red'],
                'chi2_null': fit_null['chi2'],
                'rch2_null': fit_null['red_chi2'],
                'bic_null': fit_null['bic'],
                'model_null': fit_null['model_flux'],
                'Cb_geff': fit_geff['C_blue'],
                'Cr_geff': fit_geff['C_red'],
                'chi2_geff': fit_geff['chi2'],
                'rch2_geff': fit_geff['red_chi2'],
                'bic_geff': fit_geff['bic'],
                'model_geff': fit_geff['model_flux'],
                'residuals_geff': fit_geff['residuals'],
                'delta_chi2_nu': delta_chi2_nu,
                'delta_bic': delta_bic,
                'verdict': verdict,
                'rest_wave': src['rest_wave'],
                'nuFnu_norm': src['nuFnu_norm'],
                'nuFnu_err_norm': src['nuFnu_err_norm'],
            })
            
            marker = " ★" if verdict in ('Support', 'Weak') else ""
            print(f"{target_id:>7d} {z_phot:>5.2f} {src['n_bands']:>2d} | "
                  f"{fit_null['red_chi2']:>8.2f} {fit_geff['red_chi2']:>8.2f} "
                  f"{delta_chi2_nu:>+8.2f} | {delta_bic:>+7.1f} | "
                  f"{verdict:>10s}{marker}")
        else:
            results.append({
                'id': target_id, 'field': src['field'], 'z_phot': z_phot,
                'group': group_name, 'verdict': 'FAILED',
                'delta_chi2_nu': np.nan, 'delta_bic': np.nan,
                'z_dist': z_dist_val, 'n_bands': src['n_bands'],
            })
            print(f"{target_id:>7d} {z_phot:>5.2f} {src['n_bands']:>2d} | *** FIT FAILED ***")
    
    return results


# 分别跑三组
print("\n" + "=" * 70)
print("RUNNING TRIPLE-SCAN: Gold → Silver → Noise")
print("=" * 70)

results_gold = run_group_test('Gold', gold_srcs)
results_silver = run_group_test('Silver', silver_srcs)
results_noise = run_group_test('Noise', noise_srcs)

# 合并全部结果
all_results = results_gold + results_silver + results_noise


# ==========================================
# 8. 分组统计分析
# ==========================================

def analyze_group(group_name, results_list):
    """计算单组的统计摘要"""
    valid = [r for r in results_list if r['verdict'] != 'FAILED']
    if not valid:
        return None
    
    dchi = [r['delta_chi2_nu'] for r in valid]
    dbic = [r['delta_bic'] for r in valid]
    zs = [r['z_phot'] for r in valid]
    verdicts = [r['verdict'] for r in valid]
    
    n_support = sum(1 for v in verdicts if v == 'Support')
    n_weak = sum(1 for v in verdicts if v == 'Weak')
    n_neutral = sum(1 for v in verdicts if v == 'Neutral')
    n_nullwin = sum(1 for v in verdicts if v == 'Null wins')
    n_total_favor = n_support + n_weak
    
    # Sign test: Mean Δχ²_ν 是否显著不为零？
    mean_dchi = np.mean(dchi)
    median_dchi = np.median(dchi)
    std_dchi = np.std(dchi)
    sem_dchi = std_dchi / np.sqrt(len(dchi))
    
    # t-test against zero (one-sided: is mean > 0?)
    if len(dchi) > 1:
        t_stat, t_pval = sp_stats.ttest_1samp(dchi, 0)
        # one-tailed p-value
        one_tailed_p = t_pval / 2 if mean_dchi > 0 else 1 - t_pval/2
    else:
        t_stat, one_tailed_p = np.nan, np.nan
    
    # Wilcoxon signed-rank (non-parametric)
    if len(dchi) > 5:
        try:
            wilc_stat, wilc_pval = sp_stats.wilcoxon(dchi)
        except:
            wilc_stat, wilc_pval = np.nan, np.nan
    else:
        wilc_stat, wilc_pval = np.nan, np.nan
    
    return {
        'name': group_name,
        'N': len(valid),
        'support': n_support,
        'weak': n_weak,
        'neutral': n_neutral,
        'nullwin': n_nullwin,
        'favor_pct': 100 * n_total_favor / len(valid),
        'mean_dchi': mean_dchi,
        'median_dchi': median_dchi,
        'std_dchi': std_dchi,
        'sem_dchi': sem_dchi,
        't_stat': t_stat,
        't_pvalue_one_sided': one_tailed_p,
        'wilc_stat': wilc_stat,
        'wilc_pval': wilc_pval,
        'median_dbic': np.median(dbic),
        'z_range': (min(zs), max(zs)),
        'dchi_values': dchi,
        'dbic_values': dbic,
        'zs': zs,
        'verdicts': verdicts,
        'valid_results': valid,
    }


stats_gold = analyze_group('🥇 Gold', results_gold)
stats_silver = analyze_group('🥈 Silver', results_silver)
stats_noise = analyze_group('🥉 Noise', results_noise)

all_stats = [s for s in [stats_gold, stats_silver, stats_noise] if s is not None]

# 打印统计表
print(f"\n{'='*80}")
print(f"TRIPLE-SCAN STATISTICAL COMPARISON")
print(f"{'='*80}")

header = (f"{'Group':^12s} │ {'N':>4s} │ {'Sup':>4s} {'Wk':>4s} {'Neu':>4s} {'Null':>4s} │ "
          f"{'Favor%':>6s} │ {'Mean Δχ²ν':>10s} │ {'Med Δχ²ν':>10s} │ {'t-test p':>10s} │ {'Wilcox p':>10s}")
print(header)
print('─' * 120)

for st in all_stats:
    sig_marker = "***" if st['t_pvalue_one_sided'] < 0.001 else "**" if st['t_pvalue_one_sided'] < 0.01 else "*" if st['t_pvalue_one_sided'] < 0.05 else ""
    print(f"{st['name']:^12s} │ {st['N']:>4d} │ {st['support']:>4d} {st['weak']:>4d} "
          f"{st['neutral']:>4d} {st['nullwin']:>4d} │ {st['favor_pct']:>5.1f}% │ "
          f"{st['mean_dchi']:>+9.2f} │ {st['median_dchi']:>+9.2f} │ "
          f"{st['t_pvalue_one_sided']:>9.4f}{sig_marker:>3s} │ "
          f"{st['wilc_pval']:>9.4f}")

print('─' * 120)

# 关键检验：Gold vs Noise 的差异是否显著？
if stats_gold and stats_noise:
    gold_dchi = np.array(stats_gold['dchi_values'])
    noise_dchi = np.array(stats_noise['dchi_values'])
    
    u_stat, u_pval = sp_stats.mannwhitneyu(gold_dchi, noise_dchi, alternative='greater')
    print(f"\n  🔬 MANN-WHITNEY U TEST (Gold > Noise?):  U={u_stat:.1f}, p={u_pval:.4f}")
    if u_pval < 0.05:
        print(f"  ✅ Gold 样本的 Δχ²_ν 显著高于 Noise！→ G_eff 信号具有物理选择性")
    else:
        print(f"  ⚠️ Gold 与 Noise 差异不显著")


# 保存数值统计
with open('Triple_Scan_Stats.txt', 'w') as f:
    f.write("="*70 + "\n")
    f.write("TRIPLE SCAN: Gold/Silver/Noise G_eff Test Results\n")
    f.write("="*70 + "\n\n")
    f.write(f"Parameters: A_G={window_optimal['A_G']}, γ={GAMMA_OPTIMAL}, "
            f"β={window_optimal['beta']}, z_on={window_optimal['z_on']}\n\n")
    
    f.write(f"Classification thresholds:\n")
    f.write(f"  r_eff median = {R_EFF_MEDIAN:.1f} pc\n")
    f.write(f"  r_eff P90    = {R_EFF_P90:.1f} pc\n")
    f.write(f"  IR excess P25 = {IR_P25:.2f}\n")
    f.write(f"  IR excess P75 = {IR_P75:.2f}\n")
    f.write(f"  IR excess P95 = {IR_P95:.2f}\n\n")
    
    for st in all_stats:
        f.write("-"*60 + "\n")
        f.write(f"GROUP: {st['name']}\n")
        f.write(f"-"*60 + "\n")
        f.write(f"  N total           = {st['N']}\n")
        f.write(f"  Support           = {st['support']} ({100*st['support']/st['N']:.1f}%)\n")
        f.write(f"  Weak support      = {st['weak']} ({100*st['weak']/st['N']:.1f}%)\n")
        f.write(f"  Neutral           = {st['neutral']} ({100*st['neutral']/st['N']:.1f}%)\n")
        f.write(f"  Null wins         = {st['nullwin']} ({100*st['nullwin']/st['N']:.1f}%)\n")
        f.write(f"  Total favor %     = {st['favor_pct']:.1f}%\n")
        f.write(f"  Mean Δχ²_ν        = {st['mean_dchi']:+.2f} ± {st['std_dchi']:.2f}\n")
        f.write(f"  Median Δχ²_ν      = {st['median_dchi']:+.2f}\n")
        f.write(f"  SEM               = {st['sem_dchi']:.3f}\n")
        f.write(f"  One-sample t-test : t={st['t_stat']:.3f}, p(one-tail)={st['t_pvalue_one_sided']:.4f}\n")
        if not np.isnan(st['wilc_pval']):
            f.write(f"  Wilcoxon sign-rank : W={st['wilc_stat']:.1f}, p={st['wilc_pval']:.4f}\n")
        f.write(f"  Redshift range     = {st['z_range'][0]:.2f} ~ {st['z_range'][1]:.2f}\n")
        f.write(f"  Median ΔBIC        = {st['median_dbic']:+.1f}\n\n")

print("\n✅ Stats saved: Triple_Scan_Stats.txt")


# ==========================================
# 9. 保存完整结果 CSV
# ==========================================

full_df = pd.DataFrame([{
    'id': r['id'],
    'field': r['field'],
    'group': r['group'],
    'z_phot': r['z_phot'],
    'n_bands': r['n_bands'],
    'ir_excess': r.get('ir_excess', np.nan),
    'r_eff': r.get('r_eff', np.nan),
    'z_dist_pred': r['z_dist'],
    'chi2nu_null': r.get('rch2_null', np.nan),
    'chi2nu_geff': r.get('rch2_geff', np.nan),
    'delta_chi2_nu': r['delta_chi2_nu'],
    'delta_BIC': r['delta_bic'],
    'verdict': r['verdict'],
} for r in all_results])

full_df.to_csv('Triple_Scan_Results.csv', index=False, float_format='%.6f')
print("✅ Results saved: Triple_Scan_Results.csv")


# ==========================================
# 10. 大图：三组对比综合图
# ==========================================

colors_map = {
    'Gold': '#FFD700',    # 金色
    'Silver': '#C0C0C0',  # 银色
    'Noise': '#CD853F',   # Peru 棕（控制组）
}
verdict_colors = {
    'Support': '#2ca02c',
    'Weak': '#98df8a',
    'Neutral': '#ffbb33',
    'Null wins': '#d62728',
}

fig = plt.figure(figsize=(22, 18))
gs = fig.add_gridspec(3, 3, height_ratios=[1.2, 1, 1], hspace=0.32, wspace=0.30)

# ─── Panel A (左上): 三组 Δχ²_ν 分布对比（箱线图+strip）───
ax_a = fig.add_subplot(gs[0, 0])

group_labels = []; group_data = []; group_colors = []
for st in all_stats:
    group_labels.append(st['name'].replace('🥇 ', '').replace('🥈 ', '').replace('🥉 ', ''))
    group_data.append(st['dchi_values'])
    group_colors.append(colors_map[st['name'].split()[-1]])

bp = ax_a.boxplot(group_data, labels=group_labels, patch_artist=True,
                   widths=0.5, showmeans=True, meanline=True,
                   meanprops={'color':'red','linewidth':2,'linestyle':'--'},
                   medianprops={'color':'black','linewidth':2})
for patch, color in zip(bp['boxes'], group_colors):
    patch.set_facecolor(color); patch.set_alpha(0.7); patch.set_edgecolor('black')

# strip plot overlay
for i, (data, color) in enumerate(zip(group_data, group_colors)):
    jitter = np.random.uniform(-0.15, 0.15, size=len(data))
    ax_a.scatter(np.full(len(data), i+1) + jitter, data,
                 c=color, alpha=0.4, s=12, edgecolor='none', zorder=5)

ax_a.axhline(0, color='black', linewidth=1.5, linestyle='-')
ax_a.set_ylabel(r'$\Delta \chi^2_\nu$', fontsize=13)
ax_a.set_title('(A) $\Delta\\chi^2_\\nu$ Distribution by Sample Quality', fontsize=13, fontweight='bold')

# 在每组上方标注均值
for i, st in enumerate(all_stats):
    ax_a.text(i+1, max(st['dchi_values']) + 3,
             f'Mean={st["mean_dchi"]:+.1f}\nn={st["N"]}',
             ha='center', va='bottom', fontsize=9,
             bbox=dict(facecolor='white', alpha=0.8, edgecolor=color))


# ─── Panel B (中上): 三组支持率对比柱状图 ───
ax_b = fig.add_subplot(gs[0, 1])

x_pos = np.arange(len(all_stats))
bar_width = 0.20

vtype_key_map = {'Support': 'support', 'Weak': 'weak', 'Neutral': 'neutral', 'Null wins': 'nullwin'}
for bi, vtype in enumerate(['Support', 'Weak', 'Neutral', 'Null wins']):
    vals = [st[vtype_key_map[vtype]] for st in all_stats]
    pct_vals = [100*v/st['N'] for v, st in zip(vals, all_stats)]
    bars = ax_b.bar(x_pos + bi*bar_width, pct_vals, bar_width,
                    label=vtype, color=verdict_colors[vtype],
                    edgecolor='black', alpha=0.8)
    # 数值标注
    for j, (b, v) in enumerate(zip(bars, vals)):
        if v > 0:
            ax_b.text(b.get_x() + b.get_width()/2, b.get_height() + 0.5,
                     f'{v}', ha='center', va='bottom', fontsize=7)

ax_b.set_xticks(x_pos + 1.5*bar_width)
ax_b.set_xticklabels([st['name'] for st in all_stats], fontsize=11)
ax_b.set_ylabel('% of Sources', fontsize=12)
ax_b.set_title('(B) Verdict Breakdown by Group', fontsize=13, fontweight='bold')
ax_b.legend(fontsize=9, loc='upper right')


# ─── Panel C (右上): 统计显著性表格 ───
ax_c = fig.add_subplot(gs[0, 2])
ax_c.axis('off')

table_text = "STATISTICAL SIGNIFICANCE SUMMARY\n" + "="*45 + "\n\n"
for st in all_stats:
    name_short = st['name'].replace('🥇 ','').replace('🥈 ','').replace('🥉 ','')
    table_text += f"[{name_short}]  N={st['N']}\n"
    table_text += f"  Favor G_eff: {st['favor_pct']:.1f}%\n"
    table_text += f"  Mean Δχ²_ν:  {st['mean_dchi']:+.2f} ± {st['sem_dchi']:.3f}\n"
    p_str = f"{st['t_pvalue_one_sided']:.4f}"
    if st['t_pvalue_one_sided'] < 0.001:
        sig_str = " ***"
    elif st['t_pvalue_one_sided'] < 0.01:
        sig_str = " **"
    elif st['t_pvalue_one_sided'] < 0.05:
        sig_str = " *"
    else:
        sig_str = " (n.s.)"
    table_text += f"  One-tailed t:  p = {p_str}{sig_str}\n"
    if not np.isnan(st.get('wilc_pval', np.nan)):
        wp = st['wilc_pval']
        wsig = " ***" if wp < 0.001 else " **" if wp < 0.01 else " *" if wp < 0.05 else " (n.s.)"
        table_text += f"  Wilcoxon:     W={st['wilc_stat']:.0f}, p={wp:.4f}{wsig}\n"
    table_text += "\n"

# Gold vs Noise
if stats_gold and stats_noise:
    gu, gp = u_stat, u_pval  # from earlier M-W test
    gsignif = " ✅ YES" if gp < 0.05 else " ❌ No"
    table_text += f"GOLD vs NOISE (M-W U):\n"
    table_text += f"  Is Gold > Noise?  p = {gp:.4f}{gsignif}\n"

table_text += "\nSignificance: *** p<0.001, ** p<0.01, * p<0.05"

ax_c.text(0.05, 0.98, table_text, transform=ax_c.transAxes, fontsize=11,
          family='monospace', verticalalignment='top',
          bbox=dict(facecolor='#f8f8f8', edgecolor='gray', boxstyle='round,pad=0.8'))
ax_c.set_title('(C) Statistical Tests Summary', fontsize=13, fontweight='bold')


# ─── Panel D (左下): 三组散点图（按颜色区分）───
ax_d = fig.add_subplot(gs[1, :])

group_color_list = []
for st in all_stats:
    gn = st['name'].split()[-1]  # Gold/Silver/Noise
    for _ in range(st['N']):
        group_color_list.append(colors_map.get(gn, '#333333'))

flat_zs = []; flat_dchis = []; flat_verdicts = []
for st in all_stats:
    flat_zs.extend(st['zs'])
    flat_dchis.extend(st['dchi_values'])
    flat_verdicts.extend(st['verdicts'])

sc = ax_d.scatter(flat_zs, flat_dchis, c=group_color_list, s=35,
                  edgecolor='black', alpha=0.65, linewidth=0.4)

# 各组分 bin 均值线
for idx, st in enumerate(all_stats):
    gn = st['name'].split()[-1]
    zs_arr = np.array(st['zs']); dc_arr = np.array(st['dchi_values'])
    zbins = np.linspace(min(zs_arr)-0.1, max(zs_arr)+0.1, 5)
    bcs, bms = [], []
    for i in range(len(zbins)-1):
        mask = (zs_arr >= zbins[i]) & (zs_arr < zbins[i+1])
        if sum(mask) >= 2:
            bcs.append((zbins[i]+zbins[i+1])/2)
            bms.append(np.mean(dc_arr[mask]))
    if bcs:
        ax_d.plot(bcs, bms, 'o-', color=colors_map[gn], markersize=8,
                 linewidth=2.5, markerfacecolor='white', markeredgewidth=2,
                 label=f"{gn} bin-mean", zorder=15)

ax_d.axhline(0, color='black', linewidth=1.5)
ax_d.set_xlabel(r'$z_{\rm phot}$', fontsize=14)
ax_d.set_ylabel(r'$\Delta \chi^2_\nu$', fontsize=14)
ax_d.set_title('(D) All Sources: $\\Delta\\chi^2_\\nu$ vs Redshift, Color-coded by Group', 
               fontsize=13, fontweight='bold')

# 图例
from matplotlib.patches import Patch
legend_elements = [Patch(fc=c, ec='black', alpha=0.7, label=k) 
                   for k, c in [('Gold', colors_map['Gold']), 
                                ('Silver', colors_map['Silver']),
                                ('Noise', colors_map['Noise'])]]
ax_d.legend(handles=legend_elements, loc='upper right', fontsize=10)


# ─── Panel E (底部左): Gold 样本详细分布 ───
ax_e = fig.add_subplot(gs[2, 0])

if stats_gold:
    eg_dchi = stats_gold['dchi_values']
    eg_z = stats_gold['zs']
    eg_v = stats_gold['verdicts']
    
    bins_g = np.linspace(min(eg_dchi)-1, max(eg_dchi)+1, 20)
    vc = [verdict_colors.get(v, 'gray') for v in eg_v]
    ax_e.hist(eg_dchi, bins=bins_g, color=colors_map['Gold'], alpha=0.7,
              edgecolor='black', linewidth=0.5)
    ax_e.axvline(0, color='black', lw=1.5)
    ax_e.axvline(stats_gold['mean_dchi'], color='red', lw=2, ls='--',
                label=f'Mean={stats_gold["mean_dchi"]:+.1f}')
    ax_e.set_xlabel(r'$\Delta \chi^2_\nu$', fontsize=12)
    ax_e.set_ylabel('Count', fontsize=12)
    ax_e.set_title(f'(E) Gold Only (N={stats_gold["N"]})', fontsize=12, fontweight='bold')
    ax_e.legend(fontsize=9)
    
    # 统计文本
    gt = (f"Favor: {stats_gold['favor_pct']:.0f}%\n"
          f"Mean: {stats_gold['mean_dchi']:+.2f}\n"
          f"p(t-test): {stats_gold['t_pvalue_one_sided']:.4f}")
    ax_e.text(0.97, 0.97, gt, transform=ax_e.transAxes, fontsize=9,
              va='top', ha='right', family='monospace',
              bbox=dict(fc='white', alpha=0.85))


# ─── Panel F (底部中): Silver 样本详细分布 ───
ax_f = fig.add_subplot(gs[2, 1])

if stats_silver:
    es_dchi = stats_silver['dchi_values']
    es_z = stats_silver['zs']
    
    bins_s = np.linspace(min(es_dchi)-1, max(es_dchi)+1, 25)
    ax_f.hist(es_dchi, bins=bins_s, color=colors_map['Silver'], alpha=0.7,
              edgecolor='black', linewidth=0.5)
    ax_f.axvline(0, color='black', lw=1.5)
    ax_f.axvline(stats_silver['mean_dchi'], color='red', lw=2, ls='--',
                label=f'Mean={stats_silver["mean_dchi"]:+.1f}')
    ax_f.set_xlabel(r'$\Delta \chi^2_\nu$', fontsize=12)
    ax_f.set_ylabel('Count', fontsize=12)
    ax_f.set_title(f'(F) Silver Only (N={stats_silver["N"]})', fontsize=12, fontweight='bold')
    ax_f.legend(fontsize=9)
    
    stxt = (f"Favor: {stats_silver['favor_pct']:.0f}%\n"
            f"Mean: {stats_silver['mean_dchi']:+.2f}\n"
            f"p(t-test): {stats_silver['t_pvalue_one_sided']:.4f}")
    ax_f.text(0.97, 0.97, stxt, transform=ax_f.transAxes, fontsize=9,
              va='top', ha='right', family='monospace',
              bbox=dict(fc='white', alpha=0.85))


# ─── Panel G (底部右): Noise 控制组详细分布 ───
ax_g = fig.add_subplot(gs[2, 2])

if stats_noise:
    en_dchi = stats_noise['dchi_values']
    en_z = stats_noise['zs']
    
    bins_n = np.linspace(min(en_dchi)-1, max(en_dchi)+1, 25)
    ax_g.hist(en_dchi, bins=bins_n, color=colors_map['Noise'], alpha=0.7,
              edgecolor='black', linewidth=0.5)
    ax_g.axvline(0, color='black', lw=1.5)
    ax_g.axvline(stats_noise['mean_dchi'], color='red', lw=2, ls='--',
                label=f'Mean={stats_noise["mean_dchi"]:+.1f}')
    ax_g.set_xlabel(r'$\Delta \chi^2_\nu$', fontsize=12)
    ax_g.set_ylabel('Count', fontsize=12)
    ax_g.set_title(f'(G) Noise Control (N={stats_noise["N"]})', fontsize=12, fontweight='bold')
    ax_g.legend(fontsize=9)
    
    ntxt = (f"Favor: {stats_noise['favor_pct']:.0f}%\n"
            f"Mean: {stats_noise['mean_dchi']:+.2f}\n"
            f"p(t-test): {stats_noise['t_pvalue_one_sided']:.4f}")
    ax_g.text(0.97, 0.97, ntxt, transform=ax_g.transAxes, fontsize=9,
              va='top', ha='right', family='monospace',
              bbox=dict(fc='white', alpha=0.85))

# 总标题
fig.suptitle(
    r"Triple Scan: $G_{\rm eff}$ Window Test on Classified LRD Samples"
    f"\n(Gold: z≥7, compact, moderate IR; Silver: mid-z, auxiliary; Noise: low-z/extended/extreme-IR control)",
    fontsize=15, y=1.005
)

fig.savefig('Triple_Scan_Summary.png', dpi=250, bbox_inches='tight')
print("✅ Summary saved: Triple_Scan_Summary.png")


# ==========================================
# 11. Gold 样本 SED PDF（最重要的子集！）
# ==========================================

if stats_gold and stats_gold['valid_results']:
    with PdfPages('Triple_Scan_Gold_SEDs.pdf') as pdf:
        # 封面
        fig_cov = plt.figure(figsize=(11, 8.5))
        ax_cov = fig_cov.add_subplot(111)
        ax_cov.axis('off')
        
        cov_txt = (
            f"GOLD SAMPLE SED COMPARISON\n"
            f"Purist LRD Selection: X-ray Silent, Compact, High-z, Moderate IR\n\n"
            f"Selection Criteria:\n"
            f"  • z_phot ≥ 7 (EoR LRD candidates)\n"
            f"  • r_eff < median ({R_EFF_MEDIAN:.0f} pc) — point-source like\n"
            f"  • f444W/f150W < P75 ({IR_P75:.1f}) — non-AGN infrared\n"
            f"  • Nbands ≥ 6 — good photometric quality\n\n"
            f"Statistics:\n"
            f"  N = {stats_gold['N']}\n"
            f"  Support G_eff = {stats_gold['support']} ({100*stats_gold['support']/stats_gold['N']:.1f}%)\n"
            f"  Weak = {stats_gold['weak']} ({100*stats_gold['weak']/stats_gold['N']:.1f}%)\n"
            f"  Neutral = {stats_gold['neutral']} ({100*stats_gold['neutral']/stats_gold['N']:.1f}%)\n"
            f"  Null wins = {stats_gold['nullwin']} ({100*stats_gold['nullwin']/stats_gold['N']:.1f}%)\n\n"
            f"  Mean Δχ²_ν = {stats_gold['mean_dchi']:+.2f}\n"
            f"  t-test p (one-tailed) = {stats_gold['t_pvalue_one_sided']:.4f}\n\n"
            f"Parameters:\n"
            f"  A_G={window_optimal['A_G']}, γ={GAMMA_OPTIMAL}, β={window_optimal['beta']}, z_on={window_optimal['z_on']}"
        )
        ax_cov.text(0.5, 0.5, cov_txt, transform=ax_cov.transAxes, fontsize=12,
                    ha='center', va='center', family='monospace',
                    bbox=dict(fc='white', ec='gray', boxstyle='round,pad=1'))
        pdf.savefig(fig_cov, dpi=200)
        plt.close(fig_cov)
        
        # 每个 Gold 源画 SED
        def plot_gold_sed(result):
            fg = plt.figure(figsize=(10, 7))
            gsfg = fg.add_gridspec(3, 1, height_ratios=[3, 1, 0.5], hspace=0.18)
            
            axm = fg.add_subplot(gsfg[0])
            axr = fg.add_subplot(gsfg[1], sharex=axm)
            axi = fg.add_subplot(gsfg[2]); axi.axis('off')
            
            wav_sm = np.logspace(-1, 1.2, 300)
            m_ns = sed_model(wav_sm, result['Cb_null'], result['Cr_null'], z_dist=0.0)
            m_gs = sed_model(wav_sm, result['Cb_geff'], result['Cr_geff'], z_dist=result['z_dist'])
            
            axm.plot(wav_sm, m_gs, color='#2ca02c', lw=2.5,
                    label=r'$G_{\rm eff}$ (z_dist=' + f"{result['z_dist']:.4f})")
            axm.plot(wav_sm, m_ns, color='#1f77b4', lw=1.8, ls='--',
                    label=r'Null ($z_{\rm dist}=0$)')
            axm.errorbar(result['rest_wave'], result['nuFnu_norm'], yerr=result['nuFnu_err_norm'],
                        fmt='o', color='black', ms=7, capsize=2, zorder=10)
            
            axm.set_xscale('log'); axm.set_xlim(0.1, 15)
            axm.set_ylim(0, max(max(m_gs), max(m_ns))*1.3)
            plt.setp(axm.get_xticklabels(), visible=False)
            axm.set_ylabel(r'Norm. $\nu F_\nu$', fontsize=12)
            axm.legend(loc='upper right', fontsize=10)
            
            vc = verdict_colors[result['verdict']]
            info = (f"ID={result['id']}  z={result['z_phot']:.2f}  field={result['field']}\n"
                   f"$\\Delta\\chi^2_\\nu$ = {result['delta_chi2_nu']:+.2f}  |  "
                   f"$\\Delta$BIC = {result['delta_bic']:+.1f}  |  {result['verdict']}\n"
                   f"IR_excess = {result.get('ir_excess',np.nan):.2f}  r_eff = {result.get('r_eff',np.nan):.1f} pc")
            axm.text(0.03, max(max(m_gs), max(m_ns))*1.2, info, fontsize=10,
                     bbox=dict(fc='white', alpha=0.88, ec=vc, lw=1.5), va='top')
            
            resg = result['residuals_geff']
            resn = (result['nuFnu_norm'] - result['model_null']) / result['nuFnu_err_norm']
            axr.axhline(0, color='black', lw=1)
            axr.errorbar(result['rest_wave'], resg, yerr=1.0, fmt='o', color='#2ca02c', ms=6, label=r'$G_{\rm eff}$')
            axr.errorbar(result['rest_wave'], resn, yerr=1.0, fmt='s', color='#1f77b4', ms=5, alpha=0.6, label='Null')
            axr.set_ylim(-4, 4); axr.set_xticks([0.1, 0.5, 1, 5, 10])
            axr.set_xlabel(r'Rest-frame $\lambda$ ($\mu$m)', fontsize=12)
            axr.set_ylabel(r'$\Delta \chi$', fontsize=11)
            axr.legend(fontsize=9, ncol=2)
            
            fg.suptitle(f"Gold Sample: ID={result['id']}", fontsize=13, fontweight='bold')
            return fg
        
        # 按 Δχ²_ν 排序（从好到差）
        gold_valid = sorted(stats_gold['valid_results'], key=lambda x: x['delta_chi2_nu'], reverse=True)
        for i, r in enumerate(gold_valid):
            fg = plot_gold_sed(r)
            pdf.savefig(fg, dpi=180)
            plt.close(fg)
    
    print(f"✅ Gold SED PDF saved: Triple_Scan_Gold_SEDs.pdf ({len(stats_gold['valid_results'])+1} pages)")


# ==========================================
# 12. 最终输出清单
# ==========================================

print(f"\n{'='*60}")
print(f"TRIPLE SCAN OUTPUT FILES:")
print(f"{'='*60}")
print(f"  1. Triple_Scan_Stats.txt       — 统计检验详情")
print(f"  2. Triple_Scan_Results.csv     — 含分组标签的全结果")
print(f"  3. Triple_Scan_Summary.png     — 三组对比综合大图")
print(f"  4. Triple_Scan_Gold_SEDs.pdf   — Gold 样本 SED 逐一展示")
print(f"{'='*60}")
print(f"\n✅✅✅ TRIPLE SCAN COMPLETE ✅✅✅")
