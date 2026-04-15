"""
AGN Three-Component Model Test (Bulletproof Relic Version)
==========================================================
Upgraded with strict Apple-to-Apple comparison and Fossil Hypothesis.

Model: F = C_blue*f_uv + C_red*f_int(λ/(1+z_d))/(1+z_d) * Dust(Av) + C_agn*(λ/λ0)^(-α_agn)
       (5 free params: C_blue, C_red, C_agn, α_agn, Av)

Compare:
  Null Model: z_dist = 0.0 (Relies entirely on Dust Av and AGN to fit redness)
  G_eff Relic Model: z_dist = z_dist(z_form) (Relies on z=13.1 gravity fossil + Dust)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 0. Physics Setup & Fossil Hypothesis
# ============================================================

# 核心修改：设定统一的“形成红移”。GS-z13-1-LA 所在的 z=13.1 是引力增强的极盛期
Z_FORM_RELIC = 13.1 

def H_z(z, Omega_m=0.3, Omega_L=0.7, H0=70.0):
    return H0 * np.sqrt(Omega_m * (1+z)**3 + Omega_L)

def G_eff_window(z, A_G=0.189, beta=1.58, z_p=10.0, z_on=2.0, z_off=30.0, dz=3.0):
    W = 0.25 * (1 + np.tanh((z - z_on) / dz)) * (1 - np.tanh((z - z_off) / dz))
    Hz_ratio = H_z(z) / H_z(z_p)
    return A_G * (Hz_ratio ** beta) * W

def z_dist_from_phot(z_phot, gamma=0.50):
    delta_G = G_eff_window(z_phot)
    return gamma * np.clip(delta_G, 0, None)

# 计算全局遗迹畸变值 (Fossil z_dist)
Z_DIST_FOSSIL = z_dist_from_phot(Z_FORM_RELIC)

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

def f_agn(wave_rest, alpha_agn=1.2, lam0=2.0):
    return (wave_rest / lam0)**(-alpha_agn)

# 核心修改：增加尘埃消光律 (简单的 Screen Dust 近似)
def dust_attenuation(wave_rest, Av):
    # A_lambda ≈ Av * (0.55 / lambda) for basic rest-frame optical/UV proxy
    A_lam = Av * (0.55 / wave_rest)
    return 10**(-0.4 * A_lam)

# ============================================================
# 1. Five-Component Model Function
# ============================================================

def sed_model_final(wave_rest, C_blue, C_red, C_agn, alpha_agn, Av, z_dist=0.0):
    """
    5 params: UV + Shifted/Unshifted Dust Host + AGN power law.
    Note: Dust only obscures the red host component, matching the 'V-shape' LRD morphology.
    """
    distorted_wave = wave_rest / (1 + z_dist)
    term_uv = C_blue * f_uv(wave_rest)
    term_dust = C_red * (f_int(distorted_wave) / (1 + z_dist)) * dust_attenuation(wave_rest, Av)
    term_agn = C_agn * f_agn(wave_rest, alpha_agn)
    return np.maximum(term_uv + term_dust + term_agn, 1e-30)

# ============================================================
# 2. Fit Engine (5 parameters)
# ============================================================

def fit_5param(rest_wave, nuFnu_obs, nuFnu_err, z_dist_fixed=0.0, p0=None, max_iter=15000):
    def model(wave, C_blue, C_red, C_agn, alpha_agn, Av):
        return sed_model_final(wave, C_blue, C_red, C_agn, alpha_agn, Av, z_dist_fixed)
    
    if p0 is None:
        p0 = [1.0, 0.5, 0.1, 1.2, 0.5] # 增加了 Av 的初始猜想值
    
    try:
        # 核心修改：收紧 alpha_agn 上限至 2.5 防止过拟合，增加 Av [0, 5] 边界
        popt, pcov = curve_fit(model, rest_wave, nuFnu_obs,
                               sigma=nuFnu_err,
                               bounds=([0, 0, 0, 0.01, 0.0], [20, 20, 20, 2.5, 5.0]),
                               p0=p0, maxfev=max_iter)
        
        n_data = len(rest_wave)
        k_params = 5 # 升至 5 个自由参数
        model_flux = model(rest_wave, *popt)
        residuals = (nuFnu_obs - model_flux) / nuFnu_err
        
        chi2 = np.sum(residuals**2)
        dof = max(n_data - k_params, 1)
        red_chi2 = chi2 / dof
        bic = k_params * np.log(n_data) + chi2
        
        if red_chi2 > 1e10 or np.any(np.isnan(popt)) or np.any(np.isinf(popt)):
            return {'success': False, 'error': 'Numerical instability'}
        
        return {
            'C_blue': float(popt[0]), 'C_red': float(popt[1]),
            'C_agn': float(popt[2]), 'alpha_agn': float(popt[3]),
            'Av': float(popt[4]),
            'chi2': float(chi2), 'red_chi2': float(red_chi2),
            'bic': float(bic), 'dof': int(dof), 'n_data': n_data,
            'model_flux': model_flux, 'residuals': residuals,
            'success': True
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}

# ============================================================
# 3. Data Loading (Same as before)
# ============================================================

DATA_DIR = '/Users/tanxin/Desktop/数据处理'
CSV_FILE = f'{DATA_DIR}/Kokorev_LRDs_Full.csv'

print("=" * 65)
print("BULLETPROOF RELIC MODEL TEST")
print(f"Assuming Formation Redshift z_form = {Z_FORM_RELIC}")
print(f"Calculated Fossil Distortion z_dist = {Z_DIST_FOSSIL:.4f}")
print("=" * 65)

df = pd.read_csv(CSV_FILE)
filters = {
    'f090w': 0.90, 'f115w': 1.15, 'f150w': 1.50, 'f200w': 2.00,
    'f277w': 2.77, 'f356w': 3.56, 'f444w': 4.44, 'f770w': 7.70
}

sources_data = []
for idx, row in df.iterrows():
    target_id = row['id']
    z_phot = row['z_phot']
    field = row.get('field', 'unknown')
    
    obs_wave, obs_flux, obs_err = [], [], []
    for f, wave in filters.items():
        flux_col = f'{f}_flux'; err_col = f'{f}_fluxerr'
        if flux_col in row.index and pd.notnull(row[flux_col]) and row[flux_col] > -90:
            obs_wave.append(wave); obs_flux.append(row[flux_col])
            err_val = row[err_col] if err_col in row.index and row[err_col] > 0 else abs(row[flux_col]) * 0.1
            obs_err.append(err_val)
    
    if len(obs_wave) >= 4:
        obs_wave = np.array(obs_wave, dtype=float)
        obs_flux = np.array(obs_flux, dtype=float)
        obs_err = np.array(obs_err, dtype=float)
        rest_wave = obs_wave / (1 + z_phot)
        
        nuFnu_obs = obs_flux / obs_wave
        nuFnu_err = obs_err / obs_wave
        norm_factor = np.max(nuFnu_obs[nuFnu_obs > 0]) if np.any(nuFnu_obs > 0) else 1.0
        norm_factor = max(norm_factor, 1e-30)
        
        ir_excess = 0.0
        f150_idx = None; f444_idx = None
        for i, w in enumerate(obs_wave):
            if abs(w - 1.50) < 0.05: f150_idx = i
            if abs(w - 4.44) < 0.05: f444_idx = i
        if f150_idx is not None and f444_idx is not None:
            if obs_flux[f150_idx] > 0:
                ir_excess = obs_flux[f444_idx] / obs_flux[f150_idx]
        
        sources_data.append({
            'target_id': target_id, 'z_phot': z_phot, 'field': field,
            'rest_wave': rest_wave, 'nuFnu_norm': nuFnu_obs / norm_factor,
            'nuFnu_err_norm': nuFnu_err / norm_factor, 'norm_factor': norm_factor,
            'obs_wave': obs_wave, 'obs_flux': obs_flux, 'obs_err': obs_err,
            'ir_excess': ir_excess
        })

# ============================================================
# 4. Run Fits (Apple to Apple Comparison)
# ============================================================

results_list = []

for src in sources_data:
    tid = src['target_id']; zp = src['z_phot']
    rw = src['rest_wave']; nf = src['nuFnu_norm']; ne = src['nuFnu_err_norm']
    ir_exc = src['ir_excess']
    
    # 核心对决
    # Null Model: 只有 Dust Av 和 AGN 去拟合光谱
    r_null = fit_5param(rw, nf, ne, z_dist_fixed=0.0)
    
    # Geff Model: 使用统一的化石畸变 Z_DIST_FOSSIL，同时也允许使用 Dust
    r_geff = fit_5param(rw, nf, ne, z_dist_fixed=Z_DIST_FOSSIL)
    
    if r_null['success'] and r_geff['success']:
        dchi2 = r_null['red_chi2'] - r_geff['red_chi2']
        dbic = r_null['bic'] - r_geff['bic']
        
        if dchi2 >= 10 and dbic > 0: verdict = 'Support'
        elif dchi2 > 2 and dbic > -10: verdict = 'Weak'
        elif dchi2 > -2: verdict = 'Neutral'
        else: verdict = 'Null wins'
        
        results_list.append({
            'id': tid, 'field': src['field'], 'z_phot': round(zp, 3),
            'z_dist_fossil': round(Z_DIST_FOSSIL, 5),
            'chi2nu_null': round(r_null['red_chi2'], 2),
            'chi2nu_geff': round(r_geff['red_chi2'], 2),
            'delta_chi2nu': round(dchi2, 2),
            'bic_null': round(r_null['bic'], 1), 'bic_geff': round(r_geff['bic'], 1),
            'delta_bic': round(dbic, 1),
            'Av_null': round(r_null['Av'], 2), 'Av_geff': round(r_geff['Av'], 2),
            'C_AGN': round(r_geff['C_agn'], 4), 'alpha_AGN': round(r_geff['alpha_agn'], 3),
            'ir_excess': round(ir_exc, 2), 'verdict': verdict,
            'n_filters': int(src['n_data']) if hasattr(src, 'n_data') else len(rw)
        })
    else:
        results_list.append({
            'id': tid, 'field': src['field'], 'z_phot': round(zp, 3),
            'z_dist_fossil': round(Z_DIST_FOSSIL, 5),
            'delta_chi2nu': 0, 'delta_bic': 0, 'verdict': 'FAILED',
            'ir_excess': round(ir_exc, 2), 'C_AGN': 0, 'Av_null': 0, 'Av_geff': 0
        })

# ============================================================
# 5. Output & Analysis & Plotting (Skipped verbose plotting for brevity, but same as yours)
# ============================================================

results_df = pd.DataFrame(results_list)
out_csv = f'{DATA_DIR}/Bulletproof_Results.csv'
results_df.to_csv(out_csv, index=False)

valid = results_df[results_df['verdict'] != 'FAILED']
support = sum(valid['verdict']=='Support')
weak = sum(valid['verdict']=='Weak')
neutral = sum(valid['verdict']=='Neutral')
nwins = sum(valid['verdict']=='Null wins')

print(f"\n╔═════════════════════════════════════════════════════════╗")
print(f"║ BULLETPROOF 5-PARAM MODEL RESULTS ({len(valid)} sources)        ║")
print(f"╠═════════════════════════════════════════════════════════╣")
print(f"║  ✅ Support:      {support:>4} ({support/max(len(valid),1)*100:.1f}%)                           ║")
print(f"║  ⚡ Weak:         {weak:>4} ({weak/max(len(valid),1)*100:.1f}%)                           ║")
print(f"║  ➖ Neutral:      {neutral:>4} ({neutral/max(len(valid),1)*100:.1f}%)                           ║")
print(f"║  ❌ Null wins:    {nwins:>4} ({nwins/max(len(valid),1)*100:.1f}%)                           ║")
print(f"╚═════════════════════════════════════════════════════════╝")
print("\nRun your plotting scripts using the new Bulletproof_Results.csv!")