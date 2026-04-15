#!/usr/bin/env python3
"""
Gold Positive Deep SED Analysis
===============================
Focus on Gold-sample LRDs that show significant G_eff support (z >= 7).
These are our "high-redshift witnesses" for the G_eff window model.

Output: 
  - GoldPositive_SED_Deep.pdf — multi-page deep SED analysis
"""

import numpy as np
import csv
import os

# ─── Constants ───
WORK_DIR = '/Users/tanxin/Desktop/数据处理'

# Best parameters from param scan v5
AG_BEST = 0.189
GAMMA_BEST = 0.50
BETA_BEST = 1.58
ZON_BEST = 2.0

# NIRCam filters (microns)
FILTERS = ['F090W', 'F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']
LAMBDA_CENTRAL = [0.90, 1.15, 1.50, 2.00, 2.77, 3.56, 4.10, 4.44]

# Filter width for plotting (microns)
FILTER_WIDTHS = [0.12, 0.25, 0.30, 0.40, 0.45, 0.55, 0.60, 0.70]

def geff_window(z, AG=AG_BEST, BETA=BETA_BEST, ZON=ZON_BEST):
    """Tanh window: G_eff/G_N = 1 + A_G * 0.5 * [1 + tanh(beta*(z - z_on))]"""
    return AG * 0.5 * (1.0 + np.tanh(BETA * (z - ZON)))

def z_dist_prediction(z_phot, gamma=GAMMA_BEST):
    return gamma * geff_window(z_phot)

def uv_powerlaw(lam, C_b, beta_uv=-2.5):
    """UV power law: F_nu ~ nu^beta_uv ~ lam^(-beta_uv)"""
    lam = np.atleast_1d(lam)
    return C_b * (lam ** (-beta_uv))

def dust_blackbody(lam, T_dust=80.0):
    """Dust blackbody in F_nu units."""
    h = 6.626e-34
    c = 3.0e8
    k = 1.381e-23
    lam_um = np.atleast_1d(lam) * 1e-6
    x = h * c / (lam_um * k * T_dust)
    result = np.where(x > 500, 1e-99, (lam**(-3)) / (np.exp(np.minimum(x, 500)) - 1))
    return result[0] if np.ndim(lam) == 0 else result

def sed_model(lam, C_b, C_r, z_dist, beta_uv=-2.5, T_dust=80.0, z_phot=5.0):
    """Two-component SED model with redshift distortion."""
    lam = np.atleast_1d(lam)
    result = C_b * (lam ** (-beta_uv)) + C_r * dust_blackbody(np.atleast_1d(lam)/(1+z_dist), T_dust) / (1+z_dist)
    return float(result[0]) if len(result) == 1 else result

def chi_squared(observed, model, errors):
    """Reduced chi-squared."""
    obs = np.atleast_1d(observed)
    mod = np.atleast_1d(model)
    err = np.atleast_1d(errors)
    mask = err > 0
    resid = ((obs[mask] - mod[mask]) / err[mask]) ** 2
    n = int(np.sum(mask))
    if n <= 2:
        return np.inf
    return float(np.sum(resid) / (n - 2))

def bic(chi2nu, n_data, n_params=2):
    """BIC from reduced chi2."""
    dof = max(n_data - n_params, 1)
    total_chi2 = chi2nu * dof
    return float(total_chi2 + n_params * np.log(max(n_data, 1)))

def load_full_csv(filename):
    """Load Kokorev LRD catalog."""
    data = {}
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[row['id']] = row
    return data

def get_fluxes(row, filters=FILTERS):
    """Extract flux and error arrays."""
    fluxes, errors, valid_filters = [], [], []
    for f in filters:
        fl_key = f.lower() + '_flux'
        er_key = fl_key.replace('flux', 'fluxerr')
        
        if fl_key in row:
            val = row[fl_key].strip()
            if not val or val == '' or val.upper() == 'NAN':
                continue
            try:
                fv = float(val)
            except ValueError:
                continue
            if np.isnan(fv) or abs(fv) < 1e-32 or abs(fv) > 1e10 or fv > 999:
                continue
            try:
                ev = float(row.get(er_key, '').strip() or '0')
            except (ValueError, KeyError):
                ev = abs(fv) * 0.1
            if ev <= 0:
                ev = max(abs(fv)*0.05, 1e-30)
            fluxes.append(fv)
            errors.append(ev)
            valid_filters.append(f)
    
    return np.array(fluxes), np.array(errors), valid_filters

def fit_source(fluxes, errors, lambdas, z_phot, zd_val=None):
    """Fit two-component model to source."""
    try:
        from scipy.optimize import curve_fit
        
        def model_fixed_zd(lam, C_b, C_r):
            return sed_model(lam, C_b, C_r, zd_val, z_phot=z_phot)
        
        def model_free(lam, C_b, C_r):
            return sed_model(lam, C_b, C_r, 0.0, z_phot=z_phot)
        
        # Use first point for UV estimate, last for IR
        f_uv_est = max(abs(fluxes[0]), 1e-30) if len(fluxes) > 0 else 1e-28
        f_ir_est = max(abs(fluxes[-1]) * 0.3, 1e-31) if len(fluxes) > 1 else 1e-29
        p0 = [f_uv_est, f_ir_est]
        
        popt_null, _ = curve_fit(model_free, lambdas, fluxes,
                                  sigma=errors, absolute_sigma=True,
                                  p0=p0, bounds=([0,0],[np.inf,np.inf]),
                                  maxfev=20000)
        
        popt_geff, _ = curve_fit(model_fixed_zd, lambdas, fluxes,
                                  sigma=errors, absolute_sigma=True,
                                  p0=p0, bounds=([0,0],[np.inf,np.inf]),
                                  maxfev=20000)
        
        return popt_null, popt_geff
    except Exception as e:
        print(f"  Fit failed: {e}")
        return None, None

# ─── Main ───
print("="*70)
print("GOLD POSITIVE DEEP SED ANALYSIS")
print("="*70)

# Load data
full_data = load_full_csv(os.path.join(WORK_DIR, 'Kokorev_LRDs_Full.csv'))

# Target sources: top Gold positives + top Silver z=6-7 for comparison
targets = [
    ('8559',   7.94,  5.8, 'Gold+ 🏆', 53.1),
    ('67076',  7.94,  5.8, 'Gold+',     51.5),
    ('40589',  8.01,  6.7, 'Gold+',     12.1),
    ('47021',  8.02,  3.8, 'Gold+',      8.7),
    ('36162',  6.41,  4.4, 'Silver🔥', 661.1),
    ('42093',  6.45, 16.7, 'Silver🔥', 658.6),
    ('26971',  6.41,  4.3, 'Silver🔥', 657.7),
]

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available, outputting text only")

results = []

for src_id, zph, ir_excess, label, dchi2 in targets:
    if src_id not in full_data:
        print(f"  {src_id}: NOT FOUND in catalog")
        continue
    
    row = full_data[src_id]
    
    # Get photometry
    fluxes, errors, vf = get_fluxes(row)
    nlam = np.array([LAMBDA_CENTRAL[FILTERS.index(f)] for f in vf])
    
    if len(vf) < 4:
        print(f"  {src_id}: Only {len(vf)} bands, skipping")
        continue
    
    # Compute predicted z_dist
    zd_pred = z_dist_prediction(zph)
    
    # Fit both models
    popt_null, popt_geff = fit_source(fluxes, errors, nlam, zph, zd_pred)
    
    if popt_null is None or popt_geff is None:
        continue
    
    # Model predictions on fine grid
    lam_fine = np.linspace(0.85, 4.8, 300)
    mod_null = np.array([sed_model(l, popt_null[0], popt_null[1], 0.0, z_phot=zph) for l in lam_fine])
    mod_geff = np.array([sed_model(l, popt_geff[0], popt_geff[1], zd_pred, z_phot=zph) for l in lam_fine])
    
    # Statistics
    pred_null = np.array([sed_model(l, popt_null[0], popt_null[1], 0.0, z_phot=zph) for l in nlam])
    pred_geff = np.array([sed_model(l, popt_geff[0], popt_geff[1], zd_pred, z_phot=zph) for l in nlam])
    
    chisq_nu_null = chi_squared(fluxes, pred_null, errors)
    chisq_nu_geff = chi_squared(fluxes, pred_geff, errors)
    
    dc2 = chisq_nu_null - chisq_nu_geff
    dBIC = bic(chisq_nu_null, len(vf), 2) - bic(chisq_nu_geff, len(vf), 2)
    
    results.append({
        'id': src_id, 'z': zph, 'ir': ir_excess, 'label': label,
        'zd': zd_pred, 'dc2': dc2, 'dBIC': dBIC,
        'chisq_null': chisq_nu_null, 'chisq_geff': chisq_nu_geff,
        'Cb_null': popt_null[0], 'Cr_null': popt_null[1],
        'Cb_geff': popt_geff[0], 'Cr_geff': popt_geff[1],
        'fluxes': fluxes, 'errors': errors, 'lam_obs': nlam,
        'vf': vf, 'mod_null': mod_null, 'mod_geff': mod_geff,
        'lam_fine': lam_fine
    })
    
    print(f"  {src_id} ({label}): z={zph:.2f}, zd={zd_pred:.4f}, Δχ²_ν={dc2:+.1f}, "
          f"χ²null={chisq_nu_null:.1f}, χ²geff={chisq_nu_geff:.1f}")

if HAS_MATPLOTLIB and results:
    try:
        fig = plt.figure(figsize=(18, 24))
        gs = GridSpec(len(results), 2, figure=fig, height_ratios=[1]*len(results),
                      hspace=0.35, wspace=0.28)
        
        colors_gold = {'Gold+ 🏆': '#FFD700', 'Gold+': '#FFA500'}
        colors_silver = {'Silver🔥': '#C0C0C0'}
        
        for i, r in enumerate(results):
            ax_a = fig.add_subplot(gs[i, 0])
            
            color = colors_gold.get(r['label'], colors_silver.get(r['label'], '#333'))
            
            ax_a.errorbar(r['lam_obs'], r['fluxes'] * 1e6, yerr=r['errors']*1e6,
                         fmt='o', ms=9, color='#222', ecolor='gray', capsize=3,
                         label='JWST Photometry', zorder=5, alpha=0.85)
            
            ax_a.plot(r['lam_fine'], r['mod_null'] * 1e6, '--',
                     color='#3498db', lw=2.2, alpha=0.85, label=f'Null ($z_d$=0)')
            
            ax_a.plot(r['lam_fine'], r['mod_geff'] * 1e6, '-',
                     color='#e74c3c', lw=2.5, alpha=0.95,
                     label=f'G$_{{eff}}$ ($z_d$={r["zd"]:.3f})')
            
            ax_a.axvspan(3.0, 4.6, alpha=0.06, color='red')
            ax_a.set_xscale('log')
            ax_a.set_xlabel(r'Rest-frame $\lambda$ ($\mu$m)', fontsize=11)
            ax_a.set_ylabel(r'$\lambda F_\lambda$', fontsize=10)
            ax_a.set_title(
                f'{r["label"]} ID:{r["id"]}  $z_{{phot}}={r["z"]:.2f}$  '
                f'IR={r["ir"]:.1f}  N={len(r["vf"])} bands',
                fontsize=12, fontweight='bold', color=color[:7] if isinstance(color,str) else color
            )
            ax_a.legend(loc='upper right', fontsize=9, framealpha=0.85)
            ax_a.tick_params(labelsize=9)
            
            stats_str = (f'$\\Delta\\chi^2_\\nu = {r["dc2"]:+.1f}$\n'
                        f'$\\Delta$BIC = {r["dBIC"]:+.0f}\n'
                        f'$\\chi^2_{{null}} = {r["chisq_null"]:.1f}$\n'
                        f'$\\chi^2_{{Geff}} = {r["chisq_geff"]:.1f}$')
            bbox_props = dict(boxstyle='round,pad=0.35', facecolor='white', alpha=0.88, edgecolor=color[:7] if isinstance(color,str) else color)
            ax_a.text(0.04, 0.96, stats_str, transform=ax_a.transAxes, fontsize=9,
                     verticalalignment='top', fontfamily='monospace',
                     bbox=bbox_props, zorder=10)
            
            # Residuals panel
            ax_b = fig.add_subplot(gs[i, 1])
            
            resids_null = (r['fluxes'] - np.array([sed_model(l, r['Cb_null'], r['Cr_null'], 0.0, z_phot=r['z']) for l in r['lam_obs']])) / r['errors']
            resids_geff = (r['fluxes'] - np.array([sed_model(l, r['Cb_geff'], r['Cr_geff'], r['zd'], z_phot=r['z']) for l in r['lam_obs']])) / r['errors']
            
            x_pos = range(len(resids_null))
            width = 0.32
            
            ax_b.bar([p-width/2 for p in x_pos], resids_null, width,
                     color='#3498db', alpha=0.75, label='Null Resid.')
            ax_b.bar([p+width/2 for p in x_pos], resids_geff, width,
                     color='#e74c3c', alpha=0.85, label='G$_{eff}$ Resid.')
            
            ax_b.axhline(y=0, color='black', lw=0.8, ls='-')
            ax_b.axhline(y=-1, color='gray', lw=0.5, ls='--', alpha=0.5)
            ax_b.axhline(y=+1, color='gray', lw=0.5, ls='--', alpha=0.5)
            
            ax_b.set_xticks(x_pos)
            ax_b.set_xticklabels([f.replace('F','') for f in r['vf']], rotation=45, ha='right', fontsize=8)
            ax_b.set_ylabel('Residual ($\\sigma$)', fontsize=10)
            ax_b.set_title('Fit Residuals by Band', fontsize=11)
            ax_b.legend(loc='upper left', fontsize=8, framealpha=0.8)
            ax_b.tick_params(labelsize=8)
            
            rms_null = np.sqrt(np.mean(resids_null**2))
            rms_geff = np.sqrt(np.mean(resids_geff**2))
            ax_b.text(0.97, 0.96, f'RMS: {rms_null:.2f} -> {rms_geff:.2f}', transform=ax_b.transAxes,
                    fontsize=9, ha='right', va='top', bbox=dict(facecolor='white', alpha=0.8))
        
        suptitle = (
            'High-z G_eff Witnesses: Gold Positives vs Silver Comparison\n'
            f'AG={AG_BEST}, gamma={GAMMA_BEST}, beta={BETA_BEST}, zon={ZON_BEST}'
        )
        fig.suptitle(suptitle, fontsize=14, fontweight='bold', y=0.998)
        
        plt.savefig(os.path.join(WORK_DIR, 'GoldPositive_SED_Deep.png'), dpi=180, bbox_inches='tight')
        plt.savefig(os.path.join(WORK_DIR, 'GoldPositive_SED_Deep.pdf'), bbox_inches='tight')
        print(f"\nSaved: GoldPositive_SED_Deep.png & .pdf")

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Plotting error: {e}")

print("\nDone.")
