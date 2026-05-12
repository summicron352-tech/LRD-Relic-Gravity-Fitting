#!/usr/bin/env python3
"""
JADES DR5 Cross-Field Validation of G_eff(Σ) Signal — CORRECTED VERSION
========================================================================
Replicate the COSMOS 5.6σ Σ-color correlation in an independent JADES GOODS-N/S sample.

CRITICAL FIX (2026-05-07):
  The original cross-validation used WRONG definitions:
    OLD Σ = log(f444) - 2log(R)  ← f444 appears in BOTH Σ and color → math artifact
    OLD color = m444 - m150       ← magnitude difference, not flux ratio

  Now matching PathA exactly:
    CORRECT Σ = log10(L_bol_proxy) - 2*log10(R)  [or equivalently -2logR since Σ≈size]
    CORRECT color = log10(F444W/F150W)             [flux ratio, same as PathA]
    CORRECT method = Partial Spearman ctrl z_phot + log_Lbol_proxy

  L_bol proxy: JADES has no stellar mass/L_bol column. We estimate from broadband
  photometry using a multi-band flux sum calibrated against COSMOS L_bol values.
  KEY FINDING: In COSMOS, Σ correlates with -2logR at r=0.9997 because L_bol has
  only ~0.77 dex scatter vs ~1.16 dex from size. So -2logR is a near-perfect proxy.

Selection criteria (Rinaldi et al.):
  1. Photo-z > 4
  2. R_eff ≤ 0.06"
  3. F277W - F444W > 0.5 mag
  4. F115W - F200W > -0.5 mag
  5. Detection in key bands (FLAG=0)

Author: Xin Tan / 土鳖 (2026-05-07)
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, rankdata, norm
from numpy.linalg import lstsq
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import warnings
warnings.filterwarnings('ignore')
import os
import time

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR = '/Users/tanxin/Desktop/数据处理/重新整理投稿APJ'
DATA_DIR = '/Users/tanxin/Desktop/数据处理/04_JADES_DR5'
COSMOS_DIR = '/Users/tanxin/Desktop/数据处理/重新整理投稿APJ/源数据'
OUT_DIR = os.path.join(BASE_DIR, '源数据', 'JADES_CrossValidation')
os.makedirs(OUT_DIR, exist_ok=True)

CATALOGS = {
    'GOODS-S': os.path.join(DATA_DIR, 'GOODS-S',
                            'hlsp_jades_jwst_nircam_goods-s_photometry_v5.0_catalog.fits'),
    'GOODS-N': os.path.join(DATA_DIR, 'GOODS-N',
                            'hlsp_jades_jwst_nircam_goods-n_photometry_v5.0_catalog.fits'),
}

COSMO = FlatLambdaCDM(H0=70, Om0=0.3)

CUTS = {
    'z_min': 4.0,
    'reff_max_arcsec': 0.06,
    'f277w_f444w_min': 0.5,
    'f115w_f200w_min': -0.5,
}

COSMOS_BASELINE = {
    'N': 253,
    'spearman_rho': 0.341,
    'spearman_sigma': 5.6,
    'ks_D': 0.574,
    'ks_sigma': 6.2,
}


# =============================================================================
# Step 0: Calibrate L_bol proxy from COSMOS data
# =============================================================================
def calibrate_lbol_proxy():
    """
    Train an L_bol proxy from COSMOS Kokorev catalog.
    Since JADES has no stellar mass/L_bol, we estimate from multi-band fluxes.

    Key finding: Σ in COSMOS is dominated by -2log(R) (r=0.9997 with Σ).
    So -2logR alone is an excellent proxy. But for partial correlation
    (controlling log_Lbol), we need the actual L_bol estimate.
    """
    csv_path = os.path.join(COSMOS_DIR, 'Kokorev_LRDs_Full.csv')
    df = pd.read_csv(csv_path)

    mask = (df['lbol'].notna()) & (df['lbol'] > 0) & \
           df['f444w_flux'].notna() & (df['f444w_flux'] > 0) & \
           df['f356w_flux'].notna() & (df['f356w_flux'] > 0) & \
           df['f277w_flux'].notna() & (df['f277w_flux'] > 0) & \
           df['f150w_flux'].notna() & (df['f150w_flux'] > 0)
    d = df[mask].copy()

    # L_bol in Kokorev is already log10(L_bol/L_sun), range [43.3, 46.9]
    log_lbol_true = d['lbol'].values  # already log10

    # Multi-band flux features (log)
    log_f444 = np.log10(d['f444w_flux'].values)
    log_f356 = np.log10(d['f356w_flux'].values)
    log_f277 = np.log10(d['f277w_flux'].values)
    log_f150 = np.log10(d['f150w_flux'].values)

    # Feature matrix: individual bands + total flux
    X = np.column_stack([
        log_f150, log_f277, log_f356, log_f444,
        np.log10(d['f150w_flux'].values + d['f277w_flux'].values +
                 d['f356w_flux'].values + d['f444w_flux'].values)
    ])

    # Linear regression: log_Lbol = a0 + a1*x1 + ...
    X_aug = np.column_stack([np.ones(len(X)), X])
    coef, _, _, _ = lstsq(X_aug, log_lbol_true, rcond=None)

    # Predict and evaluate
    pred = X_aug @ coef
    resid = log_lbol_true - pred
    r_pred = np.corrcoef(log_lbol_true, pred)[0, 1]
    rmse = np.sqrt(np.mean(resid**2))

    print(f"\n{'='*60}")
    print(f"  Step 0: L_bol Proxy Calibration (COSMOS training)")
    print(f"{'='*60}")
    print(f"  Training N = {len(d)}")
    print(f"  Features: log(f150), log(f277), log(f356), log(f444), log(f_total)")
    print(f"  Prediction r = {r_pred:.4f}")
    print(f"  RMSE = {rmse:.3f} dex")
    print(f"  Coefficients: a0={coef[0]:.3f}, a150={coef[1]:.3f}, "
          f"a277={coef[2]:.3f}, a356={coef[3]:.3f}, a444={coef[4]:.3f}, "
          f"a_tot={coef[5]:.3f}")

    # Also verify: how much of Σ variance comes from size?
    neg2logR = -2 * np.log10(d['r_eff_50_phys'].values)
    sigma_old = np.log10(d['lbol'].values) - 2*np.log10(d['r_eff_50_phys'].values)
    r_size_sigma = np.corrcoef(neg2logR, sigma_old)[0, 1]
    print(f"\n  KEY: corr(-2logR, Σ) = {r_size_sigma:.4f}")
    print(f"  → Σ is essentially a SIZE proxy!")

    return {
        'coef': coef,
        'rmse': rmse,
        'r_pred': r_pred,
    }


def predict_log_lbol(flux_150, flux_277, flux_356, flux_444, coef):
    """Predict log10(L_bol/L_sun) from JADES fluxes using calibrated coefficients."""
    log_f150 = np.log10(np.maximum(flux_150, 1e-10))
    log_f277 = np.log10(np.maximum(flux_277, 1e-10))
    log_f356 = np.log10(np.maximum(flux_356, 1e-10))
    log_f444 = np.log10(np.maximum(flux_444, 1e-10))
    log_ftot = np.log10(np.maximum(flux_150 + flux_277 + flux_356 + flux_444, 1e-10))

    X = np.column_stack([
        np.ones(len(flux_150)),
        log_f150, log_f277, log_f356, log_f444, log_ftot
    ])
    return X @ coef


# =============================================================================
# Step 1: Load JADES catalog
# =============================================================================
def load_jades_catalog(field):
    """Load key HDUs from JADES DR5 FITS catalog."""
    fpath = CATALOGS[field]
    if not os.path.exists(fpath):
        print(f"  ❌ File not found: {fpath}")
        return None

    print(f"\n📂 Loading {field} ({os.path.getsize(fpath)/1e9:.2f} GB)")
    t0 = time.time()

    def to_native(recarr):
        result = {}
        for name in recarr.dtype.names:
            arr = recarr[name]
            if arr.ndim == 1:
                result[name] = arr.astype(arr.dtype.newbyteorder('='))
        return result

    with fits.open(fpath) as hdul:
        flag = pd.DataFrame(to_native(hdul['FLAG'].data))
        size = pd.DataFrame(to_native(hdul['SIZE'].data))
        kron = pd.DataFrame(to_native(hdul['KRON'].data))
        pz_data = to_native(hdul['PHOTOZ'].data)
        photoz = pd.DataFrame(pz_data)

    print(f"  Loaded in {time.time()-t0:.1f}s, {len(flag):,} sources")
    return {'field': field, 'flag': flag, 'size': size, 'kron': kron, 'photoz': photoz}


# =============================================================================
# Step 2: Apply Rinaldi et al. LRD selection
# =============================================================================
def select_lrd_candidates(data):
    """Apply Rinaldi et al. (2026) LRD selection criteria."""
    field = data['field']
    flag, size, kron, photoz = data['flag'], data['size'], data['kron'], data['photoz']
    N_total = len(flag)

    print(f"\n{'='*65}")
    print(f"  LRD Selection: {field} (N_total = {N_total:,})")
    print(f"{'='*65}")

    df = pd.DataFrame({
        'ID': flag['ID'].values,
        'RA': flag['RA'].values,
        'DEC': flag['DEC'].values,
        'z_phot': photoz['z_a'].values,
        'z_ml': photoz['z_ml'].values,
    })

    # Size: R_eff = sqrt(A * B)
    A = size['A'].values
    B = size['B'].values
    df['reff_arcsec'] = np.sqrt(A * B)

    # Photometry
    for band in ['F277W', 'F444W', 'F115W', 'F200W', 'F150W', 'F356W']:
        flux_col = f'{band}_KRON'
        err_col = f'{band}_KRON_e'
        if flux_col in kron.columns:
            flux = kron[flux_col].values.astype(float)
            df[f'{band}_flux'] = flux
            df[f'{band}_flag'] = flag[f'{band}_FLAG'].values
            with np.errstate(invalid='ignore', divide='ignore'):
                df[f'{band}_mag'] = np.where(flux > 0,
                    -2.5 * np.log10(flux) + 31.4, np.nan)
        else:
            print(f"  ⚠️  {flux_col} not found!")
            df[f'{band}_flux'] = np.nan
            df[f'{band}_mag'] = np.nan

    # Selection cuts
    cuts = {}

    mask_z = np.isfinite(df['z_phot']) & (df['z_phot'] > 0)
    cuts['valid_z'] = mask_z.sum()
    print(f"\n  Step 0: Valid photo-z          → N = {cuts['valid_z']:,}")

    mask_z4 = mask_z & (df['z_phot'] > CUTS['z_min'])
    cuts['z_gt_4'] = mask_z4.sum()
    print(f"  Step 1: z_phot > {CUTS['z_min']}              → N = {cuts['z_gt_4']:,}")

    mask_compact = mask_z4 & (df['reff_arcsec'] <= CUTS['reff_max_arcsec'])
    cuts['compact'] = mask_compact.sum()
    print(f"  Step 2: R_eff ≤ {CUTS['reff_max_arcsec']}\"             → N = {cuts['compact']:,}")

    mask_color = mask_compact & np.isfinite(df['F277W_mag']) & np.isfinite(df['F444W_mag'])
    df['color_277_444'] = df['F444W_mag'] - df['F277W_mag']
    mask_color = mask_color & (df['color_277_444'] > CUTS['f277w_f444w_min'])
    cuts['red_color'] = mask_color.sum()
    print(f"  Step 3: F277W-F444W > {CUTS['f277w_f444w_min']}      → N = {cuts['red_color']:,}")

    mask_bd = mask_color & np.isfinite(df['F115W_mag']) & np.isfinite(df['F200W_mag'])
    df['color_115_200'] = df['F200W_mag'] - df['F115W_mag']
    mask_bd = mask_bd & (df['color_115_200'] > CUTS['f115w_f200w_min'])
    cuts['bd_reject'] = mask_bd.sum()
    print(f"  Step 4: F115W-F200W > {CUTS['f115w_f200w_min']}    (BD cut) → N = {cuts['bd_reject']:,}")

    mask_snr = mask_bd.copy()
    for band in ['F277W', 'F444W', 'F150W', 'F356W']:
        if f'{band}_flag' in df.columns:
            mask_snr = mask_snr & (df[f'{band}_flag'] == 0)
    cuts['final'] = mask_snr.sum()
    print(f"  Step 5: Detection in key bands   → N = {cuts['final']:,}")

    mask_final = mask_snr & np.isfinite(df['F150W_mag']) & np.isfinite(df['F356W_mag'])
    cuts['with_colors'] = mask_final.sum()
    print(f"  Step 6: Valid F150W+F356W        → N = {cuts['with_colors']:,}")

    lrd = df[mask_final].copy()
    print(f"\n  ✅ {field} LRD candidates: {len(lrd):,}")
    return lrd, cuts


# =============================================================================
# Step 3: Compute Σ and color — CORRECT DEFINITIONS
# =============================================================================
def compute_sigma_and_color(df, lbol_coef):
    """
    Compute Σ and color matching PathA definitions EXACTLY.

    Σ = log10(L_bol_proxy) - 2*log10(R_eff_pc)   [L_bol from multi-band regression]
    color = log10(F444W / F150W)                   [flux ratio, same as PathA]
    """
    z = df['z_phot'].values
    da = COSMO.angular_diameter_distance(z).value  # Mpc
    reff_arcsec = df['reff_arcsec'].values
    reff_kpc = reff_arcsec * da * 1000 * np.pi / (180.0 * 3600.0)

    df['reff_kpc'] = reff_kpc

    # ── Σ Definition 1: L_bol proxy - 2log(R) ──
    # Estimate L_bol from multi-band fluxes
    log_lbol_pred = predict_log_lbol(
        df['F150W_flux'].values,
        df['F277W_flux'].values,
        df['F356W_flux'].values,
        df['F444W_flux'].values,
        lbol_coef
    )
    df['log_Lbol_proxy'] = log_lbol_pred

    mask_valid = (reff_kpc > 0) & np.isfinite(reff_kpc) & np.isfinite(log_lbol_pred)

    log_sigma_lbol = np.full(len(df), np.nan)
    log_sigma_lbol[mask_valid] = log_lbol_pred[mask_valid] - 2.0 * np.log10(reff_kpc[mask_valid])
    df['logSigma'] = log_sigma_lbol

    # ── Σ Definition 2: Pure size proxy -2log(R) ──
    # In COSMOS, this is r=0.9997 correlated with the "full" Σ
    # Use this as independent cross-check
    neg2logR = np.full(len(df), np.nan)
    neg2logR[mask_valid] = -2.0 * np.log10(reff_kpc[mask_valid])
    df['neg2logR'] = neg2logR

    # ── Σ Definition 3 (WRONG, for comparison): f444-based ──
    f444 = df['F444W_flux'].values
    log_sigma_wrong = np.full(len(df), np.nan)
    m_wrong = mask_valid & (f444 > 0)
    log_sigma_wrong[m_wrong] = np.log10(f444[m_wrong]) - 2.0 * np.log10(reff_kpc[m_wrong])
    df['logSigma_WRONG_f444'] = log_sigma_wrong

    # Center Σ
    med = np.nanmedian(log_sigma_lbol)
    df['logSigma_centered'] = log_sigma_lbol - med

    # ── Color: log(F444W/F150W) — same as PathA ──
    with np.errstate(invalid='ignore', divide='ignore'):
        df['color_444_150'] = np.log10(df['F444W_flux'].values / df['F150W_flux'].values)
        df['color_444_356'] = np.log10(df['F444W_flux'].values / df['F356W_flux'].values)

    # Also compute magnitude-based color for comparison with old (wrong) version
    df['mag_color_444_150'] = df['F444W_mag'] - df['F150W_mag']

    N_valid = np.sum(mask_valid & np.isfinite(df['color_444_150']))
    print(f"  N with valid Σ + color: {N_valid}")
    print(f"  logΣ range: [{np.nanmin(log_sigma_lbol):.2f}, {np.nanmax(log_sigma_lbol):.2f}]")
    print(f"  color range: [{np.nanmin(df['color_444_150']):.2f}, {np.nanmax(df['color_444_150']):.2f}]")

    return df


# =============================================================================
# Step 4: Partial Spearman (matching PathA)
# =============================================================================
def partial_spearman(x, y, ctrl_cols, data):
    """Spearman partial correlation: regress controls from ranks, then correlate residuals."""
    rx = rankdata(x)
    ry = rankdata(y)
    X_ctrl = np.column_stack([np.ones(len(data))] + [data[c].values for c in ctrl_cols])

    resid_x = rx - X_ctrl @ lstsq(X_ctrl, rx, rcond=None)[0]
    resid_y = ry - X_ctrl @ lstsq(X_ctrl, ry, rcond=None)[0]

    return spearmanr(resid_x, resid_y)


# =============================================================================
# Step 5: Statistical tests
# =============================================================================
def run_statistics(df, label=""):
    """
    Full statistical analysis with CORRECT definitions.
    Tests multiple Σ definitions to show the impact of the bug fix.
    """
    results = {}

    # ── Test A: CORRECT Σ (L_bol proxy - 2logR) vs CORRECT color (log FR) ──
    print(f"\n  {'─'*50}")
    print(f"  TEST A: CORRECT Σ (L_bol-2logR) vs log(F444W/F150W)")
    print(f"  {'─'*50}")

    mask = np.isfinite(df['logSigma']) & np.isfinite(df['color_444_150'])
    d = df[mask].copy()
    N = len(d)
    results['N'] = N
    results['Sigma_def'] = 'L_bol_proxy - 2logR'

    if N < 10:
        print(f"  ⚠️  Too few sources ({N})")
        return results

    # Raw Spearman
    rho_raw, p_raw = spearmanr(d['logSigma'], d['color_444_150'])
    sig_raw = abs(norm.ppf(p_raw/2)) if p_raw > 0 else np.inf
    results['A_raw_rho'] = rho_raw
    results['A_raw_sigma'] = sig_raw
    print(f"  Raw Spearman:         ρ = {rho_raw:+.4f} ({sig_raw:.1f}σ)")

    # Partial Spearman ctrl z + L_bol
    rho_p, p_p = partial_spearman(
        d['logSigma'].values, d['color_444_150'].values,
        ['z_phot', 'log_Lbol_proxy'], d
    )
    sig_p = abs(norm.ppf(p_p/2)) if p_p > 0 else np.inf
    results['A_partial_rho'] = rho_p
    results['A_partial_sigma'] = sig_p
    print(f"  Partial (ctrl z+L_bol): ρ_p = {rho_p:+.4f} ({sig_p:.1f}σ)")

    # KS test
    med_s = np.median(d['logSigma'])
    high = d[d['logSigma'] >= med_s]['color_444_150'].values
    low = d[d['logSigma'] < med_s]['color_444_150'].values
    D_ks, p_ks = stats.ks_2samp(high, low)
    sig_ks = abs(norm.ppf(p_ks/2)) if p_ks > 0 else np.inf
    results['ks_D'] = D_ks
    results['ks_sigma'] = sig_ks
    results['delta_color'] = np.mean(high) - np.mean(low)
    print(f"  KS test:              D = {D_ks:.4f} ({sig_ks:.1f}σ)")
    print(f"  Δcolor (high-Σ-low):  {results['delta_color']:+.3f} dex")

    # ── Test B: Pure size proxy -2log(R) vs color ──
    print(f"\n  {'─'*50}")
    print(f"  TEST B: Pure size -2log(R) vs log(F444W/F150W)")
    print(f"  {'─'*50}")

    mask_b = np.isfinite(df['neg2logR']) & np.isfinite(df['color_444_150'])
    db = df[mask_b]
    rho_b, p_b = spearmanr(db['neg2logR'], db['color_444_150'])
    sig_b = abs(norm.ppf(p_b/2)) if p_b > 0 else np.inf
    results['B_size_rho'] = rho_b
    results['B_size_sigma'] = sig_b
    print(f"  Raw Spearman:  ρ = {rho_b:+.4f} ({sig_b:.1f}σ)")

    rho_bp, p_bp = partial_spearman(
        db['neg2logR'].values, db['color_444_150'].values,
        ['z_phot', 'log_Lbol_proxy'], db
    )
    sig_bp = abs(norm.ppf(p_bp/2)) if p_bp > 0 else np.inf
    results['B_size_partial_rho'] = rho_bp
    results['B_size_partial_sigma'] = sig_bp
    print(f"  Partial (ctrl z+L): ρ_p = {rho_bp:+.4f} ({sig_bp:.1f}σ)")

    # ── Test C: WRONG Σ (f444-based) vs WRONG color (mag diff) — for comparison ──
    print(f"\n  {'─'*50}")
    print(f"  TEST C: WRONG Σ (f444-2logR) vs mag_color (old definitions)")
    print(f"  {'─'*50}")

    mask_c = np.isfinite(df['logSigma_WRONG_f444']) & np.isfinite(df['mag_color_444_150'])
    dc = df[mask_c]
    rho_c, p_c = spearmanr(dc['logSigma_WRONG_f444'], dc['mag_color_444_150'])
    sig_c = abs(norm.ppf(p_c/2)) if p_c > 0 else np.inf
    results['C_wrong_rho'] = rho_c
    results['C_wrong_sigma'] = sig_c
    print(f"  Raw Spearman:  ρ = {rho_c:+.4f} ({sig_c:.1f}σ)")
    print(f"  (This was the WRONG version that gave negative ρ!)")

    # ── Comparison with COSMOS ──
    print(f"\n  {'='*50}")
    print(f"  📐 COMPARISON WITH COSMOS BASELINE ({label})")
    print(f"  {'='*50}")
    cb = COSMOS_BASELINE
    print(f"  {'':35s} {'COSMOS':>10s} {'GOODS':>10s}")
    print(f"  {'N':35s} {cb['N']:>10d} {N:>10d}")
    print(f"  {'Partial ρ (ctrl z+L_bol)':35s} {cb['spearman_rho']:>+10.3f} {rho_p:>+10.3f}")
    print(f"  {'Partial σ':35s} {cb['spearman_sigma']:>10.1f} {sig_p:>10.1f}")
    print(f"  {'KS D':35s} {cb['ks_D']:>10.3f} {D_ks:>10.3f}")
    print(f"  {'KS σ':35s} {cb['ks_sigma']:>10.1f} {sig_ks:>10.1f}")

    return results


# =============================================================================
# Main pipeline
# =============================================================================
def main():
    print("=" * 70)
    print("  JADES DR5 Cross-Validation: G_eff(Σ) Signal — CORRECTED")
    print("  Σ = log10(L_bol_proxy) - 2log10(R)  |  color = log(F444/F150)")
    print("=" * 70)

    # ── Step 0: Calibrate L_bol proxy ──
    cal = calibrate_lbol_proxy()
    lbol_coef = cal['coef']

    all_lrds = []
    all_cuts = {}

    for field in ['GOODS-S', 'GOODS-N']:
        data = load_jades_catalog(field)
        if data is None:
            continue

        lrd, cuts = select_lrd_candidates(data)
        all_cuts[field] = cuts

        if len(lrd) == 0:
            print(f"  ⚠️  No LRD candidates in {field}!")
            continue

        print(f"\n  📐 Computing Σ and color for {field}...")
        lrd = compute_sigma_and_color(lrd, lbol_coef)
        lrd['field'] = field
        all_lrds.append(lrd)

    if len(all_lrds) == 0:
        print("\n  ❌ No LRD candidates!")
        return

    df = pd.concat(all_lrds, ignore_index=True)
    N = len(df)

    print(f"\n{'='*70}")
    print(f"  🎯 Combined: {N:,} LRD candidates")
    print(f"  (GOODS-S: {(df['field']=='GOODS-S').sum()}, "
          f"GOODS-N: {(df['field']=='GOODS-N').sum()})")
    print(f"{'='*70}")

    # Redshift distribution
    z = df['z_phot'].values
    print(f"\n  Redshift distribution:")
    for zlo, zhi in [(4,5),(5,6),(6,7),(7,8),(8,9),(9,12)]:
        n = np.sum((z >= zlo) & (z < zhi))
        print(f"    z ∈ [{zlo},{zhi}): {n:>5d}  ({100*n/N:.1f}%)")

    # ── Combined statistics ──
    print(f"\n{'='*70}")
    print(f"  STATISTICAL ANALYSIS: Combined GOODS-N+S")
    print(f"{'='*70}")
    stats_comb = run_statistics(df, label="GOODS Combined")

    # ── Per-field statistics ──
    print(f"\n{'='*70}")
    print(f"  PER-FIELD ANALYSIS")
    print(f"{'='*70}")
    per_field = {}
    for field in ['GOODS-S', 'GOODS-N']:
        df_f = df[df['field'] == field]
        if len(df_f) > 5:
            per_field[field] = run_statistics(df_f, label=field)

    # ── Save results ──
    print(f"\n{'='*70}")
    print(f"  💾 Saving results...")
    print(f"{'='*70}")

    out_csv = os.path.join(OUT_DIR, 'JADES_LRD_candidates_v2.csv')
    df.to_csv(out_csv, index=False)
    print(f"  → {out_csv} ({len(df)} rows)")

    # Summary CSV
    summary = pd.Series({
        'N_total': N,
        'N_GOODS_S': (df['field']=='GOODS-S').sum(),
        'N_GOODS_N': (df['field']=='GOODS-N').sum(),
        'Sigma_definition': 'log10(L_bol_proxy) - 2*log10(R)',
        'Color_definition': 'log10(F444W/F150W)',
        'Method': 'Partial Spearman ctrl z+Lbol',
        'partial_rho': stats_comb.get('A_partial_rho'),
        'partial_sigma': stats_comb.get('A_partial_sigma'),
        'raw_rho': stats_comb.get('A_raw_rho'),
        'raw_sigma': stats_comb.get('A_raw_sigma'),
        'ks_D': stats_comb.get('ks_D'),
        'ks_sigma': stats_comb.get('ks_sigma'),
        'delta_color': stats_comb.get('delta_color'),
        'size_proxy_rho': stats_comb.get('B_size_rho'),
        'size_proxy_sigma': stats_comb.get('B_size_sigma'),
        'WRONG_old_rho': stats_comb.get('C_wrong_rho'),
        'WRONG_old_sigma': stats_comb.get('C_wrong_sigma'),
        'COSMOS_rho': COSMOS_BASELINE['spearman_rho'],
        'COSMOS_sigma': COSMOS_BASELINE['spearman_sigma'],
        'lbol_proxy_RMSE': cal['rmse'],
    })
    out_sum = os.path.join(OUT_DIR, 'cross_validation_summary_v2.csv')
    summary.to_csv(out_sum, header=['value'])
    print(f"  → {out_sum}")

    cuts_df = pd.DataFrame(all_cuts)
    out_cuts = os.path.join(OUT_DIR, 'selection_cuts_v2.csv')
    cuts_df.to_csv(out_cuts)
    print(f"  → {out_cuts}")

    # ── Verdict ──
    print(f"\n{'='*70}")
    print(f"  🏁 FINAL VERDICT")
    print(f"{'='*70}")

    sig_p = stats_comb.get('A_partial_sigma', 0)
    sig_ks = stats_comb.get('ks_sigma', 0)
    sig_raw = stats_comb.get('A_raw_sigma', 0)

    print(f"\n  CORRECT definitions (L_bol-2logR vs log FR):")
    print(f"    Raw Spearman:       {sig_raw:.1f}σ")
    print(f"    Partial (ctrl z+L): {sig_p:.1f}σ")
    print(f"    KS test:            {sig_ks:.1f}σ")

    rho_wrong = stats_comb.get('C_wrong_rho', 0)
    sig_wrong = stats_comb.get('C_wrong_sigma', 0)
    print(f"\n  WRONG definitions (f444-2logR vs mag diff) for reference:")
    print(f"    Raw Spearman: ρ = {rho_wrong:+.3f} ({sig_wrong:.1f}σ)")

    if sig_p >= 4.0:
        verdict = "✅ STRONG CONFIRMATION"
    elif sig_p >= 3.0:
        verdict = "✅ CONFIRMED (consistent with smaller GOODS area)"
    elif sig_p >= 2.0:
        verdict = "⚠️  MARGINAL"
    else:
        verdict = "❌ NOT DETECTED with correct definitions"

    print(f"\n  {verdict}")
    print(f"\n  📝 Note: COSMOS baseline = {COSMOS_BASELINE['spearman_sigma']}σ (ρ={COSMOS_BASELINE['spearman_rho']})")

    if sig_p > 0:
        ratio = COSMOS_BASELINE['spearman_sigma'] / sig_p
        print(f"  GOODS signal is {ratio:.1f}× weaker/stronger than COSMOS")

    print(f"\n  DONE! 🐢")


if __name__ == '__main__':
    main()
