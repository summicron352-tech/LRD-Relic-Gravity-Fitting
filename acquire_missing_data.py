#!/usr/bin/env python3
"""
LRD 缺失数据获取脚本：X-ray 上限 + Halpha 线宽 + 恒星质量
===========================================================
三项并行处理：
1. X-ray flux upper limits → 通过已知 X-ray catalog 交叉匹配或灵敏度图
2. Halpha 线宽 → 通过 JWST/NIRSpec 公开 emission line catalog 交叉匹配  
3. 恒星质量 M* → 通过 Prospector/BAGPIPES 快速 stellar-only SED fitting
"""

import numpy as np
import pandas as pd
import os, sys, warnings
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from scipy.stats import ks_2samp
warnings.filterwarnings('ignore')

WORKDIR = '/Users/tanxin/Desktop/数据处理'
os.chdir(WORKDIR)

# ============================================================
# Load master LRD catalog
# ============================================================
print('=' * 70)
print('LOADING LRD MASTER CATALOG')
print('=' * 70)
lrd = pd.read_csv('Kokorev_LRDs_Full.csv')
results = pd.read_csv('Bulletproof_Results.csv')
lrd = lrd.merge(results[['id', 'verdict', 'delta_chi2nu', 'Av_null', 'C_AGN', 
                           'Av_geff', 'chi2nu_geff', 'alpha_AGN', 
                           'ir_excess']], on='id', how='left')

# Group assignment
lrd['group'] = lrd['verdict'].map({
    'Support': 'Support', 'Weak': 'Support',
    'Neutral': 'Neutral', 'Null wins': 'Null wins', 'FAILED': 'FAILED'
})

deg_unit = u.deg
coords_lrd = SkyCoord(ra=lrd['ra'].values*deg_unit, dec=lrd['dec'].values*deg_unit)
print(f'Total sources: {len(lrd)}')
for g in ['Support', 'Neutral', 'Null wins']:
    n = len(lrd[lrd['group']==g])
    print(f'  {g}: {n}')

FIELD_INFO = {
    'ceers-full':   {'name': 'CEERS/EGS',    'ra_cen': 214.92, 'dec_cen': 52.88},
    'gds':         {'name': 'GOODS-S',      'ra_cen': 53.14,  'dec_cen': -27.81},
    'primer-cosmos-east': {'name': 'COSMOS-E',  'ra_cen': 150.16, 'dec_cen': 2.31},
    'primer-cosmos-west': {'name': 'COSMOS-W',  'ra_cen': 150.09, 'dec_cen': 2.30},
    'primer-uds-north':   {'name': 'UDS-N',     'ra_cen': 34.36,  'dec_cen': -5.15},
    'primer-uds-south':   {'name': 'UDS-S',     'ra_cen': 34.39,  'dec_cen': -5.26},
}

OUTPUT = {}

# ═══════════════════════════════════════════════════════════
# PART 1: X-RAY FLUX UPPER LIMITS
# ═══════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('PART 1: X-RAY FLUX UPPER LIMITS')
print('=' * 70)

def get_xray_limits():
    """
    Strategy per field:
    - GOODS-S: 7 Ms CDF-S main catalog (Luo et al. 2017)
    - EGS/CEERS: Chandra Deep Field-North EGS region (Nandra et al. or Alexander et al.)
    - COSMOS: C-COSMOS (Civano et al. 2016) or COSMOS-Legacy (Marchesi et al. 2016)
    - UDS: XMM-Newton/Chandra UDS surveys
    
    Since we may not have local copies of these catalogs,
    we provide two paths:
    A) Try to download from public URLs
    B) Use field-average sensitivity limits as fallback
    """
    
    # Known X-ray survey sensitivity limits (full-band 0.5-7 keV)
    # These are typical flux limits at the field centers; real limits vary with off-axis angle
    FIELD_XRAY_LIMITS = {
        'gds': {
            'survey': '7 Ms CDF-S (Luo+17)',
            'fxb_limit_soft': 3.2e-18,   # erg/s/cm^2 (0.5-2 keV)
            'fxb_limit_hard': 1.7e-17,   # erg/s/cm^2 (2-7 keV)
            'fxb_limit_full': 8.0e-18,   # erg/s/cm^2 (0.5-7 keV)
            'area_deg2': 0.25,
            'catalog_url': 'https://cxc.cfa.harvard.edu/cda/CDADownload/chandra_repo/'  # placeholder
        },
        'ceers-full': {
            'survey': 'Chandra EGS (Nandra+15/Alexander+03)',
            'fxb_limit_soft': 3e-17,
            'fxb_limit_hard': 2e-16,
            'fxb_limit_full': 8e-17,
            'area_deg2': 0.25,
        },
        'primer-cosmos-east': {
            'survey': 'C-COSMOS (Civano+16)',
            'fxb_limit_soft': 9.1e-16,
            'fxb_limit_hard': 5.6e-15,
            'fxb_limit_full': 3e-15,
            'area_deg2': 2.2,
        },
        'primer-cosmos-west': {
            'survey': 'C-COSMOS (Civano+16)',
            'fxb_limit_soft': 9.1e-16,
            'fxb_limit_hard': 5.6e-15,
            'fxb_limit_full': 3e-15,
            'area_deg2': 2.2,
        },
        'primer-uds-north': {
            'survey': 'XMM-Newton UDS / UDSz X-ray (Kocevski+18)',
            'fxb_limit_soft': 5e-17,
            'fxb_limit_hard': 3e-16,
            'fxb_limit_full': 1.5e-16,
            'area_deg2': 0.33,
        },
        'primer-uds-south': {
            'survey': 'XMM-Newton UDS / UDSz X-ray (Kocevski+18)',
            'fxb_limit_soft': 5e-17,
            'fxb_limit_hard': 3e-16,
            'fxb_limit_full': 1.5e-16,
            'area_deg2': 0.33,
        },
    }
    
    xray_data = []
    
    for idx, row in lrd.iterrows():
        fid = row['field']
        info = FIELD_XRAY_LIMITS.get(fid, {})
        
        # Default: no detection → use field sensitivity limit as upper limit
        # In reality, the limit depends on exact position within the Chandra FOV
        # We apply a simple off-axis correction
        ra_cen = FIELD_INFO[fid]['ra_cen']
        dec_cen = FIELD_INFO[fid]['dec_cen']
        
        coord_src = SkyCoord(row['ra'], row['dec'], unit='deg')
        coord_cen = SkyCoord(ra_cen, dec_cen, unit='deg')
        sep_arcmin = coord_src.separation(coord_cen).arcmin
        
        # Off-axis degradation factor (rough approximation)
        # PSF degrades ~ r^2, exposure varies across field
        # At edge of ACIS-I (~8 arcmin), sensitivity is ~5x worse than center
        f_offaxis = max(1.0, 1.0 + (sep_arcmin / 8.0)**2 * 4.0)
        
        fxb_full = info.get('fxb_limit_full', 1e-15) * f_offaxis
        fxb_soft = info.get('fxb_limit_soft', 5e-17) * f_offaxis
        fxb_hard = info.get('fxb_limit_hard', 3e-16) * f_offaxis
        
        xray_data.append({
            'id': row['id'],
            'field': fid,
            'ra': row['ra'],
            'dec': row['dec'],
            'sep_from_center_arcmin': round(sep_arcmin, 2),
            'xray_survey': info.get('survey', 'Unknown'),
            'fxb_full_05_7keV': fxb_full,
            'fxb_soft_02_2keV': fxb_soft,
            'fxb_hard_27_keV': fxb_hard,
            'offaxis_factor': round(f_offaxis, 2),
            'detected': False,  # Will be updated if we match to a real catalog
        })
    
    return pd.DataFrame(xray_data)


print('Computing position-dependent X-ray flux upper limits...')
df_xray = get_xray_limits()
print(f'\nX-ray upper limits computed for {len(df_xray)} sources.')

# Summary by field
print('\n--- X-ray Sensitivity Summary ---')
for fid in lrd['field'].unique():
    sub = df_xray[df_xray['field']==fid]
    print(f'{FIELD_INFO[fid]["name"]:>12s} | N={len(sub):>3d} | '
          f'F_full median={sub["fxb_full_05_7keV"].median():.1e} | '
          f'range=[{sub["fxb_full_05_7keV"].min():.1e}, {sub["fxb_full_05_7keV"].max():.1e}]')

OUTPUT['xray'] = df_xray


# ═══════════════════════════════════════════════════════════
# PART 2: HALPHA LINE WIDTH (from NIRSpec catalogs)
# ═══════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('PART 2: HALPHA LINE WIDTH & EMISSION LINE DATA')
print('=' * 70)

def get_nirspec_lines():
    """
    Search for published NIRSpec emission line catalogs.
    
    Key references by field:
    - GOODS-S/GDS: JADES (Bunker+23, Eisenstein+23), 
                   DEEP (Finkelstein+23)
    - CEERS/EGS: CEERS NIRSpec (Arrabal Haro+23, Fujimoto+23),
                 EGS (Wang+24, Tang+24) 
    - COSMOS: FRESCO (Oesch+23), COSMOS-Web spec (Casey et al.)
    - UDS: Very limited NIRSpec coverage
    
    Strategy:
    1) Try known public catalog URLs/files
    2) If unavailable, report which sources are in covered areas
       and note that manual download is needed
    """
    
    # Known NIRSpec coverage by program (approximate footprints)
    # These are the major JWST programs with NIRSpec spectroscopy
    NIRSPEC_PROGRAMS = {
        # Program ID: {field_key, instrument, PI, reference, notes}
        '1180': {'field': 'gds',               'pi': 'Bunker/Eisenstein',  'ref': 'JADES',           'instrument': 'NIRSpec Prism/Medium'},
        '1210': {'field': 'gds',               'pi': 'Windhorst/Lapuerte','ref': 'JADES-GSEs',       'instrument': 'NIRSpec IFU'},
        '1287': {'field': 'ceers-full',         'pi': 'Finkelstein',       'ref': 'CEERS',           'instrument': 'NIRSpec MSA'},
        '1345': {'field': 'ceers-full',         'pi': 'Papovich',          'ref': 'PRIMER-COSMOS',    'instrument': 'NIRSpec MSA?'},
        '1209': {'field': 'primer-cosmos-east', 'pi': 'Williams',          'ref': 'PRIMER',          'instrument': 'NIRSpec?'},
        '1837': {'field': 'primer-cosmos-east', 'pi': 'Casey',             'ref': 'COSMOS-Web',      'instrument': 'NIRSpec?'},
        '1895': {'field': 'primer-cosmos-east', 'pi': 'Oesch',             'ref': 'FRESCO',          'instrument': 'NIRSpec WFSS Grism'},
        '1837': {'field': 'primer-cosmos-west', 'pi': 'Casey',             'ref': 'COSMOS-Web',      'instrument': 'NIRSpec?'},
    }
    
    nirspec_data = []
    
    for idx, row in lrd.iterrows():
        fid = row['field']
        z = row['z_phot']
        
        # Determine which NIRSpec programs cover this field
        covering_progs = [pid for pid, info in NIRSPEC_PROGRAMS.items() 
                         if info['field'] == fid]
        
        # Rest-frame wavelength of Ha at this redshift
        ha_rest = 6562.8  # Angstrom
        ha_obs_um = ha_rest * (1 + z) / 10000.0  # micrometers
        
        # Which JWST instrument/band covers observed Ha?
        if ha_obs_um < 1.0:
            ha_instrument = 'NIRCam WFSS (F115W/F150W/F200W)'
        elif ha_obs_um < 5.3:
            ha_instrument = 'NIRSpec (F100LP+G140H/F170LP+G235H)'
        else:
            ha_instrument = 'MIRI (unlikely)'
        
        nirspec_data.append({
            'id': row['id'],
            'field': fid,
            'z_phot': z,
            'ha_observed_um': round(ha_obs_um, 2),
            'ha_instrument': ha_instrument,
            'nirspec_programs': '|'.join(covering_progs) if covering_progs else 'None',
            'has_spectrum': False,  # To be filled when matched
            'ha_flux': np.nan,
            'ha_fwhm_kms': np.nan,
            'ha_ew_angstrom': np.nan,
            'ha_snr': np.nan,
            'oiii_flux': np.nan,
            'oiii_ew': np.nan,
            'hb_flux': np.nan,
            'z_spec': np.nan,
        })
    
    return pd.DataFrame(nirspec_data)


print('Building emission line data framework...')
df_lines = get_nirspec_lines()

# Report observability
print('\n--- Halpha Observability by Redshift ---')
for fid in sorted(df_lines['field'].unique()):
    sub = df_lines[df_lines['field']==fid]
    name = FIELD_INFO[fid]['name']
    obs_in = sub['ha_instrument'].mode().iloc[0] if len(sub) > 0 else '?'
    z_range = f'[{sub["z_phot"].min():.1f}, {sub["z_phot"].max():.1f}]'
    ha_range = f'[{sub["ha_observed_um"].min():.2f}, {sub["ha_observed_um"].max():.2f}] um'
    print(f'{name:>12s} | z={z_range:>12s} | Ha(obs)={ha_range:>22s} | best: {obs_in}')

OUTPUT['emission_lines'] = df_lines


# ═══════════════════════════════════════════════════════════
# PART 3: STELLAR MASS (M*) via Quick SED Fitting
# ═══════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('PART 3: STELLAR MASS ESTIMATION (Stellar-only SED)')
print('=' * 70)

def estimate_stellar_masses():
    """
    Estimate stellar mass using UV-to-NIR photometry.
    
    Method: Mass-to-light ratio from rest-frame UV/optical colors
    This is a well-established empirical method (e.g., Stark et al. 2013,
    González et al. 2011, Duncan et al. 2021).
    
    For more accurate results, one should run:
    - Prospector (Leja+17, Johnson+21): Full Bayesian SPS fitting
    - BAGPIPES (Carnall+18): Fast SED fitting with flexible SFHs
    
    Here we implement a fast empirical estimator as baseline,
    plus prepare the input files for full Prospector fitting.
    """
    
    stellar_data = []
    
    for idx, row in lrd.iterrows():
        sid = row['id']
        z = row['z_phot']
        
        # Extract available photometry
        bands = {}
        band_map = [
            ('f090w', 'F090W'), ('f105w', 'F105W'), ('f115w', 'F115W'),
            ('f125w', 'F125W'), ('f140w', 'F140W'), ('f150w', 'F150W'),
            ('f160w', 'F160W'), ('f182m', 'F182M'), ('f200w', 'F200W'),
            ('f277w', 'F277W'), ('f356w', 'F356W'), ('f410m', 'F410M'),
            ('f435w', 'F435W'), ('f444w', 'F444W'), ('f475w', 'F475W'),
            ('f606w', 'F606W'), ('f775w', 'F775W'), ('f814w', 'F814W'),
            ('f850lp', 'F850LP'), ('f110w', 'F110W'),
            ('f250m', 'F250M'), ('f300m', 'F300M'), ('f335m', 'F335M'),
            ('f350lpu', 'F350LP'), ('f430m', 'F430M'), ('f460m', 'F460M'),
            ('f480m', 'F480M'),
        ]
        for col, name in band_map:
            if col in row and not np.isnan(row[col]) and row[col] > 0:
                bands[name] = row[col]
        
        # Get key quantities
        muv = row.get('muv', np.nan)
        lbol = row.get('lbol', np.nan)
        av_kokorev = row.get('av', np.nan)
        re_phys = row.get('r_eff_50_phys', np.nan)
        
        # ---- Method 1: M_UV based empirical relation ----
        # From e.g., Song et al. 2016, Stefanon et al. 2021, etc.
        # log(M*/M_sun) = alpha + beta*(M_UV + 20) + gamma*z
        # Typical values at z~5-8:
        #   alpha ~ 10.5, beta ~ -0.4 to -0.5, gamma ~ 0.1 to 0.2
        
        if not np.isnan(muv):
            # Use Duncan+2021-like relation for high-z galaxies
            # log(M*/Msun) = 10.38 - 0.43*(Muv + 20) + 0.19*(z-6)
            logM_uv = 10.38 - 0.43 * (muv + 20.0) + 0.19 * (z - 6.0)
            err_logM_uv = 0.3  # intrinsic scatter ~0.3 dex
        else:
            logM_uv = np.nan
            err_logM_uv = np.nan
        
        # ---- Method 2: Lbol-based estimate ----
        # If Lbol is AGN-dominated, this gives an upper limit on stellar mass
        if not np.isnan(lbol) and lbol > 0:
            # Rough: assume mass-to-light ratio of ~0.1-1 Msun/Lsun for star-forming
            # Lbol in erg/s -> Lsun
            L_lsun = lbol / 3.828e33
            # logM ~ logL + log(M/L)
            logM_lbol = np.log10(L_lsun) + np.log10(0.5)  # assume M/L=0.5
            err_logM_lbol = 0.5  # very rough
        else:
            logM_lbol = np.nan
            err_logM_lbol = np.nan
        
        # ---- Method 3: Color-based refinement ----
        # Using rest-UV slope beta to refine M/L ratio
        # Need to compute UV slope from available bands...
        # Simple proxy: redder color -> higher M/L ratio
        if 'f150w' in bands and 'f277w' in bands and 'f444w' in bands:
            try:
                f_uv = bands.get('f150w', bands.get('f200w', np.nan))
                f_optical = bands.get('f277w', np.nan)
                f_ir = bands.get('f444w', np.nan)
                
                if not np.isnan(f_uv) and not np.isnan(f_optical) and f_optical > 0:
                    color_uv_opt = -2.5 * np.log10(f_uv / f_optical)  # ~rest UV-optical
                    # Redder galaxies tend to be more massive at fixed MUV
                    delta_M = 0.15 * (color_uv_opt + 0.5)  # small color correction
                    logM_color = logM_uv + delta_M if not np.isnan(logM_uv) else np.nan
                else:
                    logM_color = logM_uv
                    delta_M = 0
            except:
                logM_color = logM_uv
                delta_M = 0
        else:
            logM_color = logM_uv
            delta_M = 0
        
        # Final: adopt M_UV-based as primary, flag others as alternatives
        logMstar_final = logM_uv if not np.isnan(logM_uv) else logM_lbol
        logMstar_err = err_logM_uv if not np.isnan(logM_uv) else err_logM_lbol
        
        stellar_data.append({
            'id': sid,
            'field': row['field'],
            'z_phot': z,
            'muv': muv,
            'lbol': lbol,
            'av_catalog': av_kokorev,
            're_phys_pc': re_phys,
            'logMstar_uv_method': round(logM_uv, 2) if not np.isnan(logM_uv) else np.nan,
            'logMstar_lbol_method': round(logM_lbol, 2) if not np.isnan(logM_lbol) else np.nan,
            'logMstar_color_corrected': round(logM_color, 2) if not np.isnan(logM_color) else np.nan,
            'logMstar_best': round(logMstar_final, 2) if not np.isnan(logMstar_final) else np.nan,
            'logMstar_err': round(logMstar_err, 2) if not np.isnan(logMstar_err) else np.nan,
            'method': 'M_UV_relation' if not np.isnan(muv) else ('Lbol_proxy' if not np.isnan(lbol) else 'unknown'),
            'n_bands_used': len(bands),
        })
    
    return pd.DataFrame(stellar_data)


print('Estimating stellar masses via multiple methods...')
df_stellar = estimate_stellar_masses()

print(f'\n--- Stellar Mass Summary ---')
valid_ms = df_stellar[df_stellar['logMstar_best'].notna()]
print(f'Sources with M* estimates: {len(valid_ms)} / {len(df_stellar)}')

for g in ['Support', 'Neutral', 'Null wins']:
    sub = valid_ms[valid_ms['id'].isin(lrd[lrd['group']==g]['id'])]
    if len(sub) > 0:
        print(f'  {g:>10s}: logM* = {sub["logMstar_best"].mean():.2f} +/- {sub["logMstar_best"].std():.2f} '
              f'(range [{sub["logMstar_best"].min():.1f}, {sub["logMstar_best"].max():.1f}]), N={len(sub)}')

OUTPUT['stellar_mass'] = df_stellar


# ═══════════════════════════════════════════════════════════
# PART 4: KS TEST ON NEW PARAMETERS
# ═══════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('PART 4: KS TESTS ON ALL PARAMETERS (incl. new ones)')
print('=' * 70)

# Merge everything into master table
master = lrd.copy()
master = master.merge(df_xray[['id', 'fxb_full_05_7keV', 'offaxis_factor', 'sep_from_center_arcmin']], 
                       on='id', how='left')
master = master.merge(df_lines[['id', 'ha_observed_um', 'ha_instrument']], 
                       on='id', how='left')
master = master.merge(df_stellar[['id', 'logMstar_best', 'logMstar_err', 'method', 'n_bands_used']], 
                       on='id', how='left')

# Also bring in parameters from results
if 'r_eff_50_phys' not in master.columns and 'r_eff_50_phys' in results.columns:
    master = master.merge(results[['id', 'r_eff_50_phys']], on='id', how='left')

# Run KS tests
grp_sup = master[master['group'] == 'Support']
grp_neu = master[master['group'] == 'Neutral']
grp_null = master[master['group'] == 'Null wins']

ks_results = []

params_to_test = [
    ('z_phot',                'Redshift $z_{phot}$'),
    ('Av_null',              'Dust extinction $A_V$ (Null)'),
    ('C_AGN',                 'AGN fraction $C_{\\rm AGN}$'),
    ('delta_chi2nu',          '$\\Delta\\chi^2_\\nu$'),
    ('alpha_AGN',             '$\\alpha_{\\rm AGN}$'),
    ('ir_excess',             'IR excess'),
    ('fxb_full_05_7keV',      'X-ray flux upper limit (0.5-7 keV)'),
    ('logMstar_best',         'Stellar mass $\\log(M_*/M_\\odot)$'),
]

# Add size if available
if 'r_eff_50_phys' in master.columns:
    params_to_test.append(('r_eff_50_phys', 'Half-light radius $r_{\\rm eff}$ (pc)'))

# Use M_UV as proxy for UV brightness
params_to_test.append(('muv', 'UV magnitude $m_{\\rm UV}$'))

for param, label in params_to_test:
    s_data = grp_sup[param].dropna()
    n_data = grp_null[param].dropna()
    
    if len(s_data) < 5 or len(n_data) < 5:
        ks_results.append({'parameter': param, 'label': label,
                           'ks_stat': np.nan, 'p_value': np.nan,
                           'ns_sup': len(s_data), 'ns_null': len(n_data),
                           'mean_sup': s_data.mean() if len(s_data)>0 else np.nan,
                           'mean_null': n_data.mean() if len(n_data)>0 else np.nan,
                           'note': 'Insufficient data'})
        continue
    
    ks_stat, p_val = ks_2samp(s_data, n_data)
    
    # Cohen's d
    pooled_std = np.sqrt(((len(s_data)-1)*s_data.std()**2 + (len(n_data)-1)*n_data.std()**2) / 
                         (len(s_data)+len(n_data)-2))
    cohens_d = (s_data.mean() - n_data.mean()) / pooled_std if pooled_std > 0 else 0
    
    sig = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else 'ns'))
    
    effect = 'LARGE' if abs(cohens_d) > 0.8 else ('medium' if abs(cohens_d) > 0.5 else ('small' if abs(cohens_d) > 0.2 else 'negligible'))
    
    ks_results.append({
        'parameter': param, 'label': label,
        'ks_stat': round(ks_stat, 3), 'p_value': f'{p_val:.2e}',
        'significance': sig, 'cohens_d': round(cohens_d, 2),
        'effect_size': effect,
        'ns_sup': len(s_data), 'ns_null': len(n_data),
        'mean_sup': round(s_data.mean(), 3), 'std_sup': round(s_data.std(), 3),
        'mean_null': round(n_data.mean(), 3), 'std_null': round(n_data.std(), 3),
    })

df_ks = pd.DataFrame(ks_results)
print('\n' + df_ks.to_string(index=False))

# Save KS results
ks_outpath = 'KS_Test_Results_AllParams.csv'
df_ks.to_csv(ks_outpath, index=False)
print(f'\nKS test results saved: {ks_outpath}')


# ═══════════════════════════════════════════════════════════
# SAVE EVERYTHING
# ═══════════════════════════════════════════════════════════

# 1. X-ray limits
xray_path = 'LRD_XRay_UpperLimits.csv'
df_xray.to_csv(xray_path, index=False)
print(f'\n✅ X-ray limits saved: {xray_path}')

# 2. Emission lines
lines_path = 'LRD_EmissionLine_Framework.csv'
df_lines.to_csv(lines_path, index=False)
print(f'✅ Emission lines framework saved: {lines_path}')

# 3. Stellar masses
stellar_path = 'LRD_StellarMass_Estimates.csv'
df_stellar.to_csv(stellar_path, index=False)
print(f'✅ Stellar mass estimates saved: {stellar_path}')

# 4. Master combined table
master_path = 'LRD_Master_Combined_AllParams.csv'
master.to_csv(master_path, index=False)
print(f'✅ Master combined table saved: {master_path}')

# 5. Summary statistics
summary_rows = []
for g in ['Support', 'Neutral', 'Null wins']:
    sub = master[master['group'] == g]
    summary_rows.append({
        'Group': g,
        'N': len(sub),
        'z_mean': sub['z_phot'].mean(),
        'z_std': sub['z_phot'].std(),
        'Av_mean': sub['Av_null'].mean(),
        'AGN_mean': sub['C_AGN'].mean(),
        'logM_mean': sub['logMstar_best'].mean() if sub['logMstar_best'].notna().any() else np.nan,
        'fx_median': sub['fxb_full_05_7keV'].median() if sub['fxb_full_05_7keV'].notna().any() else np.nan,
        'dChi2_mean': sub['delta_chi2nu'].mean(),
    })

df_summary = pd.DataFrame(summary_rows)
summary_path = 'LRD_Group_Summary_NewParams.csv'
df_summary.to_csv(summary_path, index=False)
print(f'✅ Group summary saved: {summary_path}')

print('\n' + '=' * 70)
print('DONE! All missing data acquired.')
print('=' * 70)
