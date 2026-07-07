#!/usr/bin/env python3
"""
Reproduce Mazzolari+ (2026) three LRD radio models and cross-compare
with our LRD sample data.

References:
  - Mazzolari+ 2026a (A&A 706, A372): Radio properties of JWST AGN
  - Mazzolari+ 2026b (AASKAII): SKAO predictions for LRD

Author: Xin Tan
Date: 2026-07-06
"""

import numpy as np
from scipy import stats, interpolate
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import json
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Cosmology & Constants
# ============================================================================
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
c_kms = 299792.458  # km/s

# ============================================================================
# Section 1: Reproduce Mazzolari+ 3-Model Radio Predictions
# ============================================================================

def luminosity_distance(z):
    """Luminosity distance in cm."""
    return cosmo.luminosity_distance(z).to(u.cm).value

def rest_to_obs_flux(L_nu_rest, nu_rest, z, alpha=-0.5):
    """
    Convert rest-frame luminosity L_nu (erg/s/Hz) to observed flux S_nu (uJy).
    S_nu = L_nu_rest * (1+z)^(1+alpha) / (4π dL^2)
    Returns flux in microJansky.
    """
    dL_cm = luminosity_distance(z)
    # L_nu in erg/s/Hz -> convert to W/Hz (×1e-7), then to Jy (×1e26), then to uJy (×1e6)
    # S = L_nu * (1+z)^(1+alpha) / (4π dL^2)
    flux_cgs = L_nu_rest * (1 + z)**(1 + alpha) / (4 * np.pi * dL_cm**2)  # erg/s/cm^2/Hz
    flux_jy = flux_cgs * 1e23   # convert: 1 Jy = 1e-23 erg/s/cm^2/Hz, so multiply by 1e23
    flux_ujy = flux_jy * 1e6    # uJy
    return flux_ujy

def fundamental_plane_L5GHz(L_X_2_10keV, M_BH):
    """
    Wang+2024 Fundamental Plane for RQ AGN.

    The FP gives νL_ν at 5 GHz in erg/s.
    We divide by ν=5e9 Hz to get L_ν(5 GHz) in erg/s/Hz.
    """
    nu_5GHz = 5e9
    nuLnu_5GHz = 10**(0.47 * np.log10(L_X_2_10keV) + 0.29 * np.log10(M_BH) + 17.06)
    return nuLnu_5GHz / nu_5GHz

def d2020_bol_to_X(L_bol):
    """
    Bolometric to X-ray (2-10 keV) luminosity for Type 1 AGN.

    Uses standard bolometric correction K_X = L_bol / L_X(2-10 keV).
    For Type 1 AGN: K_X ~ 10-50 (Duras+2020, Lusso+2012).
    This simplified fit avoids the numerical issues of the published
    polynomial form at typical high-z BLAGN luminosities.
    """
    log_Lbol = np.log10(np.atleast_1d(L_bol))
    log_KX = 1.2 + 0.1 * (log_Lbol - 44)  # K_X ~ 16 at L_bol=1e45
    K_X = np.clip(10**log_KX, 8, 100)  # reasonable bounds for Type 1
    return L_bol / K_X

def free_free_tau(nu_rest, Ne=1e23, Te=5e4, r=0.01):
    """
    Free-free optical depth (Mazzolari Eq. 5):
    τ_ff = 3 × (Ne/1e22)^2 × (Te/1e4)^-1.35 × (ν_rest/GHz)^-2.1 × (r/1pc)^-1

    Parameters:
    - Ne: free electron column density [cm^-2]
    - Te: electron temperature [K]
    - r: thickness of ionized medium [pc]
    """
    tau = 3.0 * (Ne/1e22)**2 * (Te/1e4)**-1.35 * (nu_rest/1.0)**-2.1 * (r/1.0)**-1
    return tau

def star_formation_radio_flux(z, M_star, frequencies_obs):
    """
    Compute star-formation radio flux following Mazzolari Sect 3.1.4.
    """
    # MS SFR from Popesso+2023
    age_Gyr = cosmo.age(z).value  # age of universe at z in Gyr
    a0, a1, a2, a3, a4 = 0.17, 0.12, 2.44, -0.24, 0.43
    t = age_Gyr
    SFR_MS = 10**(a0 + a1*t - np.log10(1 + (M_star / 10**(a2 + a3*t))**-a4))

    # Radio luminosity from SFR via Novak+2017 / Delvecchio+2021
    q_IR = 2.646 * (1 + z)**(-0.023) - 0.148 * (np.log10(M_star) - 10)
    # L_1.4GHz [W/Hz]
    f_IMF = 1.0  # Chabrier IMF
    L_1p4_WHz = SFR_MS * 1e24 / (f_IMF * 10**q_IR)
    # erg/s/Hz
    L_1p4 = L_1p4_WHz * 1e7

    # Scale to other frequencies using alpha_SF = -0.7
    alpha_SF = -0.7
    fluxes = {}
    for freq in frequencies_obs:
        nu_rest = freq * (1 + z)
        L_nu = L_1p4 * (nu_rest / 1.4)**alpha_SF
        fluxes[freq] = rest_to_obs_flux(L_nu, nu_rest, z, alpha=alpha_SF)

    return fluxes, SFR_MS

# ============================================================================
# Section 2: Reproduce Figure 2 from Mazzolari+ SKAO paper
# ============================================================================

def compute_model_predictions(z, L_bol, M_BH_values, frequencies_obs, scenarios):
    """
    Compute radio flux predictions for all three models.
    Returns dict[scenario][M_BH] = {freq: flux_uJy}
    """
    results = {}

    for scenario in scenarios:
        results[scenario] = {}
        for M_BH in M_BH_values:
            results[scenario][M_BH] = {}

            # Compute X-ray luminosity
            L_X_expected = d2020_bol_to_X(L_bol)

            if scenario == 'compton_thick':
                # Standard AGN, X-ray obscured but radio unaffected
                L_X = L_X_expected

            elif scenario == 'intrinsically_weak':
                # Intrinsically X-ray weak (2 dex fainter)
                L_X = L_X_expected * 0.01

            elif scenario == 'free_free':
                # Standard radio but suppressed by free-free absorption
                L_X = L_X_expected  # standard X-ray production

            # Compute rest-frame 5 GHz luminosity via FP
            L_5GHz = fundamental_plane_L5GHz(L_X, M_BH)

            for freq in frequencies_obs:
                nu_rest_5GHz = 5.0 * (1 + z)
                nu_obs_GHz = freq
                nu_rest_GHz = nu_obs_GHz * (1 + z)

                # Scale from 5 GHz to observed frequency using alpha=-0.5
                L_nu_rest = L_5GHz * (nu_rest_GHz / 5.0)**(-0.5)

                flux = rest_to_obs_flux(L_nu_rest, nu_rest_GHz, z, alpha=-0.5)

                if scenario == 'free_free':
                    tau = free_free_tau(nu_rest_GHz)
                    flux *= np.exp(-tau)

                results[scenario][M_BH][freq] = flux

    return results

# ============================================================================
# Section 3: SKAO Sensitivity limits (from Mazzolari+ Table 2)
# ============================================================================

SKAO_SENSITIVITIES = {
    'SKA-Low':    {'freq': 0.2,  'rms_1hr_uJy': 9.3,  'res_arcsec': 8.0},
    'Band-1':     {'freq': 0.8,  'rms_1hr_uJy': 2.3,  'res_arcsec': 1.15},
    'Band-2':     {'freq': 1.4,  'rms_1hr_uJy': 1.1,  'res_arcsec': 0.7},
    'Band-5a':    {'freq': 6.5,  'rms_1hr_uJy': 0.7,  'res_arcsec': 0.12},
    'Band-5b':    {'freq': 11.5, 'rms_1hr_uJy': 0.84, 'res_arcsec': 0.07},
}

def skao_5sigma_limit(rms_1hr_uJy, t_obs_hr):
    """SKAO 5-sigma detection limit for integration time t_obs (hours)."""
    return 5 * rms_1hr_uJy / np.sqrt(t_obs_hr)

# ============================================================================
# Section 4: OBSERVED radio upper limits (from Mazzolari+ 2026a Table 2)
# ============================================================================

OBSERVED_LIMITS = {
    # From stack of 37 BLAGN on GOODS-N
    'freqs_GHz': np.array([0.144, 1.5, 3.0, 5.5, 10.0]),
    'rms_mean_uJy': np.array([4.63, 0.22, 0.15, 0.39, 0.12]),
    'rms_median_uJy': np.array([7.34, 0.27, 0.18, 0.44, 0.15]),
    'peak_mean_uJy': np.array([-5.92, 0.18, 0.15, 0.11, 0.05]),
    'peak_median_uJy': np.array([5.14, 0.23, 0.18, 0.56, 0.13]),
}

# ---- Literature stacking results (from SKAO paper Table 1) ----
LITERATURE_STACKS = {
    'Akins+2025 (P, N=434, 3GHz)': {'freq': 3.0, 'flux_3sig': 0.36},
    'Akins+2025 (P, N=434, 1.28GHz)': {'freq': 1.28, 'flux_3sig': 1.08},
    'Perger+2024 (P, N=919, 1.4GHz)': {'freq': 1.4, 'flux_3sig': 18.0},
    'Perger+2024 (P, N=919, 3GHz)': {'freq': 3.0, 'flux_3sig': 11.0},
    'Gloudemans+2025 (P, N=593, 3GHz)': {'freq': 3.0, 'flux_3sig': 0.41},
    'Gloudemans+2025 (P, N=593, 1.3GHz)': {'freq': 1.3, 'flux_3sig': 0.44},
    'Gloudemans+2025 (P, N=593, 144MHz)': {'freq': 0.144, 'flux_3sig': 17.0},
    'Mazzolari+2026a (S, N=37, 144MHz)': {'freq': 0.144, 'flux_3sig': 14.0},
    'Mazzolari+2026a (S, N=37, 1.4GHz)': {'freq': 1.4, 'flux_3sig': 0.65},
    'Mazzolari+2026a (S, N=37, 3GHz)': {'freq': 3.0, 'flux_3sig': 0.45},
    'Mazzolari+2026a (S, N=37, 5GHz)': {'freq': 5.0, 'flux_3sig': 1.16},
    'Mazzolari+2026a (S, N=37, 10GHz)': {'freq': 10.0, 'flux_3sig': 0.37},
}

# ============================================================================
# Section 5: Main Reproduction — Figure 2 of Mazzolari+ SKAO Paper
# ============================================================================

def reproduce_mazzolari_fig2(savepath='figures/mazzolari_reproduction_fig2.png'):
    """Reproduce Mazzolari+ SKAO paper Figure 2."""

    # Parameters from the paper
    L_bol = 1e45  # erg/s (median for JWST AGN)
    M_BH_values = [1e7, 3.16e7, 1e8]  # 7.0, 7.5, 8.0 in log
    M_star = 10**8.5  # M_sun
    z_values = [5.5, 3.0]

    frequencies = np.array([0.15, 0.2, 0.3, 0.5, 0.7, 0.8, 1.0, 1.4, 2.0, 3.0, 5.0, 6.5, 8.0, 10.0, 11.5])
    scenarios = ['compton_thick', 'intrinsically_weak', 'free_free']
    colors = {'compton_thick': '#d62728', 'intrinsically_weak': '#1f77b4', 'free_free': '#2ca02c'}
    labels_cn = {'compton_thick': 'Compton-thick AGN',
                 'intrinsically_weak': 'Intrinsically Weak',
                 'free_free': 'Free-Free Absorption'}

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    for row, z in enumerate(z_values):
        results = compute_model_predictions(z, L_bol, M_BH_values, frequencies, scenarios)
        sf_fluxes, sfr = star_formation_radio_flux(z, M_star, frequencies)

        for col, scenario in enumerate(scenarios):
            ax = axes[row, col]

            # Plot AGN predictions
            for i, M_BH in enumerate(M_BH_values):
                fluxes = [results[scenario][M_BH].get(f, np.nan) for f in frequencies]
                if i == 1:  # Main prediction (M_BH = 10^7.5)
                    ax.plot(frequencies, fluxes, '-', color=colors[scenario], linewidth=2.5,
                           label=f'AGN ($M_{{\\rm BH}}=10^{{{np.log10(M_BH):.1f}}}\\,M_\\odot$)')
                elif i == 0:
                    ax.plot(frequencies, fluxes, '--', color=colors[scenario], linewidth=1.0, alpha=0.6,
                           label=f'$M_{{\\rm BH}}=10^{{7.0}}\\,M_\\odot$')
                else:
                    ax.plot(frequencies, fluxes, '-.', color=colors[scenario], linewidth=1.0, alpha=0.6,
                           label=f'$M_{{\\rm BH}}=10^{{8.0}}\\,M_\\odot$')

            # Plot SF contribution
            sf_low = [sf_fluxes[f] / 2 for f in frequencies]  # rough scatter
            sf_high = [sf_fluxes[f] * 2 for f in frequencies]
            ax.fill_between(frequencies, sf_low, sf_high, alpha=0.15, color='gray',
                           label='SF (MS, ±0.3 dex)')

            # SKAO 5-sigma limits
            t_obs_dict = {
                'SKA-Low': 1,
                'Band-1': 30,
                'Band-2': [1, 30, 1000],
                'Band-5a': 30,
                'Band-5b': 100,
            }

            band_colors = {
                'SKA-Low': 'gray',
                'Band-1': 'blue',
                'Band-2': 'green',
                'Band-5a': 'gold',
                'Band-5b': 'orange',
            }

            for band_name, band_info in SKAO_SENSITIVITIES.items():
                freq = band_info['freq']
                rms = band_info['rms_1hr_uJy']

                if band_name == 'Band-2':
                    for t_obs, ls in zip([1, 30, 1000], ['--', '-', ':']):
                        lim_5sig = skao_5sigma_limit(rms, t_obs)
                        ax.axhline(y=lim_5sig, color='green', linestyle=ls, alpha=0.5, linewidth=1.2)
                        if row == 0:
                            ax.text(0.17, lim_5sig * 1.1, f'Band-2 ({t_obs}h)', fontsize=6, color='green', alpha=0.7)
                elif band_name == 'Band-2':
                    pass
                else:
                    t_obs = t_obs_dict.get(band_name, 30)
                    lim_5sig = skao_5sigma_limit(rms, t_obs)
                    ax.axhline(y=lim_5sig, color=band_colors.get(band_name, 'black'),
                              linestyle='-', alpha=0.5, linewidth=1.2)
                    if row == 0:
                        ax.text(0.17, lim_5sig * 1.1, f'{band_name} ({t_obs}h)',
                               fontsize=5.5, color=band_colors.get(band_name, 'black'), alpha=0.7)

            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(0.12, 13)
            ax.set_ylim(1e-3, 50)
            ax.set_xlabel('Observed Frequency [GHz]', fontsize=11)
            ax.set_ylabel('Flux Density [μJy]', fontsize=11)
            ax.set_title(f'{labels_cn[scenario]}\nz={z}, $L_{{\\rm bol}}=10^{{{45}}}$ erg/s, '
                        f'$M_* = 10^{{8.5}}\\,M_\\odot$', fontsize=11)
            ax.legend(fontsize=7, loc='lower left')
            ax.grid(True, alpha=0.3)

    plt.suptitle('Mazzolari+ (2026) SKAO Predictions — Reproduction', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(savepath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {savepath}")
    return fig

# ============================================================================
# Section 6: Cross-comparison with our LRD data
# ============================================================================

def load_lrd_data():
    """Load our LRD catalogs and compute expected radio luminosities."""
    import pandas as pd
    from astropy.io import fits

    data = {}

    # Try loading the Kokorev 260 catalog
    try:
        kokorev = pd.read_csv('/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/data/kokorev_260_sb.csv')
        # Use correct column names: z_phot, logMstar_best, logSigma_Mstar
        kokorev = kokorev.rename(columns={'z_phot': 'z', 'logMstar_best': 'logMstar',
                                          'logSigma_Mstar': 'logSigma'})
        data['kokorev'] = kokorev
        print(f"Loaded Kokorev 260: {len(kokorev)} sources")
    except Exception as e:
        print(f"Could not load Kokorev catalog: {e}")

    # Try loading deGraaff cross-matched sample
    try:
        path1 = pd.read_csv('/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/data/path1_merged_38sources.csv')
        # Use correct columns: z_phot, logMstar_best, lbol, Balmer_dec_total
        path1 = path1.rename(columns={'z_phot': 'z', 'logMstar_best': 'logMstar',
                                      'lbol': 'Lbol', 'Balmer_dec_total': 'BD',
                                      'logLHa_total': 'logLHa'})
        data['path1'] = path1
        print(f"Loaded Path1: {len(path1)} sources")
    except Exception as e:
        print(f"Could not load Path1: {e}")

    # Try Weibel+ catalog
    try:
        weibel = fits.open('/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/data/external/weibel2026/W26_bhstar_dominated_sample_phot_catalog.fits')
        wdata = weibel[1].data
        # Convert to simple dict of arrays
        z = wdata['z_phot_eazy']
        Lbol = wdata['Lbol']
        bhs_contrib = wdata['bhs_template_contribution']
        data['weibel'] = {'z': z, 'Lbol': Lbol, 'bhs_contrib': bhs_contrib,
                         'N': len(z), 'catalog': wdata}
        print(f"Loaded Weibel+: {len(z)} sources")
    except Exception as e:
        print(f"Could not load Weibel catalog: {e}")

    return data

def estimate_radio_for_lrds(lrd_data, z_col='z', Mstar_col='logMstar',
                             Lbol_col='Lbol', MBH_col=None):
    """
    Estimate expected radio luminosities for LRDs under each model.
    Handles both DataFrames and dict-based catalogs.
    """
    results = {}

    for name, catalog in lrd_data.items():
        try:
            # Handle dict-based catalog (e.g., Weibel)
            if isinstance(catalog, dict):
                z = np.asarray(catalog['z']).flatten()
                n_src = catalog['N']
                M_BH_est = np.full(n_src, 1e7)  # default

                if 'Lbol' in catalog:
                    # Weibel+ Lbol is already linear [erg/s]
                    L_bol = np.asarray(catalog['Lbol']).flatten()
                    # Filter invalid/missing/zero values
                    valid = np.isfinite(L_bol) & (L_bol > 1e40)
                    L_bol = np.where(valid, L_bol, 1e45)
                else:
                    L_bol = np.full(n_src, 1e45)

            elif isinstance(catalog, np.ndarray):  # FITS
                z = catalog['z_phot']
                n_src = len(z)
                M_BH_est = np.full(n_src, 1e7)
                L_bol = np.full(n_src, 1e45)

            else:  # DataFrame
                z = catalog[z_col].values if z_col in catalog.columns else np.full(len(catalog), 5.5)
                n_src = len(z)

                # If we have stellar masses, estimate M_BH from local relations
                if Mstar_col in catalog.columns:
                    log_Mstar = catalog[Mstar_col].values
                    # Reines & Volonteri 2015: log M_BH = 7.45 + 1.05 * log(M*/1e11)
                    log_MBH_est = 7.45 + 1.05 * (log_Mstar - 11.0)
                    # Clip to physically reasonable range
                    log_MBH_est = np.clip(log_MBH_est, 6.0, 9.0)
                    M_BH_est = 10**log_MBH_est
                else:
                    M_BH_est = np.full(n_src, 1e7)  # default

                # Bolometric luminosity
                if Lbol_col and Lbol_col in catalog.columns:
                    L_bol = 10**catalog[Lbol_col].values
                    # Handle bad values
                    L_bol = np.where(np.isfinite(L_bol) & (L_bol > 0), L_bol, 1e45)
                else:
                    L_bol = np.full(n_src, 1e45)  # median from literature

            # Compute for three models at 1.4 GHz and 5 GHz
            results[name] = {'N': n_src, 'z': z, 'M_BH': M_BH_est, 'L_bol': L_bol}

            # Model 1: Compton-thick
            L_X_std = np.array([d2020_bol_to_X(lb) for lb in L_bol])
            L_5GHz_m1 = fundamental_plane_L5GHz(L_X_std, M_BH_est)

            # Model 2: Intrinsically weak
            L_X_weak = L_X_std * 0.01
            L_5GHz_m2 = fundamental_plane_L5GHz(L_X_weak, M_BH_est)

            # Model 3: Free-free
            L_5GHz_m3 = L_5GHz_m1.copy()  # same intrinsic, absorbed

            for freq_obs in [1.4, 5.0]:
                nu_rest = freq_obs * (1 + np.mean(z))
                tau = free_free_tau(nu_rest)

                flux_m1 = np.array([rest_to_obs_flux(L, 5.0*(1+zi), zi) for L, zi in zip(L_5GHz_m1, z)])
                flux_m2 = np.array([rest_to_obs_flux(L, 5.0*(1+zi), zi) for L, zi in zip(L_5GHz_m2, z)])
                flux_m3 = np.array([rest_to_obs_flux(L * np.exp(-free_free_tau(5.0*(1+zi))),
                                                      5.0*(1+zi), zi) for L, zi in zip(L_5GHz_m1, z)])

                results[name][f'flux_m1_{freq_obs}GHz'] = flux_m1
                results[name][f'flux_m2_{freq_obs}GHz'] = flux_m2
                results[name][f'flux_m3_{freq_obs}GHz'] = flux_m3

            print(f"\n{name} (N={n_src}):")
            print(f"  Model 1 (Compton-thick): median flux @1.4GHz = {np.median(flux_m1):.3f} μJy")
            print(f"  Model 2 (Intr. weak):     median flux @1.4GHz = {np.median(flux_m2):.3f} μJy")
            print(f"  Model 3 (Free-free):      median flux @1.4GHz = {np.median(flux_m3):.6f} μJy")

        except Exception as e:
            print(f"Error processing {name}: {e}")

    return results

# ============================================================================
# Section 7: Dyson model radio predictions (simple estimate)
# ============================================================================

def dyson_radio_estimate():
    """
    First-principles estimate of radio emission from Dyson Entrainment model.

    Unlike AGN models (accretion disk + corona → synchrotron), the Dyson model
    has NO accretion disk and NO magnetized corona. Radio emission can only
    come from:

    1. Stellar wind synchrotron: weak, from stellar magnetic fields
    2. Collision-induced shocks: transient, non-thermal electrons
    3. Free-free emission from ionized gas (thermal bremsstrahlung)

    We estimate upper limits for each.
    """
    print("\n" + "="*70)
    print("DYSON MODEL RADIO EMISSION ESTIMATES")
    print("="*70)

    # Stellar parameters
    N_stars = 1e6  # ~10^6 stars in the dense core
    B_surf = 100  # Gauss, typical stellar surface field
    R_star = 1e11  # cm, ~1.4 R_sun for massive stars

    # Synchrotron from stellar winds
    # Each massive star: wind mass loss ~1e-6 M_sun/yr, v_wind ~2000 km/s
    Mdot_per_star = 1e-6 * 2e33 / 3.15e7  # g/s
    v_wind = 2000e5  # cm/s
    B_ISM = 1e-4  # G, typical ISM field compressed by winds

    # Total wind kinetic power
    L_wind_total = 0.5 * N_stars * Mdot_per_star * v_wind**2  # erg/s
    print(f"Total stellar wind kinetic power: {L_wind_total:.2e} erg/s")

    # Synchrotron efficiency (very low for wind-ISM interaction)
    # Typically ~1e-5 to 1e-4 of wind power
    epsilon_sync = 1e-5
    L_sync = epsilon_sync * L_wind_total
    print(f"Estimated synchrotron luminosity: {L_sync:.2e} erg/s")

    # Convert to 1.4 GHz radio luminosity (assuming flat-ish spectrum)
    nu_1p4 = 1.4e9  # Hz
    L_nu_sync = L_sync / nu_1p4  # rough, assuming power distributed over ~1 dex in freq
    print(f"L_nu(1.4 GHz) synchrotron: {L_nu_sync:.2e} erg/s/Hz")

    # At z=5.5, convert to observed flux
    z = 5.5
    flux_sync_uJy = rest_to_obs_flux(L_nu_sync, nu_1p4*(1+z), z, alpha=-0.7)
    print(f"Observed flux @1.4 GHz: {flux_sync_uJy:.6f} μJy")

    # Free-free emission from ionized gas (Dyson core)
    # Core parameters: ne ~ 10^6 cm^-3 in inner 0.01 pc
    Te = 1e4  # K
    ne = 1e6  # cm^-3 (dense nuclear core, ~10^6 e-/cm^3)
    r_core_pc = 0.03  # pc — inner core radius
    r_core_cm = r_core_pc * 3.086e18  # cm
    V = (4/3) * np.pi * r_core_cm**3  # cm^3

    # Free-free volume emissivity (Rybicki & Lightman Eq. 5.14)
    # j_ff = 6.8e-38 * ne^2 * Te^-0.5 * exp(-hν/kT) * g_ff  erg/s/cm^3/Hz
    g_ff = 5.0  # Gaunt factor at ~1 GHz, Te~1e4
    hnu_over_kT = 4.8e-11 * nu_1p4 / (1.38e-16 * Te)  # hν/kT
    j_ff_nu = 6.8e-38 * ne**2 * Te**-0.5 * np.exp(-hnu_over_kT) * g_ff  # erg/s/cm^3/Hz
    L_nu_ff = j_ff_nu * V  # total free-free luminosity
    flux_ff_uJy = rest_to_obs_flux(L_nu_ff, nu_1p4*(1+z), z, alpha=0.0)  # free-free is optically thin, flat
    print(f"L_nu(1.4 GHz) free-free: {L_nu_ff:.2e} erg/s/Hz")
    print(f"Observed flux @1.4 GHz (free-free): {flux_ff_uJy:.6f} μJy")

    # SKAO Band-2 sensitivity
    skao_band2_5sig_1000h = skao_5sigma_limit(1.1, 1000)
    skao_band5a_5sig_100h = skao_5sigma_limit(0.7, 100)

    print(f"\nSKAO Detection Limits:")
    print(f"  Band-2 (1.4 GHz, 1000h, 5σ): {skao_band2_5sig_1000h:.4f} μJy")
    print(f"  Band-5a (6.5 GHz, 100h, 5σ): {skao_band5a_5sig_100h:.4f} μJy")
    print(f"\nDyson Model Predictions @1.4 GHz:")
    print(f"  Synchrotron: {flux_sync_uJy:.6f} μJy — {'DETECTABLE' if flux_sync_uJy > skao_band2_5sig_1000h else 'UNDETECTABLE'} even at 1000h")
    print(f"  Free-free:   {flux_ff_uJy:.6f} μJy — {'DETECTABLE' if flux_ff_uJy > skao_band2_5sig_1000h else 'UNDETECTABLE'} even at 1000h")

    return {
        'L_sync': L_sync,
        'L_nu_sync': L_nu_sync,
        'flux_sync_uJy': flux_sync_uJy,
        'L_nu_ff': L_nu_ff,
        'flux_ff_uJy': flux_ff_uJy,
        'skao_band2_5sig_1000h': skao_band2_5sig_1000h,
        'skao_band5a_5sig_100h': skao_band5a_5sig_100h,
    }

# ============================================================================
# Section 8: Observability Diagnostic Plot
# ============================================================================

def plot_mazzolari_vs_dyson(savepath='figures/mazzolari_dyson_comparison.png'):
    """Plot Mazzolari's three models alongside Dyson model predictions."""

    z = 5.5
    L_bol = 1e45
    M_BH = 3.16e7  # 10^7.5
    M_star = 10**8.5

    frequencies = np.logspace(np.log10(0.1), np.log10(15), 100)
    scenarios = ['compton_thick', 'intrinsically_weak', 'free_free']

    results = compute_model_predictions(z, L_bol, [M_BH], frequencies, scenarios)
    sf_fluxes, sfr = star_formation_radio_flux(z, M_star, frequencies)
    dyson = dyson_radio_estimate()

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    colors = {'compton_thick': '#d62728', 'intrinsically_weak': '#1f77b4', 'free_free': '#2ca02c'}
    linestyles = {'compton_thick': '-', 'intrinsically_weak': '--', 'free_free': ':'}

    for scenario in scenarios:
        fluxes = [results[scenario][M_BH][f] for f in frequencies]
        ax.plot(frequencies, fluxes, color=colors[scenario], linestyle=linestyles[scenario],
                linewidth=2.5, label=f'Mazzolari: {scenario.replace("_", " ").title()}')

    # SF
    ax.fill_between(frequencies, [sf_fluxes[f]/1.5 for f in frequencies],
                   [sf_fluxes[f]*1.5 for f in frequencies],
                   alpha=0.12, color='gray', label='Star Formation (±0.18 dex)')

    # Dyson model limits
    ax.axhline(y=dyson['flux_sync_uJy'], color='purple', linestyle='--', linewidth=2,
              label=f'Dyson: Synchrotron ({dyson["flux_sync_uJy"]:.1e} μJy)')
    ax.axhline(y=dyson['flux_ff_uJy'], color='darkviolet', linestyle=':', linewidth=2,
              label=f'Dyson: Free-Free ({dyson["flux_ff_uJy"]:.1e} μJy)')

    # SKAO limits
    for band, t_obs, color, ls in [('Band-2', 1000, 'green', '-'),
                                     ('Band-5a', 100, 'gold', '-'),
                                     ('Band-5b', 100, 'orange', '-')]:
        info = SKAO_SENSITIVITIES[band]
        lim = skao_5sigma_limit(info['rms_1hr_uJy'], t_obs)
        ax.axhline(y=lim, color=color, linestyle=ls, alpha=0.5, linewidth=1.5)
        ax.text(12, lim * 1.15, f'{band} ({t_obs}h, 5σ)', fontsize=8, color=color, ha='right')

    # Current limits from stacking
    stack_freqs = OBSERVED_LIMITS['freqs_GHz']
    stack_rms = OBSERVED_LIMITS['rms_mean_uJy']
    ax.plot(stack_freqs, 3 * stack_rms, 'ko', markersize=8, markerfacecolor='none',
           label='Current Stack Limits (3σ, 37 BLAGN)')

    # Local LRD analog Rodriguez & Mirabel 2026
    ax.axhline(y=0.03, color='cyan', linestyle='--', alpha=0.6, linewidth=1.5)
    ax.text(0.12, 0.035, 'Local LRD (J1047+0739, z=0.168)', fontsize=8, color='cyan')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Observed Frequency [GHz]', fontsize=13)
    ax.set_ylabel('Flux Density [μJy]', fontsize=13)
    ax.set_title(f'Mazzolari+ Models vs Dyson Entrainment (z={z}, $L_{{\\rm bol}}=10^{{45}}$ erg/s)',
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='lower left', ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.08, 15)
    ax.set_ylim(1e-8, 50)

    plt.tight_layout()
    plt.savefig(savepath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {savepath}")

# ============================================================================
# Section 9: Statistical test — probability of non-detection
# ============================================================================

def monte_carlo_nondetection_probability():
    """
    Reproduce Mazzolari Table 3: probability of non-detecting all 37 BLAGN
    given standard RQ AGN relations.
    """
    print("\n" + "="*70)
    print("MONTE CARLO: Probability of non-detecting all 37 BLAGN")
    print("="*70)

    np.random.seed(42)
    N_sources = 37
    N_trials = 10000

    # Reference 3σ limit at 3 GHz (stack): 0.45 uJy
    # Convert to rest-frame 5 GHz luminosity at z=5.2
    z_mean = 5.2
    flux_limit_3sig = 0.45  # uJy

    # Relations with their scatters
    relations = {
        'L_X - L_rad (D\'Amato+2022)': {'scatter': 0.5, 'extra_scatter': 0.27},
        'Fundamental Plane (Bariuan+2022)': {'scatter': 0.39, 'extra_scatter': 0.27},
        'L_opt - L_rad (log R=1)': {'scatter': 0.42, 'extra_scatter': 0.0},
        'L_opt - L_rad (log R=0)': {'scatter': 0.42, 'extra_scatter': 0.0},
        'L_[OIII] - L_rad (de Vries+2007)': {'scatter': 0.8, 'extra_scatter': 0.0},
    }

    # Use median expected luminosity from the paper
    # From stacked 3σ: L_5GHz < 10^39.4 erg/s (stack)
    # Expected RQ: L_5GHz ~ 10^39.5-10^40.5 erg/s

    probs = {}
    for name, params in relations.items():
        total_scatter = np.sqrt(params['scatter']**2 + params['extra_scatter']**2)

        # Expected median log L_5GHz ~ 39.8 (from Fig 3, stack position)
        not_detected_count = 0
        for trial in range(N_trials):
            # Draw expected L_5GHz from Gaussian
            log_L_exp = np.random.normal(39.8, total_scatter, N_sources)
            # Observed limit for each source ~39.7 (individual 3σ)
            # If all expected < observed limit → all undetected
            if np.all(log_L_exp < 39.7):
                not_detected_count += 1

        prob = not_detected_count / N_trials
        probs[name] = prob
        print(f"  {name}: P = {prob:.2e} (scatter={total_scatter:.2f} dex)")

    return probs

# ============================================================================
# Main
# ============================================================================

if __name__ == '__main__':
    print("="*70)
    print("MAZZOLARI+ LRD RADIO ANALYSIS — REPRODUCTION & CROSS-COMPARISON")
    print("="*70)

    # 1. Reproduce Figure 2
    print("\n[1/4] Reproducing Mazzolari+ SKAO Figure 2...")
    reproduce_mazzolari_fig2()

    # 2. Reproduce statistical test
    print("\n[2/4] Monte Carlo non-detection probability...")
    probs = monte_carlo_nondetection_probability()

    # 3. Dyson model radio
    print("\n[3/4] Dyson model radio estimates...")
    dyson_results = dyson_radio_estimate()

    # 4. Cross-comparison plot
    print("\n[4/4] Generating Mazzolari vs Dyson comparison plot...")
    plot_mazzolari_vs_dyson()

    # 5. Load our LRD data and compute expected radio fluxes
    print("\n[Bonus] Computing expected radio fluxes for our LRD samples...")
    try:
        lrd_data = load_lrd_data()
        if lrd_data:
            lrd_radio = estimate_radio_for_lrds(lrd_data)

            # Save results
            output = {
                'monte_carlo_probs': probs,
                'dyson_model': dyson_results,
                'skao_limits': {band: {'rms_1hr': info['rms_1hr_uJy'],
                                        '5sig_1h': skao_5sigma_limit(info['rms_1hr_uJy'], 1),
                                        '5sig_30h': skao_5sigma_limit(info['rms_1hr_uJy'], 30),
                                        '5sig_100h': skao_5sigma_limit(info['rms_1hr_uJy'], 100),
                                        '5sig_1000h': skao_5sigma_limit(info['rms_1hr_uJy'], 1000)}
                               for band, info in SKAO_SENSITIVITIES.items()},
            }

            with open('results/mazzolari_analysis_results.json', 'w') as f:
                json.dump(output, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else x)
            print("\nResults saved to results/mazzolari_analysis_results.json")
    except Exception as e:
        print(f"LRD data load failed: {e}")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
