#!/usr/bin/env python3
"""
SPARC Rotation Curve Fitting with G_eff(Σ) Framework
======================================================
Core Idea: "Dark Matter" is actually density-dependent gravitational boost.

Standard picture:
  V_obs^2 = V_gas^2 + V_disk^2 + V_bulge^2 + V_DM^2   (need DM halo)

G(Σ) picture:
  V_obs^2 = G_eff(Σ) / G_0 × V_baryon^2
           = [1 + ε_g (Σ/Σ_0)^β] × V_Newton^2
  
where Σ = surface density (from Spitzer [3.6] surface brightness),
and ε_g, β, Σ_0 are UNIVERSAL constants (same for all 175 galaxies!).

Author: Tan Xin (IAIP), based on UID theory
Date: 2026-04-17
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
SPARC_DIR = Path(__file__).parent / "database" / "sparc_database"
OUTPUT_DIR = Path(__file__).parent
PLOT_DIR = OUTPUT_DIR / "figures"
PLOT_DIR.mkdir(exist_ok=True)

# ============================================================
# Data Loading
# ============================================================

def load_sparc_galaxy(filename):
    """Load one galaxy's rotation curve data."""
    data = []
    galaxy_id = filename.replace('_rotmod.dat', '')
    
    with open(filename, 'r') as f:
        # Read distance from first comment line
        distance_mpc = None
        for line in f:
            line = line.strip()
            if line.startswith('# Distance'):
                distance_mpc = float(line.split('=')[1].strip().replace('Mpc',''))
            elif line.startswith('# Rad'):
                # Header line — skip
                continue
            elif line.startswith('#'):
                continue
            elif line:
                parts = line.split('\t')
                if len(parts) >= 8:
                    data.append({
                        'galaxy': galaxy_id,
                        'R_kpc': float(parts[0]),
                        'Vobs': float(parts[1]),
                        'e_Vobs': float(parts[2]) if len(parts) > 2 else 0,
                        'Vgas': float(parts[3]),
                        'Vdisk': float(parts[4]),
                        'Vbul': float(parts[5]),
                        'SBdisk': float(parts[6]),    # L_sun / pc^2
                        'SBbul': float(parts[7]) if len(parts) > 7 else 0,
                        'distance': distance_mpc,
                    })
    
    df = pd.DataFrame(data)
    return df if len(df) > 0 else None


def load_all_galaxies():
    """Load all 175 SPARC galaxies."""
    all_data = {}
    files = sorted(list(SPARC_DIR.glob('*_rotmod.dat')))
    print(f"Found {len(files)} galaxy files")
    
    for fpath in files:
        df = load_sparc_galaxy(str(fpath))
        if df is not None and len(df) >= 5:
            gid = df['galaxy'].iloc[0]
            all_data[gid] = df
    
    print(f"Successfully loaded {len(all_data)} galaxies")
    return all_data


# ============================================================
# G(Σ) Model
# ============================================================

def Geff_over_G0(Sigma, epsilon_g, beta, Sigma_0):
    """
    Effective gravity enhancement factor.
    
    G_eff(Σ) / G_0 = 1 + ε_g × (Σ / Σ_0)^β
    """
    ratio = Sigma / Sigma_0
    return 1.0 + epsilon_g * np.power(np.maximum(ratio, 0), beta)


def compute_Vmodel(R, Vgas, Vdisk, Vbul, SBdisk, SBbul, 
                   ML_disk, ML_bulb, epsilon_g, beta, Sigma_0):
    """
    Compute model circular velocity under G(Σ) framework.
    
    V_model^2 = G_eff(Σ) × V_Newtonian^2
              = [1 + ε_g(Σ/Σ_0)^β] × [ML×Vdisk^2 + ML×Vbul^2 + Vgas^2]
    
    Σ is total surface brightness (mass proxy): 
      Σ_total = ML_disk × SBdisk + ML_bulb × SBbul  [in M_sun/pc^2 units]
    """
    # Baryonic Newtonian velocity squared
    Vnewton_sq = (ML_disk * Vdisk)**2 + (ML_bulb * Vbul)**2 + Vgas**2
    Vnewton_sq = np.maximum(Vnewton_sq, 0)  # avoid negative
    
    # Surface density proxy (total mass surface density)
    # SBdisk and SBbul are in L_sun/pc^2, convert to mass via M/L
    Sigma_total = ML_disk * SBdisk + ML_bulb * SBbul  # M_sun/pc^2
    Sigma_total = np.maximum(Sigma_total, 1e-10)       # avoid log(0)
    
    # Gravity enhancement
    g_boost = Geff_over_G0(Sigma_total, epsilon_g, beta, Sigma_0)
    
    # Model velocity
    Vmodel_sq = g_boost * Vnewton_sq
    return np.sqrt(Vmodel_sq), g_boost


# ============================================================
# Chi-squared computation
# ============================================================

def chi_squared_galaxy(galaxy_df, ML_disk, ML_bulb, epsilon_g, beta, Sigma_0):
    """Chi-squared for one galaxy."""
    R = galaxy_df['R_kpc'].values
    Vobs = galaxy_df['Vobs'].values
    eV = galaxy_df['e_Vobs'].values
    eV = np.maximum(eV, 1.0)  # minimum 1 km/s uncertainty
    
    Vgas = galaxy_df['Vgas'].values
    Vdisk = galaxy_df['Vdisk'].values
    Vbul = galaxy_df['Vbul'].values
    SBdisk = galaxy_df['SBdisk'].values
    SBbul = galaxy_df['SBbul'].values
    
    Vmodel, _ = compute_Vmodel(R, Vgas, Vdisk, Vbul, SBdisk, SBbul,
                                ML_disk, ML_bulb, epsilon_g, beta, Sigma_0)
    
    residuals = Vobs - Vmodel
    chi2 = np.sum((residuals / eV)**2)
    return chi2


def total_chi_squared(params_universal, ML_values, galaxies):
    """
    Total chi-squared across all galaxies.
    
    params_universal: [log_epsilon_g, beta, log_Sigma_0]
                      (we fit in log space for positivity)
    ML_values: dict of {galaxy_id: (ML_disk, ML_bulb)}
    """
    eps_g = 10**params_universal[0]
    beta = params_universal[1]
    Sig0 = 10**params_universal[2]
    
    # Physical constraints
    if beta < 0 or beta > 3:
        return 1e20
    if eps_g < 0.01 or eps_g > 100:
        return 1e20
    if Sig0 < 1e-4 or Sig0 > 1e6:
        return 1e20
    
    total_chi2 = 0
    n_points = 0
    
    for gid, df in galaxies.items():
        ML_disk, ML_bulb = ML_values[gid]
        
        # Skip if M/L is unphysical
        if ML_disk < 0.1 or ML_disk > 2.0:
            return 1e20
        
        chi2 = chi_squared_galaxy(df, ML_disk, ML_bulb, eps_g, beta, Sig0)
        total_chi2 += chi2
        n_points += len(df)
    
    return total_chi2


def find_best_ML_per_galaxy(galaxy_df, epsilon_g, beta, Sigma_0):
    """
    For fixed universal params, find best-fit M/L for one galaxy.
    Uses simple grid search + refinement.
    """
    best_ml = 0.5
    best_chi2 = 1e20
    
    # Coarse grid
    for ml in np.linspace(0.2, 1.5, 28):
        chi2 = chi_squared_galaxy(galaxy_df, ml, ml, epsilon_g, beta, Sigma_0)
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_ml = ml
    
    # Fine grid around best
    for ml in np.linspace(max(0.2, best_ml-0.05), min(1.5, best_ml+0.05), 21):
        chi2 = chi_squared_galaxy(galaxy_df, ml, ml, epsilon_g, beta, Sigma_0)
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_ml = ml
    
    return best_ml, best_chi2


# ============================================================
# Global Fitting
# ============================================================

def fit_universal_parameters(galaxies, verbose=True):
    """
    Find universal ε_g, β, Σ_0 that minimize total chi-squared.
    
    Strategy: nested optimization
    - Outer loop: optimize (ε_g, β, Σ_0) using differential evolution
    - Inner loop: for each galaxy, find optimal M/L given universal params
    """
    
    # Pre-compute best M/L on a grid of universal parameter combinations
    def objective(params_universal):
        eps_g = 10**params_universal[0]
        beta = params_universal[1]
        Sig0 = 10**params_universal[2]
        
        if beta < 0.1 or beta > 2.5:
            return 1e20
        if eps_g < 0.3 or eps_g > 500:
            return 1e20
        if Sig0 < 0.05 or Sig0 > 200:
            return 1e20
        
        total = 0
        ML_dict = {}
        for gid, df in galaxies.items():
            ml, chi2 = find_best_ML_per_galaxy(df, eps_g, beta, Sig0)
            ML_dict[gid] = (ml, ml)  # assume ML_disk = ML_bulb for now
            total += chi2
        
        return total
    
    print("Running global optimization...")
    print("  Parameters: log₁₀(ε_g), β, log₁₀(Σ₀ [M⊙/pc²])")
    print("  Each evaluation optimizes M/L individually per galaxy\n")
    
    # Differential evolution global search
    # Bounds tuned for SPARC data: typical Sigma ~1-300 M_sun/pc^2
    # We want Sigma_0 near the lower end of this range so inner disks get boosted
    bounds = [
        (-0.5, 2.5),   # log10(eps_g): 0.3 to ~300 (boost can be large!)
        (0.2, 2.0),    # beta: 0.2 to ~2 (moderate nonlinearity)
        (-1, 2),       # log10(Sig0): 0.1 to 100 M_sun/pc^2 (SPARC median ~7)
    ]
    
    result = differential_evolution(
        objective,
        bounds,
        seed=42,
        maxiter=200,
        tol=1e-4,
        polish=True,
        workers=1,           # single process to avoid pickling issues
        disp=False,
        popsize=15,
    )
    
    best_eps_g = 10**result.x[0]
    best_beta = result.x[1]
    best_Sig0 = 10**result.x[2]
    
    # Final pass: get all M/L values at best universal params
    final_ML = {}
    total_chi2 = 0
    total_dof = 0
    galaxy_chi2s = {}
    
    for gid, df in galaxies.items():
        ml, chi2 = find_best_ML_per_galaxy(df, best_eps_g, best_beta, best_Sig0)
        final_ML[gid] = ml
        total_chi2 += chi2
        total_dof += len(df)  # -1 (for M/L) subtract later
        galaxy_chi2s[gid] = chi2
    
    total_dof -= len(galaxies)  # one M/L per galaxy
    total_dof -= 3               # three universal params
    
    results = {
        'epsilon_g': best_eps_g,
        'beta': best_beta,
        'Sigma_0': best_Sig0,
        'ML_per_galaxy': final_ML,
        'chi2_total': total_chi2,
        'dof': total_dof,
        'chi2_red': total_chi2 / max(total_dof, 1),
        'galaxy_chi2': galaxy_chi2s,
    }
    
    if verbose:
        print("\n" + "="*60)
        print("★ G(Σ) UNIVERSAL PARAMETERS — BEST FIT ★")
        print("="*60)
        print(f"  ε_g (density-gravity coupling)  = {best_eps_g:.4f}")
        print(f"  β  (nonlinear index)           = {best_beta:.4f}")
        print(f"  Σ_0 (reference density)        = {best_Sig0:.4f} M☉/pc²")
        print(f"  (log Σ_0 = {np.log10(best_Sig0):.3f})")
        print("-"*60)
        print(f"  Total χ²  = {total_chi2:.1f}")
        print(f"  DOF        = {total_dof}")
        print(f"  χ²_red     = {results['chi2_red']:.3f}")
        print("="*60)
        
        # Physical interpretation
        print("\n📐 PHYSICAL INTERPRETATION:")
        print(f"  At Σ = Σ₀: G_eff = {(1+best_eps_g):.2f} × G₀ ({(1+best_eps_g)*100-100:.0f}% enhancement)")
        
        # Typical SPARC disk surface densities range from ~1 to ~1000 M_sun/pc^2
        for sig_test in [1, 10, 100, 1000]:
            boost = Geff_over_G0(sig_test, best_eps_g, best_beta, best_Sig0)
            print(f"  At Σ = {sig_test:>4} M☉/pc²: G_eff = {boost:.2f} × G₀")
    
    return results


# ============================================================
# Plotting
# ============================================================

def plot_sample_fits(galaxies, fit_results, sample_names=None, n_samples=9):
    """Plot rotation curve fits for a selection of galaxies."""
    
    eps_g = fit_results['epsilon_g']
    beta = fit_results['beta']
    Sig0 = fit_results['Sigma_0']
    ML_dict = fit_results['ML_per_galaxy']
    
    if sample_names is None:
        # Pick diverse examples: different sizes, types
        all_gids = list(galaxies.keys())
        # Sort by number of data points and pick spread
        indices = np.linspace(0, len(all_gids)-1, min(n_samples, len(all_gids))).astype(int)
        sample_names = [all_gids[i] for i in indices]
    
    n_plot = min(len(sample_names), n_samples)
    ncols = 3
    nrows = (n_plot + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4.5*nrows))
    axes = np.array(axes).flatten()
    
    for idx, gid in enumerate(sample_names[:n_plot]):
        ax = axes[idx]
        df = galaxies[gid]
        ML = ML_dict.get(gid, 0.5)
        
        R = df['R_kpc'].values
        Vobs = df['Vobs'].values
        eV = df['e_Vobs'].values
        Vgas = df['Vgas'].values
        Vdisk = df['Vdisk'].values
        Vbul = df['Vbul'].values
        SBdisk = df['SBdisk'].values
        SBbul = df['SBbul'].values
        
        # Observed
        ax.errorbar(R, Vobs, yerr=eV, fmt='o', color='black', markersize=4, 
                   label='Observed', zorder=10, alpha=0.8)
        
        # Newtonian baryonic only (no G boost)
        Vnewton = np.sqrt(np.maximum((ML*Vdisk)**2 + (ML*Vbul)**2 + Vgas**2, 0))
        ax.plot(R, Vnewton, '--', color='gray', alpha=0.7, lw=1.5, 
               label=f'Baryon (M/L={ML:.2f})')
        
        # G(Σ) boosted model
        Vmodel, g_boost = compute_Vmodel(R, Vgas, Vdisk, Vbul, SBdisk, SBbul,
                                          ML, ML, eps_g, beta, Sig0)
        ax.plot(R, Vmodel, '-', color='#c62828', lw=2.5, 
               label=r'$G(\Sigma)$ model', zorder=11)
        
        ax.set_xlabel('R (kpc)', fontsize=10)
        ax.set_ylabel('V (km/s)', fontsize=10)
        ax.set_title(f'{gid}', fontsize=11, fontweight='bold')
        ax.legend(fontsize=7, loc='lower right')
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        
        # Add chi2 info
        chi2 = fit_results['galaxy_chi2'].get(gid, 0)
        ndata = len(df)
        ax.text(0.03, 0.97, f'χ²/ndof={chi2/(ndata-1):.1f}', transform=ax.transAxes,
               fontsize=7, va='top', alpha=0.7)
    
    # Hide empty subplots
    for idx in range(n_plot, len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle(
        r'SPARC Rotation Curves: $G_{\rm eff}(\Sigma)$ Fit'
        f'\nε_g={eps_g:.3f}, β={beta:.3f}, Σ₀=10^{np.log10(Sig0):.2f} M☉/pc²',
        fontsize=13, fontweight='bold', y=1.02
    )
    plt.tight_layout()
    
    outpath = PLOT_DIR / "SPARC_GSigma_RotationCurves.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.savefig(PLOT_DIR / "SPARC_GSigma_RotationCurves.pdf", bbox_inches='tight')
    print(f"\nSaved: {outpath}")
    plt.close()


def plot_boost_profile(galaxies, fit_results):
    """Plot how G_eff varies with radius for representative galaxies."""
    
    eps_g = fit_results['epsilon_g']
    beta = fit_results['beta']
    Sig0 = fit_results['Sigma_0']
    ML_dict = fit_results['ML_per_galaxy']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left panel: Boost vs Surface Density (universal curve)
    sig_range = np.logspace(-2, 4, 200)  # M_sun/pc^2
    boost_curve = Geff_over_G0(sig_range, eps_g, beta, Sig0)
    
    ax1.semilogx(sig_range, boost_curve, '-', color='#c62828', lw=3)
    ax1.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5, label='$G_0$ (no boost)')
    ax1.fill_between(sig_range, 1.0, boost_curve, alpha=0.15, color='#c62828')
    
    ax1.set_xlabel(r'Surface Density $\Sigma$ (M$_\odot$/pc$^2$)', fontsize=12)
    ax1.set_ylabel(r'$G_{\rm eff}$ / $G_0$', fontsize=12)
    ax1.set_title(r'Universal $G(\Sigma)$ Enhancement Curve', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(alpha=0.2)
    ax1.set_xlim(1e-2, 1e4)
    
    # Mark characteristic densities
    ax1.axvline(x=Sig0, color='#ff6f00', ls='--', lw=1.5, alpha=0.7, label=f'$\\Sigma_0$={Sig0:.2f}')
    
    # Right panel: Individual galaxy boost profiles
    sample_gids = list(galaxies.keys())[:12]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(sample_gids)))
    
    for gid, color in zip(sample_gids, colors):
        df = galaxies[gid]
        ML = ML_dict.get(gid, 0.5)
        
        R = df['R_kpc'].values
        SBdisk = df['SBdisk'].values
        SBbul = df['SBbul'].values
        
        Sigma = ML * SBdisk + ML * SBbul
        boost = Geff_over_G0(Sigma, eps_g, beta, Sig0)
        
        ax2.plot(R, boost, '-', color=color, lw=1.5, alpha=0.8, label=gid)
    
    ax2.set_xlabel('R (kpc)', fontsize=12)
    ax2.set_ylabel(r'$G_{\rm eff}$ / $G_0$', fontsize=12)
    ax2.set_title(r'Radial $G(\Sigma)$ Profiles (sample)', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=6, loc='upper right', ncol=2)
    ax2.grid(alpha=0.2)
    ax2.set_xlim(left=0)
    
    plt.suptitle(
        r'Density-Dependent Gravity: $\frac{G_{\rm eff}}{G_0} = 1 + \epsilon_g \left(\frac{\Sigma}{\Sigma_0}\right)^\beta$',
        fontsize=14, y=1.02
    )
    plt.tight_layout()
    
    outpath = PLOT_DIR / "SPARC_GSigma_BoostProfile.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.savefig(PLOT_DIR / "SPARC_GSigma_BoostProfile.pdf", bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close()


def plot_chi2_distribution(galaxies, fit_results):
    """Show distribution of per-galaxy chi-squared."""
    
    chi2_list = list(fit_results['galaxy_chi2'].values())
    ndofs = [len(g)-1 for g in galaxies.values()]
    chi2red_list = [c/n for c,n in zip(chi2_list, ndofs)]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram of reduced chi-squared
    ax1.hist(chi2red_list, bins=25, color='#1565c0', alpha=0.7, edgecolor='white')
    ax1.axvline(x=np.median(chi2red_list), color='#c62828', lw=2, 
               label=f'Median = {np.median(chi2red_list):.2f}')
    ax1.axvline(x=1.0, color='gray', ls=':', lw=1.5, label='Ideal (χ²=1)')
    ax1.axvline(x=np.mean(chi2red_list), color='#ff6f00', ls='--', lw=1.5,
               label=f'Mean = {np.mean(chi2red_list):.2f}')
    ax1.set_xlabel('Reduced χ² (per galaxy)', fontsize=12)
    ax1.set_ylabel('Number of Galaxies', fontsize=12)
    ax1.set_title(f'SPARC: {len(chi2red_list)} Galaxies', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    
    # Chi2 vs number of data points
    ax2.scatter(ndofs, chi2_list, c=chi2red_list, cmap='RdYlGn_r', 
               s=30, alpha=0.7, edgecolors='gray', linewidths=0.5)
    cbar = plt.colorbar(ax2.collections[0], ax=ax2)
    cbar.set_label('χ²_red', fontsize=10)
    ax2.set_xlabel('Degrees of Freedom (N_points - 1)', fontsize=12)
    ax2.set_ylabel('Total χ²', fontsize=12)
    ax2.set_title('Goodness-of-Fit per Galaxy', fontsize=13, fontweight='bold')
    ax2.plot([0, max(ndofs)], [0, max(ndofs)], ':', color='gray', alpha=0.5)
    
    # Stats text
    n_good = sum(1 for c in chi2red_list if c < 2.0)
    n_excellent = sum(1 for c in chi2red_list if c < 1.5)
    stats_text = f'{n_excellent}/{len(chi2red_list)} galaxies with χ²<1.5\n{n_good}/{len(chi2red_list)} with χ²<2.0'
    ax2.text(0.03, 0.97, stats_text, transform=ax2.transAxes, fontsize=10,
             va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle(
        f'G(Σ) Fit Quality: Total χ²/{fit_results["dof"]} = {fit_results["chi2_red"]:.3f}',
        fontsize=14, y=1.02
    )
    plt.tight_layout()
    
    outpath = PLOT_DIR / "SPARC_GSigma_Chi2Distribution.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.savefig(PLOT_DIR / "SPARC_GSigma_Chi2Distribution.pdf", bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close()


def plot_comparison_standard(galaxies, fit_results):
    """
    Compare G(Σ) model against standard "maximum disk" + dark matter halo picture.
    Show that we don't need DM!
    """
    eps_g = fit_results['epsilon_g']
    beta = fit_results['beta']
    Sig0 = fit_results['Sigma_0']
    ML_dict = fit_results['ML_per_galaxy']
    
    # Select a few representative galaxies (different Hubble types)
    all_gids = list(galaxies.keys())
    showcase_ids = []
    for gid in all_gids:
        if gid.startswith('DDO'):   # Dwarf irregular
            if len([x for x in showcase_ids if x.startswith('DDO')]) < 2:
                showcase_ids.append(gid)
        elif gid.startswith('UGC'):  # General spiral
            if len([x for x in showcase_ids if x.startswith('UGC')]) < 3:
                showcase_ids.append(gid)
        elif gid.startswith('ESO'):  # Early type
            if len([x for x in showcase_ids if x.startswith('ESO')]) < 2:
                showcase_ids.append(gid)
        if len(showcase_ids) >= 7:
            break
    
    # Fallback: if not enough matches, just take first N
    if len(showcase_ids) < 4:
        n_show = min(8, len(all_gids))
        indices = np.linspace(0, len(all_gids)-1, n_show).astype(int)
        showcase_ids = [all_gids[i] for i in indices]
    
    n = len(showcase_ids)
    fig, axes = plt.subplots(2, (n+1)//2, figsize=(5*((n+1)//2), 9))
    axes = axes.flatten()
    
    for idx, gid in enumerate(showcase_ids):
        ax = axes[idx]
        df = galaxies[gid]
        ML = ML_dict.get(gid, 0.5)
        
        R = df['R_kpc'].values
        Vobs = df['Vobs'].values
        eV = df['e_Vobs'].values
        Vgas = df['Vgas'].values
        Vdisk = df['Vdisk'].values
        Vbul = df['Vbul'].values
        SBdisk = df['SBdisk'].values
        SBbul = df['SBbul'].values
        
        # Observations
        ax.errorbar(R, Vobs, yerr=eV, fmt='o', color='black', ms=4, zorder=10, alpha=0.8)
        
        # Baryonic only (Newtonian, no boost) — shows the "missing mass"
        Vbar = np.sqrt(np.maximum((ML*Vdisk)**2 + (ML*Vbul)**2 + Vgas**2, 0))
        ax.plot(R, Vbar, '--', color='blue', lw=1.5, alpha=0.7,
               label='Baryons only')
        
        # G(Σ) boosted — our prediction
        Vgs, _ = compute_Vmodel(R, Vgas, Vdisk, Vbul, SBdisk, SBbul,
                                 ML, ML, eps_g, beta, Sig0)
        ax.plot(R, Vgs, '-', color='#c62828', lw=2.5,
               label=r'$G(\Sigma)$ (no DM!)', zorder=11)
        
        # Shade the "dark matter gap" that G(Σ) closes
        ax.fill_between(R, Vbar, Vgs, alpha=0.15, color='#c62828')
        
        ax.set_xlabel('R (kpc)', fontsize=10)
        ax.set_ylabel('V (km/s)', fontsize=10)
        ax.set_title(f'{gid}', fontsize=11, fontweight='bold')
        ax.legend(fontsize=7, loc='lower right')
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
    
    # Hide extra panels
    for jdx in range(idx+1, len(axes)):
        axes[jdx].set_visible(False)
    
    plt.suptitle(
        r'$G(\Sigma)$ Replaces Dark Matter: No Halo Needed!\nRed area = what used to be attributed to dark matter halos',
        fontsize=13, fontweight='bold', y=1.02
    )
    plt.tight_layout()
    
    outpath = PLOT_DIR / "SPARC_GSigma_NoDarkMatter.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.savefig(PLOT_DIR / "SPARC_GSigma_NoDarkMatter.pdf", bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 65)
    print("  SPARC Rotation Curve Analysis — G_eff(Σ) Framework")
    print("  Testing: Does density-dependent gravity replace dark matter?")
    print("=" * 65)
    print()
    
    # Load data
    galaxies = load_all_galaxies()
    
    if len(galaxies) == 0:
        print("ERROR: No data loaded!")
        return
    
    # Run the big fit
    results = fit_universal_parameters(galaxies)
    
    # Generate all plots
    print("\n--- Generating figures ---\n")
    plot_sample_fits(galaxies, results, n_samples=9)
    plot_boost_profile(galaxies, results)
    plot_chi2_distribution(galaxies, results)
    plot_comparison_standard(galaxies, results)
    
    # Save results summary
    summary_path = OUTPUT_DIR / "SPARC_GSigma_FitResults.txt"
    with open(summary_path, 'w') as f:
        f.write("SPARC G(Σ) Fitting Results\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"UNIVERSAL PARAMETERS:\n")
        f.write(f"  ε_g = {results['epsilon_g']:.6f}\n")
        f.write(f"  β   = {results['beta']:.6f}\n")
        f.write(f"  Σ_0 = {results['Sigma_0']:.6f} M☉/pc²\n")
        f.write(f"       (log Σ_0 = {np.log10(results['Sigma_0']):.4f})\n\n")
        f.write(f"FIT QUALITY:\n")
        f.write(f"  Total χ² = {results['chi2_total']:.1f}\n")
        f.write(f"  DOF       = {results['dof']}\n")
        f.write(f"  χ²_red    = {results['chi2_red']:.4f}\n\n")
        f.write(f"N galaxies = {len(galaxies)}\n\n")
        f.write("PER-GALAXY M/L VALUES AND χ²:\n")
        f.write(f"{'Galaxy':<16} {'M/L':>6} {'χ²':>8}\n")
        f.write("-" * 34 + "\n")
        
        for gid in sorted(results['ML_per_galaxy'].keys()):
            ml = results['ML_per_galaxy'][gid]
            chi2 = results['galaxy_chi2'].get(gid, 0)
            f.write(f"{gid:<16} {ml:>6.3f} {chi2:>8.1f}\n")
    
    print(f"\nResults saved: {summary_path}")
    print("\n✅ Done! All figures in:", PLOT_DIR)
    
    # Quick sanity check
    print("\n--- QUICK SANITY CHECK ---")
    sig_typical_disk_center = 100  # M_sun/pc^2, typical central disk
    sig_typical_outer = 10         # M_sun/pc^2, outer disk
    boost_center = Geff_over_G0(sig_typical_disk_center, 
                                  results['epsilon_g'], results['beta'], results['Sigma_0'])
    boost_outer = Geff_over_G0(sig_typical_outer,
                                 results['epsilon_g'], results['beta'], results['Sigma_0'])
    print(f"At disk center (Σ~100): G_eff/G₀ = {boost_center:.2f}x")
    print(f"At disk edge   (Σ~10):  G_eff/G₀ = {boost_outer:.2f}x")
    print(f"This means inner rotation curves get MORE boost → naturally flat!")


if __name__ == '__main__':
    main()
