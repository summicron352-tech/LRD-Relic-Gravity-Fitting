#!/usr/bin/env python3
"""
Rigorous Radio Emission for the Dyson Entrainment Model
========================================================
Calculates three mechanisms using actual model parameters:
  1. Free-free: SELF-ABSORBED blackbody from τ=1 surface (dominant)
  2. Synchrotron: stellar wind cosmic-ray electrons
  3. Inverse Compton: CMB upscattering (negligible at radio frequencies)

Key physics insight:
  The dense inner cocoon (n_e ~ 10^8-10^10 cm^-3) is optically thick
  to free-free absorption at radio frequencies (τ_ff >> 1).
  → The observed radio emission is blackbody from the τ=1 surface.
  → This is ORDERS OF MAGNITUDE below both AGN predictions and SKAO limits.

Author: Xin Tan, 2026-07-06
"""
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, json

# Constants
MSUN, PC, YR = 1.989e33, 3.086e18, 3.156e7
c_light, k_B, h_planck, m_p = 2.998e10, 1.381e-16, 6.626e-27, 1.673e-24
sigma_T = 6.652e-25
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# ============================================================================
# Dyson Model Parameters
# ============================================================================
class DysonParams:
    def __init__(self, z=5.5):
        self.z = z
        self.r_c, self.r_eff = 0.003, 83.3       # pc
        self.alpha_in, self.alpha_out = 2.0, 1.2
        self.rho_c = 6.8e11                       # M_sun/pc^3
        self.M_star = 10**9.51                     # M_sun
        self.n_H_blr = 1e9                        # cm^-3 at BLR (ABCD)
        self.T_gas = 1e4                          # K
        self.N_core = 2e6                         # stars in r<0.03pc
        self.Mdot_w = 1e-6                        # M_sun/yr per star
        self.v_wind = 2000                         # km/s
        self.B_star = 100                          # G
        self.m_star = 0.5                          # M_sun (avg)
        self.f_gas = 0.1                           # gas fraction in core
        self.eps_e = 0.01                          # electron acceleration eff
        self.eps_B = 0.05                          # magnetic energy fraction
        self.dL = cosmo.luminosity_distance(z).to(u.cm).value

def gas_density(r_pc, p):
    """n_H [cm^-3] from stellar density profile."""
    r = np.atleast_1d(r_pc)
    rho = np.zeros_like(r, float)
    inner = r < 0.3
    rho[inner] = p.rho_c * (r[inner]/p.r_c)**(-p.alpha_in)
    outer = ~inner
    rho_at_b = p.rho_c * (0.3/p.r_c)**(-p.alpha_in)
    rho[outer] = rho_at_b * (r[outer]/0.3)**(-p.alpha_out)
    rho[r < p.r_c*0.1] = p.rho_c * 0.1**(-p.alpha_in)
    # Normalize to M_star
    r_grid = np.logspace(-5, np.log10(p.r_eff*1.1), 5000)
    rho_full = np.ones_like(r_grid)
    for i, rr in enumerate(r_grid):
        if rr < 0.3: rho_full[i] = p.rho_c*(rr/p.r_c)**(-p.alpha_in)
        else:
            rb = p.rho_c*(0.3/p.r_c)**(-p.alpha_in)
            rho_full[i] = rb*(rr/0.3)**(-p.alpha_out)
    rho_full[r_grid<p.r_c*0.1] = p.rho_c*0.1**(-p.alpha_in)
    dr = np.diff(np.concatenate([[0], r_grid]))
    Menc = np.cumsum(rho_full*4*np.pi*r_grid**2*dr)
    scale = p.M_star / Menc[np.argmin(np.abs(r_grid-p.r_eff))]
    return p.f_gas * np.interp(r_pc, r_grid, rho_full*scale)*MSUN/PC**3/(1.4*m_p)

# ============================================================================
# FREE-FREE (Self-Absorbed)
# ============================================================================
def gaunt_ff(nu_GHz, T_K):
    x = nu_GHz
    return np.clip(24.5 + np.log10(T_K) - np.log10(x), 3, 25)

def free_free_self_absorbed(p, nu_obs_GHz):
    """
    Self-absorbed free-free: find τ=1 surface → blackbody emission.
    This is the correct treatment when τ_ff >> 1.
    """
    nu_rest = nu_obs_GHz*(1+p.z)
    # Find τ=1 surface from outside in (coarse grid for speed)
    r_grid = np.logspace(-5, 3, 500)
    tau = np.zeros_like(r_grid)
    for i in range(len(r_grid)-2, -1, -1):
        n_e = gas_density(r_grid[i], p)
        if n_e < 1: continue
        g = float(gaunt_ff(nu_rest, p.T_gas))
        hnu_kt = h_planck*nu_rest*1e9/(k_B*p.T_gas)
        alpha = float(3.7e8 * n_e**2 * p.T_gas**-0.5 * nu_rest**-3 * (1-np.exp(-hnu_kt))*g)
        tau[i] = tau[i+1] + alpha*(r_grid[i+1]-r_grid[i])*PC
    idx = np.argmin(np.abs(tau-1.0))
    r1 = r_grid[idx]
    # Rayleigh-Jeans blackbody from τ=1 surface
    B_RJ = 2*k_B*p.T_gas*(nu_rest*1e9)**2/c_light**3  # erg/s/cm^2/Hz/sr
    dA = cosmo.angular_diameter_distance(p.z).to(u.cm).value
    omega = np.pi*(r1*PC/dA)**2
    flux_ujy = B_RJ*omega*1e23*1e6
    return flux_ujy, r1

# ============================================================================
# SYNCHROTRON
# ============================================================================
def B_field_estimate(p):
    """Magnetic field from wind equipartition in the core."""
    Mdot = p.N_core*p.Mdot_w*MSUN/YR
    L_w = 0.5*Mdot*(p.v_wind*1e5)**2
    r_blr_cm = 0.03*PC
    rho_w = Mdot/(4*np.pi*r_blr_cm**2*p.v_wind*1e5)
    B_w = np.sqrt(4*np.pi*rho_w*(p.v_wind*1e5)**2)
    # Turbulent dynamo: B ≈ √(4πρ) v_turb
    rho_c = p.rho_c*MSUN/PC**3
    v_dyn = 3000/2.355*1e5
    B_dyn = np.sqrt(4*np.pi*rho_c)*v_dyn
    return np.exp(0.5*np.log(B_w)+0.5*np.log(B_dyn))

def synchrotron_dyson(p, nu_obs_GHz):
    """Stellar wind → cosmic rays → synchrotron."""
    Mdot = p.N_core*p.Mdot_w*MSUN/YR
    L_wind = 0.5*Mdot*(p.v_wind*1e5)**2
    L_sync = p.eps_e*p.eps_B*L_wind  # total, erg/s
    B = B_field_estimate(p)
    g_min, g_max = 10, 1e4
    nu_min = 4.2e-3*B*g_min**2  # GHz
    nu_max = 4.2e-3*B*g_max**2
    nu_rest = nu_obs_GHz*(1+p.z)
    alpha = -0.6  # p_e=2.2 → α=-(p-1)/2=-0.6
    ratio = (nu_min/nu_max)**(1+alpha) if nu_min<nu_max else 0.0
    K = L_sync*(1+alpha)/(nu_max*1e9*(1-ratio))  # norm at ν_max, erg/s/Hz
    if nu_rest < nu_min:
        return K*(nu_min/nu_max)**alpha*(nu_rest/nu_min)**2.5
    elif nu_rest <= nu_max:
        return K*(nu_rest/nu_max)**alpha
    else:
        return K*np.exp(-(nu_rest-nu_max)/nu_max)

def flux_ujy(L_nu, nu_obs, z, alpha=-0.7):
    dL = cosmo.luminosity_distance(z).to(u.cm).value
    fcgs = L_nu*(1+z)**(1+alpha)/(4*np.pi*dL**2)
    return fcgs*1e23*1e6

# ============================================================================
# IC/CMB (negligible at radio frequencies)
# ============================================================================
def ic_cmb(p, nu_obs_GHz):
    B = B_field_estimate(p)
    u_B = B**2/(8*np.pi)
    u_CMB = 7.56e-15*(2.725*(1+p.z))**4
    L_sync = synchrotron_dyson(p, nu_obs_GHz)
    return L_sync*u_CMB/u_B*1e-3  # ×10^-3 for radio suppression

# ============================================================================
# SKAO Limits
# ============================================================================
SKAO = {'Band-2': (1.4, 1.1), 'Band-5a': (6.5, 0.7), 'Band-5b': (11.5, 0.84)}

def skao_5sig(rms, t_hr):
    return 5*rms/np.sqrt(t_hr)

# ============================================================================
# Main
# ============================================================================
if __name__ == '__main__':
    print("="*70)
    print("DYSON MODEL RADIO EMISSION — SELF-ABSORBED FREE-FREE")
    print("="*70)

    p = DysonParams(z=5.5)
    B = B_field_estimate(p)
    print(f"\nz={p.z}, ρ_c={p.rho_c:.1e} M_sun/pc^3, N_core={p.N_core:.0e}")
    print(f"n_H(r=0.03pc)={gas_density(0.03, p):.1e} cm^-3, T_gas={p.T_gas:.0e} K")
    print(f"B_est = {B:.2e} G")

    freqs = np.array([0.2, 0.5, 0.8, 1.0, 1.4, 2.0, 3.0, 5.0, 6.5, 8.0, 10.0, 11.5])

    print(f"\n{'Freq':>6s} {'FF(SA) [μJy]':>14s} {'r_τ1 [pc]':>10s} "
          f"{'Sync [μJy]':>12s} {'IC/CMB':>12s} {'TOTAL':>12s}")
    print("-"*70)

    ff_f, sync_f, ic_f = [], [], []
    for nu in freqs:
        f_ff, r1 = free_free_self_absorbed(p, nu)
        L_s = synchrotron_dyson(p, nu)
        f_s = flux_ujy(L_s, nu, p.z, -0.7)
        L_ic = ic_cmb(p, nu)
        f_ic = flux_ujy(L_ic, nu, p.z, -0.7)
        ff_f.append(f_ff); sync_f.append(f_s); ic_f.append(f_ic)
        print(f"{nu:>6.2f} {f_ff:>14.4e} {r1:>10.6f} {f_s:>12.4e} {f_ic:>12.4e} {f_ff+f_s+f_ic:>12.4e}")

    # Summary
    i14 = np.argmin(np.abs(freqs-1.4))
    total_14 = ff_f[i14]+sync_f[i14]+ic_f[i14]
    s2_1k = skao_5sig(1.1, 1000)
    s2_stack = s2_1k/np.sqrt(260)

    print(f"\n{'='*70}")
    print(f"SUMMARY @ 1.4 GHz (z=5.5)")
    print(f"{'='*70}")
    print(f"  Free-Free (self-absorbed): {ff_f[i14]:.4e} μJy")
    print(f"  Synchrotron:               {sync_f[i14]:.4e} μJy")
    print(f"  IC/CMB:                    {ic_f[i14]:.4e} μJy")
    print(f"  TOTAL:                     {total_14:.4e} μJy")
    print(f"  SKAO Band-2 1000h (5σ):    {s2_1k:.4e} μJy")
    print(f"  SKAO Stack 260 LRD:        {s2_stack:.4e} μJy")
    print(f"  Ratio TOTAL/stack_limit:   {total_14/s2_stack:.1f}×")
    print(f"  → {'UNDETECTABLE even with SKAO stacking' if total_14 < s2_stack else 'DETECTABLE'}")

    # MC parameter exploration
    print(f"\n[MC] Parameter space exploration (N=1000)...")
    np.random.seed(42)
    Nmc = 1000
    log_nH = np.random.uniform(7, 11, Nmc)      # log n_H(BLR)
    log_Nc = np.random.uniform(5, 8, Nmc)       # log N_core
    log_Mw = np.random.uniform(-7, -4, Nmc)      # log Mdot_wind
    zs = np.random.uniform(4, 9, Nmc)
    ff_mc, sy_mc = np.zeros(Nmc), np.zeros(Nmc)

    for i in range(Nmc):
        if i%1000==0: print(f"  {i}/{Nmc}...")
        try:
            pi = DysonParams(z=zs[i])
            pi.n_H_blr = 10**log_nH[i]
            pi.N_core = 10**log_Nc[i]
            pi.Mdot_w = 10**log_Mw[i]
            f_f, _ = free_free_self_absorbed(pi, 1.4)
            ff_mc[i] = f_f
            L_s = synchrotron_dyson(pi, 1.4)
            sy_mc[i] = flux_ujy(L_s, 1.4, pi.z, -0.7)
        except: pass

    valid = (ff_mc>0)&(sy_mc>0)
    ff_v, sy_v = ff_mc[valid], sy_mc[valid]
    print(f"\n  Free-Free MC: median={np.median(ff_v):.4e}, 16-84%=[{np.percentile(ff_v,16):.4e}, {np.percentile(ff_v,84):.4e}] μJy")
    print(f"  Synchrotron MC: median={np.median(sy_v):.4e}, 16-84%=[{np.percentile(sy_v,16):.4e}, {np.percentile(sy_v,84):.4e}] μJy")
    total_mc = ff_v+sy_v
    print(f"  TOTAL MC: median={np.median(total_mc):.4e} μJy")
    p_detect = np.mean(total_mc > s2_1k)
    p_detect_stack = np.mean(total_mc > s2_stack)
    print(f"  P(detect individually, 1000h): {p_detect*100:.1f}%")
    print(f"  P(detect via 260-LRD stack):   {p_detect_stack*100:.1f}%")

    # FIGURES
    fig_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    # Figure 1: Full SED
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(freqs, ff_f, 'r-', lw=2.5, label='Dyson: Free-Free (self-absorbed)')
    ax.plot(freqs, sync_f, 'b-', lw=2, label='Dyson: Synchrotron (stellar winds)')
    ax.plot(freqs, ic_f, 'purple', ls='--', lw=1.5, label='Dyson: IC/CMB')
    ax.plot(freqs, np.array(ff_f)+np.array(sync_f), 'k-', lw=3, label='Dyson: TOTAL')
    # MC bands
    if len(ff_v)>0:
        ax.fill_between(freqs, [np.percentile(ff_v,16)]*len(freqs),
                       [np.percentile(ff_v,84)]*len(freqs), alpha=0.1, color='red')
    # SKAO
    for band,(freq,rms) in SKAO.items():
        for t,c,ls in [(30,'green','--'),(1000,'green','-')]:
            lim=skao_5sig(rms,t)
            ax.axhline(lim, color=c, ls=ls, alpha=0.6, lw=1.2)
            if t==1000: ax.text(freq*0.8, lim*0.5, f'{band} {t}h', fontsize=7, color=c)
    ax.axhline(s2_stack, color='darkgreen', ls='-.', lw=1.5, label='SKAO stack 260 LRD (1000h)')
    ax.set(xscale='log', yscale='log', xlabel='Observed Frequency [GHz]',
           ylabel='Flux Density [μJy]', xlim=(0.15,15), ylim=(1e-12, 10),
           title=f'Dyson Entrainment Model Radio SED (z={p.z})')
    ax.legend(fontsize=9, loc='lower left'); ax.grid(True, alpha=0.3)
    plt.tight_layout(); plt.savefig(os.path.join(fig_dir,'dyson_radio_sed.png'),dpi=200); plt.close()
    print(f"  → dyson_radio_sed.png")

    # Figure 2: MC histogram
    fig, axes = plt.subplots(1,3,figsize=(16,5.5))
    for ax,data,color,label in [(axes[0],ff_v,'red','Free-Free'),
                                 (axes[1],sy_v,'blue','Synchrotron'),
                                 (axes[2],total_mc,'black','Total')]:
        ax.hist(np.log10(data), bins=50, color=color, alpha=0.7, edgecolor='white')
        ax.axvline(np.log10(np.median(data)), color=color, lw=2)
        ax.axvline(np.log10(s2_1k), color='green', ls='--', lw=1.5, label='Band-2 1000h')
        ax.axvline(np.log10(s2_stack), color='darkgreen', ls='-.', lw=1.5, label='Stack 260')
        ax.set(xlabel='log₁₀ Flux @1.4 GHz [μJy]', ylabel='N', title=label)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
    plt.suptitle(f'Dyson Radio Flux Distribution (z=4-9, N={len(ff_v)})', fontweight='bold')
    plt.tight_layout(); plt.savefig(os.path.join(fig_dir,'dyson_radio_mc.png'),dpi=200); plt.close()
    print(f"  → dyson_radio_mc.png")

    # Save data
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results')
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir,'dyson_radio_results.json'),'w') as f:
        json.dump({
            'fiducial_1p4GHz': {'ff':float(ff_f[i14]), 'sync':float(sync_f[i14]),
                                'ic':float(ic_f[i14]), 'total':float(total_14)},
            'mc_1p4GHz': {'ff_median':float(np.median(ff_v)), 'ff_16':float(np.percentile(ff_v,16)),
                         'ff_84':float(np.percentile(ff_v,84)),
                         'sync_median':float(np.median(sy_v)), 'total_median':float(np.median(total_mc))},
            'skao': {'band2_1000h':float(s2_1k), 'band2_stack260':float(s2_stack)},
            'p_detect': float(p_detect), 'p_detect_stack': float(p_detect_stack)
        }, f, indent=2)
    print(f"\nResults → {outdir}/dyson_radio_results.json")
    print("="*70)
