#!/usr/bin/env python3
"""
MCMC Bayesian inference for G_eff(Sigma) framework — v2 (reparameterized, flux-based Sigma)

Model: color = c0 + c1 * (Sigma/Sigma0)^beta
- c0: baseline color offset
- c1: amplitude of Sigma dependence (= A_c * eps_g in the original parameterization)
- beta: power-law index (fixed to literature value 0.097)

Sigma definition (matches fit_GSigma_params.py Path-A):
  logSigma = log10(f444w_flux) - 2*log10(r_eff_phys)
  Sigma = 10^(logSigma - median(logSigma))  [centered]

This reparameterization eliminates the A_c × eps_g degeneracy present in v1.
"""

import numpy as np
import pandas as pd
import emcee
import corner
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.optimize import minimize
import json
import time
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──
BASE = "/Users/tanxin/Desktop/数据处理/重新整理投稿APJ"
KOKOREV = f"{BASE}/源数据/Kokorev_LRDs_Full.csv"
GPARAMS = f"{BASE}/源数据/GSigma_Parameters.csv"
OUTDIR = f"{BASE}/图表/robustness_analysis"

# ── Load data (same as fit_GSigma_params.py) ──
print("=" * 60)
print("MCMC v2: Bayesian G_eff(Sigma) parameter inference")
print("=" * 60)

dfA = pd.read_csv(KOKOREV)
f444 = dfA['f444w_flux'].values.astype(float)
f150 = dfA['f150w_flux'].values.astype(float)
reff = dfA['r_eff_50_phys'].values.astype(float)

# Compute magnitudes and color
m444 = -2.5 * np.log10(np.clip(f444, 1e-10, None)) + 23.9
m150 = -2.5 * np.log10(np.clip(f150, 1e-10, None)) + 23.9
color = m444 - m150

# Compute Sigma (flux-based, centered) — identical to fit_GSigma_params.py
lS_raw = np.log10(np.clip(f444, 1e-10, None)) - 2*np.log10(np.clip(np.abs(reff), 1e-5, None))
lS_centered = lS_raw - np.nanmedian(lS_raw)
Sigma_centered = 10 ** lS_centered

# Quality mask
mask = (np.isfinite(color) & np.isfinite(lS_centered) & (f444 > 0) & (f150 > 0) & (reff > 0))
c = color[mask]
S = Sigma_centered[mask]
lS = lS_centered[mask]
N = len(c)

print(f"Sources: {len(dfA)} total, {N} valid")
print(f"Color range: [{c.min():.2f}, {c.max():.2f}] mag")
print(f"log(Sigma/Sigma0) range: [{lS.min():.2f}, {lS.max():.2f}]")

# ── Fixed beta from prior fit ──
gp = pd.read_csv(GPARAMS)
beta_fixed = float(gp.loc[gp['method'] == 'Path-A(Color)', 'beta'].values[0])
eps_g_lit = float(gp.loc[gp['method'] == 'Path-A(Color)', 'eps_g'].values[0])
print(f"beta_fixed = {beta_fixed:.4f} (from Path-A fitting)")
print(f"eps_g (literature) = {eps_g_lit:.4f}")

# ── Raw Spearman (reproduce the 5.6σ signal) ──
rho_raw, p_raw = spearmanr(S, c)
print(f"\nRaw Spearman rho(Sigma, color) = {rho_raw:.4f} (p = {p_raw:.2e})")
sigma_z = abs(rho_raw) * np.sqrt(N)
print(f"Significance (Gaussian approx): {sigma_z:.1f} sigma")

# ── Model definition ──
def model(params, S_arr, beta):
    c0, c1 = params
    return c0 + c1 * S_arr ** beta

def log_likelihood(params, c_obs, S_arr, beta, sigma):
    c0, c1 = params
    if c1 > 50 or c1 < -50:  # soft bound
        return -np.inf
    mu = model(params, S_arr, beta)
    return -0.5 * np.sum(((c_obs - mu) / sigma) ** 2)

def log_prior(params):
    c0, c1 = params
    if not (-10 < c0 < 10):
        return -np.inf
    if not (-20 < c1 < 20):
        return -np.inf
    return 0.0

def log_posterior(params, c_obs, S_arr, beta, sigma):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, c_obs, S_arr, beta, sigma)

# ── Noise estimate ──
sigma_data = 1.4826 * np.median(np.abs(c - np.median(c)))
print(f"MAD-based sigma_color = {sigma_data:.4f}")

# ── MLE via Nelder-Mead ──
def neg_log_post(params):
    return -log_posterior(params, c, S, beta_fixed, sigma_data)

best_result = None
best_val = np.inf
for c0_init in [-2.0, -1.5, -2.5, -1.8]:
    for c1_init in [0.05, 0.1, 0.2, 0.5, 1.0]:
        try:
            res = minimize(neg_log_post, [c0_init, c1_init], method='Nelder-Mead',
                          options={'maxiter': 20000, 'xatol': 1e-8, 'fatol': 1e-8})
            if res.fun < best_val:
                best_val = res.fun
                best_result = res
        except:
            pass

p0_opt = best_result.x
print(f"\nMLE (Nelder-Mead): c0 = {p0_opt[0]:.4f}, c1 = {p0_opt[1]:.4f}")

# ── MCMC ──
ndim = 2
nwalkers = 32
nsteps = 5000
nburn = 1000

np.random.seed(42)
p0_walkers = p0_opt + 1e-3 * np.random.randn(nwalkers, ndim)
for i in range(nwalkers):
    while not (-10 < p0_walkers[i, 0] < 10) or not (-20 < p0_walkers[i, 1] < 20):
        p0_walkers[i] = p0_opt + 1e-3 * np.random.randn(ndim)

print(f"\nMCMC: {nwalkers} walkers × {nsteps} steps, {nburn} burn-in")
print(f"Posterior samples: {nwalkers * (nsteps - nburn)}")
print("Running...")

t0 = time.time()
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,
                                args=(c, S, beta_fixed, sigma_data))
sampler.run_mcmc(p0_walkers, nsteps, progress=True)
elapsed = time.time() - t0
print(f"MCMC completed in {elapsed:.1f}s")

# ── Discard burn-in ──
samples = sampler.get_chain(discard=nburn, thin=1, flat=True)
print(f"Posterior shape: {samples.shape}")

# ── Constraints ──
labels = ['c0', 'c1']
results_summary = {}
for i, label in enumerate(labels):
    q = np.percentile(samples[:, i], [16, 50, 84])
    results_summary[label] = {
        'median': float(q[1]),
        'lo': float(q[0]),
        'hi': float(q[2]),
        'err_lo': float(q[1] - q[0]),
        'err_hi': float(q[2] - q[1])
    }
    print(f"  {label} = {q[1]:.4f} (+{q[2]-q[1]:.4f}, -{q[1]-q[0]:.4f})")

# ── Derived: implied eps_g = c1 / A_c (where A_c is the original amplitude) ──
# Original model: color = A_c * eps_g * S^beta → c1 = A_c * eps_g
# We can't separate them without additional info, but c1 itself is the physically relevant quantity
c1_med = results_summary['c1']['median']

# ── Model-observation correlation ──
model_best = model(p0_opt, S, beta_fixed)
rho_model, p_model = spearmanr(model_best, c)
print(f"\nSpearman rho(model, observed) = {rho_model:.4f} (p = {p_model:.2e})")

# Posterior predictive: Spearman for 1000 random samples
np.random.seed(123)
idx = np.random.choice(len(samples), min(2000, len(samples)), replace=False)
spearman_posterior = []
for i in idx:
    m_i = model(samples[i], S, beta_fixed)
    r, _ = spearmanr(m_i, c)
    spearman_posterior.append(r)
spearman_posterior = np.array(spearman_posterior)
print(f"Posterior Spearman rho: {np.median(spearman_posterior):.4f} "
      f"[{np.percentile(spearman_posterior, 5):.4f}, {np.percentile(spearman_posterior, 95):.4f}]")

# ── Goodness of fit ──
chi2 = np.sum(((c - model_best) / sigma_data) ** 2)
chi2_red = chi2 / (N - ndim)
print(f"Chi-squared / dof = {chi2_red:.2f}")

# Residual trend check
residuals = c - model_best
rho_res, p_res = spearmanr(S, residuals)
print(f"Spearman rho(Sigma, residual) = {rho_res:.4f} (p = {p_res:.2e})")

# ── Save results ──
results = {
    'version': 'v2_reparameterized_fluxSigma',
    'status': 'completed',
    'model': 'color = c0 + c1 * (Sigma/Sigma0)^beta',
    'sigma_definition': 'log10(f444w_flux) - 2*log10(reff_phys), centered by median',
    'n_sources': N,
    'beta_fixed': float(beta_fixed),
    'eps_g_literature': float(eps_g_lit),
    'sigma_data': float(sigma_data),
    'nwalkers': nwalkers,
    'nsteps': nsteps,
    'nburn': nburn,
    'mle': {'c0': float(p0_opt[0]), 'c1': float(p0_opt[1])},
    'best_fit': results_summary,
    'validation': {
        'spearman_raw': {'rho': float(rho_raw), 'p': float(p_raw), 'sigma_significance': float(sigma_z)},
        'spearman_model_vs_obs': {'rho': float(rho_model), 'p': float(p_model)},
        'posterior_spearman': {
            'median': float(np.median(spearman_posterior)),
            '5th': float(np.percentile(spearman_posterior, 5)),
            '95th': float(np.percentile(spearman_posterior, 95))
        },
        'chi2_red': float(chi2_red),
        'residual_trend': {'rho': float(rho_res), 'p': float(p_res)}
    },
    'summary': {
        'c0': f"{results_summary['c0']['median']:.4f} (+{results_summary['c0']['err_hi']:.4f}, -{results_summary['c0']['err_lo']:.4f})",
        'c1': f"{results_summary['c1']['median']:.4f} (+{results_summary['c1']['err_hi']:.4f}, -{results_summary['c1']['err_lo']:.4f})",
        'spearman_raw': f"{rho_raw:.3f} ({sigma_z:.1f}σ)",
        'chi2_red': f"{chi2_red:.2f}",
        'elapsed': f"{elapsed:.1f}s"
    }
}

with open(f"{OUTDIR}/mcmc_results.json", 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults → {OUTDIR}/mcmc_results.json")

# ── Figure 1: Corner plot ──
fig = corner.corner(samples, labels=[r'$c_0$', r'$c_1$'],
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_fmt='.4f',
                    title_kwargs={'fontsize': 13},
                    label_kwargs={'fontsize': 14},
                    hist_kwargs={'density': True, 'color': '#2196F3', 'alpha': 0.7},
                    color='#2196F3')

# MLE markers
for i_ax in [0, 2]:
    fig.axes[i_ax].axvline(p0_opt[i_ax // 2], color='red', ls='--', lw=1.5, alpha=0.8)

fig.suptitle(r'MCMC Posterior: $\mathrm{color} = c_0 + c_1\,(\Sigma/\Sigma_0)^{\,\beta}$'
             f'\n$\\beta$ = {beta_fixed:.3f} (fixed), N = {N}',
             fontsize=14, y=1.06)

fig.savefig(f"{OUTDIR}/Fig_mcmc_corner_GSigma.png", dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close(fig)
print(f"Corner plot → {OUTDIR}/Fig_mcmc_corner_GSigma.png")

# ── Figure 2: Model fit + residuals ──
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Data + model
ax1 = axes[0]
si = np.argsort(lS)
lSs = lS[si]

# Posterior draw band
idx_draw = np.random.choice(len(samples), 200, replace=False)
for i in idx_draw:
    y_draw = model(samples[i], S[si], beta_fixed)
    ax1.plot(lSs, y_draw, color='#2196F3', alpha=0.03, lw=0.5)

# Best fit
y_best = model(p0_opt, S[si], beta_fixed)
ax1.plot(lSs, y_best, 'r-', lw=2.5, label='MLE fit', zorder=5)

# Data
ax1.scatter(lS, c, c='#2166ac', alpha=0.55, s=22, edgecolors='none', zorder=3,
            label=f'LRDs (N={N})')

ax1.axhline(y=np.median(c), color='gray', ls='--', lw=1, alpha=0.5)
ax1.axvline(x=0, color='gray', ls=':', lw=1, alpha=0.7)
ax1.set_xlabel(r'$\log(\Sigma / \Sigma_0)$', fontsize=13)
ax1.set_ylabel(r'$m_{\rm F444W} - m_{\rm F150W}$ [mag]', fontsize=13)
ax1.set_title(f'Model Fit to {N} LRDs', fontsize=13, fontweight='bold')
ax1.legend(fontsize=10, loc='lower left')

# Formula + stats box
c0s = results_summary['c0']
c1s = results_summary['c1']
box_text = (
    f"$c_0$ = {c0s['median']:.3f} $^{{+{c0s['err_hi']:.3f}}}_{{-{c0s['err_lo']:.3f}}}$\n"
    f"$c_1$ = {c1s['median']:.3f} $^{{+{c1s['err_hi']:.3f}}}_{{-{c1s['err_lo']:.3f}}}$\n"
    f"$\\beta$ = {beta_fixed:.3f} (fixed)\n"
    f"$\\rho$ = {rho_raw:.3f} ({sigma_z:.1f}$\\sigma$)"
)
ax1.text(0.03, 0.97, box_text, transform=ax1.transAxes, fontsize=10, va='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))

# Panel B: Residuals
ax2 = axes[1]
ax2.scatter(model_best, residuals, c='#2166ac', alpha=0.5, s=18, edgecolors='none')
ax2.axhline(0, color='red', ls='-', lw=1.5)
ax2.axhline(-sigma_data, color='gray', ls=':', lw=1, alpha=0.5)
ax2.axhline(sigma_data, color='gray', ls=':', lw=1, alpha=0.5)
ax2.set_xlabel('Predicted color [mag]', fontsize=13)
ax2.set_ylabel('Residual (obs − pred) [mag]', fontsize=13)
ax2.set_title(f'Residuals ($\\chi^2_{{\\rm red}}$ = {chi2_red:.2f})', fontsize=13, fontweight='bold')

res_text = f"$\\rho$(resid, $\\Sigma$) = {rho_res:.3f}\np = {p_res:.2e}"
ax2.text(0.03, 0.03, res_text, transform=ax2.transAxes, fontsize=10, va='bottom',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.85))

plt.tight_layout()
fig.savefig(f"{OUTDIR}/Fig_mcmc_model_fit.png", dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close(fig)
print(f"Model fit → {OUTDIR}/Fig_mcmc_model_fit.png")

# ── Figure 3: Posterior predictive Spearman distribution ──
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(spearman_posterior, bins=50, color='#2196F3', alpha=0.7, density=True,
        edgecolor='white')
ax.axvline(np.median(spearman_posterior), color='red', lw=2, ls='-',
           label=f'Posterior median = {np.median(spearman_posterior):.3f}')
ax.axvline(rho_raw, color='orange', lw=2, ls='--',
           label=f'Direct Spearman = {rho_raw:.3f}')
ax.set_xlabel(r'Spearman $\rho$(model, observed)', fontsize=13)
ax.set_ylabel('Density', fontsize=13)
ax.set_title('Posterior Predictive: Model-Observation Correlation', fontsize=13)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(f"{OUTDIR}/Fig_mcmc_spearman_posterior.png", dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close(fig)
print(f"Spearman posterior → {OUTDIR}/Fig_mcmc_spearman_posterior.png")

print("\n" + "=" * 60)
print("MCMC v2 COMPLETE — ALL OUTPUTS SAVED")
print("=" * 60)
