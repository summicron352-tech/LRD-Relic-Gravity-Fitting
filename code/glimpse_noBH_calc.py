"""
GLIMPSE-17775: Can G_eff(Σ) explain the broad-line FWHM without a central BH?
===========================================================
Physical idea:
  v^2(r) = G_eff(Σ) × M(<r) / r
  where M(<r) is enclosed stellar mass within radius r,
  G_eff = G_N × [1 + ε_g(Σ/Σ₀)^β]

  At small r (BLR scale ~0.01-0.1 pc), if stars are densely packed,
  the orbital velocity can be significant. We test whether G_eff enhancement
  can produce v ~ 1000-3000 km/s without a central point mass.
"""

import numpy as np

# ── Constants ──
G_N = 4.30091e-3          # [pc·(km/s)^2 / Msun]
MSUN = 1.989e33            # [g]
PC = 3.086e18              # [cm/pc]
YR = 3.156e7               # [s/yr]
kpc = 1e3                  # pc in kpc

# ── GLIMPSE-17775 parameters ──
logMstar = 9.51
Mstar = 10**logMstar       # [Msun]
r_eff = 83.3               # [pc]
z = 5.66

# Surface density
Sigma = Mstar / (np.pi * r_eff**2)
logSigma = np.log10(Sigma)
print("=" * 65)
print("GLIMPSE-17775: Basic Parameters")
print("=" * 65)
print(f"  M*     = {Mstar:.2e} Msun")
print(f"  r_eff  = {r_eff:.1f} pc")
print(f"  Σ      = {Sigma:.2e} Msun/pc^2")
print(f"  log Σ  = {logSigma:.2f}")

# ── G_eff(Σ) framework parameters (from COSMOSWeb calibration) ──
# From the UID framework: G_eff(Σ) = G_N [1 + ε_g (Σ/Σ₀)^β]
# Calibrated from Σ-color correlation strength
# Using COSMOSWeb results: at z=7-9, ρ~+0.516; Δmag_max~2.47 mag
# ε_g ~ 0.1-0.3 (depends on calibration), β ~ 0.5-1.0
# Σ₀ ~ 10^3 (typical pivot)

# Let me use a conservative calibration:
# At log Σ ~ 5.17 and z ~ 5.66, the effective G enhancement
# should be substantial based on the Δmag_max ~ 2.47 mag

# Parameter choices (conservative to aggressive)
epsilon_g = 0.15    # coupling strength
beta = 0.7          # density exponent
Sigma0 = 1e3        # pivot density [Msun/pc^2]

G_ratio = 1.0 + epsilon_g * (Sigma / Sigma0)**beta
G_eff = G_N * G_ratio

print(f"\n  G_eff/G_N = {G_ratio:.2f}")
print(f"  (ε_g={epsilon_g}, β={beta}, Σ₀=10³)")
print()

# ── Velocity from extended mass distribution (no BH) ──
# For a Sersic profile, the enclosed mass at radius r is:
# M(<r) = M* × γ(2n, b_n (r/r_eff)^(1/n)) / Γ(2n)
# where γ is the lower incomplete gamma function

# Simplified: for a Plummer/King profile, at small r:
# M(<r) ~ M* × (r/r_eff)^(3-n)  for n~1-4 Sersic
# At BLR scales (r << r_eff), the mass enclosed is very small

# Let's use a more physically motivated approach:
# The stellar mass density at small radius for a Sersic n=4 (de Vaucouleurs):
# ρ(r) ∝ r^(-p) where p ≈ 1 - 0.6097/n + 0.05463/n^2 ≈ 0.85 for n=4
# So at r << r_eff, M(<r) ∝ r^(3-p) ≈ r^(2.15)

# In nuclear star clusters and compact galaxies, the central density
# can be very high. Let's use the core density approach.

# Core density estimate from M* and r_eff (assuming Sersic n~4):
# ρ_0 ≈ M* / (8 × r_eff³)  (approximate central density)
rho0 = Mstar / (8 * r_eff**3)    # [Msun/pc³]
rho0_cgs = rho0 * MSUN / PC**3   # [g/cm³]

print(f"  Central stellar density: {rho0:.2e} Msun/pc³")
print(f"                          {rho0_cgs:.2e} g/cm³")
print()

# ── Test 1: Velocity at r_eff (macroscopic) ──
# At ~r_eff, a significant fraction of M* is enclosed
M_enc_eff = 0.5 * Mstar    # ~50% within r_eff for typical Sersic
v_reff = np.sqrt(G_eff * M_enc_eff / r_eff)
FWHM_reff = 2.355 * v_reff    # converting σ to FWHM

print("─" * 65)
print("Test 1: Velocity at r_eff (global scale)")
print("─" * 65)
print(f"  v(r=r_eff)     = {v_reff:.0f} km/s")
print(f"  FWHM (2.355σ)  = {FWHM_reff:.0f} km/s")
print()

# ── Test 2: Gas cloud at BLR-scale radii ──
# Broad lines come from gas at ~0.01-0.1 pc from center
# In the UID/no-BH picture, the gas is in the nuclear star cluster core
# The enclosed mass at r << r_eff depends on central density profile

# For a Sersic n=4, the 3D deprojected density at small r is:
# ρ(r) ≈ ρ_0 × (r/r_eff)^(-p) where p ≈ 0.85 for n=4
# M(<r) = 4π ∫ ρ(r') r'^2 dr' = 4π ρ_0 r_eff^p / (3-p) × r^(3-p)

# Critical question: at BLR radii (0.01-0.1 pc), what is M(<r)?
p_sersic = 0.85
M_enc_at_r = lambda r: 4 * np.pi * rho0 * r_eff**p_sersic / (3 - p_sersic) * r**(3 - p_sersic)

print("─" * 65)
print("Test 2: Velocity at BLR-scale radii (no BH)")
print("─" * 65)
print(f"  (Sersic n=4, ρ∝r^(-{p_sersic:.2f}) inner profile)")
print()

radii_pc = np.array([0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0])
print(f"  {'r [pc]':>8s}  {'M(<r) [Msun]':>14s}  {'v [km/s]':>10s}  {'FWHM [km/s]':>14s}  {'G_eff/G_N':>10s}")
print("  " + "-" * 60)
for r in radii_pc:
    Menc = M_enc_at_r(r)
    v_local = np.sqrt(G_eff * Menc / r)
    fwhm_local = 2.355 * v_local
    # Note: at these small radii, the density is high → high Σ_local → G_eff even larger
    # Local Σ at radius r: Σ_local ≈ M(<r)/(π r^2)
    Sigma_local = Menc / (np.pi * r**2)
    G_ratio_local = 1.0 + epsilon_g * (Sigma_local / Sigma0)**beta
    G_eff_local = G_N * G_ratio_local
    v_local_geff = np.sqrt(G_eff_local * Menc / r)
    fwhm_local_geff = 2.355 * v_local_geff
    print(f"  {r:>8.2f}  {Menc:>14.2e}  {v_local_geff:>8.0f}  {fwhm_local_geff:>12.0f}  {G_ratio_local:>8.2f}")

# ── Test 3: What BH mass would G_eff(Σ) framework predict as a "false" virial estimate? ──
# If an observer assumes G=G_N and derives M_BH from:
#   M_BH, virial = f × R_BLR × v^2 / G_N
# But the actual v^2 is from G_eff × M(<r)/r, then:
#   M_BH, virial = f × R_BLR × [G_eff/G_N × M(<r)/r] / G_N
# The ratio M_BH, virial / M_true = G_eff/G_N × [R_BLR × f × M(<r) / (r × M_true)]
# Since the observer assumes R_BLR ~ L^0.5 (R-L relation), they estimate:

# For a typical AGN with L_bol ~ 10^45.6:
# R_BLR ≈ 0.1 × (L_bol/10^45)^0.5 ≈ 0.025-0.05 pc (standard R-L relation)
L_bol = 10**45.6   # erg/s
R_BLR = 0.1 * (L_bol / 1e45)**0.5   # [pc]
# But if there's no central BH and the line width comes from extended potential,
# the actual R_BLR is ill-defined. The "virial" product gives a spurious BH mass.

print()
print("─" * 65)
print("Test 3: Spurious BH mass from assuming virial equilibrium")
print("─" * 65)
print(f"  Assuming R_L relation → R_BLR ≈ {R_BLR:.3f} pc")
print(f"  f_virial ≈ 3 (typical)")
f_vir = 3.0

# At R_BLR, compute M_enc and velocity
Menc_BLR = M_enc_at_r(R_BLR)
# But wait - at 0.03 pc, the stellar mass is tiny. The gas at BLR radii
# would be responding to the central potential well. Without a BH,
# the gas velocity at r ~ R_BLR is:
v_at_RBLR = np.sqrt(G_eff * Menc_BLR / R_BLR)
M_BH_spurious = f_vir * R_BLR * v_at_RBLR**2 / G_N
M_BH_spurious_geff = f_vir * R_BLR * (G_eff * Menc_BLR / R_BLR) / G_N
# = f_vir × Menc_BLR × G_eff/G_N

# Actually this is trivially:
M_BH_spurious_geff = f_vir * Menc_BLR * G_ratio

# Kokorev+ measured M_BH ~ 10^6.7
M_BH_kokorev = 10**6.7
print(f"  M_enc at R_BLR = {Menc_BLR:.2e} Msun")
print(f"  v_at_RBLR (G_eff) = {v_at_RBLR:.0f} km/s")
print(f"  Spurious M_BH = {M_BH_spurious_geff:.2e} Msun = 10^{np.log10(M_BH_spurious_geff):.1f}")
print(f"  Kokorev+ M_BH = {M_BH_kokorev:.2e} Msun")
print(f"  Ratio = {M_BH_spurious_geff/M_BH_kokorev:.2f}")
print()

# ── Test 4: Central concentration needed ──
# For the stellar mass alone to produce v ~ FWHM/2.355 ~ 3000/2.355 ~ 1270 km/s
# at the observed BLR radius, we need much higher central concentration.
# This tests: what central density is needed?

target_FWHM = 3000  # km/s (typical broad line FWHM)
target_v = target_FWHM / 2.355
print("─" * 65)
print("Test 4: What central concentration is needed?")
print("─" * 65)
print(f"  Target FWHM ≈ {target_FWHM} km/s → v = {target_v:.0f} km/s")
print()

# At radius r, v^2 = G_eff × M(<r) / r
# M(<r) = v^2 × r / G_eff = v^2 × r / (G_N × G_ratio)
# For target_v = 1274 km/s at r = R_BLR = 0.03 pc:
r_test = R_BLR
Menc_needed = target_v**2 * r_test / (G_N * G_ratio)
print(f"  At r = {r_test:.3f} pc:")
print(f"    M_enc needed = {Menc_needed:.2e} Msun")
print(f"    Available (Sersic) = {Menc_BLR:.2e} Msun")
print(f"    Shortfall factor = {Menc_needed/Menc_BLR:.0f}x")

# What about at 0.1 pc?
r_test2 = 0.1
Menc_at_01 = M_enc_at_r(r_test2)
Menc_needed_01 = target_v**2 * r_test2 / (G_N * G_ratio)
print(f"\n  At r = {r_test2:.2f} pc:")
print(f"    M_enc needed = {Menc_needed_01:.2e} Msun")
print(f"    Available (Sersic) = {Menc_at_01:.2e} Msun")
print(f"    Shortfall factor = {Menc_needed_01/Menc_at_01:.0f}x")

# ── Test 5: Minimum central BH mass under G_eff ──
# If the core contains a compact cluster (not a single BH), 
# what mass must be packed into the central ~pc to produce the lines?
print()
print("─" * 65)
print("Test 5: Minimum central mass needed (G_eff enhanced)")
print("─" * 65)
# For a point mass M_c at center, v^2(r) = G_eff × M_c / r
# At r = R_BLR = 0.03 pc, v_target = 1274 km/s:
M_central_needed = target_v**2 * R_BLR / (G_N * G_ratio)
M_central_solar = M_central_needed / MSUN * PC / PC  # wrong units...

# Let me redo this properly
# G_N = 4.30091e-3 pc·(km/s)^2 / Msun
# v^2 = G_eff × M_c / r
# M_c = v^2 × r / G_eff
M_central_needed = target_v**2 * R_BLR / (G_eff)
print(f"  R_BLR = {R_BLR:.3f} pc")
print(f"  v_target = {target_v:.0f} km/s")
print(f"  G_eff/G_N = {G_ratio:.2f}")
print(f"  M_central needed = {M_central_needed:.2e} Msun")
print(f"  = 10^{np.log10(M_central_needed):.1f} Msun")
print(f"  Kokorev+ M_BH = {M_BH_kokorev:.2e} = 10^{np.log10(M_BH_kokorev):.1f} Msun")
print(f"  Ratio to Kokorev+ = {M_central_needed/M_BH_kokorev:.2f}")

# ── Test 6: What if the BLR radius is much smaller? ──
# Some models suggest LRD broad lines come from even smaller radii
print()
print("─" * 65)
print("Test 6: Sensitivity to R_BLR")
print("─" * 65)
test_rs = np.logspace(-3, 0, 8)  # 0.001 to 1 pc
print(f"  {'r [pc]':>10s}  {'M_central [Msun]':>16s}  {'log M_central':>14s}  {'FWHM [km/s]':>12s}")
print("  " + "-" * 55)
for r in test_rs:
    # At this radius, what central mass produces v_target?
    Mc = target_v**2 * r / (G_eff)
    # Conversely, if we assume M_central = 10^6.7 * G_ratio (mock virial error),
    # what velocity does it produce?
    v_at_r = np.sqrt(G_eff * M_BH_kokorev / r)
    fwhm_at_r = 2.355 * v_at_r
    print(f"  {r:>10.4f}  {Mc:>16.2e}  {np.log10(Mc):>13.1f}  {fwhm_at_r:>10.0f}")

print()
print("=" * 65)
print("Summary")
print("=" * 65)
print(f"  G_eff/G_N at GLIMPSE-17775's Σ: {G_ratio:.2f}")
print(f"  Stellar mass alone at r_eff: v ~ {v_reff:.0f} km/s (FWHM ~ {FWHM_reff:.0f} km/s)")
print(f"  For BLR-scale (~{R_BLR:.2f} pc):")
print(f"    Enclosed stellar mass too small by ~{max(1, int(Menc_needed/max(Menc_BLR,1e-10))):.0f}x")
print(f"    → A central point mass of ~10^{np.log10(M_central_needed):.1f} Msun still needed")
print(f"    → But this is {M_central_needed/M_BH_kokorev:.1f}x smaller than Kokorev+ virial estimate")
print()
print("Conclusion: G_eff(Σ) alone cannot eliminate the need for a central mass,")
print("but it reduces the required mass by G_eff/G_N ~ {:.1f}x compared to standard virial.".format(G_ratio))
print("The virial BH mass is overestimated by roughly G_eff/G_N factor.")

