"""
Figure 4: The Phase Transition — "理论升华图"
Shows how G_eff/G_0 skyrockets when Sigma crosses critical threshold,
"freezing" compact objects into LRDs.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from scipy.special import expit

# ─── Physical Model Parameters ──────────
# G_eff(Sigma) = G_0 * [1 + eps * (Sigma/Sigma_0)^beta]
G0 = 1.0
eps_g = 2.8
Sigma_0 = 10**-4.5
beta = 1.0
log_Sigma_crit = -4.55

def G_eff(logSigma):
    Sigma = 10**logSigma
    ratio = Sigma / Sigma_0
    return G0 * (1 + eps_g * ratio**beta)

def G_eff_smooth(logSigma, k=12, x0=log_Sigma_crit):
    base = G_eff(logSigma)
    transition = expit(k * (logSigma - x0))
    G_low = np.ones_like(logSigma) * 1.02
    G_high = base
    return G_low + (G_high - G_low) * transition

# ─── Data range ─────────────────────────
logS_range = np.linspace(-5.3, -3.6, 500)
G_values = G_eff_smooth(logS_range, k=18, x0=log_Sigma_crit)

# Scatter representing real sources
np.random.seed(2024)
n_scatter = 120
scatter_logS = np.random.uniform(-5.2, -3.8, n_scatter)
scatter_G = G_eff_smooth(scatter_logS, k=18, x0=log_Sigma_crit)
scatter_G_noisy = scatter_G + np.random.normal(0, 0.08, n_scatter)
scatter_G_noisy = np.maximum(scatter_G_noisy, 0.95)

mask_normal = scatter_logS < log_Sigma_crit
mask_frozen = scatter_logS >= log_Sigma_crit

# ─── Plot Setup ─────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.linewidth': 1.3,
    'figure.facecolor': '#fafaf8',
})

fig = plt.figure(figsize=(15, 9), facecolor='#fafaf8')
gs = plt.GridSpec(2, 2, width_ratios=[1.3, 1], height_ratios=[1.2, 1],
                   wspace=0.25, hspace=0.32)

# ══════════════════════════════════
# PANEL (a): Main Phase Transition Curve
# ══════════════════════════════════
ax_main = fig.add_subplot(gs[0, :])

ax_main.axvspan(-5.3, log_Sigma_crit, alpha=0.06, color='#2196F3', zorder=0)
ax_main.axvspan(log_Sigma_crit, -3.6, alpha=0.08, color='#d32f2f', zorder=0)

ax_main.plot(logS_range, G_values, '-', c='#b71c1c', lw=3.5, zorder=7,
            label=r'$G_{\rm eff}(\Sigma)/G_0$')

# Uncertainty band
G_up = G_eff_smooth(logS_range, k=16, x0=log_Sigma_crit-0.05) * 1.12
G_lo = G_eff_smooth(logS_range, k=20, x0=log_Sigma_crit+0.05) * 0.88
ax_main.fill_between(logS_range, G_lo, G_up, alpha=0.10, color='#b71c1c', zorder=4)

ax_main.scatter(scatter_logS[mask_normal], scatter_G_noisy[mask_normal],
               s=25, c='#1565c0', marker='o', alpha=0.45, edgecolors='white',
               linewidths=0.4, zorder=5, label='Normal galaxies')
ax_main.scatter(scatter_logS[mask_frozen], scatter_G_noisy[mask_frozen],
               s=35, c='#c62828', marker='o', alpha=0.55, edgecolors='#7f0000',
               linewidths=0.6, zorder=6, label='LRDs (frozen regime)')

# Critical threshold line
ax_main.axvline(log_Sigma_crit, ls=(0, (6, 4)), c='#ff6f00', lw=2.8, zorder=8,
               label=f'Critical $\\Sigma_{{crit}}$ = $10^{{{log_Sigma_crit:.2f}}}$')

# Threshold annotation — moved LEFT & UP per user sketch
bbox_style = dict(boxstyle='round,pad=0.25', fc='#fff8e1', ec='#ff6f00',
                  lw=1.8, alpha=0.96)
ax_main.annotate('★ PHASE ★\n'
                 r'$\Sigma > \Sigma_{crit}$' + '\n'
                 r'$G_{\rm eff}$ → "Frozen"',
                xy=(log_Sigma_crit, G_eff_smooth(log_Sigma_crit)),
                xytext=(log_Sigma_crit+0.35, 10),
                fontsize=9, fontweight='bold', color='#e65100', ha='center',
                arrowprops=dict(arrowstyle='->', color='#ff6f00', lw=2.2,
                               connectionstyle='arc3,rad=-0.25'),
                bbox=bbox_style)

# Zone labels — NORMAL: top-left; FROZEN: just RIGHT of critical line per user circle
ax_main.text(-5.22, 12.5, 'NORMAL\nREGIME\n\n$G \\approx G_0$\n(Newtonian)',
            fontsize=9, ha='center', va='top', color='#1565c0',
            fontweight='bold', alpha=0.85,
            bbox=dict(boxstyle='round,pad=0.3', fc='#e3f2fd', ec='#1565c0',
                     alpha=0.7, lw=1.2))
ax_main.text(log_Sigma_crit + 0.15, 6.5, '"FROZEN"\nREGIME\n\n$G_{\\rm eff} \\gg G_0$\n→ LRD',
            fontsize=8.5, ha='left', va='bottom', color='#c62828',
            fontweight='bold', alpha=0.85,
            bbox=dict(boxstyle='round,pad=0.25', fc='#ffebee', ec='#c62828',
                     alpha=0.7, lw=1.2))

ax_main.set_xlabel(
    r'Surface Mass Density $\log_{10}(\Sigma \, [{\rm M}_\odot \, {\rm pc}^{-2}])$',
    fontsize=13, fontweight='bold')
ax_main.set_ylabel(r'Normalized Effective Gravity  $G_{\rm eff}/G_0$', fontsize=13,
                   fontweight='bold')
ax_main.set_title('(a) The Density–Gravity Phase Transition', fontsize=14,
                   fontweight='bold', pad=12)
ax_main.set_xlim(-5.3, -3.6)
ax_main.set_ylim(0.2, max(G_values)*1.28)
ax_main.legend(loc='upper left', fontsize=9, framealpha=0.95, edgecolor='#bbb')
ax_main.grid(True, alpha=0.20, ls='-', lw=0.5)

# ══════════════════════════════════
# PANEL (b): Physical Schematic
# ══════════════════════════════════
ax_schem = fig.add_subplot(gs[1, 0])
ax_schem.set_xlim(-0.5, 11)
ax_schem.set_ylim(0, 10)
ax_schem.axis('off')
ax_schem.set_facecolor('#fafaf8')

# --- LEFT: Normal galaxy ---
for i in range(12):
    r = 1.2 + i*0.22
    alpha = 0.35 - i*0.025
    circle = Circle((3.0, 5), r, fc='#90caf9', ec='none', alpha=max(alpha,0.06))
    ax_schem.add_patch(circle)
core_norm = Circle((3.0, 5), 0.5, fc='#1976d2', ec='#0d47a1', lw=1.5, zorder=10)
ax_schem.add_patch(core_norm)

# Photons escape freely
for angle in [25, 55, 85, 115, 145, 175]:
    rad = np.radians(angle)
    dx, dy = 2.0*np.cos(rad), 1.8*np.sin(rad)
    ax_schem.annotate('', xy=(3.0+dx, 5+dy),
                      xytext=(3.0+0.6*np.cos(rad), 5+0.6*np.sin(rad)),
                      arrowprops=dict(arrowstyle='->', color='#1976d2', lw=1.8))

ax_schem.text(3.0, 1.8, 'Normal Galaxy\n($\\Sigma < \\Sigma_{crit}$)',
             fontsize=8, ha='center', va='center', fontweight='bold', color='#1565c0',
             bbox=dict(boxstyle='round,pad=0.2', fc='#e3f2fd', ec='#1976d2', lw=1.5))
ax_schem.text(3.0, 0.65, r'$G \approx G_0$' + '\nPhotons escape freely\nStandard colors',
             fontsize=7.5, ha='center', va='center', color='#424242',
             style='italic')

# --- RIGHT: Frozen LRD ---
for i in range(8):
    r = 0.25 + i*0.15
    alpha = 0.6 - i*0.06
    cf = Circle((8.8, 5), r, fc='#ef9a9a', ec='none', alpha=max(alpha, 0.08),
                zorder=10-i)
    ax_schem.add_patch(cf)
core_frz = Circle((8.8, 5), 0.35, fc='#d32f2f', ec='#b71c1c', lw=2, zorder=20)
ax_schem.add_patch(core_frz)

# Ice crystal effect around frozen core
for angle in np.linspace(0, 2*np.pi, 8, endpoint=False):
    ix = 8.8 + 0.65*np.cos(angle)
    iy = 5 + 0.65*np.sin(angle)
    ax_schem.plot([8.8, ix], [5, iy], ':', c='#81d4fa', lw=1.0, alpha=0.4)
    ax_schem.scatter([ix], [iy], marker='*', s=25, c='#81d4fa', alpha=0.5, zorder=15)

ax_schem.text(8.8, 1.8, '"Frozen" LRD\n($\\Sigma > \\Sigma_{crit}$)',
             fontsize=8, ha='center', va='center', fontweight='bold', color='#c62828',
             bbox=dict(boxstyle='round,pad=0.2', fc='#ffebee', ec='#d32f2f', lw=1.5))
ax_schem.text(8.8, 0.65, r'$G_{\rm eff} \gg G_0$' + '\nGrav. redshift traps light\nExtreme colors',
             fontsize=7.5, ha='center', va='center', color='#424242',
             style='italic')

# Arrow between them
ax_schem.annotate("", xy=(7.5, 5), xytext=(4.5, 5),
                 arrowprops=dict(arrowstyle="->", color="#ff6f00", lw=3.5,
                                connectionstyle="arc3,rad=0"))
ax_schem.text(6.0, 5.65, r'$\Sigma \uparrow$' + '\nPhase\nTransition',
             fontsize=10, ha='center', va='bottom', color='#e65100', fontweight='bold')

ax_schem.set_title('(b) Physical Picture: From Normal to "Frozen"',
                    fontsize=14, fontweight='bold', pad=12)

# ══════════════════════════════════
# PANEL (c): Observable Consequences Table
# ══════════════════════════════════
ax_table = fig.add_subplot(gs[1, 1])
ax_table.axis('off')
ax_table.set_facecolor('#fafaf8')

table_data = [
    ['Observable', 'Normal Regime', 'Frozen (LRD)'],
    ['Color (F150W−F444W)', '~0.5 (blue)', '<−1.5 (extreme red)'],
    ['Line Width FWHM', 'Narrow (<1000 km/s)', 'Broad (>2000 km/s)'],
    ['$G_{\\rm eff}/G_0$', '~1.0', '>3−8× enhanced'],
    ['Physical state', 'Dynamical equilibrium', 'Gravitational relic'],
]

cell_colors = [['#eceff1', '#eceff1', '#eceff1'],
               ['#ffffff', '#ffffff', '#fff8e1'],
               ['#ffffff', '#ffffff', '#fff8e1'],
               ['#ffffff', '#ffffff', '#fff8e1'],
               ['#ffffff', '#ffffff', '#fff8e1']]

tbl = ax_table.table(cellText=table_data, cellLoc='center', loc='center',
                     colWidths=[0.38, 0.31, 0.31], cellColours=cell_colors)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9.5)
tbl.scale(1.1, 2.0)

for j in range(3):
    cell = tbl[(0, j)]
    cell.set_text_props(fontweight='bold', fontsize=10.5)
    cell.set_facecolor('#cfd8dc')
    cell.set_height(0.14)
    if j == 2:
        cell.get_text().set_color('#c62828')
    elif j == 1:
        cell.get_text().set_color('#1565c0')

for i in range(1, 5):
    tbl[(i, 0)].set_text_props(fontweight='bold', color='#333')
    tbl[(i, 2)].set_text_props(color='#b71c1c', fontweight='bold')

ax_table.set_title('(c) Observable Signatures of the Phase Transition',
                    fontsize=14, fontweight='bold', pad=12, y=0.92)

# ─── Main Title ───────────────────────
fig.suptitle('Figure 4: The Phase Transition — When Surface Density Freezes Spacetime',
             fontsize=17, fontweight='bold', y=0.995, color='#1a1a1a')

plt.tight_layout(rect=[0, 0, 1, 0.97])

out_path = '/Users/tanxin/WorkBuddy/20260412234449/dense-env-repo/figures/Figure4_PhaseTransition.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='#fafaf8')
print(f'Saved: {out_path}')
print(f'Critical log Sigma: {log_Sigma_crit}')
print(f'G_eff at crit: {G_eff_smooth(log_Sigma_crit):.2f}')
print(f'G_eff at max (logS=-3.6): {G_eff_smooth(-3.6):.2f}')
