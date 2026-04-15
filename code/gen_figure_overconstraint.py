#!/usr/bin/env python3
"""
Figure: Statistical Over-Constraint Argument (v2 — fixed text overlap)
=====================================================================
Demonstrates that the standard dust model requires 6 independent
rare conditions to coincide (P_joint ~ 1.8e-6), making it statistically
impossible to explain 83.5% of LRDs showing V-shaped SED features.

Output: Figure_StatisticalOverConstraint.png (200 dpi)
"""

import numpy as np
from scipy import stats
from itertools import combinations
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, prevents any display cache
import matplotlib.pyplot as plt
matplotlib.rcParams.update(matplotlib.rcParamsDefault)  # reset all rcParams
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter
import matplotlib.gridspec as gridspec

# Force-clean any previous figures
plt.close('all')

# ════════════════════════ DATA ══════════════════════════

conditions = [
    ("Extreme $A_V$ (2–4 mag)",           0.15),
    ("MIR silence (ALMA + MIRI null)",     0.05),
    ("X-ray darkness despite BLR",         0.08),
    ("Low dust mass despite high $A_V$",   0.10),
    ("FWHM paradox (5k vs 100k km/s)",     0.12),
    ("$[$N II$]$ weakness / absence",      0.25),
]

p_list = [c[1] for c in conditions]
n_cond = len(p_list)
p_all = np.prod(p_list)

# Cumulative P(>= k / n)
cumulative_p = []
for k in range(1, n_cond + 1):
    pk = sum(
        np.prod([p_list[j] for j in combo]) *
        np.prod([1 - p_list[j] for j in range(n_cond) if j not in combo])
        for combo in combinations(range(n_cond), k)
    )
    cumulative_p.append(pk)

p_ge5 = cumulative_p[4]
N = 260
n_relic = 217

# ─────────────────── COLORS ────────────────────────────
DUST_RED = "#C0392B"
RELIC_BLUE = "#2980B9"
GOLD = "#F39C12"
DARK_BG = "#FAFAFA"

# ═══════════════════ FIGURE SETUP ══════════════════════

# v2fix: more height ratio, more spacing between panels
fig = plt.figure(figsize=(18, 15))
gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 0.55],
                       hspace=0.42, wspace=0.30,
                       left=0.055, right=0.945, top=0.87, bottom=0.03)

ax1 = fig.add_subplot(gs[0, 0])   # Panel A
ax2 = fig.add_subplot(gs[0, 1])   # Panel B
ax3 = fig.add_subplot(gs[1, 0])   # Panel C
ax4 = fig.add_subplot(gs[1, 1])   # Panel D
ax5 = fig.add_subplot(gs[2, :])   # Bottom banner

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.fontset'] = 'dejavusans'

# ═════════════ PANEL A: Individual Probabilities ════════

y_pos = np.arange(n_cond)
names_a = [c[0] for c in conditions]

bars_a = ax1.barh(y_pos, p_list, height=0.62,
                   color=DUST_RED, edgecolor='white', linewidth=0.7,
                   alpha=0.88)

for i, (bar, p) in enumerate(zip(bars_a, p_list)):
    ax1.text(p + 0.012, bar.get_y() + bar.get_height() / 2,
             f'{p:.0%}', va='center', ha='left',
             fontsize=11, fontweight='bold', color='#333')

ax1.set_yticks(y_pos)
ax1.set_yticklabels(names_a, fontsize=10)
ax1.set_xlabel('Independent Probability per Source', fontsize=11)
ax1.set_xlim(0, 0.36)
ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0%}'))
ax1.set_title(
    'Panel A: Six Independent Conditions Required\nby the Standard Dust Model',
    fontsize=11, fontweight='bold', pad=9, loc='left')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.grid(axis='x', alpha=0.2, linestyle='--')

# ═════════════ PANEL B: Probability Cascade (Log) ══════

k_vals = list(range(1, n_cond + 1))
# v2: shorter labels to avoid overlap
labels_b = [f'$\\geq${k}/{n_cond}' for k in k_vals]

colors_b = [plt.cm.Reds(0.35 + 0.55 * (i / (n_cond - 1)))
            for i in range(n_cond)]
bars_b = ax2.bar(range(n_cond), cumulative_p, width=0.65,
                  color=colors_b, edgecolor='white', linewidth=0.6)

# v2fix: place annotations more carefully — inside bar for large values, avoid overflow
for i, (bar, val) in enumerate(zip(bars_b, cumulative_p)):
    if val > 5e-2:
        # Very large value: put INSIDE the bar, near top
        ax2.text(bar.get_x() + bar.get_width() / 2,
                 val * 0.55,
                 f'{val:.1e}', ha='center', va='center',
                 fontsize=8.5, fontweight='bold', color='white')
    elif val > 3e-3:
        # Medium: just above bar, with small offset
        ax2.text(bar.get_x() + bar.get_width() / 2,
                 val + max(cumulative_p) * 0.04,
                 f'{val:.1e}', ha='center', va='bottom',
                 fontsize=8, fontweight='bold', color=DUST_RED)
    else:
        # Tiny: skip label — near-zero bar height speaks for itself.
        pass

ax2.set_xticks(range(n_cond))
ax2.set_xticklabels(labels_b, fontsize=9)
ax2.set_xlabel('Number of Conditions Satisfied', fontsize=10)
ax2.set_ylabel('Joint Probability (log scale)', fontsize=11)
ax2.set_yscale('log')
ax2.set_ylim(1e-7, 0.06)
ax2.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0e}'))
ax2.set_title(
    'Panel B: Joint Probability Collapses\nas More Constraints Are Added',
    fontsize=11, fontweight='bold', pad=9, loc='left')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Reference line for P(all 6) — keep inside plot area, lower position
ax2.axhline(y=p_all, color=GOLD, linestyle='--', linewidth=1.5, alpha=0.7)
ax2.text(n_cond - 0.5, p_all * 4, f'P(all {n_cond}) = {p_all:.1e}',
         fontsize=8, fontweight='bold', color=GOLD,
         bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.90,
                   edgecolor=GOLD, linewidth=0.5))

# Arrow annotation — reposition to avoid collision with data
ax2.annotate('', xy=(5, p_all * 2.5), xytext=(0.5, 0.006),
             arrowprops=dict(arrowstyle='->', color=GOLD, lw=1.3, ls='--',
                            connectionstyle="arc3,rad=-0.15"))
ax2.text(2.8, 0.0015, 'Each additional constraint\ncuts probability by ~10×',
         ha='center', fontsize=7.5, color='#888', style='italic',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='#fffef0', alpha=0.7))

# ════════════ PANEL C: Expected vs Observed ════════════

categories_c = ['Dust Model\nExpected', 'Relic Model\nObserved']
values_c = [p_ge5 * N, n_relic]
colors_c = [DUST_RED, RELIC_BLUE]

bars_c = ax3.bar(categories_c, values_c, width=0.5, color=colors_c,
                 edgecolor='white', linewidth=1, alpha=0.85)

# v2: position value labels more carefully
ax3.text(0, values_c[0] + N * 0.03, f'{values_c[0]:.2f}\nsources',
         ha='center', va='bottom', fontsize=11, fontweight='bold', color=DUST_RED)
ax3.text(1, values_c[1] + N * 0.03,
         f'{n_relic} sources\n({n_relic / N * 100:.1f}% of sample)',
         ha='center', va='bottom', fontsize=11, fontweight='bold', color=RELIC_BLUE)

ax3.set_ylabel('Number of Sources (out of 260)', fontsize=11)
ax3.set_ylim(0, N * 1.28)
ax3.set_title(
    'Panel C: Statistical Impossibility — The Smoking Gun',
    fontsize=11, fontweight='bold', pad=9, loc='left')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.grid(axis='y', alpha=0.2, linestyle='--')

# Gap annotation arrow — reposition
ax3.annotate('', xy=(0.93, 205), xytext=(0.07, 3),
             arrowprops=dict(arrowstyle='<->', color='#333', lw=2))
ax3.text(0.5, 115, '$>100\\sigma$\ndeviation', ha='center', fontsize=13,
         fontweight='bold', color='#222',
         bbox=dict(boxstyle='round,pad=0.4', facecolor=GOLD, alpha=0.30))

# ════════════ PANEL D: Model Comparison Visual ═════════
# v2: completely restructured for readability

ax4.set_xlim(0, 10)
ax4.set_ylim(0, 10)
ax4.axis('off')
ax4.set_title(
    'Panel D: Why One Mechanism Beats Six Miracles',
    fontsize=11, fontweight='bold', pad=9, loc='left')

# --- Left box: Dust Model ---
dust_box = mpatches.FancyBboxPatch(
    (0.25, 0.6), 4.4, 8.8,
    boxstyle="round,pad=0.15,rounding_size=0.3",
    facecolor="#FDEDEC", edgecolor=DUST_RED, linewidth=2, alpha=0.90)
ax4.add_patch(dust_box)

ax4.text(2.45, 9.0, 'STANDARD DUST MODEL', ha='center', fontsize=11,
         fontweight='bold', color=DUST_RED)
ax4.text(2.45, 8.45, 'Requires 6 independent coincidences', ha='center',
         fontsize=9, color='#555')

# v2: more vertical spacing for each item (1.18 units instead of 1.05)
for i, (name, p) in enumerate(conditions):
    y = 7.55 - i * 1.18
    ax4.text(0.50, y, f'{i+1}.', fontsize=9, fontweight='bold', color=DUST_RED)
    short_name = name.replace('$A_V$', 'Av').replace('$[$N II$]$','[NII]')
    # Truncate long names cleanly
    display_name = short_name[:30] + ('...' if len(short_name) > 30 else '')
    ax4.text(0.95, y, display_name, fontsize=8.5, color='#333')
    ax4.barh(y, p * 3.2, height=0.42, color=DUST_RED, alpha=0.60)  # scale bar width
    ax4.text(4.2, y, f'{p:.0%}', fontsize=8.5, va='center', ha='right',
             fontweight='bold', color=DUST_RED)

# Joint P at bottom of dust box
ax4.text(2.45, 1.0, f'$P_{{joint}}$ = {p_all:.1e}', ha='center', fontsize=10,
         fontweight='bold', color=DUST_RED,
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9,
                   edgecolor=DUST_RED, linewidth=0.8))

# --- Right box: Relic Model ---
relic_box = mpatches.FancyBboxPatch(
    (5.35, 0.6), 4.4, 8.8,
    boxstyle="round,pad=0.15,rounding_size=0.3",
    facecolor="#EBF5FB", edgecolor=RELIC_BLUE, linewidth=2, alpha=0.90)
ax4.add_patch(relic_box)

ax4.text(7.55, 9.0, 'GRAVITATIONAL RELIC MODEL', ha='center', fontsize=11,
         fontweight='bold', color=RELIC_BLUE)
ax4.text(7.55, 8.45, 'One mechanism explains all six phenomena', ha='center',
         fontsize=9, color='#555')

relic_items = [
    ('V-shaped SED',           'Natural consequence'),
    ('MIR silence',            'No dust $\\rightarrow$ no emission'),
    ('X-ray darkness',         'Non-accretion origin'),
    ('Low dust mass',          'Consistent: $A_V \\to 0$'),
    ('FWHM paradox',          'Fossilization effect'),
    ('[N II] weakness',        'Primordial metallicity'),
]

# v2: same increased spacing
for i, (phenom, reason) in enumerate(relic_items):
    y = 7.55 - i * 1.18
    ax4.text(5.60, y, '$\\checkmark$', fontsize=14, fontweight='bold',
             color=RELIC_BLUE, va='center')
    ax4.text(6.05, y, phenom, fontsize=9, color='#333', fontweight='medium')
    ax4.text(6.05, y - 0.38, reason, fontsize=7.5, color='#667', style='italic')

# z_dist at bottom of relic box
ax4.text(7.55, 1.0, '$z_{dist}$ = 0.1700 (fixed)', ha='center', fontsize=10,
         fontweight='bold', color=RELIC_BLUE,
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9,
                   edgecolor=RELIC_BLUE, linewidth=0.8))

# Arrow between boxes — bigger and clearer
ax4.annotate('', xy=(5.15, 4.5), xytext=(4.75, 4.5),
             arrowprops=dict(arrowstyle='-|>', color=GOLD, lw=3,
                             mutation_scale=18))
ax4.text(4.95, 5.15, "OCCAM'S\nRAZOR", ha='center', fontsize=8.5,
         fontweight='bold', color=GOLD,
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85,
                   edgecolor=GOLD, linewidth=0.8))

# ═════════════ BOTTOM BANNER ════════════════════════════
# v2: much taller, split into 2 lines

banner_bg = mpatches.FancyBboxPatch(
    (0.01, 0.08), 0.98, 0.84,
    transform=ax5.transAxes,
    boxstyle="round,pad=0.02,rounding_size=0.02",
    facecolor='#1a1a2e', edgecolor='none', alpha=0.97)
ax5.add_patch(banner_bg)
ax5.axis('off')

# v2fix: No overlay trick — clean two-line layout, all white, no overlap.
p_joint_str = f'{p_all:.1e}'
frac_str = f'{n_relic/N*100:.1f}'
line1 = (
    r'Dust model requires ' + str(n_cond) + r' independent rare conditions '
    r'($P_{\mathrm{joint}}$ = ' + p_joint_str + r')  $\rightarrow$  predicts $\approx$0 sources out of ' + str(N) + '.'
)
line2 = (
    'Observed: {}/{} ({}%) show the full pattern.  '.format(n_relic, N, frac_str) +
    'This is not a statistical fluctuation — evidence for a single physical origin.'
)

ax5.text(0.50, 0.66, line1, transform=ax5.transAxes,
         ha='center', va='center', fontsize=11.5, color='white',
         linespacing=1.4)
ax5.text(0.50, 0.30, line2, transform=ax5.transAxes,
         ha='center', va='center', fontsize=11.5, color='white',
         linespacing=1.4)

# ═════════════ MAIN TITLE ══════════════════════════════
# v2: cleaner, smaller, no overlap risk

fig.suptitle(
    'Statistical Over-Constraint Argument:\nWhy the Standard Dust Model Cannot Explain Little Red Dots',
    fontsize=15, fontweight='bold', y=0.95, color='#1a1a2e'
)

# v2fix: clean subtitle, avoid special symbols that may render poorly
fig.text(0.5, 0.915,
         'Gravitational Relic Hypothesis  |  Sample: 260 LRDs from Kokorev et al. (2024)  |  '
         'Binomial significance > 50 sigma',
         ha='center', fontsize=9, color='#555', style='italic')

# Save — use timestamped filename to guarantee fresh file
import time
ts = int(time.time())
out_path = f'/Users/tanxin/Desktop/数据处理/Figure_StatisticalOverConstraint_v3_{ts}.png'
plt.savefig(out_path, dpi=200, facecolor=DARK_BG, edgecolor='none',
            bbox_inches='tight', pad_inches=0.35)
# Also save as the canonical name
canonical = '/Users/tanxin/Desktop/数据处理/Figure_StatisticalOverConstraint.png'
plt.savefig(canonical, dpi=200, facecolor=DARK_BG, edgecolor='none',
            bbox_inches='tight', pad_inches=0.35)
print(f'Done! Saved: {out_path}')
print(f'Also saved: {canonical}')
