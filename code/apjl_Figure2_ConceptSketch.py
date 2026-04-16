#!/usr/bin/env python3
"""
Figure 2 — Concept Sketch: Density-Dependent Gravitational Redshift
================================================================
Corrected physics v3 — fixed all text overlaps.

ALL LRDs have intrinsically blue spectra. What we observe as "red" is
the result of gravitational redshift from photons climbing out of
the potential well. The key difference is HOW MUCH shift.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle

# ════════════════════════════════════════════════════════
# SETTINGS
# ════════════════════════════════════════════════════════
FIG_WIDTH_INCH = 7.3
FIG_HEIGHT_INCH = 5.8   # taller for bottom text spacing
DPI = 300

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 8,
    'axes.labelsize': 8.5,
    'axes.titlesize': 9.5,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
})

OUTPUT_DIR = '/Users/tanxin/Desktop/数据处理/LRD_GitHub_Release_v2/APJL_Letter_DensityDependentRedshift'

C = {
    'intrinsic_blue': '#1565c0',
    'observed_red':  '#c62828',
    'observed_mild': '#e65100',
    'well_shallow': '#42a5f5',
    'well_deep':    '#e65100',
    'core_color':   '#ffca28',
    'accent_green': '#2e7d32',
    'text_dark':    '#212121',
    'text_muted':   '#616161',
    'arrow_gray':   '#757575',
}

# ════════════════════════════════════════════════════════
# FIGURE SETUP — generous margins
# ════════════════════════════════════════════════════════
fig = plt.figure(figsize=(FIG_WIDTH_INCH, FIG_HEIGHT_INCH), dpi=DPI)

gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.26, wspace=0.18,
                      left=0.05, right=0.95, top=0.82, bottom=0.06)

ax_LT = fig.add_subplot(gs[0, 0])
ax_RT = fig.add_subplot(gs[0, 1])
ax_LB = fig.add_subplot(gs[1, 0])
ax_RB = fig.add_subplot(gs[1, 1])

for ax in [ax_LT, ax_RT]:
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.set_aspect('equal'); ax.axis('off')

# ════════════════════════════════════════════════════════
# HELPERS
# ════════════════════════════════════════════════════════
def draw_well(ax, cx, cy, depth, width, color, af=0.08):
    x = np.linspace(cx-width/2, cx+width/2, 150)
    y = cy+0.35*depth - depth*(1-((x-cx)/(width/2))**2)
    ax.fill_between(x, cy-0.4, y, color=color, alpha=af, zorder=1)
    ax.plot(x, y, color=color, lw=2, zorder=2)

def draw_core(ax, cx, cy, r):
    ax.add_patch(Circle((cx,cy), r, fc=C['core_color'], ec='#333', lw=1, alpha=0.9, zorder=10))
    for i in range(2):
        gr = r*(1.18+i*0.1)
        ax.add_patch(Circle((cx,cy), gr, fc='none', ec=C['core_color'],
                   lw=0.65-i*0.22, alpha=0.28-i*0.08, zorder=9-i))

def photon(ax, s, e, col):
    t=np.linspace(0,1,50); my=(s[1]+e[1])/2+0.45
    x=s[0]+t*(e[0]-s[0]); y=s[1]*(1-t)**2+2*my*t*(1-t)+e[1]*t**2
    ax.plot(x,y,color=col,lw=1.3,zorder=15,alpha=0.85)
    ax.scatter([e[0]],[e[1]],c=col,s=22,zorder=16,marker='*',ec='white',lw=0.3)

def sed_shift_arrow(ax, x, y1, y2, lbl):
    ax.annotate('', xy=(x,y2), xytext=(x,y1),
                arrowprops=dict(arrowstyle='->',color=C['arrow_gray'],lw=1.2,shrinkA=0,shrinkB=0),zorder=20)
    if lbl:
        ax.text(x+0.22,(y1+y2)/2,lbl,fontsize=6.5,color=C['arrow_gray'],va='center',fontstyle='italic')

def make_sed(waves, pp, pv, bs, rs, sh=0):
    w=np.array(waves); fi=np.zeros_like(w)
    mb=w<pp; mr=w>=pp
    fi[mb]=pv-bs*(pp-w[mb]); fi[mr]=pv+rs*(w[mr]-pp); fi=np.maximum(fi,0.02)
    if sh>0:
        ws=w+sh; fs=np.interp(w,ws,fi,left=0.01,right=fi[-1]*1.04)
    else: fs=fi.copy()
    return w,fs,fi

# ════════════════════════════════════════════════════════
# TOP HEADER — key message
# ════════════════════════════════════════════════════════
header = (
    'All LRDs are intrinsically blue. Deeper wells '
    r'$\rightarrow$ larger redshift-like shift $\rightarrow$ higher F$_{444W}$/F$_{150W}$ at higher $\Sigma$.'
)
fig.text(0.50, 0.97, header, ha='center', va='top', fontsize=7.5,
         color=C['text_dark'], fontstyle='italic',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#f5f5f5', ec='#aaa', alpha=0.95, lw=0.4))

# ════════════════════════════════════════════════════════
# (a) LEFT TOP — Low Σ SHALLOW WELL
# ════════════════════════════════════════════════════════
ax_LT.text(5, 9.7, '(a)  Low $\\Sigma$ (Diffuse)', ha='center', va='top',
           fontsize=9, fontweight='bold', color=C['text_dark'])
ax_LT.text(5, 9.15, 'Shallow Potential Well', ha='center', va='top',
           fontsize=7.5, color=C['text_muted'], style='italic')

draw_well(ax_LT, 5, 4.5, depth=1.8, width=6.5, color=C['well_shallow'])
draw_core(ax_LT, 5, 4.5, r=0.95)

photon(ax_LT, (5, 5.7), (8.0, 7.6), C['observed_mild'])

# Labels — carefully spaced, NO overlaps
ax_LT.text(5, 4.5, 'Core', ha='center', va='center', fontsize=6, color='white', fontweight='bold', zorder=11)
ax_LT.text(1.3, 7.4, r'$z_{\rm grav}$ small', fontsize=8, color=C['observed_mild'], fontweight='bold')
ax_LT.text(7.8, 2.6, 'Shallow well', ha='center', fontsize=7, color=C['well_shallow'], fontweight='bold')
# Blue source label — lower left, compact single line
ax_LT.text(0.6, 2.8, 'Intrinsically\nblue', fontsize=6.5, color=C['intrinsic_blue'], ha='center',
           bbox=dict(boxstyle='round,pad=0.2', facecolor='#e3f2fd', ec='#90caf9', lw=0.4))
# Delta-lambda annotation — moved away from photon path
ax_LT.annotate('', xy=(7.3, 6.2), xytext=(6.6, 5.5),
               arrowprops=dict(arrowstyle='->', color=C['arrow_gray'], lw=0.9, connectionstyle='arc3,rad=0.15'))
ax_LT.text(7.5, 6.55, r'$\Delta\lambda$ small', fontsize=6, color=C['arrow_gray'], ha='center')

# ════════════════════════════════════════════════════════
# (b) RIGHT TOP — High Σ DEEP WELL
# ════════════════════════════════════════════════════════
ax_RT.text(5, 9.7, '(b)  High $\\Sigma$ (Compact)', ha='center', va='top',
           fontsize=9, fontweight='bold', color=C['text_dark'])
ax_RT.text(5, 9.15, 'Deep Potential Well', ha='center', va='top',
           fontsize=7.5, color=C['text_muted'], style='italic')

draw_well(ax_RT, 5, 4.5, depth=3.8, width=6.0, color=C['well_deep'], af=0.10)
draw_core(ax_RT, 5, 4.5, r=0.5)

photon(ax_RT, (5, 5.4), (8.0, 7.6), C['observed_red'])

# Labels — carefully spaced
ax_RT.text(5, 4.5, 'Core', ha='center', va='center', fontsize=5.5, color='#333', fontweight='bold', zorder=11)
ax_RT.text(1.0, 7.6, r'$z_{\rm grav}$ large', fontsize=8, color=C['observed_red'], fontweight='bold')
ax_RT.text(7.8, 2.3, 'Deep well', ha='center', fontsize=7, color=C['well_deep'], fontweight='bold')
# Blue source label — same position, slightly wider for 3 lines
ax_RT.text(0.6, 2.8, 'Same intrinsic\nblue source', fontsize=6, color=C['intrinsic_blue'], ha='center',
           bbox=dict(boxstyle='round,pad=0.2', facecolor='#e3f2fd', ec='#90caf9', lw=0.4))
# Delta-lambda annotation
ax_RT.annotate('', xy=(7.3, 6.2), xytext=(6.0, 5.0),
               arrowprops=dict(arrowstyle='->', color=C['arrow_gray'], lw=1.0, connectionstyle='arc3,rad=0.15'))
ax_RT.text(7.5, 6.55, r'$\Delta\lambda$ large', fontsize=6, color=C['arrow_gray'], ha='center', fontweight='bold')

# Σ↑ arrow between top panels — raised to avoid collision
fig.patches.append(mpatches.FancyArrowPatch(
    (0.48, 0.68), (0.52, 0.68),
    transform=fig.transFigure, arrowstyle='->', mutation_scale=16,
    lw=1.8, color=C['accent_green'], zorder=30))
fig.text(0.50, 0.72, r'$\Sigma \uparrow$', ha='center', va='bottom',
         fontsize=10, fontweight='bold', color=C['accent_green'])

# ════════════════════════════════════════════════════════
# BOTTOM ROW: SED COMPARISON
# ════════════════════════════════════════════════════════
wavs = np.array([0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.5, 4.4])

wi, fi_in, _ = make_sed(wavs, pp=1.5, pv=0.8, bs=0.8, rs=0.45, sh=0)
wo_l, fo_l, _ = make_sed(wavs, pp=1.5, pv=0.8, bs=0.8, rs=0.45, sh=0.32)
wo_r, fo_r, _ = make_sed(wavs, pp=1.5, pv=0.8, bs=0.8, rs=0.45, sh=0.80)

# ── (c) Low-Σ SED ──
ax_LB.set_xlim(0.2, 5.0); ax_LB.set_ylim(0, 1.12)
ax_LB.set_xlabel('Rest-frame Wavelength [$\\mu$m]', fontsize=8)
ax_LB.set_ylabel('Normalized Flux', fontsize=8)
ax_LB.set_title(r'(c)  Low $\Sigma$ LRD: Small Shift', fontsize=8.5, fontweight='bold', pad=4)

ax_LB.plot(wi, fi_in, color=C['intrinsic_blue'], lw=1.2, alpha=0.30, ls='--', zorder=8)
ax_LB.plot(wo_l, fo_l, color=C['observed_mild'], lw=2.2, zorder=10)

# Band markers — placed below x-axis area to avoid curve overlap
for wl, nm in [(1.5,'F150W'), (4.4,'F444W')]:
    ax_LB.axvline(x=wl, ymin=0.02, ymax=0.06, color='#888', lw=1.2, zorder=5)
    ax_LB.text(wl, 0.085, nm, fontsize=6, ha='center', color='#666')

sed_shift_arrow(ax_LB, 1.15, 0.22, 0.44, r'$z_{\rm grav}$')

ax_LB.tick_params(labelsize=6.5)
for s in ['top','right']: ax_LB.spines[s].set_visible(False)

# ── (d) High-Σ SED ──
ax_RB.set_xlim(0.2, 5.0); ax_RB.set_ylim(0, 1.22)
ax_RB.set_xlabel('Rest-frame Wavelength [$\\mu$m]', fontsize=8)
ax_RB.set_ylabel('Normalized Flux', fontsize=8)
ax_RB.set_title(r'(d)  High $\Sigma$ LRD: Large Shift', fontsize=8.5, fontweight='bold', pad=4)

ax_RB.plot(wi, fi_in, color=C['intrinsic_blue'], lw=1.2, alpha=0.22, ls='--', zorder=8)
ax_RB.plot(wo_r, fo_r, color=C['observed_red'], lw=2.2, zorder=10)

for wl, nm in [(1.5,'F150W'), (4.4,'F444W')]:
    ax_RB.axvline(x=wl, ymin=0.02, ymax=0.055, color='#888', lw=1.2, zorder=5)
    ax_RB.text(wl, 0.08, nm, fontsize=6, ha='center', color='#666')

sed_shift_arrow(ax_RB, 1.15, 0.18, 0.58, r'$z_{\rm grav}\uparrow$')

ax_RB.tick_params(labelsize=6.5)
for s in ['top','right']: ax_RB.spines[s].set_visible(False)

# Insight box — inside plot area but clear of data
ax_RB.text(0.96, 0.88, 'Larger $\\rightarrow$\nhigher F444/F150',
           transform=ax_RB.transAxes, fontsize=7, color=C['observed_red'],
           ha='right', va='top', fontweight='bold',
           bbox=dict(boxstyle='round,pad=0.25', facecolor='#fff8e1', ec='#ffb74d', lw=0.5, alpha=0.95))

# Legend / shared note between SED panels
ax_LB.text(0.03, 0.88, '--- Intrinsic (same)\n— Observed (shifted)',
           transform=ax_LB.transAxes, fontsize=6, color='#555', va='top',
           bbox=dict(boxstyle='round,pad=0.2', facecolor='white', ec='#ccc', lw=0.3, alpha=0.9))

# Save
for fmt in ['png','pdf']:
    fp = f'{OUTPUT_DIR}/Figure2_ConceptSketch_DensityDependentReddening.{fmt}'
    fig.savefig(fp, dpi=DPI, facecolor='white', bbox_inches='tight', pad_inches=0.08)
    print(f'Saved: {fp}')

plt.close(fig)
print('Done.')
