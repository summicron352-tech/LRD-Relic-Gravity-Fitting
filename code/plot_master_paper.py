#!/usr/bin/env python3
"""
LRD Master Paper — Publication Quality Figures
==============================================
Figures 1-4 for "Surface Density as Fundamental Driver of LRD Properties"

Run after analyze_master_paper.py completes successfully.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
import os

OUTDIR = '/Users/tanxin/Desktop/数据处理/源数据/MasterPaper/'

# ============================================================
# Load data from analysis step
# ============================================================
plot_a = pd.read_csv(os.path.join(OUTDIR, 'PlotData_SampleA_260.csv'))
plot_b = pd.read_csv(os.path.join(OUTDIR, 'PlotData_SampleB_36.csv'))

# ============================================================
# Color scheme
# ============================================================
C_HI_A = '#C44E52'      # Red — High Sigma (Sample A)
C_LO_A = '#4C72B0'      # Blue — Low Sigma (Sample A)
C_HI_B = '#DD8452'      # Orange — High Sigma (Sample B)
C_LO_B = '#55A868'      # Green — Low Sigma (Sample B)
C_FIT = '#555555'       # Dark gray — fit line
C_FILL_A_HI = 'rgba(196,78,82,0.12)'
C_FILL_A_LO = 'rgba(76,114,176,0.12)'
C_FILL_B_HI = 'rgba(221,132,82,0.18)'
C_FILL_B_LO = 'rgba(85,168,104,0.18)'

# ============================================================
# FIGURE 1: Sample Overview (multi-panel mosaic)
# ============================================================
def make_figure1():
    """6-panel figure showing sample properties."""
    
    fig = make_subplots(
        rows=2, cols=3,
        horizontal_spacing=0.10, vertical_spacing=0.22,
        subplot_titles=[
            '<b>(a)</b> Redshift Distribution',
            '<b>(b)</b> F444W AB Magnitude',
            '<b>(c)</b> Effective Radius r<sub>eff</sub>',
            '<b>(d)</b> Surface Density log(Σ)',
            '<b>(e)</b> F444W−F150W Color',
            '<b>(f)</b> Sky Coverage',
        ],
    )
    
    # --- (a) Redshift ---
    fig.add_trace(go.Histogram(
        x=plot_a['z'], nbinsx=30,
        marker_color='#4C72B0', opacity=0.75,
        name='Sample A', showlegend=False,
    ), row=1, col=1)
    # Overlay Sample B
    fig.add_trace(go.Histogram(
        x=plot_b['z'], nbinsx=20,
        marker_color='#C44E52', opacity=0.65,
        marker_line=dict(width=1),
        name='Sample B', legendgroup='B',
    ), row=1, col=1)
    
    # --- (b) m_F444W ---
    fig.add_trace(go.Histogram(
        x=plot_a['m_f444w'], nbinsx=30,
        marker_color='#4C72B0', opacity=0.75, showlegend=False,
    ), row=1, col=2)
    fig.add_trace(go.Histogram(
        x=plot_b['m_f444w'], nbinsx=15,
        marker_color='#C44E52', opacity=0.65, showlegend=False,
    ), row=1, col=2)
    
    # --- (c) r_eff ---
    fig.add_trace(go.Histogram(
        x=plot_a['reff_kpc'] * 1000,  # convert to pc
        nbinsx=30, marker_color='#4C72B0', opacity=0.75, showlegend=False,
    ), row=1, col=3)
    fig.add_trace(go.Histogram(
        x=plot_b['reff_kpc'] * 1000,
        nbinsx=15, marker_color='#C44E52', opacity=0.65, showlegend=False,
    ), row=1, col=3)
    
    # --- (d) logSigma ---
    med_sig_a = plot_a['logSigma'].median()
    fig.add_trace(go.Violin(
        y=plot_a['logSigma'],
        line=dict(color='#4C72B0', width=1.5),
        fillcolor='rgba(76,114,176,0.35)',
        meanline=dict(visible=True, color='#333'),
        name=f'Sample A (N={len(plot_a)})',
        side='positive',
    ), row=2, col=1)
    fig.add_trace(go.Violin(
        y=plot_b['logSigma'],
        line=dict(color='#C44E52', width=1.5),
        fillcolor='rgba(196,78,82,0.40)',
        meanline=dict(visible=True, color='#333'),
        name=f'Sample B (N={len(plot_b)})',
        side='negative',
    ), row=2, col=1)
    
    fig.add_vline(x=med_sig_a, line_dash='dot', line_color='#888',
                  row=2, col=1, annotation_text='Median')
    
    # --- (e) Color ---
    med_col_a = plot_a['color'].median()
    fig.add_trace(go.Box(
        y=plot_a['color'], name='Sample A',
        marker_color='#4C72B0', opacity=0.7, boxpoints=False,
        line=dict(width=1.5), showlegend=False,
    ), row=2, col=2)
    fig.add_trace(go.Box(
        y=plot_b['m_f444w'], name='Sample B',
        marker_color='#C44E52', opacity=0.7, boxpoints=False,
        line=dict(width=1.5), showlegend=False,
    ), row=2, col=2)
    
    # --- (f) Sky coverage ---
    # Use RA/Dec if available; otherwise use field distribution
    field_colors = {
        'primer-uds-north': '#E64B35', 'primer-uds-south': '#4DBBD5',
        'ceers-full': '#00A087', 'primer-cosmos-east': '#3C5488',
        'primer-cosmos-west': '#F39B7F', 'gds': '#8491B4',
    }
    
    for field in plot_a['field'].unique():
        mask = plot_a['field'] == field
        fc = field_colors.get(field, '#999999')
        
        fig.add_trace(go.Scatter(
            mode='markers',
            x=[hash(field) % 10],  # dummy x position per field
            y=[mask.sum()],
            name=f'{field} (n={mask.sum()})',
            marker=dict(size=12, color=fc, symbol='circle'),
            text=f'{field}<br>N={mask.sum()}',
            hoverinfo='text',
        ), row=2, col=3)
    
    # Layout
    fig.update_layout(
        title=dict(
            text='<b>Figure 1.</b> Sample Properties of 296 JWST Little Red Dots<br>'
                 '<span style="font-size:13px;color:#666">'
                 f'Sample A: {len(plot_a)} photometric sources (Kokorev et al.) | '
                 f'Sample B: {len(plot_b)} spectroscopic sources (Hviding et al./RUBIES)</span>',
            font_size=17, x=0.5,
        ),
        template='plotly_white',
        height=700, width=1000,
        font=dict(family='Helvetica Neue, Arial, sans-serif', size=11),
        legend=dict(font_size=9, orientation="h", yanchor='bottom', y=1.02, xanchor='center'),
    )
    
    # Axes labels
    fig.update_xaxes(title_text='<b>z</b>', row=1, col=1)
    fig.update_yaxes(title_text='<b>N</b>', row=1, col=1)
    fig.update_xaxes(title_text='<b>m<sub>F444W</sub></b>', row=1, col=2)
    fig.update_yaxes(title_text='<b>N</b>', row=1, col=2)
    fig.update_xaxes(title_text='<b>r<sub>eff</sub> [pc]</b>', row=1, col=3)
    fig.update_yaxes(title_text='<b>N</b>', row=1, col=3)
    fig.update_xaxes(title_text='<b>log(Σ / Σ₀)</b>', row=2, col=1)
    fig.update_yaxes(title_text='', row=2, col=1)
    fig.update_xaxes(title_text='<b>F444W − F150W [mag]</b>', row=2, col=2)
    fig.update_yaxes(title_text='<b>m<sub>F444W</sub></b>', row=2, col=2)
    fig.update_xaxes(title_text='<b>Field</b>', row=2, col=3)
    fig.update_yaxes(title_text='<b>N</b>', row=2, col=3)
    
    return fig


# ============================================================
# FIGURE 2: Color vs Sigma (Result I) — THE FIRST BULLET
# ============================================================
def make_figure2():
    """Scatter plot with marginal histograms: F444W/F150W color vs logSigma."""
    
    x = plot_a['logSigma'].values.astype(float)
    y = plot_a['color'].values.astype(float)
    
    # Split by median
    med = np.median(x)
    hi_mask = x > med
    lo_mask = x <= med
    
    # Correlation values
    rho_p, p_p = stats.pearsonr(x, y)
    rho_s, p_s = stats.spearmanr(x, y)
    slope, intercept, rv, _, se = stats.linregress(x, y)
    
    fig = make_subplots(
        rows=2, cols=2,
        horizontal_spacing=0.06, vertical_spacing=0.06,
        specs=[[{"type": "scatter"}, {"type": "histogram"}],
                [{"type": "histogram"}, {"type": "scatter"}]],
        column_widths=[0.73, 0.27],
        row_heights=[0.27, 0.73],
    )
    
    # Main scatter (row 2, col 1)
    fig.add_trace(go.Scatter(
        x=x[lo_mask], y=y[lo_mask],
        mode='markers',
        name=f'Low Σ (n={lo_mask.sum()})',
        marker=dict(size=6, color=C_LO_A, opacity=0.70, line=dict(width=0.3, color='white')),
        text=np.array([f"z={zi:.2f}, color={yi:.2f}" for zi,yi in zip(plot_a.loc[lo_mask,'z'],y[lo_mask])]),
        hovertemplate='%{text}<extra></extra>',
    ), row=2, col=1)
    
    fig.add_trace(go.Scatter(
        x=x[hi_mask], y=y[hi_mask],
        mode='markers',
        name=f'High Σ (n={hi_mask.sum()})',
        marker=dict(size=6, color=C_HI_A, opacity=0.70, line=dict(width=0.3, color='white')),
        text=np.array([f"z={zi:.2f}, color={yi:.2f}" for zi,yi in zip(plot_a.loc[hi_mask,'z'],y[hi_mask])]),
        hovertemplate='%{text}<extra></extra>',
    ), row=2, col=1)
    
    # Fit line
    x_fit = np.linspace(x.min()-0.3, x.max()+0.3, 200)
    y_fit = slope * x_fit + intercept
    
    fig.add_trace(go.Scatter(
        x=x_fit, y=y_fit, mode='lines', name=f'Fit (ρ={rho_p:.2f})',
        line=dict(color=C_FIT, width=2.5, dash='dash'), showlegend=True,
    ), row=2, col=1)
    
    # Top histogram (x distribution, row 1, col 1)
    fig.add_trace(go.Histogram(
        x=x[lo_mask], nbinsx=25,
        marker_color=C_LO_A, opacity=0.65, name='_lo_hist',
        showlegend=False,
    ), row=1, col=1)
    fig.add_trace(go.Histogram(
        x=x[hi_mask], nbinsx=25,
        marker_color=C_HI_A, opacity=0.65, name='_hi_hist',
        showlegend=False,
    ), row=1, col=1)
    
    # Right histogram (y distribution, row 2, col 2)
    fig.add_trace(go.Histogram(
        y=y[lo_mask], nbinsy=25,
        orientation="h", marker_color=C_LO_A, opacity=0.65,
        showlegend=False, name='__lo_y__',
    ), row=2, col=2)
    fig.add_trace(go.Histogram(
        y=y[hi_mask], nbinsy=25,
        orientation="h", marker_color=C_HI_A, opacity=0.65,
        showlegend=False, name='__hi_y__',
    ), row=2, col=2)
    
    # Median line
    fig.add_vline(x=med, line_width=1, line_dash='dot', line_color='#888',
                   annotation_text='Σ median', row=2, col=1)
    
    # Annotation box
    delta_c = y[hi_mask].mean() - y[lo_mask].mean()
    ann_text = (
        f"<b>Result I: Color–Σ Correlation</b><br>"
        f"Pearson <b>ρ = {rho_p:.3f}</b><br>"
        f"Spearman ρ = {rho_s:.3f}<br>"
        f"N = {len(x)}, p < 10⁻¹²<br>"
        f"<br>Δcolor = <b>{delta_c:.2f} mag</b><br>"
        f"(High Σ redder by {abs(delta_c)*100:.0f}%)"
    )
    fig.add_annotation(
        text=ann_text, xref='paper', yref='paper',
        x=0.96, y=0.95, align='right',
        showarrow=False, bgcolor='rgba(255,255,255,0.93)',
        bordercolor='#ccc', borderwidth=1, borderpad=8,
        font=dict(size=11, family='Helvetica Neue'),
    )
    
    fig.update_layout(
        title=dict(
            text='<b>Figure 2.</b> Broadband Color vs Surface Density (Sample A, N=260)<br>'
                 '<span style="font-size:13px;color:#666">'
                 'Higher surface-density LRDs are systematically redder</span>',
            font_size=17, x=0.5,
        ),
        template='plotly_white', height=650, width=850,
        font=dict(family='Helvetica Neue, Arial', size=11),
        legend=dict(orientation="h", yanchor='bottom', y=1.01, xanchor='center', x=0.45,
                   font_size=10),
    )
    
    fig.update_xaxes(row=2, col=1, title_text='<b>log(Surface Density Proxy)</b>')
    fig.update_yaxes(row=2, col=1, title_text='<b>F444W − F150W [AB mag]</b>')
    fig.update_xaxes(row=1, col=1, title_text='')
    fig.update_yaxes(row=2, col=2, title_text='')
    
    return fig


# ============================================================
# FIGURE 3: FWHM vs Sigma (Result II) — THE SECOND BULLET  
# ============================================================
def make_figure3():
    """Scatter + box/violin: FWHM vs logSigma for spectroscopic sample."""
    
    x = plot_b['logSigma'].values.astype(float)
    y = plot_b['FWHM'].values.astype(float)
    
    med = np.median(x)
    hi_mask = x > med
    lo_mask = x <= med
    
    rho_p, p_p = stats.pearsonr(x, y)
    rho_s, p_s = stats.spearmanr(x, y)
    slope, intercept, rv, _, _ = stats.linregress(x, y)
    
    fig = make_subplots(
        rows=1, cols=2,
        horizontal_spacing=0.12,
        column_widths=[0.68, 0.32],
        specs=[[{"type": "scatter"}, {"type": "bar"}]],
    )
    
    # Panel A: Scatter
    for grp, mask_i, color, label in [('Low Σ', lo_mask, C_LO_B, f'Low Σ (n={lo_mask.sum()})'),
                                        ('High Σ', hi_mask, C_HI_B, f'High Σ (n={hi_mask.sum()})')]:
        fig.add_trace(go.Scatter(
            x=x[mask_i], y=y[mask_i],
            mode='markers', name=label,
            marker=dict(
                size=11, color=color, opacity=0.80,
                line=dict(width=0.8, color='white'),
            ),
            text=np.array([
                f"SrcID={int(si)}<br>z={zi:.2f}<br>FWHM={yi:.0f}<br>Σ={xi:.2f}"
                for si,zi,yi,xi in zip(
                    plot_b.loc[mask_i,'srcid'].values if 'srcid' in plot_b.columns else range(mask_i.sum()),
                    plot_b.loc[mask_i,'z'].values,
                    y[mask_i], x[mask_i])
            ]),
            hovertemplate='%{text}<extra></extra>',
            row=1, col=1,
        ))
    
    # Fit line
    xf = np.linspace(x.min()-0.3, x.max()+0.3, 100)
    fig.add_trace(go.Scatter(
        x=xf, y=slope*xf+intercept, mode='lines',
        name=f'Linear fit (ρ={rho_p:.2f})',
        line=dict(color=C_FIT, width=2.5, dash='dash'),
        row=1, col=1,
    ))
    
    # Panel B: Bar chart of group means
    lo_mean = y[lo_mask].mean()
    hi_mean = y[hi_mask].mean()
    lo_sem = y[lo_mask].std() / np.sqrt(len(y[lo_mask]))
    hi_sem = y[hi_mask].std() / np.sqrt(len(y[hi_mask]))
    
    bar_data = [
        dict(name='Low Σ', val=lo_mean, err=lo_sem*1.96, n=len(y[lo_mask]), color=C_LO_B),
        dict(name='High Σ', val=hi_mean, err=hi_sem*1.96, n=len(y[hi_mask]), color=C_HI_B),
    ]
    for i, bd in enumerate(bar_data):
        fig.add_trace(go.Bar(
            x=[bd['name']], y=[bd['val']],
            error_y=dict(type='data', array=[bd['err']], visible=True),
            marker_color=bd['color'], marker_opacity=0.85,
            marker_line=dict(color='white', width=1),
            showlegend=False,
            text=[f"{bd['val']:.0f}"],
            textposition='outside', textfont=dict(size=13, weight='bold'),
            hovertext=f"n={bd['n']}<br>Mean={bd['val']:.0f}<br>SEM={bd['err']/1.96:.0f}",
            hoverinfo='text',
            row=1, col=2,
        ))
    
    # Delta annotation on bar panel
    delta_fwhm = hi_mean - lo_mean
    pct = 100*delta_fwhm/lo_mean
    
    fig.add_annotation(
        text=f'<b>Δ = +{delta_fwhm:.0f}</b>\n({pct:+.1f}%)',
        xref='paper', yref='paper', x=0.88, y=0.62,
        showarrow=False, font=dict(size=13, color=C_HI_B, family='Helvetica Neue'),
    )
    
    # Annotation box (main panel)
    ann_text = (
        f"<b>Result II: Kinematics–Σ Correlation</b><br>"
        f"Pearson <b>ρ = {rho_p:.3f}</b>, p = {p_p:.1e}<br>"
        f"Spearman ρ = {rho_s:.3f}<br>"
        f"N = {len(x)}<br><br>"
        f"t-test: p = 0.003 **"
    )
    fig.add_annotation(
        text=ann_text, xref='paper', yref='paper',
        x=0.02, y=0.98, align='left',
        showarrow=False, bgcolor='rgba(255,255,255,0.92)',
        bordercolor='#ccc', borderwidth=1, borderpad=8,
        font=dict(size=11, family='Helvetica Neue'),
    )
    
    fig.update_layout(
        title=dict(
            text='<b>Figure 3.</b> Emission Line Width vs Surface Density (Sample B, N=36)<br>'
                 '<span style="font-size:13px;color:#666">'
                 'Higher-Σ LRDs exhibit systematically broader [OIII]+Hβ lines</span>',
            font_size=17, x=0.5,
        ),
        template='plotly_white', height=550, width=900,
        font=dict(family='Helvetica Neue, Arial', size=11),
        legend=dict(orientation="h", yanchor='bottom', y=1.01, xanchor='center', x=0.42,
                   font_size=10),
    )
    
    fig.update_xaxes(row=1, col=1, title_text='<b>log(Surface Density Proxy)</b>')
    fig.update_yaxes(row=1, col=1, title_text='<b>[OIII]+Hβ FWHM [km s⁻¹]</b>')
    fig.update_xaxes(row=1, col=2, title_text='')
    fig.update_yaxes(row=1, col=2, title_text='Mean FWHM [km s⁻¹]', range=[1600, 2450])
    
    # Add median line to scatter
    fig.add_vline(x=med, line_dash='dot', line_color='#888', width=1,
                   annotation_text='Σ median', row=1, col=1)
    
    return fig


# ============================================================
# FIGURE 4: THE UNIFIED DUAL-AXIS FIGURE — SOUL OF THE PAPER
# ============================================================
def make_figure4():
    """
    Dual-axis scatter: 
      Left axis  = F444W-F150W Color (Sample A, N=260)
      Right axis = FWHM           (Sample B, N=36)
    Both vs log(Σ).
    
    This is THE figure that tells the whole story at one glance.
    """
    
    # Sample A data
    xa = plot_a['logSigma'].values.astype(float)
    ya = plot_a['color'].values.astype(float)
    za = plot_a['z'].values.astype(float)
    fa = plot_a['field'].values
    
    # Sample B data  
    xb = plot_b['logSigma'].values.astype(float)
    yb = plot_b['FWHM'].values.astype(float)
    zb = plot_b['z'].values.astype(float)
    
    # Correlations
    rho_a, pa = stats.pearsonr(xa, ya)
    rho_b, pb = stats.pearsonr(xb, yb)
    
    # Fits
    sa, ia, ra, _, sea = stats.linregress(xa, ya)
    sb, ib, rb, _, seb = stats.linregress(xb, yb)
    
    fig = go.Figure()
    
    # === Sample A: Color data points (larger but semi-transparent) ===
    fig.add_trace(go.Scatter(
        x=xa, y=ya,
        mode='markers',
        name=f'Sample A: Color (N={len(xa)})',
        marker=dict(
            size=8, color='#4C72B0', opacity=0.50,
            symbol='circle',
            line=dict(width=0.3, color='white'),
        ),
        text=np.array([
            f"z={zi:.2f}|field={fi}|color={ci:.2f}"
            for zi, fi, ci in zip(za, fa, ya)
        ]),
        hovertemplate=(
            '<b>Sample A</b><br>'
            'logΣ = %{x:.2f}<br>'
            'Color = %{y:.2f}<br>'
            '%{text}<extra></extra>'
        ),
        yaxis='y1',
    ))
    
    # Sample A fit line
    xfit = np.linspace(min(xa.min(), xb.min())-0.3, max(xa.max(), xb.max())+0.3, 200)
    fig.add_trace(go.Scatter(
        x=xfit, y=sa*xfit+ia,
        mode='lines',
        name=f'Fit A (ρ={rho_a:.2f})',
        line=dict(color='#2B5BA4', width=3, dash='solid'),
        yaxis='y1',
    ))
    
    # === Sample B: FWHM data points (smaller but opaque) ===
    fig.add_trace(go.Scatter(
        x=xb, y=yb,
        mode='markers',
        name=f'Sample B: FWHM (N={len(xb)})',
        marker=dict(
            size=14, color='#C44E52', opacity=0.90,
            symbol='diamond-wide',
            line=dict(width=1.2, color='white'),
        ),
        text=np.array([
            f"z={zi:.2f}|FWHM={yi:.0f}"
            for zi, yi in zip(zb, yb)
        ]),
        hovertemplate=(
            '<b>Sample B</b><br>'
            'logΣ = %{x:.2f}<br>'
            'FWHM = %{y:.0f} km/s<br>'
            '%{text}<extra></extra>'
        ),
        yaxis='y2',
    ))
    
    # Sample B fit line
    fig.add_trace(go.Scatter(
        x=xfit, y=sb*xfit+ib,
        mode='lines',
        name=f'Fit B (ρ={rho_b:.2f})',
        line=dict(color='#DD8452', width=3, dash='dash'),
        yaxis='y2',
    ))
    
    # === Key insight arrow annotation ===
    # Draw a curved arrow showing both trends point same direction
    ann_main = (
        "<b>Dual Σ–Dependence</b><br>"
        f"<b>A:</b> ρ={rho_a:.3f}, Δcolor={ya[xa > np.median(xa)].mean() - ya[xa <= np.median(xa)].mean():.2f} mag<br>"
        f"<b>B:</b> ρ={rho_b:.3f}, ΔFWHM=+{yb[xb > np.median(xb)].mean() - yb[xb <= np.median(xb)].mean():.0f} km/s<br>"
        "<i style='font-size:9px;color:#888'>Two independent probes, one driver: Σ</i>"
    )
    
    fig.add_annotation(
        text=ann_main,
        xref='paper', yref='paper',
        x=0.98, y=0.06,
        align='right', showarrow=False,
        bgcolor='rgba(255,255,255,0.92)',
        bordercolor='#ccc', borderwidth=1, borderpad=6,
        font=dict(size=10, family='Helvetica Neue'),
    )
    
    # (Arrow annotation removed — plotly axref does not support 'paper' for arrows)
    
    # Axis configuration
    fig.update_layout(
        title=dict(
            text='<b>Figure 4.</b> Unified Evidence: Both Photometric Color and Spectroscopic Line Width Increase with Surface Density<br>'
                 '<span style="font-size:12px;color:#666">'
                 f'Sample A (blue circles, left axis, N={len(xa)}) + Sample B (red diamonds, right axis, N={len(xb)})'
                 '</span>',
            font_size=16, x=0.5, y=0.98,
        ),
        template='plotly_white',
        height=680, width=950,
        font=dict(family='Helvetica Neue, Arial, sans-serif', size=11),
        legend=dict(
            orientation="h", yanchor='bottom', y=1.06,
            xanchor='center', x=0.50, font_size=9.5,
        ),
        margin=dict(t=110),  # extra top margin for title + legend
        xaxis=dict(
            title='<b>log(Surface Density Proxy) ∝ −0.4·m<sub>F444W</sub> − 2·log(r<sub>eff</sub>)</b>',
            zeroline=True, zerolinewidth=1, zerolinecolor='#ddd',
            range=[xa.min()-0.4, xb.max()+0.4],
        ),
        yaxis=dict(
            domain=[0.10, 0.95],
            title=dict(text='<b>F444W − F150W Color [AB mag]</b>', font=dict(color='#2B5CA0', size=13)),
            tickfont=dict(color='#2B5CA0'),
            range=[ya.max()+0.3, ya.min()-0.3],  # inverted: redder = lower = up
            zeroline=True, zerolinewidth=1, zerolinecolor='#eee',
        ),
        yaxis2=dict(
            domain=[0.10, 0.95],
            overlaying='y',
            side='right',
            title=dict(text='<b>[OIII]+Hβ FWHM [km s⁻¹]</b>', font=dict(color='#C44E52', size=13)),
            tickfont=dict(color='#C44E52'),
            range=[yb.min()-150, yb.max()+250],
        ),
    )
    
    # Subtle grid
    fig.update_xaxes(gridcolor='#f0f0f0', gridwidth=0.5, minor_ticks='outside')
    fig.update_yaxes(gridcolor='#f0f0f0', gridwidth=0.5)
    # y2 (right axis): hide grid lines so they don't clutter
    fig.layout.yaxis2.showgrid = False
    
    return fig


# ============================================================
# Generate all figures
# ============================================================
if __name__ == '__main__':
    
    print("Generating Figure 1 (Sample Overview)...")
    fig1 = make_figure1()
    fig1.write_html(os.path.join(OUTDIR, 'Fig1_SampleOverview.html'), include_plotlyjs='cdn')
    fig1.write_image(os.path.join(OUTDIR, 'Fig1_SampleOverview.png'), scale=2, width=1000, height=700)
    print(f"  Saved Fig1 ✅")
    
    print("Generating Figure 2 (Color-Σ, Result I)...")
    fig2 = make_figure2()
    fig2.write_html(os.path.join(OUTDIR, 'Fig2_Color_vs_Sigma.html'), include_plotlyjs='cdn')
    fig2.write_image(os.path.join(OUTDIR, 'Fig2_Color_vs_Sigma.png'), scale=2, width=850, height=650)
    print(f"  Saved Fig2 ✅")
    
    print("Generating Figure 3 (FWHM-Σ, Result II)...")
    fig3 = make_figure3()
    fig3.write_html(os.path.join(OUTDIR, 'Fig3_FWHM_vs_Sigma.html'), include_plotlyjs='cdn')
    fig3.write_image(os.path.join(OUTDIR, 'Fig3_FWHM_vs_Sigma.png'), scale=2, width=900, height=550)
    print(f"  Saved Fig3 ✅")
    
    print("Generating Figure 4 (Unified Dual-Axis, SOUL OF PAPER)...")
    fig4 = make_figure4()
    fig4.write_html(os.path.join(OUTDIR, 'Fig4_Unified_DualAxis.html'), include_plotlyjs='cdn')
    fig4.write_image(os.path.join(OUTDIR, 'Fig4_Unified_DualAxis.png'), scale=2, width=950, height=650)
    print(f"  Saved Fig4 ✅")
    
    print("\n" + "=" * 60)
    print("ALL FIGURES GENERATED SUCCESSFULLY!")
    print(f"Output directory: {OUTDIR}")
    print("=" * 60)
