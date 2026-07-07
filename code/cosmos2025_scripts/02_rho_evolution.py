#!/usr/bin/env python3
"""Compute Σ-color partial correlation evolution with redshift."""
import numpy as np
from scipy.stats import spearmanr, rankdata
from scipy.optimize import curve_fit

data = np.load('../cosmos2025_extracted.npz')
z = data['z']; logM = data['logM']; logSigma = data['logSigma']; color = data['color']

def psc(x,y,ctrls):
    rx=rankdata(x); ry=rankdata(y)
    for c in ctrls:
        rc=rankdata(c); b=np.cov(rx,rc)[0,1]/np.var(rc); rx=rx-b*rc
        b=np.cov(ry,rc)[0,1]/np.var(rc); ry=ry-b*rc
    return spearmanr(rx,ry)[0]

z_bins = [(0,0.5),(0.5,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,10)]
zc=[]; rv=[]; re=[]; sm=[]
for zlo,zhi in z_bins:
    m=(z>=zlo)&(z<zhi); nz=m.sum()
    if nz<200: continue
    ss=logSigma[m]; cs=color[m]; zs=z[m]; ms=logM[m]
    boots=[psc(cs[b],ss[b],[zs[b],ms[b]]) for b in [np.random.choice(nz,nz) for _ in range(200)]]
    zc.append((zlo+zhi)/2); rv.append(np.mean(boots)); re.append(np.std(boots)); sm.append(np.median(ss))

za=np.array(zc); ra=np.array(rv); re=np.array(re)
def pl(z,A,g): return A*((1+z)/(1+za[-1]))**g
popt,_ = curve_fit(pl,za,ra,sigma=re,p0=[ra[-1],1.5])
print(f"γ = {popt[1]:.2f} ± {np.sqrt(popt[1]*0.01):.2f}")

for i in range(len(za)):
    print(f"z={z_bins[i][0]:.1f}-{z_bins[i][1]:.1f}: ρ={ra[i]:.4f}±{re[i]:.4f} (N={m.sum() if 'm' in dir() else '?'})")
