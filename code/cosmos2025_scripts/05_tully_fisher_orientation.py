#!/usr/bin/env python3
"""Tully-Fisher + Orientation test for Σ-color correlation."""
import numpy as np
from scipy.stats import spearmanr, rankdata

data = np.load('../cosmos2025_with_env.npz')
z = data['z']; logM = data['logM']; logSigma = data['logSigma']; color = data['color']

def psc(x,y,c):
    rx=rankdata(x); ry=rankdata(y)
    for cc in c:
        rc=rankdata(cc); b=np.cov(rx,rc)[0,1]/np.var(rc); rx=rx-b*rc
        b=np.cov(ry,rc)[0,1]/np.var(rc); ry=ry-b*rc
    return spearmanr(rx,ry)[0]

alpha_tf=3.5; logM_ref=10.2; logV_ref=np.log10(200)
logV_rot=logV_ref+(logM-logM_ref)/alpha_tf; V_rot=10**logV_rot

print("=== Within V_rot bins ===")
for vlo,vhi in [(0,100),(100,150),(150,200),(200,300)]:
    m=(V_rot>=vlo)&(V_rot<vhi)
    if m.sum()<100: continue
    r=psc(color[m],logSigma[m],[z[m],logM[m]])
    print(f"V={vlo}-{vhi}: ρ={r:.4f} (n={m.sum()})")

print("\n=== Orientation test (needs axratio from FITS) ===")
from astropy.io import fits
hdul = fits.open('../COSMOSWeb_mastercatalog_v1.1.fits', memmap=True)
phot=hdul[1].data; lephare=hdul[2].data
axratio=phot['axratio_sersic']

zfinal=lephare['zfinal']; mass_med=lephare['mass_med']
obj_type=lephare['Type'] if 'Type' in lephare.columns.names else np.zeros(len(zfinal))
r_deg=phot['radius_sersic']; warn=phot['warn_flag']
m4=phot['mag_auto_f444w']; m1=phot['mag_auto_f150w']

mask=((obj_type==0)&(zfinal>0)&(zfinal<99)&(mass_med>-90)&(r_deg>0)&(warn==0)&
      np.isfinite(axratio)&(axratio>0)&(axratio<=1))
z2=zfinal[mask]; logM2=mass_med[mask]; ax2=axratio[mask]
r_deg2=r_deg[mask]; m42=m4[mask]; m12=m1[mask]

from astropy.cosmology import Planck18; import astropy.units as u
r_rad2=np.deg2rad(r_deg2); DA2=Planck18.angular_diameter_distance(z2).to(u.pc).value
r_pc2=r_rad2*DA2; Sigma2=10**logM2/(np.pi*r_pc2**2)
logSigma2=np.log10(Sigma2); color2=0.4*(m12-m42)
logV2=np.log10(200)+(logM2-10.2)/3.5; V2=10**logV2

good=np.isfinite(logSigma2)&np.isfinite(color2)
z2=z2[good]; logM2=logM2[good]; logSigma2=logSigma2[good]; color2=color2[good]
ax2=ax2[good]; V2=V2[good]

mv=(V2>=150)&(V2<250); edge=mv&(ax2<0.3); face=mv&(ax2>0.7)
re=psc(color2[edge],logSigma2[edge],[z2[edge],logM2[edge]])
rf=psc(color2[face],logSigma2[face],[z2[face],logM2[face]])
print(f"V=150-250: ρ_edge={re:.4f} (n={edge.sum()}), ρ_face={rf:.4f} (n={face.sum()})")
print(f"Dust prediction: ρ_edge >> ρ_face. Observed: {'FAILS dust' if re<rf else 'consistent'}")

hdul.close()
