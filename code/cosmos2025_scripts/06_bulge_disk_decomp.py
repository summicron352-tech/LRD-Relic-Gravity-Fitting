#!/usr/bin/env python3
"""B+D decomposition: bulge vs disk Σ driver analysis."""
from astropy.io import fits
from astropy.cosmology import Planck18
import astropy.units as u
import numpy as np
from scipy.stats import spearmanr, rankdata

hdul = fits.open('../COSMOSWeb_mastercatalog_v1.1.fits', memmap=True)
phot=hdul[1].data; lephare=hdul[2].data; bd=hdul[6].data

zfinal=lephare['zfinal']; mass_med=lephare['mass_med']
obj_type=lephare['Type'] if 'Type' in lephare.columns.names else np.zeros(len(zfinal))
r_s=phot['radius_sersic']; warn=phot['warn_flag']
m444=phot['mag_auto_f444w']; m150=phot['mag_auto_f150w']

b_flux=bd['flux_model_bulge_f444w']; d_flux=bd['flux_model_disk_f444w']
b_r=bd['bulge_radius_deg']; d_r=bd['disk_radius_deg']

bt=np.full(len(b_flux),np.nan)
v=(b_flux>0)&(d_flux>0)&np.isfinite(b_flux)&np.isfinite(d_flux)
bt[v]=b_flux[v]/(b_flux[v]+d_flux[v])

mask=((obj_type==0)&(zfinal>0)&(zfinal<99)&(mass_med>-90)&(r_s>0)&(warn==0)&
      np.isfinite(bt)&np.isfinite(b_r)&np.isfinite(d_r))
z=zfinal[mask]; logM=mass_med[mask]; r_deg=r_s[mask]
bt_vals=bt[mask]; b_rad=b_r[mask]; d_rad=d_r[mask]
m4=m444[mask]; m1=m150[mask]

r_rad=np.deg2rad(r_deg); DA=Planck18.angular_diameter_distance(z).to(u.pc).value
r_pc=r_rad*DA; Sigma=10**logM/(np.pi*r_pc**2)
logSigma=np.log10(Sigma); color=0.4*(m1-m4)

b_pc=np.deg2rad(b_rad)*DA; d_pc=np.deg2rad(d_rad)*DA
M_b=10**logM*bt_vals; M_d=10**logM*(1-bt_vals)
logSigma_b=np.log10(np.clip(M_b/(np.pi*b_pc**2),1e-2,1e15))
logSigma_d=np.log10(np.clip(M_d/(np.pi*d_pc**2),1e-2,1e15))

good=np.isfinite(logSigma)&np.isfinite(color)&np.isfinite(logSigma_b)&np.isfinite(logSigma_d)
z=z[good]; logM=logM[good]; logSigma=logSigma[good]; color=color[good]
bt_vals=bt_vals[good]; logSigma_b=logSigma_b[good]; logSigma_d=logSigma_d[good]

def psc(x,y,c):
    rx=rankdata(x); ry=rankdata(y)
    for cc in c:
        rc=rankdata(cc); b=np.cov(rx,rc)[0,1]/np.var(rc); rx=rx-b*rc
        b=np.cov(ry,rc)[0,1]/np.var(rc); ry=ry-b*rc
    return spearmanr(rx,ry)[0]

for blo,bhi,label in [(0,0.2,'PureDisk'),(0.2,0.4,'DiskDom'),(0.4,0.6,'Interm'),
                       (0.6,0.8,'BulgeDom'),(0.8,1.0,'PureBulge')]:
    m=(bt_vals>=blo)&(bt_vals<bhi)
    if m.sum()<100: continue
    r=psc(color[m],logSigma[m],[z[m],logM[m]])
    print(f"{label}: ρ_ctrl={r:+.4f} (n={m.sum()})")

for blo,bhi,label in [(0,0.3,'DiskDom'),(0.3,0.7,'Interm'),(0.7,1.0,'BulgeDom')]:
    m=(bt_vals>=blo)&(bt_vals<bhi)
    rb=psc(color[m],logSigma_b[m],[z[m],logM[m]])
    rd=psc(color[m],logSigma_d[m],[z[m],logM[m]])
    print(f"{label} head-to-head: ρ_b={rb:+.4f} vs ρ_d={rd:+.4f} → {'BULGE' if rb>rd else 'DISK'}")

hdul.close()
