#!/usr/bin/env python3
"""CIGALE SFH: sSFR vs Σ — Jeans mass suppression test."""
from astropy.io import fits
from astropy.cosmology import Planck18
import astropy.units as u
import numpy as np
from scipy.stats import spearmanr, rankdata

hdul = fits.open('../COSMOSWeb_mastercatalog_v1.1.fits', memmap=True)
phot=hdul[1].data; lephare=hdul[2].data; cg=hdul[4].data

zfinal=lephare['zfinal']; mass_med=lephare['mass_med']
obj_type=lephare['Type'] if 'Type' in lephare.columns.names else np.zeros(len(zfinal))
r_deg=phot['radius_sersic']; warn=phot['warn_flag']
m444=phot['mag_auto_f444w']; m150=phot['mag_auto_f150w']
sfr1=cg['sfh_sfr_bin1']; sfr100=cg['sfr_100myr']; age=cg['age_form']

mask=((obj_type==0)&(zfinal>0)&(zfinal<99)&(mass_med>-90)&(r_deg>0)&(warn==0)&
      np.isfinite(sfr1)&(sfr1>0)&np.isfinite(sfr100)&(sfr100>0)&np.isfinite(age)&(age>0))
z=zfinal[mask]; logM=mass_med[mask]; r=r_deg[mask]
m4=m444[mask]; m1=m150[mask]; sfr1=sfr1[mask]; sfr100=sfr100[mask]; age=age[mask]

r_rad=np.deg2rad(r); DA=Planck18.angular_diameter_distance(z).to(u.pc).value
r_pc=r_rad*DA; Sigma=10**logM/(np.pi*r_pc**2); logSigma=np.log10(Sigma)
ssfr1=np.log10(sfr1/(10**logM)); ssfr100=np.log10(sfr100/(10**logM))

good=np.isfinite(logSigma)&np.isfinite(ssfr1)&np.isfinite(ssfr100)
z=z[good]; logM=logM[good]; logSigma=logSigma[good]
ssfr1=ssfr1[good]; ssfr100=ssfr100[good]; age=age[good]

def psc(x,y,c):
    rx=rankdata(x); ry=rankdata(y)
    for cc in c:
        rc=rankdata(cc); b=np.cov(rx,rc)[0,1]/np.var(rc); rx=rx-b*rc
        b=np.cov(ry,rc)[0,1]/np.var(rc); ry=ry-b*rc
    return spearmanr(rx,ry)[0]

print(f"ρ(Σ, sSFR_bin1 | z,M) = {psc(logSigma, ssfr1, [z,logM]):+.4f}")
print(f"ρ(Σ, sSFR_100  | z,M) = {psc(logSigma, ssfr100, [z,logM]):+.4f}")
print(f"ρ(Σ, age       | z,M) = {psc(logSigma, age, [z,logM]):+.4f}")

young=age<200
print(f"\nYoung (<200Myr): ρ(Σ,sSFR|z,M) = {psc(logSigma[young], ssfr1[young], [z[young],logM[young]]):+.4f} (n={young.sum()})")

for slo,shi in [(0,1),(1,2),(2,3),(3,4),(4,8)]:
    m=(logSigma>=slo)&(logSigma<shi)
    if m.sum()>50: print(f"logΣ {slo}-{shi}: med sSFR={np.median(ssfr1[m]):.2f}, med age={np.median(age[m]):.0f}Myr")

hdul.close()
