#!/usr/bin/env python3
"""Environmental density modulation via 10th nearest neighbor."""
import numpy as np
from scipy.stats import spearmanr, rankdata
from sklearn.neighbors import NearestNeighbors

data = np.load('../cosmos2025_extracted.npz')
z = data['z']; logM = data['logM']; logSigma = data['logSigma']; color = data['color']

# Need to recompute with RA/Dec from FITS
from astropy.io import fits
hdul = fits.open('../COSMOSWeb_mastercatalog_v1.1.fits', memmap=True)
phot = hdul[1].data; lephare = hdul[2].data

zfinal=lephare['zfinal']; mass_med=lephare['mass_med']
obj_type=lephare['Type'] if 'Type' in lephare.columns.names else np.zeros(len(zfinal))
r_deg=phot['radius_sersic']; warn=phot['warn_flag']
m444=phot['mag_auto_f444w']; m150=phot['mag_auto_f150w']
ra_arr=phot['ra']; dec_arr=phot['dec']

mask=((obj_type==0)&(zfinal>0)&(zfinal<99)&(mass_med>-90)&(r_deg>0)&(warn==0)&
      np.isfinite(ra_arr)&np.isfinite(dec_arr))
z2=zfinal[mask]; logM2=mass_med[mask]; r_deg2=r_deg[mask]
m4442=m444[mask]; m1502=m150[mask]; ra=ra_arr[mask]; dec=dec_arr[mask]

from astropy.cosmology import Planck18; import astropy.units as u
r_rad2=np.deg2rad(r_deg2); DA2=Planck18.angular_diameter_distance(z2).to(u.pc).value
r_pc2=r_rad2*DA2; Sigma2=10**logM2/(np.pi*r_pc2**2)
logSigma2=np.log10(Sigma2); color2=0.4*(m1502-m4442)

valid=np.isfinite(logSigma2)&np.isfinite(color2)
z2=z2[valid]; logM2=logM2[valid]; logSigma2=logSigma2[valid]; color2=color2[valid]
ra=ra[valid]; dec=dec[valid]

ra_med=np.median(ra); dec_med=np.median(dec)
dRA=(ra-ra_med)*np.cos(np.deg2rad(dec_med)); dDec=dec-dec_med
coords=np.column_stack([dRA,dDec])

nn=NearestNeighbors(n_neighbors=11,algorithm='kd_tree'); nn.fit(coords)
distances,_=nn.kneighbors(coords); d10=distances[:,10]
Sigma_env=10/(np.pi*d10**2); logSigma_env=np.log10(Sigma_env)

print(f"logΣ_env median: {np.median(logSigma_env):.2f}")

def psc(x,y,c):
    rx=rankdata(x); ry=rankdata(y)
    for cc in c:
        rc=rankdata(cc); b=np.cov(rx,rc)[0,1]/np.var(rc); rx=rx-b*rc
        b=np.cov(ry,rc)[0,1]/np.var(rc); ry=ry-b*rc
    return spearmanr(rx,ry)[0]

for zlo,zhi,label in [(7,9,'z=7-9'),(0,1,'z=0-1')]:
    mz=(z2>=zlo)&(z2<zhi); senv=logSigma_env[mz]; p50=np.median(senv)
    r_low=psc(color2[mz][senv<=p50],logSigma2[mz][senv<=p50],[z2[mz][senv<=p50],logM2[mz][senv<=p50]])
    r_hi=psc(color2[mz][senv>p50],logSigma2[mz][senv>p50],[z2[mz][senv>p50],logM2[mz][senv>p50]])
    print(f"{label}: low-env ρ={r_low:.4f}, high-env ρ={r_hi:.4f}, Δ={r_hi-r_low:+.4f}")

np.savez('../cosmos2025_with_env.npz', z=z2, logM=logM2, logSigma=logSigma2,
         color=color2, ra=ra, dec=dec, logSigma_env=logSigma_env, n=len(z2))
print("Saved with env density")
hdul.close()
