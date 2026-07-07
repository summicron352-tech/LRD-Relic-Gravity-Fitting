#!/usr/bin/env python3
"""Matched-pair reversal rate analysis with BallTree spatial matching."""
import numpy as np
from scipy.stats import spearmanr
from sklearn.neighbors import BallTree

data = np.load('../cosmos2025_with_env.npz')
z = data['z']; logSigma = data['logSigma']; color = data['color']
ra = data['ra']; dec = data['dec']

paired_z=[(0,1),(1,3),(3,5),(5,7),(7,9)]
for zlo,zhi in paired_z:
    m=(z>=zlo)&(z<zhi); idx=np.where(m)[0]; n=len(idx)
    if n<1000: continue
    ns=min(3000,n); np.random.seed(42); sidx=np.random.choice(idx,ns,replace=False)
    all_ra=ra[m]; all_dec=dec[m]; all_z=z[m]; all_ls=logSigma[m]; all_cl=color[m]
    idx_map={orig:loc for loc,orig in enumerate(idx)}
    t=BallTree(np.column_stack([all_ra,all_dec]),metric='haversine')
    npairs=0; nrev=0
    for i in range(ns):
        si_orig=sidx[i]; si_loc=idx_map[si_orig]
        _,neigh=t.query([np.column_stack([all_ra,all_dec])[si_loc]],k=4)
        for ji in neigh[0][1:]:
            j_orig=idx[ji]
            if abs(all_z[si_loc]-all_z[ji])>0.3: continue
            if abs(all_ls[si_loc]-all_ls[ji])<0.3: continue
            npairs+=1
            if (all_ls[si_loc]>all_ls[ji])!=(all_cl[si_loc]>all_cl[ji]): nrev+=1
            if npairs>=3000: break
        if npairs>=3000: break
    rate=nrev/npairs*100 if npairs>0 else 0
    print(f"z={zlo}-{zhi}: {npairs} pairs, {rate:.1f}% reversal")
