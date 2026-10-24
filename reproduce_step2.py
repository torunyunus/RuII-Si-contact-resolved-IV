#!/usr/bin/env python3
"""Recompute finite-bias photocurrent and dark-current controls from a deposited I–V archive.
Usage: python reproduce_step2.py /path/to/deposit --out ./step2_outputs
Archive must contain PAIRING_MAP.csv and raw_iv folders; requires numpy, pandas.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

def load(path):
    a=pd.read_csv(path,header=None,comment='#').iloc[:,:3].apply(pd.to_numeric,errors='coerce').dropna().to_numpy()
    v,i=a[:,0],a[:,2]; ix=np.argsort(v);v=v[ix];i=i[ix]
    uv,inverse=np.unique(v,return_inverse=True)
    accum=np.bincount(inverse,weights=i)/np.bincount(inverse)
    return uv,accum

def calculate(deposit):
    grid=np.linspace(-1,0,1001)
    target={'1a','1b','1e','6a','6b','6e'}
    out=[];traces=[]
    pairs=pd.read_csv(deposit/'PAIRING_MAP.csv')
    for p in pairs.itertuples(index=False):
        # Each record uses its own interpolants. No contact files are pooled.
        vd,dark=load(deposit/p.dark_file);vl,light=load(deposit/p.light_file)
        dl=np.interp(grid,vd,dark);ll=np.interp(grid,vl,light);ph=ll-dl
        norm=np.abs(ph)/abs(ph[0]);lnorm=np.abs(ll)/abs(ll[0])
        crossed=np.where(norm<=0.5)[0]
        v50=np.nan
        if len(crossed):
            j=int(crossed[0]);v50=-1 if j==0 else float(np.interp(.5,[norm[j],norm[j-1]],[grid[j],grid[j-1]]))
        delta_ph=ph[-1]-ph[0];delta_dark=dl[-1]-dl[0]
        six=(p.rpm==1000 and p.contact in target)
        out.append(dict(rpm=p.rpm,contact=p.contact,group=('six' if six else 'other24') if p.rpm==1000 else '2000',
            R_minus075=np.interp(-.75,grid,norm),R_minus050=np.interp(-.5,grid,norm),R_minus025=np.interp(-.25,grid,norm),R_zero=norm[-1],
            light_R_minus05=np.interp(-.5,grid,lnorm),V50_V=v50,delta_I_ph_uA=delta_ph*1e6,delta_I_dark_uA=delta_dark*1e6,
            dark_to_ph_delta_fraction=abs(delta_dark)/abs(delta_ph)))
        if p.rpm==1000:
            for v,r in zip(grid,norm): traces.append((p.contact, 'six' if six else 'other24',v,r))
    return pd.DataFrame(out),pd.DataFrame(traces,columns=['contact','group','voltage_V','normalized_photocurrent'])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('deposit',type=Path);ap.add_argument('--out',type=Path,default=Path('step2_outputs'));args=ap.parse_args()
    m,c=calculate(args.deposit);args.out.mkdir(parents=True,exist_ok=True)
    m.to_csv(args.out/'finite_bias_metrics.csv',index=False);c.to_csv(args.out/'normalized_curves_1000.csv',index=False)
    for g in ['six','other24']:
        a=m[m.group==g]
        print(g,'n',len(a),'R(-0.5) median',a.R_minus050.median(),'range',a.R_minus050.min(),a.R_minus050.max(),
              'V50 found',a.V50_V.notna().sum())
    s=m[m.group=='six']
    print('Six-contact max dark-current fraction:',s.dark_to_ph_delta_fraction.max())

if __name__=='__main__':main()
