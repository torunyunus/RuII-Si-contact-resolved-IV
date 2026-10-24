#!/usr/bin/env python3
"""Reproduce contact-level reverse-bias/Pmax ranks from the deposited raw I-V files."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
BASE=Path(__file__).resolve().parent

def parse(path):
    rows=[]
    for line in path.read_bytes().decode("latin1").splitlines():
        parts=line.strip().split(",")
        if len(parts)>=3:
            try: rows.append(tuple(map(float,parts[:3])))
            except ValueError: pass
    a=np.asarray(rows); return a[:,0],a[:,2]

def prepare(v,i):
    order=np.argsort(v); v,i=v[order],i[order]; uv=np.unique(v)
    return uv,np.array([i[v==x].mean() for x in uv])

pairs=pd.read_csv(BASE/"PAIRING_MAP.csv")
records=[]
for r in pairs.itertuples(index=False):
    vd,dark=prepare(*parse(BASE/r.dark_file)); vl,light=prepare(*parse(BASE/r.light_file))
    iph=abs(np.interp(-1.0,vl,light)-np.interp(-1.0,vd,dark))
    grid=np.linspace(0.0,0.5,3001); pmax=max(0.0,np.max(-grid*np.interp(grid,vl,light)))
    records.append({"rpm":r.rpm,"contact":r.contact,"Iph_minus1_uA":iph*1e6,"Pmax_uW":pmax*1e6,"dark_file":r.dark_file,"light_file":r.light_file})
data=pd.DataFrame(records); all_contacts=[]; summary=[]
for rpm,g in data.groupby("rpm",sort=True):
    g=g.copy().reset_index(drop=True); n=len(g)
    g["reverse_bias_rank"]=g.Iph_minus1_uA.rank(ascending=False,method="average")
    g["power_rank"]=g.Pmax_uW.rank(ascending=False,method="average")
    g["absolute_rank_displacement"]=(g.reverse_bias_rank-g.power_rank).abs()
    g["normalized_rank_displacement_pct"]=g.absolute_rank_displacement/(n-1)*100
    k=int(np.ceil(0.25*n)); tr=set(g.nsmallest(k,"reverse_bias_rank").contact); tp=set(g.nsmallest(k,"power_rank").contact)
    f=g.absolute_rank_displacement.sum(); fmax=n*n/2 if n%2==0 else (n*n-1)/2; rho,p=spearmanr(g.Iph_minus1_uA,g.Pmax_uW)
    summary.append({"rpm":rpm,"n":n,"spearman_rho":rho,"nominal_p":p,"median_normalized_rank_displacement_pct":g.normalized_rank_displacement_pct.median(),"q1_pct":g.normalized_rank_displacement_pct.quantile(.25),"q3_pct":g.normalized_rank_displacement_pct.quantile(.75),"normalized_footrule":f/fmax,"top_quartile_k":k,"top_quartile_overlap":len(tr&tp),"overlapping_contacts":";".join(sorted(tr&tp))})
    all_contacts.append(g)
pd.concat(all_contacts,ignore_index=True).to_csv(BASE/"contact_metrics.csv",index=False)
pd.DataFrame(summary).to_csv(BASE/"rank_summary.csv",index=False)
print(pd.DataFrame(summary).to_string(index=False))

# Contact-deletion sensitivity of the finite-set rank relationship
loco=[]
for rpm,g in data.groupby("rpm",sort=True):
    for omitted in g.contact:
        subset=g.loc[g.contact!=omitted]
        rho,_=spearmanr(subset.Iph_minus1_uA,subset.Pmax_uW)
        loco.append({"rpm":rpm,"removed_contact":omitted,"n_remaining":len(subset),"spearman_rho":rho})
(Path(__file__).resolve().parent/"derived").mkdir(exist_ok=True)
pd.DataFrame(loco).to_csv(Path(__file__).resolve().parent/"derived"/"LOCO_spearman.csv",index=False)
