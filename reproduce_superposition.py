#!/usr/bin/env python3
"""Reproducible fixed-voltage superposition check from 50 paired signed I–V exports."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parent
OUT=BASE/'derived'
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(BASE))
from reproduce_rank_analysis import parse,prepare
SIX={'1a','1b','1e','6a','6b','6e'}
VPOINTS=[-0.5,0.,0.08,0.10,0.12]
entries=[]; curves=[]
V=np.linspace(-1,0.2,1201)
for rec in pd.read_csv(BASE/'PAIRING_MAP.csv').itertuples(index=False):
 vd,dark=prepare(*parse(BASE/rec.dark_file)); vl,light=prepare(*parse(BASE/rec.light_file))
 if vd.min()>-1 or vl.min()>-1 or vd.max()<0.2 or vl.max()<0.2:
  raise ValueError('Insufficient voltage domain '+str(rec.contact))
 iph=lambda vv: np.interp(vv,vl,light)-np.interp(vv,vd,dark)
 ref=float(iph(-1))
 g='six' if rec.rpm==1000 and rec.contact in SIX else ('other24' if rec.rpm==1000 else '2000')
 dr=float(np.interp(0,vd,dark)-np.interp(-1,vd,dark))
 phr=float(iph(0)-ref)
 row={'rpm':rec.rpm,'contact':rec.contact,'group':g,'Iph_minus1_uA':ref*1e6,'dark_delta_minus1_to_0_uA':dr*1e6,'dark_delta_fraction':abs(dr/phr) if phr else np.nan}
 for v in VPOINTS:
  key=f'R_{v:+.2f}V'
  row[key]=float(iph(v)/ref) if abs(ref)>1e-12 else np.nan
 entries.append(row)
 if rec.rpm==1000:
  y=iph(V)/ref
  for vv,yy in zip(V,y): curves.append((rec.contact,g,float(vv),float(yy)))
metrics=pd.DataFrame(entries)
metrics.to_csv(OUT/'V15_15_Contact_Superposition_Check.csv',index=False)
curves=pd.DataFrame(curves,columns=['contact','group','V','signed_ratio'])
sm=metrics[metrics.rpm==1000].groupby('group').agg(n=('contact','size'),**{p+'_median':(p,'median') for p in ['R_-0.50V','R_+0.00V','R_+0.08V','R_+0.10V','R_+0.12V']},**{p+'_min':(p,'min') for p in ['R_+0.08V','R_+0.10V','R_+0.12V']},**{p+'_max':(p,'max') for p in ['R_+0.08V','R_+0.10V','R_+0.12V']})
sm.to_csv(OUT/'V15_15_Group_Summary.csv')
print(sm.T.to_string())
fig,ax=plt.subplots(figsize=(7.0,4.2))
for grp,alpha,ls in [('other24',0.24,'-'),('six',0.65,'--')]:
 for label,gp in curves[curves.group==grp].groupby('contact'):
  ax.plot(gp.V.values,gp.signed_ratio.values,linewidth=0.75 if grp=='other24' else 1.05,alpha=alpha,ls=ls,color='0.48' if grp=='other24' else '0.07')
 x=curves[curves.group==grp].groupby('V',sort=True).signed_ratio.median()
 ax.plot(x.index,x.values,linewidth=2.2,color='0.15' if grp=='six' else '0.38',ls='--' if grp=='six' else '-',label=('Six low-response contacts' if grp=='six' else 'Remaining 24 contacts')+' (median)')
ax.axvline(0,color='0.7',linewidth=0.9)
ax.axhline(1,color='0.65',linewidth=0.7,linestyle=':')
ax.axvline(.1,color='0.6',linewidth=0.8,linestyle=':')
ax.set_xlim(-1,0.2); ax.set_ylim(-0.08,1.2)
ax.set_xlabel('Applied voltage (V)')
ax.set_ylabel(r'Signed normalized photocurrent  $I_{ph}(V)/I_{ph}(-1\,\mathrm{V})$')
ax.legend(loc='upper right',fontsize=8)
ax.grid(alpha=.12)
fig.tight_layout()
fig.savefig(OUT/'Figure_S12_Superposition_1000rpm.png',dpi=250)
plt.close(fig)
