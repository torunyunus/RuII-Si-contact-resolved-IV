#!/usr/bin/env python3
"""Region-specific QC and signed-photocurrent sensitivity analysis.
Run: python reproduce_step3.py <deposit_directory> --out <output_directory>
The deposit must contain PAIRING_MAP.csv and matched raw_iv files.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_dat(filename):
    f = pd.read_csv(filename, header=None, names=['V', 'abs_current_A', 'signed_current_A'])
    return f.apply(pd.to_numeric, errors='coerce').dropna()


def aligned(f):
    t = f.groupby('V', sort=True)['signed_current_A'].mean()
    return t.index.to_numpy(), t.to_numpy()


def roots_linear(grid, y):
    j = np.where(y[:-1] * y[1:] < 0)[0]
    return [float(grid[k] - y[k] * (grid[k+1]-grid[k]) / (y[k+1]-y[k])) for k in j]


def reproduce(deposit: Path, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(deposit / 'PAIRING_MAP.csv')
    details = []
    alerts = []
    curves = {}
    vrev = np.linspace(-1, 0, 10001)
    for p in pairs.itertuples(index=False):
        d = read_dat(deposit / p.dark_file)
        l = read_dat(deposit / p.light_file)
        vd, di = aligned(d)
        vl, li = aligned(l)
        dark = lambda v: np.interp(v, vd, di)
        light = lambda v: np.interp(v, vl, li)
        ip = lambda v: light(v) - dark(v)
        names=[]
        for kind, frame in [('dark', d), ('illuminated', l)]:
            u = frame.V.to_numpy(); delta=np.diff(u)
            negatives = np.where(delta < 0)[0]
            if len(negatives):
                names.append('nonmonotonic_voltage')
                alerts.append(dict(rpm=p.rpm,contact=p.contact,trace=kind,flag='nonmonotonic_voltage',
                                   count=len(negatives),min_flag_voltage_V=float(np.min(u[negatives])),
                                   max_flag_voltage_V=float(np.max(u[negatives]))))
            near = frame[frame.signed_current_A.abs() >= .099]
            if len(near):
                names.append('near_100mA_compliance')
                alerts.append(dict(rpm=p.rpm,contact=p.contact,trace=kind,flag='near_100mA_compliance',
                                   count=len(near),min_flag_voltage_V=float(near.V.min()),
                                   max_flag_voltage_V=float(near.V.max())))
        iph_m1=ip(-1.0)
        if iph_m1 > 0:
            names.append('positive_Iph_minus1V')
        grid = np.linspace(0,.5,3001)
        power = -grid * light(grid)
        pmax = max(0, float(np.max(power))) * 1e6
        trace_ph = ip(vrev)
        rcross = roots_linear(vrev, trace_ph)
        row = dict(rpm=int(p.rpm), contact=p.contact, QA_flags=';'.join(sorted(set(names))),
                   QC_clear=(len(names)==0),n_dark=len(d),n_light=len(l),
                   I_dark_minus1V_uA=float(dark(-1.0)*1e6),
                   I_light_minus1V_uA=float(light(-1.0)*1e6),
                   I_ph_minus1V_uA=float(iph_m1*1e6),
                   I_ph_zeroV_uA=float(ip(0)*1e6),
                   I_light_zeroV_uA=float(light(0)*1e6),
                   reverse_range_sign_changes=len(rcross),
                   first_Iph_zero_crossing_V=rcross[0] if rcross else np.nan,
                   Pmax_uW=pmax)
        details.append(row)
        if p.rpm==2000:
            curves[p.contact] = trace_ph * 1e6
    df=pd.DataFrame(details)
    df.to_csv(out/'step3_contact_QC_and_photocurrent.csv',index=False,float_format='%.10g')
    alarms=pd.DataFrame(alerts)
    alarms.to_csv(out/'step3_voltage_alerts.csv',index=False,float_format='%.10g')
    a=df[df.rpm==2000]
    subsets=[('All paired contacts',a,'abs_Iph'),
             ('Excluding positive Iph at -1 V',a[a.I_ph_minus1V_uA < 0],'abs_Iph'),
             ('Excluding high-forward-current alerts',a[~a.QA_flags.str.contains('near_100mA_compliance')],'abs_Iph'),
             ('Excluding voltage-order alerts',a[~a.QA_flags.str.contains('nonmonotonic_voltage')],'abs_Iph'),
             ('Excluding all three QC categories',a[a.QC_clear],'abs_Iph'),
             ('All: illuminated current magnitude at -1 V',a,'abs_Ilight')]
    rows=[]
    for name,g,metric in subsets:
        x=g.I_ph_minus1V_uA.abs() if metric=='abs_Iph' else g.I_light_minus1V_uA.abs()
        rows.append(dict(subset=name,n=len(g),reverse_bias_descriptor='|Iph(-1 V)|' if metric=='abs_Iph' else '|Ilight(-1 V)|',
                         spearman_rho=spearmanr(x,g.Pmax_uW).statistic,
                         median_Pmax_uW=g.Pmax_uW.median(),contacts=';'.join(g.contact)))
    s=pd.DataFrame(rows)
    s.to_csv(out/'step3_2000rpm_rank_sensitivity.csv',index=False,float_format='%.9g')
    # Publication figure: signed dark-subtracted photocurrent, no absolute-value transformation.
    fig,ax=plt.subplots(figsize=(7.4,4.4))
    for contact,y in curves.items():
        if contact not in ('1d','1f'):
            ax.plot(vrev,y,lw=.85,alpha=.42,color='#9e9e9e')
    ax.plot(vrev,curves['1d'],lw=2.2,label='Contact 1d',color='#00598a')
    ax.plot(vrev,curves['1f'],lw=2.2,label='Contact 1f',color='#ba4b31',linestyle='--')
    ax.axhline(0,lw=0.8,color='black',alpha=.65)
    ax.axvline(-1,lw=.7,color='black',ls=':',alpha=.5)
    ax.set(xlabel='Applied voltage (V)',ylabel=r'Signed $I_{\rm ph}=I_{\rm light}-I_{\rm dark}$ ($\mu$A)',
           xlim=(-1,0))
    ax.legend(loc='upper right',frameon=False,fontsize=9)
    ax.grid(alpha=.13)
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(out/'Figure_S11_Signed_Photocurrent_2000rpm.png',dpi=260)
    plt.close(fig)
    print(s[['subset','n','spearman_rho','median_Pmax_uW']].to_string(index=False))
    print('2000rpm flagged',len(a[~a.QC_clear]),'high-current alerts',alarms[(alarms.rpm==2000)&(alarms.flag=='near_100mA_compliance')]['count'].sum())
    print('1d,1f crossings',a[a.contact.isin(['1d','1f'])][['contact','I_ph_minus1V_uA','I_ph_zeroV_uA','first_Iph_zero_crossing_V']].to_string(index=False))
    return df,alarms,s

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('deposit',type=Path);ap.add_argument('--out',type=Path,default=Path('step3_output'))
    ar=ap.parse_args();reproduce(ar.deposit,ar.out)
