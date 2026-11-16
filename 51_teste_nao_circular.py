#!/usr/bin/env python3
"""Noncircular wall contact: protocol 50, one fixed benchmark."""
from pathlib import Path
import csv, json, math, statistics, time
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

P=Path(__file__).resolve().parent
RX=.36; RY=.28209479177387814; Q=.045; L=.75; K=20.; W=.002; A=.1; TAU=1e-8
WIDTHS=(.02,.055,.13,.31)

def base(theta):return RY*math.cos(theta)+Q*math.cos(2*theta)

def cases():
    rows=[]
    for b in (0.,A):
        for side in ('zero','pi'):
            for e in WIDTHS:
                angle=e if side=='zero' else math.pi-e
                rows.append(dict(id=f'b{b}_{side}_e{e}',h=b-base(angle),boundary=b,side=side,width=e))
    for h in (.05,RY-Q+.12):rows.append(dict(id=f'control_h{h}',h=h,boundary='',side='control',width=''))
    return rows

def calculate(c,method):
    start=time.perf_counter_ns(); h=c['h']; root_evals=0; law_evals=0
    def raw_gap(t):return h+RY*math.cos(t)+Q*math.cos(2*t)
    def gap(t,b):
        nonlocal root_evals
        root_evals+=1
        return raw_gap(t)-b
    roots=[];flagged=False
    if method!='gk':
        for b in (0.,A):
            f0,f1=gap(0.,b),gap(math.pi,b)
            if f0*f1<0:roots.append(brentq(lambda t:gap(t,b),0.,math.pi,xtol=1e-13,rtol=1e-14))
        roots.sort();flagged=any(min(t,math.pi-t)<math.pi/16 for t in roots)
    points=roots if method=='all' or (method=='guarded' and flagged) else None
    def traction(t):
        nonlocal law_evals
        law_evals+=1;g=raw_gap(t)
        ds=math.hypot(RX*math.cos(t),-RY*math.sin(t)-2*Q*math.sin(2*t))
        v=K*g if g<=0 else 0. if g>=A else 4*W*g/A**2*(1-(g/A)**2)
        return v*ds
    result=quad(traction,0.,math.pi,epsabs=TAU/(2*L),epsrel=0.,points=points,limit=100,full_output=1)
    elapsed=time.perf_counter_ns()-start
    value,error,info=result[:3]
    return dict(force=-2*L*value,estimate=2*L*error,law_evals=law_evals,root_geometry_evals=root_evals,total_geometry_evals=root_evals+law_evals,ns=elapsed,flagged=flagged,root_count=len(roots),warning=len(result)>3)

def reference(c,dps):
    with mp.workdps(dps):
        # Convert actual binary64 parameters and heights, not rounded display values.
        h=mp.mpf(c['h']);rx=mp.mpf(RX);ry=mp.mpf(RY);q=mp.mpf(Q);a=mp.mpf(A);w=mp.mpf(W);k=mp.mpf(K);ell=mp.mpf(L)
        def gap(t):return h+ry*mp.cos(t)+q*mp.cos(2*t)
        bounds=[mp.mpf(0),mp.pi]
        for b in (mp.mpf(0),a):
            if (gap(0)-b)*(gap(mp.pi)-b)<0:
                lo=mp.mpf(0);hi=mp.pi
                for _ in range(250):
                    mid=(lo+hi)/2
                    if gap(mid)>b:lo=mid
                    else:hi=mid
                bounds.append((lo+hi)/2)
        bounds.sort()
        def f(t):
            g=gap(t);ds=mp.sqrt((rx*mp.cos(t))**2+(-ry*mp.sin(t)-2*q*mp.sin(2*t))**2)
            p=k*g if g<=0 else mp.mpf(0) if g>=a else 4*w*g/a**2*(1-(g/a)**2)
            return p*ds
        return -2*ell*sum(mp.quad(f,[a0,b0]) for a0,b0 in zip(bounds,bounds[1:]))

def main():
    cc=cases(); refs={c['id']:reference(c,45) for c in cc}
    selected=[c for c in cc if c['id'] in ('b0.0_pi_e0.055','b0.1_zero_e0.13')]
    checks=[{'id':c['id'],'abs_difference':float(abs(reference(c,65)-refs[c['id']]))} for c in selected]
    assert len(checks)==2
    records=[]
    for c in cc:
        for m in ('gk','guarded','all'):
            calculate(c,m)
            trials=[calculate(c,m) for _ in range(3)]
            r=trials[0];err=abs(r['force']-float(refs[c['id']]))
            records.append({**c,'method':m,**{k:v for k,v in r.items() if k!='ns'},'median_ns':statistics.median(t['ns'] for t in trials),'reference':float(refs[c['id']]),'abs_error':err,'failed':err>TAU,'false_safe':err>TAU and r['estimate']<=TAU and not r['warning']})
    with (P/'51_nao_circular_resultados.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    summary={'n_cases':len(cc),'tolerance_pN':TAU,'reference_checks':checks,'methods':{}}
    for method in ('gk','guarded','all'):
        r=[x for x in records if x['method']==method]
        summary['methods'][method]={'n':len(r),'failed':sum(x['failed'] for x in r),'false_safe':sum(x['false_safe'] for x in r),'warnings':sum(x['warning'] for x in r),'flagged':sum(x['flagged'] for x in r),'max_error':max(x['abs_error'] for x in r),'law_evals':sum(x['law_evals'] for x in r),'root_geometry_evals':sum(x['root_geometry_evals'] for x in r),'total_geometry_evals':sum(x['total_geometry_evals'] for x in r),'sum_median_ms':sum(x['median_ns'] for x in r)/1e6}
    (P/'51_nao_circular_resumo.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    for r in records:
        if r['failed']:print('FAILED',r['id'],r['method'],r['abs_error'],r['estimate'],r['warning'])

if __name__=='__main__':main()
