#!/usr/bin/env python3
"""Fixed geometric indicator versus adaptive GK, protocol 44."""
from pathlib import Path
import csv, json, math, statistics, time
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
RADII=['0.28209479177387814','0.36']; WIDTHS=['0.01','0.03','0.08','0.15','0.23','0.38']
L=.75;K=20.;W=.002;A=.1;TAU=1e-8;THRESHOLD=math.pi/16

def make_cases():
    out=[]
    for rs in RADII:
        r=float(rs)
        for b in (0.,.1):
            for side in ('zero','pi'):
                for es in WIDTHS:
                    e=float(es);h=b+(-1 if side=='zero' else 1)*r*math.cos(e)
                    out.append(dict(id=f'r{rs}_b{b}_{side}_e{es}',r=r,h=h,b=b,side=side,e=e))
        for h in [.05,r+.12]:
            out.append(dict(id=f'r{rs}_control_{h}',r=r,h=h,b=None,side='control',e=None))
    return out

def roots(h,r):
    return sorted(math.acos(c) for b in [0.,A] if -1<(c:=(b-h)/r)<1)

def run(case,method):
    start=time.perf_counter_ns();r=case['r'];h=case['h'];fac=2*r*L
    crossings=[]; flagged=False
    if method!='gk':
        crossings=roots(h,r)
        flagged=any(min(t,math.pi-t)<THRESHOLD for t in crossings)
    points=(crossings if method=='all' or (method=='guarded' and flagged) else None)
    count=0
    def integrand(t):
        nonlocal count
        count+=1;g=h+r*math.cos(t)
        if g<=0:return K*g
        if g>=A:return 0.
        return 4*W*g/A**2*(1-(g/A)**2)
    result=quad(integrand,0.,math.pi,epsabs=TAU/fac,epsrel=0.,points=points,limit=100,full_output=1)
    elapsed=time.perf_counter_ns()-start
    val,err,info=result[:3]
    return {'force':-fac*val,'estimate':fac*err,'evaluations':count,'ns':elapsed,'flagged':flagged,'root_count':len(crossings),'segments_used':len(points) if points else 0,'warning':len(result)>3}

def ref(case, direct=False):
    with mp.workdps(60):
        # The benchmark parameters are their actual binary64 input values.
        h=mp.mpf(case['h']);r=mp.mpf(case['r']);a=mp.mpf(A);w=mp.mpf(W);k=mp.mpf(K);ell=mp.mpf(L)
        cuts=sorted([mp.mpf(0),mp.pi]+[mp.acos((b-h)/r) for b in (mp.mpf(0),a) if -1<(b-h)/r<1])
        def deriv(t):
            g=h+r*mp.cos(t)
            return k*g if g<=0 else 4*w*g/a**2*(1-(g/a)**2) if g<a else mp.mpf(0)
        if direct:
            value=sum(mp.quad(deriv,[l,u]) for l,u in zip(cuts,cuts[1:]))
        else:
            def primitive(t):
                s=mp.sin(t); p1=h*t+r*s
                p3=h**3*t+3*h*h*r*s+3*h*r*r*(t/2+mp.sin(2*t)/4)+r**3*(s-s**3/3)
                return p1,p3
            value=mp.mpf(0)
            for lo,hi in zip(cuts,cuts[1:]):
                g=h+r*mp.cos((lo+hi)/2); l1,l3=primitive(lo);u1,u3=primitive(hi)
                if g<=0:value+=k*(u1-l1)
                elif g<a:value+=4*w/a**2*((u1-l1)-(u3-l3)/a**2)
        return -2*r*ell*value

def main():
    cases=make_cases();records=[];refs={}
    for c in cases:
        refs[c['id']]=ref(c)
    # Two preselected direct, high-precision integrals are independent of the primitive.
    selected=[c for c in cases if (c['r']==float(RADII[0]) and c['b']==0 and c['side']=='pi' and c['e']==.03) or (c['r']==float(RADII[1]) and c['b']==.1 and c['side']=='zero' and c['e']==.15)]
    validation=[{'id':c['id'],'direct_vs_primitive_abs':float(abs(ref(c,True)-refs[c['id']]))} for c in selected]
    for c in cases:
        for method in ['gk','guarded','all']:
            run(c,method)
            trials=[run(c,method) for _ in range(3)]
            r=trials[0];err=abs(r['force']-float(refs[c['id']]))
            records.append({**c,'method':method,**{k:v for k,v in r.items() if k!='ns'},'median_ns':statistics.median(t['ns'] for t in trials),'true_reference':float(refs[c['id']]),'abs_error':err,'failed':err>TAU,'false_safe':err>TAU and r['estimate']<=TAU and not r['warning']})
    with (HERE/'45_indicador_gk_resultados.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    summary={'tolerance_pN':TAU,'n_cases':len(cases),'python_numpy':np.__version__,'scipy':scipy.__version__,'reference_checks':validation,'methods':{}}
    for method in ['gk','guarded','all']:
        rows=[r for r in records if r['method']==method]
        summary['methods'][method]={'failed':sum(r['failed'] for r in rows),'false_safe':sum(r['false_safe'] for r in rows),'flagged':sum(r['flagged'] for r in rows),'warnings':sum(r['warning'] for r in rows),'max_error':max(r['abs_error'] for r in rows),'total_law_evaluations':sum(r['evaluations'] for r in rows),'sum_median_time_ms':sum(r['median_ns'] for r in rows)/1e6}
    (HERE/'45_indicador_gk_resumo.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
    for row in records:
        if row['failed']:print('FAILED',row['method'],row['id'],'error',row['abs_error'],'estimate',row['estimate'],'flagged',row['flagged'])

if __name__=='__main__':main()
