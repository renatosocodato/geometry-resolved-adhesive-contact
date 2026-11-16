#!/usr/bin/env python3
"""Frozen wall-force diagnostic; no end caps, see protocol 41."""
import csv,math,json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
import mpmath as mp

HERE=Path(__file__).resolve().parent
R=0.28209479177387814; L=.75; K=20.; W=.002; A=.1
HEIGHTS=[.235,.255,.275,.295,.315,.335,.365,.385,.405]
BETAS=[0.,.08]; ORDERS=[8,16,32,64];TAU=1.11347e-4

def angular(h,r):
    breaks=[0.,math.pi]
    for b in [0.,A]:
        v=(b-h)/r
        if -1.<v<1.:breaks.append(math.acos(v))
    breaks.sort()
    def primitives(t):
        s=math.sin(t)
        i1=h*t+r*s
        i3=h**3*t+3*h*h*r*s+3*h*r*r*(t/2+math.sin(2*t)/4)+r**3*(s-s**3/3)
        return i1,i3
    ans=0.
    for l,u in zip(breaks,breaks[1:]):
        g=h+r*math.cos((l+u)/2); a1,a3=primitives(l);b1,b3=primitives(u)
        if g<=0:ans+=K*(b1-a1)
        elif g<A:ans+=4*W/A**2*((b1-a1)-(b3-a3)/A**2)
    return ans

def geometry(h,beta):
    length=math.hypot(L,beta);rb=R*L/length; factor=-2*R*length
    cuts=[0.,1.]
    if beta:
        for b in [0.,A]:
            for sign in [-1,1]:
                u=.5+(b-sign*rb-h)/beta
                if 0.<u<1.:cuts.append(u)
    return rb,factor,sorted(cuts)

def force(h,beta,n,split):
    rb,factor,cuts=geometry(h,beta)
    if not split:cuts=[0.,1.]
    q,rem=divmod(n,len(cuts)-1);ans=0.
    for i,(l,u) in enumerate(zip(cuts,cuts[1:])):
        x,w=leggauss(q+(i<rem))
        for node,weight in zip(x,w):
            s=(l+u)/2+node*(u-l)/2
            ans+=weight*(u-l)/2*angular(h+beta*(s-.5),rb)
    return factor*ans

def reference(h,beta):
    rb,factor,cuts=geometry(h,beta)
    val=err=0.
    for l,u in zip(cuts,cuts[1:]):
        v,e=quad(lambda s:angular(h+beta*(s-.5),rb),l,u,epsabs=1e-12,epsrel=1e-12,limit=100)
        val+=v;err+=e
    return factor*val,abs(factor)*err

def mp_reference(hstr,betastr):
    with mp.workdps(40):
        r=mp.mpf(str(R));l=mp.mpf('.75');beta=mp.mpf(betastr);h=mp.mpf(hstr);a=mp.mpf('.1');w=mp.mpf('.002');k=mp.mpf(20)
        length=mp.sqrt(l*l+beta*beta);rb=r*l/length;factor=-2*r*length
        def angular_mp(z):
            bounds=sorted([mp.mpf(0),mp.pi]+[mp.acos((b-z)/rb) for b in (mp.mpf(0),a) if -1<(b-z)/rb<1])
            def integral(x):
                st=mp.sin(x);i1=z*x+rb*st;i3=z**3*x+3*z*z*rb*st+3*z*rb**2*(x/2+mp.sin(2*x)/4)+rb**3*(st-st**3/3)
                return i1,i3
            total=mp.mpf(0)
            for left,right in zip(bounds,bounds[1:]):
                g=z+rb*mp.cos((left+right)/2);x1,x3=integral(left);y1,y3=integral(right)
                if g<=0:total+=k*(y1-x1)
                elif g<a:total+=4*w/a**2*((y1-x1)-(y3-x3)/a**2)
            return total
        cuts=[mp.mpf(0),mp.mpf(1)]
        if beta:
            for b in (mp.mpf(0),a):
                for sign in (-1,1):
                    u=mp.mpf('.5')+(b-sign*rb-h)/beta
                    if 0<u<1:cuts.append(u)
        cuts.sort()
        return factor*sum(mp.quad(lambda x: angular_mp(h+beta*(x-mp.mpf('.5'))),[left,right]) for left,right in zip(cuts,cuts[1:]))

def main():
    rows=[]
    for beta in BETAS:
        for h in HEIGHTS:
            ref,re=reference(h,beta); rb,factor,cuts=geometry(h,beta)
            for n in ORDERS:
                for method in ['global','split']:
                    value=force(h,beta,n,method=='split')
                    rows.append(dict(beta=beta,h=h,n=n,method=method,force=value,reference=ref,reference_scipy_error_estimate=re,force_error=abs(value-ref),interior_cuts=len(cuts)-2))
    with (HERE/'42_geometria_resultados.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    lookup={(r['beta'],r['h'],r['n'],r['method']):r for r in rows};checks=[]
    for beta in BETAS:
        for h in HEIGHTS:
            for method in ['global','split']:
                for n in [8,16,32]:
                    a=lookup[beta,h,n,method];b=lookup[beta,h,2*n,method]
                    d=abs(a['force']-b['force']);e=b['force_error']
                    checks.append(dict(beta=beta,h=h,method=method,n=n,difference=d,fine_error=e,false_safe=(d<=TAU and e>TAU),underestimate=(d<e and e>1e-12)))
    with (HERE/'42_geometria_indicador.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=checks[0]);w.writeheader();w.writerows(checks)
    summary={}
    for beta in BETAS:
        for method in ['global','split']:
            data=[x for x in checks if x['beta']==beta and x['method']==method]
            raw=[x for x in rows if x['beta']==beta and x['method']==method]
            summary[f'beta{beta}_{method}']={'false_safe':sum(x['false_safe'] for x in data),'underestimate':sum(x['underestimate'] for x in data),'cases':len(data),'maximum_force_error':max(x['force_error'] for x in raw)}
    summary['reference_checks']={}
    for beta in BETAS:
        d=mp_reference('0.275',str(beta)); scipy_ref,e=reference(.275,beta)
        summary['reference_checks'][str(beta)]={'scipy':scipy_ref,'mp_40':str(d),'difference':abs(scipy_ref-float(d)),'scipy_error_estimate':e}
    (HERE/'42_geometria_resumo.json').write_text(json.dumps(summary,indent=2,default=lambda x:x.item())+'\n')
    print(json.dumps(summary,indent=2,default=lambda x:x.item()))
    for x in checks:
        if x['false_safe']: print('FALSE SAFE',x)

if __name__=='__main__':main()
