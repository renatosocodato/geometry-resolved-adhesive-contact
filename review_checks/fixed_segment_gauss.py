"""Exploratory comparator after review: fixed Gauss-Legendre on known branches."""
import csv,math,importlib.util,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
p=Path(__file__).resolve().parent;root=p.parent
def load(name):
 spec=importlib.util.spec_from_file_location(name,root/name);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
c=load('45_indicador_gk.py');n=load('51_teste_nao_circular.py')
oldc={r['id']:r for r in csv.DictReader((root/'45_indicador_gk_resultados.csv').open()) if r['method']=='gk'}
oldn={r['id']:r for r in csv.DictReader((root/'51_nao_circular_resultados.csv').open()) if r['method']=='gk'}
results=[]
for geometry,items,old in [('circular',c.make_cases(),oldc),('noncircular',n.cases(),oldn)]:
 for case in items:
  h=case['h'];root_calls=[0]
  if geometry=='circular':
   r=case['r']; roots=c.roots(h,r)
   def fun(t):
    g=h+r*math.cos(t);return c.K*g if g<=0 else 0. if g>=c.A else 4*c.W*g/c.A**2*(1-(g/c.A)**2)
   factor=-2*r*c.L;ref=float(old[case['id']]['true_reference'])
  else:
   def gap(t,b):root_calls[0]+=1;return h+n.RY*math.cos(t)+n.Q*math.cos(2*t)-b
   roots=[]
   for b in (0.,n.A):
    f0,f1=gap(0.,b),gap(math.pi,b)
    if f0*f1<0:roots.append(brentq(lambda t:gap(t,b),0,math.pi,xtol=1e-13,rtol=1e-14))
   roots.sort()
   def fun(t):
    g=h+n.RY*math.cos(t)+n.Q*math.cos(2*t);ds=math.hypot(n.RX*math.cos(t),-n.RY*math.sin(t)-2*n.Q*math.sin(2*t));return (n.K*g if g<=0 else 0. if g>=n.A else 4*n.W*g/n.A**2*(1-(g/n.A)**2))*ds
   factor=-2*n.L;ref=float(old[case['id']]['reference'])
  bounds=[0.]+roots+[math.pi]
  for order in (4,8):
   x,w=np.polynomial.legendre.leggauss(order);total=0.;calls=0
   for lo,hi in zip(bounds,bounds[1:]):
    for z,weight in zip(x,w):total+=weight*(hi-lo)/2*fun((hi+lo)/2+z*(hi-lo)/2);calls+=1
   error=abs(factor*total-ref)
   results.append({'geometry':geometry,'id':case['id'],'order_per_segment':order,'abs_error_pN':error,'force_evaluations':calls,'root_gap_evaluations':root_calls[0]})
with (p/'fixed_segment_gauss.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
summary={}
for g in ['circular','noncircular']:
 for order in (4,8):
  a=[x for x in results if x['geometry']==g and x['order_per_segment']==order]
  summary[f'{g}_n{order}']={'cases':len(a),'max_error_pN':max(x['abs_error_pN'] for x in a),'tolerance_exceedances':int(sum(x['abs_error_pN']>1e-8 for x in a)),'force_evaluations':sum(x['force_evaluations'] for x in a),'root_gap_evaluations':sum(x['root_gap_evaluations'] for x in a)}
(p/'fixed_segment_gauss_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
