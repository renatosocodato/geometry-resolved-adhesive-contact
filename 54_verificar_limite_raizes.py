#!/usr/bin/env python3
"""Two circular-limit checks and numerical versus algebraic roots, protocol 53."""
import csv,importlib.util,json,math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

p=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('noncircular',p/'51_teste_nao_circular.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
root_differences=[];mismatches=[];root_count=0
for c in m.cases():
 for b in (0.,m.A):
  f=lambda t:c['h']+m.RY*math.cos(t)+m.Q*math.cos(2*t)-b
  numeric=[brentq(f,0,math.pi,xtol=1e-13,rtol=1e-14)] if f(0)*f(math.pi)<0 else []
  polys=np.roots([2*m.Q,m.RY,c['h']-m.Q-b])
  analytic=sorted(math.acos(float(z.real)) for z in polys if abs(z.imag)<1e-12 and -1<float(z.real)<1)
  if len(numeric)!=len(analytic):mismatches.append((c['id'],b,len(numeric),len(analytic)))
  else:
   root_count+=len(numeric)
   root_differences.extend(abs(x-y) for x,y in zip(numeric,analytic))
prior=list(csv.DictReader((p/'45_indicador_gk_resultados.csv').open()))
circular=[r for r in prior if r['method']=='gk' and r['failed']=='True']
assert len(circular)==2
limits=[]
for r in circular:
 radius=float(r['r']);m.RX=radius;m.RY=radius;m.Q=0.
 calc=float(m.reference({'h':float(r['h'])},45));base=float(r['true_reference']);diff=abs(calc-base)
 limits.append({'id':r['id'],'abs_difference_pN':diff})
out={'n_root_comparisons':root_count,'root_count_mismatches':mismatches,'max_root_difference_rad':max(root_differences,default=0.),'circular_limits':limits}
(p/'54_validacao_geometria.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
assert not mismatches and all(x['abs_difference_pN']<1e-12 for x in limits)
