"""Exploratory diagnosis, two observed failures; no changes to original benchmark."""
import csv,math,importlib.util,json
from pathlib import Path
from scipy.integrate import quad
p=Path(__file__).resolve().parent;root=p.parent
spec=importlib.util.spec_from_file_location('mrc',root/'45_indicador_gk.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
records=list(csv.DictReader((root/'45_indicador_gk_resultados.csv').open()))
failed=[r for r in records if r['method']=='gk' and r['failed']=='True'];cases={c['id']:c for c in m.make_cases()}
out=[]
for row in failed:
 c=cases[row['id']];h=c['h'];r=c['r'];fac=2*r*m.L
 near=math.acos((m.A-h)/r);interior=math.acos(-h/r);assert near<interior
 ref=float(m.ref(c,True))
 for name,points in [('near_adhesive_only',[near]),('interior_compression_only',[interior]),('both',[near,interior])]:
  calls=[0]
  def f(t):
   calls[0]+=1;g=h+r*math.cos(t)
   return m.K*g if g<=0 else 0. if g>=m.A else 4*m.W*g/m.A**2*(1-(g/m.A)**2)
  v,e,info,*msg=quad(f,0,math.pi,epsabs=m.TAU/fac,epsrel=0,points=points,limit=100,full_output=1)
  force=-fac*v;out.append({'id':c['id'],'strategy':name,'points_rad':points,'abs_error_pN':abs(force-ref),'estimated_error_pN':fac*e,'evaluations':calls[0],'subintervals':info['last'],'warning':bool(msg)})
(p/'isolated_breakpoints.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
