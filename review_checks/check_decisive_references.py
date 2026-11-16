"""Re-evaluate only the two observed tolerance exceedances with direct mp quadrature."""
from pathlib import Path
import importlib.util,csv,json,mpmath as mp
p=Path(__file__).resolve().parent; root=p.parent
spec=importlib.util.spec_from_file_location('circular',root/'45_indicador_gk.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rows=list(csv.DictReader((root/'45_indicador_gk_resultados.csv').open()))
failed=[r for r in rows if r['method']=='gk' and r['failed']=='True']
cases={c['id']:c for c in m.make_cases()};assert len(failed)==2
out=[]
for row in failed:
 c=cases[row['id']]
 # mp workdps is internal to ref; evaluate a second independent direct implementation at 80 digits.
 with mp.workdps(80):
  h,r,a,k,w,L=map(mp.mpf,[c['h'],c['r'],m.A,m.K,m.W,m.L])
  cuts=sorted([mp.mpf('0'),mp.pi]+[mp.acos((b-h)/r) for b in [mp.mpf('0'),a] if abs((b-h)/r)<1])
  def tr(t):
   g=h+r*mp.cos(t)
   return k*g if g<=0 else 4*w*g/a**2*(1-(g/a)**2) if g<a else mp.mpf('0')
  direct=-2*r*L*sum(mp.quad(tr,[x,y]) for x,y in zip(cuts,cuts[1:]))
  primitive=m.ref(c); direct60=m.ref(c,True); numerical=mp.mpf(float(row['force']))
  out.append({'id':row['id'],'force_float64':row['force'],'direct_80_digit':mp.nstr(direct,70),'direct80_vs_primitive60_abs_pN':mp.nstr(abs(direct-primitive),12),'direct80_vs_direct60_abs_pN':mp.nstr(abs(direct-direct60),12),'error_against_direct80_pN':mp.nstr(abs(numerical-direct),12),'saved_binary64_error_pN':row['abs_error'],'saved_error_estimate_pN':row['estimate'],'transition_angles_rad':[float(x) for x in cuts[1:-1]]})
(p/'decisive_reference_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
