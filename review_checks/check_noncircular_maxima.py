"""Post hoc precision convergence at each method's largest saved error."""
import csv,importlib.util,json,mpmath as mp
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parent
spec=importlib.util.spec_from_file_location('m',root/'51_teste_nao_circular.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rows=list(csv.DictReader((root/'51_nao_circular_resultados.csv').open()));cases={x['id']:x for x in m.cases()};result=[]
for method in ('gk','guarded','all'):
 row=max((x for x in rows if x['method']==method),key=lambda x:float(x['abs_error']));c=cases[row['id']]
 a=m.reference(c,45);b=m.reference(c,65)
 with mp.workdps(65):result.append({'method':method,'id':row['id'],'saved_binary64_error_pN':row['abs_error'],'reference_45_vs_65_pN':mp.nstr(abs(a-b),12),'force_error_vs_65_pN':mp.nstr(abs(mp.mpf(float(row['force']))-b),12)})
(p/'noncircular_maxima_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
