"""Post hoc descriptive subset of saved benchmark records; no new simulations."""
import csv,json
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parent
c=list(csv.DictReader((root/'45_indicador_gk_resultados.csv').open()));n=list(csv.DictReader((root/'51_nao_circular_resultados.csv').open()))
result={}
for geometry,rows,inside in [('circular',c,lambda x:float(x['h'])>=float(x['r'])),('noncircular',n,lambda x:float(x['h'])>=.28209479177387814-.045)]:
 result[geometry]={}
 for method in ('gk','guarded','all'):
  sub=[r for r in rows if r['method']==method and inside(r)]
  result[geometry][method]={'n_cases':len(sub),'tolerance_exceedances':sum(r['failed'].lower()=='true' for r in sub),'max_abs_error_pN':max(float(r['abs_error']) for r in sub)}
(p/'nonpenetrating_subset.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
