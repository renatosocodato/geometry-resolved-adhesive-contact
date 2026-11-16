#!/usr/bin/env python3
"""Independent direct binary64 primitive comparison on frozen benchmark cases."""
import csv,importlib.util,json,statistics,time
from pathlib import Path
p=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('primitive',p/'42_indicador_geometria.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
src=list(csv.DictReader((p/'45_indicador_gk_resultados.csv').open()))
rows=[]
for r in (x for x in src if x['method']=='gk'):
    radius=float(r['r']);h=float(r['h']);factor=-2*radius*.75
    def evaluate():
        start=time.perf_counter_ns();f=factor*m.angular(h,radius);return f,time.perf_counter_ns()-start
    evaluate();data=[evaluate() for _ in range(3)]
    rows.append({'id':r['id'],'force':data[0][0],'reference':float(r['true_reference']),'abs_error':abs(data[0][0]-float(r['true_reference'])),'median_ns':statistics.median(x[1] for x in data)})
with (p/'48_forma_fechada_resultados.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
summary={'n_cases':len(rows),'max_error':max(x['abs_error'] for x in rows),'failed_tau':sum(x['abs_error']>1e-8 for x in rows),'sum_median_time_ms':sum(x['median_ns'] for x in rows)/1e6}
(p/'48_forma_fechada_resumo.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
