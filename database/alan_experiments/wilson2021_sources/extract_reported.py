"""Preserve reported coefficients and CI levels, not raw light/dark effects."""
from pathlib import Path
import re,csv,json,hashlib

OUT=Path(__file__).resolve().parent
CACHE=OUT.parents[3]/'gbif_annual/alan_experiment_sources/wilson2021'
text=(CACHE/'primary.txt').read_text(encoding='utf-8')
results=text.split('Results\nAdult moth abundance',1)[1].split('Discussion',1)[0]
values=re.findall(r'β\s*=\s*([-\d.]+);\s*(85|95)% CI:\s*([-\d.]+),\s*([-\d.]+)',results)
assert len(values)==7,len(values)
labels=[('arthropod_trap_count','lighting_treatment','experimental_crossover'),
        ('yucca_moth_trap_count','lighting_treatment','experimental_crossover'),
        ('mature_fruit_count','direct_illuminance','gradient_within_manipulation'),
        ('mature_fruit_count','skyglow','observational'),
        ('larvae_per_retained_fruit','aspect','environmental_covariate'),
        ('larvae_per_retained_fruit','yucca_density','environmental_covariate'),
        ('larvae_per_retained_fruit','skyglow','observational')]
rows=[]
for (endpoint,predictor,design),(b,level,lo,hi) in zip(labels,values):
    rows.append(dict(endpoint=endpoint,predictor=predictor,design=design,beta=float(b),
                     ci_level=int(level)/100,ci_lower=float(lo),ci_upper=float(hi),
                     source_location='Results pp13-14 author-hosted PDF',
                     scale='author Poisson GLMM; fixed predictors centered/scaled per Methods',
                     status='reported_coefficient_only_scaling_raw_data_and_covariance_hold'))
with (OUT/'reported_coefficients.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
audit=dict(doi='10.1093/icb/icab010',source_url='https://www.brettseymoure.com/s/2021_Wilson-et-al-ICB.pdf',
           pdf_sha256=hashlib.sha256((CACHE/'primary.pdf').read_bytes()).hexdigest(),
           text_sha256=hashlib.sha256((CACHE/'primary.txt').read_bytes()).hexdigest(),
           author_hosted_version='advance article/manuscript pagination; final-version equivalence not checked',
           coefficient_rows=7,meta_effects_added=0,supplement_status='Publisher page403; declared online raw supplement not recovered',
           reproductive_support='5sites;15transects:8lit7dark;58plants;42fruiting plants191fruits for larvae',
           trap_support='Separate sixth site;4towers switching light treatment over7nights')
(OUT/'source_audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print('Seven reported coefficients preserved with original CI levels; no lnRR or SE inferred')
