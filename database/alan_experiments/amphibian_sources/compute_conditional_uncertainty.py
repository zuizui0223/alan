"""Conditional group-level uncertainty; see UNCERTAINTY_PLAN.md."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import t
ROOT=Path(__file__).resolve().parent
d=pd.read_csv(ROOT/'descriptive_group_means.csv')
wide=d.pivot(index=['block','light'],columns='endpoint',values='mean')
assert wide.shape==(9,4) and not wide.isna().any().any()
names=wide.columns.tolist()
a=wide.xs('Present',level='light').to_numpy();b=wide.xs('Absent',level='light').to_numpy()
assert len(a)==4 and len(b)==5
ma,mb=a.mean(0),b.mean(0)
assert (ma>0).all() and (mb>0).all()
va=np.cov(a,rowvar=False,ddof=1)/(len(a)*np.outer(ma,ma))
vb=np.cov(b,rowvar=False,ddof=1)/(len(b)*np.outer(mb,mb))
v=va+vb
y=np.log(ma/mb)
assert np.allclose(v,v.T) and np.linalg.eigvalsh(v).min()>-1e-12
rows=[]
for j,name in enumerate(names):
    s1=np.var(a[:,j],ddof=1)/(len(a)*ma[j]**2)
    s2=np.var(b[:,j],ddof=1)/(len(b)*mb[j]**2)
    assert np.isclose(s1+s2,v[j,j])
    df=(s1+s2)**2/(s1*s1/(len(a)-1)+s2*s2/(len(b)-1))
    half=t.ppf(.975,df)*np.sqrt(v[j,j])
    rows.append(dict(endpoint=name,lnRR=y[j],variance=v[j,j],welch_df=df,
                     ratio=np.exp(y[j]),lower=np.exp(y[j]-half),upper=np.exp(y[j]+half),
                     lit_groups=len(a),dark_groups=len(b),meta_adopted=False))
result=pd.DataFrame(rows)
prior=pd.read_csv(ROOT/'descriptive_arm_means.csv').set_index('endpoint')
assert np.allclose(result.set_index('endpoint').ratio,prior.loc[names,'descriptive_ratio'])
result.to_csv(ROOT/'conditional_effects.csv',index=False)
pd.DataFrame(v,index=names,columns=names).rename_axis('endpoint').to_csv(ROOT/'conditional_covariance.csv')
(ROOT/'conditional_audit.json').write_text(json.dumps(dict(clusters=9,endpoints=4,
    covariance_entries=int(v.size),minimum_eigenvalue=float(np.linalg.eigvalsh(v).min()),
    main_catalogue_modified=False,final_eligibility=False,
    interpretation='Conditional approximate delta intervals; independent source-labelled groups assumed'),indent=2)+'\n')
print(result.to_string(index=False))
