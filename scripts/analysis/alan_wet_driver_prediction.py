"""Predict candidate membership in regions excluded from fitting."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score

repo=Path(__file__).resolve().parents[2]
base=repo/"results/all_species_forest";out=base/"wet_drivers"
traits=pd.read_csv(repo.parent/"gbif_annual/pladias_driver_accounts/compiled.csv.gz",index_col="sp")
w=pd.read_csv(base/"wet_reproduction/species.csv").set_index("sp").join(pd.read_csv(base/"decline_traits/moisture_evidence.csv.gz").set_index("sp")[["site_region","site_band","island_band"]])
w=w[traits.notna().any(axis=1) & ~w.accepted.duplicated(keep=False)]
traits=traits.loc[w.index]
# Availability is included in both models so documenting more traits adds no inherent advantage.
context=pd.concat([w[["indicator","site_band","island_band"]],traits.isna().add_prefix("unknown_").astype(int)],axis=1)
p=pd.DataFrame(index=w.index)
for reg in sorted(w.site_region.unique()):
    test=w.site_region.eq(reg);train=~test
    for name,features in [("baseline",context),("plus_traits",pd.concat([context,traits],axis=1))]:
        model=make_pipeline(SimpleImputer(strategy="median",keep_empty_features=True),StandardScaler(),LogisticRegression(C=.05,max_iter=1000))
        model.fit(features[train],w.candidate[train])
        p.loc[test,name]=model.predict_proba(features[test])[:,1]
p=p.join(w[["site_region","candidate"]])
assert p.notna().all().all()
def scores(df,col):
    y=df.candidate.astype(int);pred=df[col]
    auc=[];pairs=[]
    for _,s in df.groupby("site_region"):
        a=s.candidate.sum();b=len(s)-a
        if a*b:
            auc.append(roc_auc_score(s.candidate,s[col]));pairs.append(a*b)
    return [np.mean(-y*np.log(pred)-(1-y)*np.log1p(-pred)),np.average(auc,weights=pairs)]
rows=[]
for region in ["all",*sorted(p.site_region.unique())]:
    f=p if region=="all" else p[p.site_region!=region]
    for col in ["baseline","plus_traits"]:
        loss,auc=scores(f,col)
        rows.append(dict(omitted=region,model=col,n=len(f),cases=f.candidate.sum(),regions=f.site_region.nunique(),logloss=loss,within_region_auc=auc))
pd.DataFrame(rows).to_csv(out/"regional_prediction.csv",index=False)
p.to_csv(out/"held_out_predictions.csv")
print(pd.DataFrame(rows).head(2))
