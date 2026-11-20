"""Join fulltext design reviews without treating them as extracted inclusions."""
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import pandas as pd

ROOT=Path(__file__).resolve().parent
CACHE=ROOT.parents[2]/'gbif_annual/alan_experiment_sources/fulltexts_20261120'
parts=[]
for i in [1,2,3]:
    p=ROOT/f'search/fulltext_batch{i}.csv'
    d=pd.read_csv(p).fillna('')
    assert d.pmcid.is_unique
    for r in d.itertuples():
        article=ET.parse(CACHE/(r.pmcid+'.xml')).getroot()
        doi=article.find("./front/article-meta/article-id[@pub-id-type='doi']")
        assert doi is not None and doi.text.lower().strip()==r.doi.lower().strip(),r.pmcid
    d['review_file']=p.name
    parts.append(d)
d=pd.concat(parts,ignore_index=True)
assert d.pmcid.is_unique
d['effect_extraction_complete_by_this_review']=False
d.to_csv(ROOT/'search/fulltext_design_reviews.csv',index=False)
summary=dict(reviewed_reports=len(d),decisions=d.decision.value_counts().to_dict(),
 quantitative_inclusions_established=0,
 interpretation='Methods/Results design review only; no effect extraction by this review. Additional earlier fulltext audits remain in incremental table and source notes; this is not final PRISMA.',
 reviewer='AI-assisted source-linked review, not human dual screening')
(ROOT/'search/fulltext_design_review_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
