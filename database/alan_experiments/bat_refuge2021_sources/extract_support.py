"""Reconcile published bat positions with independent lighting support; no meta effect."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('source_dir', type=Path)
args = parser.parse_args()
out = Path(__file__).resolve().parent
record = json.loads((args.source_dir / 'record.json').read_text(encoding='utf-8'))
manifest = []
for name in ['03_Analysis.R', 'TrajFinalFiltered_20190531correctedlightheight.csv']:
    data = (args.source_dir / name).read_bytes()
    metadata = next(f for f in record['files'] if f['key'] == name)
    assert 'md5:' + hashlib.md5(data).hexdigest() == metadata['checksum']
    assert len(data) == metadata['size']
    manifest.append(dict(file=name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
                         checksum=metadata['checksum'], url=metadata['links']['self']))
raw = pd.read_csv(args.source_dir / manifest[1]['file'], encoding='cp1252', low_memory=False)
# Matches author script lines 26-51: retain identifiers represented by >1 position.
sizes = raw.groupby('IdTrajNewSpecies')['ID'].transform('size')
d = raw.loc[sizes > 1].copy()
assert len(raw) == 35067 and len(d) == 28646
assert d.Imp.notna().all() and (d.Imp <= 1).all()
assert d.EdgeSide.isin(['Forest', 'Open']).all()
assert d.groupby('IdTrajNewSpecies')[['Site', 'date', 'spetrum', 'SpeciesGroup']].nunique().max().max() == 1
d['inside'] = d.EdgeSide.eq('Forest')
night = d.groupby(['Site', 'date', 'ID', 'spetrum'], as_index=False).agg(
    positions=('inside', 'size'), trajectories=('IdTrajNewSpecies', 'nunique'))
night.to_csv(out / 'sampled_nights.csv', index=False)
support = night.groupby(['Site', 'spetrum'], as_index=False).agg(
    nights=('date', 'nunique'), positions=('positions', 'sum'), trajectories=('trajectories', 'sum'))
support.to_csv(out / 'site_treatment_support.csv', index=False)
group = d.groupby(['Site', 'date', 'ID', 'spetrum', 'SpeciesGroup'], as_index=False).agg(
    positions=('inside', 'size'), trajectories=('IdTrajNewSpecies', 'nunique'),
    forest_positions=('inside', 'sum'), mean_distance_m=('distLight', 'mean'))
group['forest_position_fraction'] = group.forest_positions / group.positions
group.to_csv(out / 'night_group_descriptive.csv', index=False)
audit = dict(article_doi='10.1186/s40462-020-00238-2', dataset_doi='10.5281/zenodo.4036279',
             source_files=manifest, raw_positions=len(raw), retained_positions=len(d),
             excluded_single_positions=int((sizes == 1).sum()), sites=int(d.Site.nunique()),
             nights=len(night), site_treatment_cells=len(support),
             sites_per_treatment=support.groupby('spetrum').Site.nunique().to_dict(),
             nights_per_treatment=night.groupby('spetrum').date.nunique().to_dict(),
             positions_by_group=d.groupby('SpeciesGroup').size().to_dict(),
             trajectories=int(d.IdTrajNewSpecies.nunique()),
             license=record['metadata'].get('license'),
             new_meta_effects=0, observation='flight location, not attacks or unique animals')
(out / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n', encoding='utf-8')
print(json.dumps(audit, indent=2))
