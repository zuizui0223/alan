"""Independent source verification using Poppler line layout, no repository extractor."""
from pathlib import Path
from decimal import Decimal
import csv, hashlib, json, subprocess, datetime, argparse

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pdf', required=True, type=Path)
parser.add_argument('--selected', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
args=parser.parse_args()
PDF=args.pdf; SELECTED=args.selected; OUT=args.out
OUT.mkdir(parents=True, exist_ok=True)
EXPECTED='baa5a7e2a3e0aaa56d1d810c5bbe6e9c93414f488e4fb12dcd8f7be03d96e2b0'
assert hashlib.sha256(PDF.read_bytes()).hexdigest()==EXPECTED
text=subprocess.check_output(['pdftotext','-layout',str(PDF),'-'],text=True)
(OUT/'anic2022_poppler_layout.txt').write_text(text)
names={'Agrostis tenuis','Anthoxanthum odoratum','Hypochaeris radicata','Lotus pedunculatus'}
parsed=[]; malformed=[]
for page_number,page in enumerate(text.split('\f')[:14],start=1):
    for line_number,line in enumerate(page.splitlines(),start=1):
        tokens=line.split()
        if ' '.join(tokens[:2]) not in names:
            continue
        if len(tokens)!=8:
            malformed.append({'page':page_number,'line':line_number,'text':line})
            continue
        genus,epithet,mesocosm,light,trophic,total,root,shoot=tokens
        if light not in ['L1','L2','L3'] or trophic not in ['T1','T2','T3']:
            malformed.append({'page':page_number,'line':line_number,'text':line})
            continue
        for value in [total,root,shoot]: assert Decimal(value)>0
        parsed.append(dict(species=f'{genus} {epithet}',mesocosm=mesocosm,light=light,trophic=trophic,total_source=total,root_source=root,shoot_source=shoot,source_page=page_number))
assert not malformed
assert len(parsed)==216
assert len({(r['species'],r['mesocosm']) for r in parsed})==216
expected_rows=list(csv.DictReader(SELECTED.open()))
assert len(expected_rows)==216
reference={(r['species'],r['mesocosm']):r for r in expected_rows}
actual={(r['species'],r['mesocosm']):r for r in parsed}
missing_in_selected=sorted(set(actual)-set(reference))
missing_in_pdf=sorted(set(reference)-set(actual))
differences=[]; lexical_differences=[]
for key in sorted(set(actual)&set(reference)):
    source=actual[key];target=reference[key]
    for col in ['light','trophic']:
        if source[col]!=target[col]: differences.append({'key':key,'column':col,'source':source[col],'selected':target[col]})
    if source['source_page']!=int(target['source_page']): differences.append({'key':key,'column':'source_page','source':source['source_page'],'selected':target['source_page']})
    for organ in ['total','root','shoot']:
        original=Decimal(source[f'{organ}_source'])
        for column in [f'{organ}_source',f'{organ}_g']:
            if original!=Decimal(target[column]): differences.append({'key':key,'column':column,'source':str(original),'selected':target[column]})
        if source[f'{organ}_source']!=target[f'{organ}_source']:
            lexical_differences.append({'key':key,'column':f'{organ}_source','source':source[f'{organ}_source'],'selected':target[f'{organ}_source']})
assert not missing_in_selected and not missing_in_pdf and not differences
by_species={name:{'rows':sum(r['species']==name for r in parsed),'source_pages':sorted({r['source_page'] for r in parsed if r['species']==name})} for name in sorted(names)}
with (OUT/'independent_selected_source_rows_poppler.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(parsed[0]));writer.writeheader();writer.writerows(parsed)
result={
    'verification_date_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'source_pdf':str(PDF),
    'source_pdf_sha256':EXPECTED,
    'selected_csv':str(SELECTED),
    'selected_csv_sha256':hashlib.sha256(SELECTED.read_bytes()).hexdigest(),
    'independent_method':'Poppler pdftotext -layout; page-local line tokenization; Decimal comparison of six numeric source/value fields; no import or invocation of repository extractor',
    'pdf_pages_used':list(range(1,15)),
    'parsed_selected_species_rows':len(parsed),
    'unique_species_mesocosm_keys':len(actual),
    'species':by_species,
    'light_and_trophic_values_compared':len(actual)*2,
    'biomass_fields_compared':len(actual)*6,
    'source_pages_compared':len(actual),
    'missing_in_selected':missing_in_selected,
    'missing_in_pdf':missing_in_pdf,
    'numeric_treatment_or_page_differences':differences,
    'numeric_treatment_or_page_differences_count':len(differences),
    'format_only_source_lexeme_differences_count':len(lexical_differences),
    'format_only_source_lexeme_examples':lexical_differences[:6],
    'malformed_selected_species_source_lines':malformed,
    'visual_review':{'status':'pending','rendered_pages':[2,3,4,8,9,10,11]},
    'human_duplicate_extraction_complete':False,
    'independent_AI_source_check':True,
    'no_missing_values_imputed':True,
    'no_original_source_values_modified':True,
}
(OUT/'selected216_original_pdf_independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
