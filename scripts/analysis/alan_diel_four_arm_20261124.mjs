// Source-arm reanalysis, descriptive study-cluster resampling. Run from ALAN root:
// node scripts/analysis/alan_diel_four_arm_20261124.mjs /path/to/DvN_LiteratureDataClean_2024-09-13.csv
import {readFileSync} from 'node:fs';
import {execFileSync} from 'node:child_process';
function csv(s){let rows=[],r=[],f="",q=false;for(let i=0;i<s.length;i++){let c=s[i];if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}else if(c==='"'){q=true}else if(c===','){r.push(f);f=""}else if(c==='\n'){r.push(f.replace(/\r$/,""));rows.push(r);r=[];f=""}else f+=c}if(f||r.length){r.push(f.replace(/\r$/,""));rows.push(r)}let h=rows.shift();if(rows.some(r=>r.length!==h.length))throw Error("CSV shape");return rows.map(r=>Object.fromEntries(h.map((k,i)=>[k,r[i]])))}
function analyze(raw,wetRaw){
const a=csv(raw),wet=csv(wetRaw);
const uniq=a=>[...new Set(a)];const group=(a,f)=>{let m=new Map;for(let x of a){let k=f(x);if(!m.has(k))m.set(k,[]);m.get(k).push(x)}return m};
const keys=["study_ID","country","locality","Site_period","coordinates_lat","coordinates_lon","year_start","year_end","month_start","month_end","plant_species","plant_cultivar","additional_treatment","treatment_effectiveness_metric"];
const arms=["open pollination","day pollination","night pollination","complete exclusion"];
const g=group(a.filter(x=>["fruit set","seed set"].includes(x.treatment_effectiveness_metric)&&arms.includes(x.treatment_condition)),x=>JSON.stringify(keys.map(k=>x[k])));
const valid=[...g.values()].filter(r=>r.length===4&&new Set(r.map(x=>x.treatment_condition)).size===4);
const excluded={},d=[];
for(const z of valid){let x=z[0],metric=x.treatment_effectiveness_metric,unit=x.effectiveness_units;
let why=null;if(z.some(v=>v.exclude!=="n"||!Number.isFinite(+v.effectiveness_value)||v.effectiveness_value===""||!Number.isFinite(+v.sample_size)||+v.sample_size<=1||+v.effectiveness_value<0))why="missing_or_n_le1";else if(z.some(v=>v.effectiveness_units!==unit))why="mixed_units";else if(metric==="fruit set"&&(!/^(proportion|fruit\/flower)/i.test(unit)||z.some(v=>+v.effectiveness_value>1)))why="fruit_nonproportion";else if(metric==="seed set"&&!(/seeds?\/flower|seeds?\/ovule|seeds? per ovule|proportion.*(seeds\/ovule|seed set)/i.test(unit)))why="seed_other_denominator";
if(why){excluded[why]=(excluded[why]||0)+1;continue}
let [O,D,N,B]=arms.map(k=>+z.find(v=>v.treatment_condition===k).effectiveness_value);
if(O<=0){excluded.open_zero=(excluded.open_zero||0)+1;continue}
d.push({study:x.study_ID,plant:x.plant_species,metric,unit,O,D,N,B,n_min:Math.min(...z.map(v=>+v.sample_size)),night_loss:1-D/O,night_only_net:(N-B)/O,interaction:(O-D-N+B)/O,baseline:B/O,day_loss:1-N/O,night_loss_abs:O-D,night_only_net_abs:N-B,interaction_abs:O-D-N+B,doi:x.DOI,extra:x.additional_treatment,site:x.locality,period:x.Site_period,year:x.year_start,source_key:keys.map(k=>x[k]),source_row_ids:z.map(v=>v[""])});
}
const mean=x=>x.reduce((s,v)=>s+v,0)/x.length;const fields=["night_loss","night_only_net","interaction","baseline","day_loss","night_loss_abs","night_only_net_abs","interaction_abs"];
function agg(arr,key){return [...group(arr,key).values()].map(z=>({study:z[0].study,plant:z[0].plant,n:z.length,...Object.fromEntries(fields.map(f=>[f,mean(z.map(v=>v[f]))]))}))}
function describe(z){let s=agg(z,v=>v.study+"|"+v.plant),byStudy=agg(s,v=>v.study),byPlant=agg(s,v=>v.plant);let state=7125;const rand=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296};const boot={};let ff=z[0].metric==="seed set"?fields.slice(0,5):fields;for(const f of ff){let v=[];for(let k=0;k<9999;k++)v.push(mean(Array.from({length:byStudy.length},()=>byStudy[Math.floor(rand()*byStudy.length)][f])));v.sort((a,b)=>a-b);boot[f]=[v[249],v[9749]]}
return {cells:z.length,studyPlant:s.length,studies:byStudy.length,plants:byPlant.length,mean_study:Object.fromEntries(ff.map(f=>[f,mean(byStudy.map(v=>v[f]))])),CI:boot,mean_species:Object.fromEntries(ff.map(f=>[f,mean(byPlant.map(v=>v[f]))])),study_plant_nonadditive:{negative:s.filter(x=>x.interaction<0).length,positive:s.filter(x=>x.interaction>0).length,zero:s.filter(x=>x.interaction===0).length,night_only_20pct_but_night_loss_le10pct:s.filter(x=>x.night_only_net>=.2&&x.night_loss<=.1).length,night_loss_gt10pct:s.filter(x=>x.night_loss>.1).length,night_loss_lt0:s.filter(x=>x.night_loss<0).length}}}
const fruits=d.filter(x=>x.metric==="fruit set"),seeds=d.filter(x=>x.metric==="seed set");
const wetMatches=d.filter(x=>wet.some(y=>x.plant===y.accepted||x.plant===y.species)).map(x=>({...x,cohort:wet.filter(y=>x.plant===y.accepted||x.plant===y.species).map(y=>({candidate:y.candidate,indicator:y.indicator,accepted:y.accepted}))}));
return {source_blob:"39a82d966afcdde670603add1c881ace8791ffe9",source_ref:"1945b7c114ad87000608ce6018d0280001c04ef1",wet_blob:"388a950643cda5769a41fbcecc8b4d4194448f7c",rawRows:a.length,rawSpeciesLabels:uniq(a.map(x=>x.plant_species)).length,rawStudyLabels:uniq(a.map(x=>x.study_ID)).length,matched_four_arm_before_filters:valid.length,excluded,fruit:describe(fruits),fruit_excluding_macgregor2019:describe(fruits.filter(x=>x.study!=="Macgregor_2019")),fruit_open_ge_point2:describe(fruits.filter(x=>x.O>=.2)),seed:describe(seeds),wetMatches,cells:d};
}
const raw=readFileSync(process.argv[2],'utf8');
const wet=execFileSync('git',['show','b780370:results/all_species_forest/wet_network/species.csv'],{encoding:'utf8',maxBuffer:10*1024*1024});
console.log(JSON.stringify(analyze(raw,wet),null,2));
