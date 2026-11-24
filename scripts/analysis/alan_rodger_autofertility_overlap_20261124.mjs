// Run from ALAN root: node scripts/analysis/alan_rodger_autofertility_overlap_20261124.mjs /path/to/rodger_2021_pollinator_contribution_snapshot.csv.gz /path/to/rodger_2021_two_backbone_synonyms_20260811.csv
import {readFileSync} from 'node:fs';
import {execFileSync} from 'node:child_process';
import {gunzipSync} from 'node:zlib';
function csv(s){let rows=[],r=[],f="",q=false;for(let i=0;i<s.length;i++){let c=s[i];if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}else if(c==='"'){q=true}else if(c===','){r.push(f);f=""}else if(c==='\n'){r.push(f.replace(/\r$/,""));rows.push(r);r=[];f=""}else f+=c}if(f||r.length){r.push(f.replace(/\r$/,""));rows.push(r)}let h=rows.shift();if(rows.some(r=>r.length!==h.length))throw Error("CSV shape");return rows.map(r=>Object.fromEntries(h.map((k,i)=>[k,r[i]])))}
function analyze(raw,synRaw,wetRaw,macroRaw){
const a=csv(raw),syn=csv(synRaw),w=csv(wetRaw),m=csv(macroRaw);
if(a.length!==1528||w.length!==363||m.length!==3901||w.filter(x=>x.candidate==="True").length!==58)throw Error("Unexpected frozen source size");
const ac=Object.keys(a[0]).filter(x=>x.startsWith("auto."));
const group=(a,f)=>{let m=new Map();for(let x of a){let k=f(x);if(!m.has(k))m.set(k,[]);m.get(k).push(x)}return m};
const p=a.map(r=>({name:r["genus.species"].replaceAll("_"," "),taxon:r.taxon.replaceAll("_"," "),study:r.study,family:r["TPL.family"],source_dataset:r["source.dataset"],row:r.source_row,pc_unique_number:r["pc.unique.number"],year:r["year.study.estimate"],vals:Object.fromEntries(ac.filter(c=>r[c]!==""&&Number.isFinite(+r[c])&&+r[c]>=0).map(c=>[c,+r[c]]))}));
const sm=new Map(syn.map(x=>[x.submitted_name,x.accepted_species]));
function classify(z){let v=z.flatMap(x=>Object.values(x.vals));return !v.length?"missing":v.some(x=>x>0)&&v.some(x=>x===0)?"mixed":v.some(x=>x>0)?"positive":"zero"}
function summarize(data,wet,alias){
let names=group(p,x=>alias?(sm.get(x.name)??x.name):x.name);let dup=group(data,x=>wet?x.accepted:x.species);let matches=[];
for(let x of data){let name=wet?x.accepted:x.species;if(!name||dup.get(name).length!==1)continue;let z=names.get(name);if(!z)continue;matches.push({name,sp:x.sp,candidate:x.candidate==="True",moist:x.indicator===""?null:+x.indicator,region:x.site_region,site_band:x.site_band,island_band:x.island_band,classification:classify(z),n:z.length,nstudy:new Set(z.map(x=>x.study)).size,rowStates:z.map(x=>classify([x])),rows:z});}
let summary=[true,false].map(c=>{let d=matches.filter(x=>x.candidate===c);return {case:c,total:d.length,known:d.filter(x=>x.classification!=="missing").length,positive:d.filter(x=>x.classification==="positive").length,zero:d.filter(x=>x.classification==="zero").length,mixed:d.filter(x=>x.classification==="mixed").length,missing:d.filter(x=>x.classification==="missing").length,knownMoist:d.filter(x=>x.moist!==null).length}});
return {summary,matches};
}
const wet=summarize(w,true,false),wetS=summarize(w,true,true),macro=summarize(m,false,false),macroS=summarize(m,false,true);
const uniqueWet=w.filter(x=>x.accepted&&w.filter(y=>y.accepted===x.accepted).length===1);
const wm=new Map(wet.matches.map(x=>[x.sp,x]));
const coverage=[true,false].map(c=>{
let z=uniqueWet.filter(x=>(x.candidate==="True")===c),named=z.filter(x=>x.named_flower==="True"),pl=z.filter(x=>x.pollen_limitation!==""&&Number.isFinite(+x.pollen_limitation)),pc=z.filter(x=>wm.has(x.sp));return {case:c,cohort:z.length,named_flower:named.length,pollen_limitation:pl.length,exclusion_snapshot:pc.length,exclusion_and_pollen_limitation:pc.filter(x=>pl.includes(x)).map(x=>({name:x.accepted,PL:+x.pollen_limitation,classification:wm.get(x.sp).classification})),exclusion_and_named_flower:pc.filter(x=>named.includes(x)).length,exclusion_and_wet_aquatic_hoverfly:pc.filter(x=>x.wet_aquatic_hoverfly==="True").length,all_three:pc.filter(x=>pl.includes(x)&&x.wet_aquatic_hoverfly==="True").length};});
return {source:{rawRows:a.length,uniqueSpecies:new Set(p.map(x=>x.name)).size,studyLabels:new Set(p.map(x=>x.study)).size,autoColumns:ac,sourceOnlyExclusion:true,naturalPollinationColumnsPresent:false,sourceBlobs:{rodger:"3ec30a6533880ebfa497aa0f0800fd94d3619bf3",synonyms:"7ed435b1c00c1f8e0ea39a44b0abfc2abeeb6ec9",wet:"388a950643cda5769a41fbcecc8b4d4194448f7c",macro:"72ee46c92794d4bb70069b043e0afd3f349c09d4"}},wet,wetSyn:wetS,macro,macroSyn:macroS,coverage,taxonomy_hold:wet.matches.filter(x=>x.candidate).map(x=>x.name),PC_estimable:false};
}
const raw=gunzipSync(readFileSync(process.argv[2])).toString('utf8'),syn=readFileSync(process.argv[3],'utf8');
const gitShow=p=>execFileSync('git',['show','b780370:'+p],{maxBuffer:10*1024*1024});
console.log(JSON.stringify(analyze(raw,syn,gitShow('results/all_species_forest/wet_network/species.csv').toString('utf8'),gunzipSync(gitShow('results/all_species_forest/decline_traits/moisture_evidence.csv.gz')).toString('utf8')),null,2));
