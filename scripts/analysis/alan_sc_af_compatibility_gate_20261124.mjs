// Node.js; execute from ALAN root with original Razanajatovo CSV as argv[2].
import {readFileSync} from 'node:fs';
import {execFileSync} from 'node:child_process';
import {gunzipSync} from 'node:zlib';
function csv(s){let rows=[],r=[],f="",q=false;for(let i=0;i<s.length;i++){let c=s[i];if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}else if(c==='"'){q=true}else if(c===','){r.push(f);f=""}else if(c==='\n'){r.push(f.replace(/\r$/,""));rows.push(r);r=[];f=""}else f+=c}if(f||r.length){r.push(f.replace(/\r$/,""));rows.push(r)}let h=rows.shift();if(rows.some(r=>r.length!==h.length))throw Error("CSV shape");return rows.map(r=>Object.fromEntries(h.map((k,i)=>[k,r[i]])))}
function analyze(raw,wetRaw,macroRaw){
 const a=csv(raw),wet=csv(wetRaw),macro=csv(macroRaw);
 const group=(a,f)=>{let m=new Map;for(let x of a){let k=f(x);if(!m.has(k))m.set(k,[]);m.get(k).push(x)}return m};
 const val=x=>x!==""&&Number.isFinite(+x)&&+x>=0&&+x<=1?+x:null;
 const names=group(a,x=>x.Species.replaceAll("_"," "));
 if(a.length!==1890||names.size!==1752||wet.length!==363||macro.length!==3901)throw Error("Frozen input size changed");
 const mean=x=>x.length?x.reduce((s,v)=>s+v,0)/x.length:null;
 function getSpecies(rows,e){
  let d=rows.map(x=>({study:x.Data_source,SC:val(x["Self-compatibility_index_"+e]),AF:val(x["Autofertility_index_"+e])}));let s=[...group(d,x=>x.study).values()].map(z=>({study:z[0].study,SC:mean(z.filter(x=>x.SC!==null).map(x=>x.SC)),AF:mean(z.filter(x=>x.AF!==null).map(x=>x.AF)),paired:z.filter(x=>x.SC!==null&&x.AF!==null).map(x=>({SC:x.SC,AF:x.AF,gap:x.SC-x.AF}))}));let p=s.filter(x=>x.paired.length).map(x=>({study:x.study,SC:mean(x.paired.map(y=>y.SC)),AF:mean(x.paired.map(y=>y.AF)),gap:mean(x.paired.map(y=>y.gap))}));return {SC:mean(s.filter(x=>x.SC!==null).map(x=>x.SC)),AF:mean(s.filter(x=>x.AF!==null).map(x=>x.AF)),pairSC:mean(p.map(x=>x.SC)),pairAF:mean(p.map(x=>x.AF)),gap:mean(p.map(x=>x.gap)),nPaired:p.length,studies:s}}
 const summarize=(arr,e)=>[true,false].map(c=>{let z=arr.filter(x=>x.c===c),p=z.filter(x=>x[e].gap!==null);return {c,n:z.length,SC_n:z.filter(x=>x[e].SC!==null).length,AF_n:z.filter(x=>x[e].AF!==null).length,paired:p.length,SC:mean(z.filter(x=>x[e].SC!==null).map(x=>x[e].SC)),AF:mean(z.filter(x=>x[e].AF!==null).map(x=>x[e].AF)),pairedSC:mean(p.map(x=>x[e].pairSC)),pairedAF:mean(p.map(x=>x[e].pairAF)),gap:mean(p.map(x=>x[e].gap)),SC_ge08_AF_le02:p.filter(x=>x[e].pairSC>=.8&&x[e].pairAF<=.2).length,low_SC:p.filter(x=>x[e].pairSC<=.2).length}});
 const global=[...names.entries()].map(([name,rows])=>({name,FS:getSpecies(rows,"FS"),SFL:getSpecies(rows,"SFL")}));
 let glob=["FS","SFL"].map(e=>{let p=global.filter(x=>x[e].gap!==null),high=p.filter(x=>x[e].pairSC>=.8);return {e,n:p.length,gap:mean(p.map(x=>x[e].gap)),highSC:high.length,highSC_lowAF:high.filter(x=>x[e].pairAF<=.2).length,thresholdSensitivity:[.7,.8,.9].flatMap(sc=>[.1,.2,.3].map(af=>({sc,af,n:p.filter(x=>x[e].pairSC>=sc).length,both:p.filter(x=>x[e].pairSC>=sc&&x[e].pairAF<=af).length})))}});
 function match(input,key){
   const dup=group(input,x=>x[key]);return input.filter(x=>x[key]&&dup.get(x[key]).length===1&&names.has(x[key])).map(x=>({name:x[key],sp:x.sp,c:x.candidate==="True",m:x.indicator===""?null:+x.indicator,reg:x.site_region,site:x.site_band,island:x.island_band,FS:getSpecies(names.get(x[key]),"FS"),SFL:getSpecies(names.get(x[key]),"SFL")}));
 }
 const z=match(wet,"accepted"),m=match(macro,"species");
 function strat(arr,e,f,strict){
   const d=arr.filter(x=>x[e][f]!==null&&x.m!==null&&Number.isFinite(x.m)&&x.reg!==""&&x.site!==""&&x.island!=="");
   const bins=[...group(d,x=>[x.reg,...(strict?[x.site,x.island]:[]),x.m>7?"wet":"other"].join("|")).values()].filter(v=>v.some(x=>x.c)&&v.some(x=>!x.c));
   const nC=bins.reduce((s,x)=>s+x.filter(y=>y.c).length,0),nR=bins.reduce((s,x)=>s+x.filter(y=>!y.c).length,0);
   const diff=nC?bins.reduce((s,x)=>s+x.filter(y=>y.c).length*(mean(x.filter(y=>y.c).map(y=>y[e][f]))-mean(x.filter(y=>!y.c).map(y=>y[e][f]))),0)/nC:null;
   return {strict,endpoint:e,metric:f,mixedStrata:bins.length,nC,nR,difference:diff,permutation:null,reason:"exploratory sparse selectively recorded traits; no inference"};
 }
 const tests=["FS","SFL"].flatMap(e=>["SC","AF","gap"].flatMap(f=>[true,false].map(s=>strat(m,e,f,s))));
 return {source:{source_rows:a.length,source_species:names.size,study_labels:new Set(a.map(x=>x.Data_source)).size,blobs:{Razan:"3654b0bcf30388d4c395e3e66972b67757ce1bad",wet:"388a950643cda5769a41fbcecc8b4d4194448f7c",macro:"72ee46c92794d4bb70069b043e0afd3f349c09d4"}},global:glob,wet:{FS:summarize(z,"FS"),SFL:summarize(z,"SFL"),matches:z},macro:{FS:summarize(m,"FS"),SFL:summarize(m,"SFL"),matches:m,matchedExploratory:tests},limits:"No light exposure in reproductive source; author indices capped at 1; observational secondary data, no inverse variance or phylogeny; exact-name matching only; no mechanism causality; wet cohort small."};
}
const gitShow=p=>execFileSync('git',['show','b780370:'+p],{maxBuffer:10*1024*1024});
console.log(JSON.stringify(analyze(readFileSync(process.argv[2],'utf8'),gitShow('results/all_species_forest/wet_network/species.csv').toString('utf8'),gunzipSync(gitShow('results/all_species_forest/decline_traits/moisture_evidence.csv.gz')).toString('utf8')),null,2));
