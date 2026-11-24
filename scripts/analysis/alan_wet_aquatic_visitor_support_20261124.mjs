// Exploratory recorded aquatic-hoverfly support, not functional dependence.
// Run from repository root: node scripts/analysis/alan_wet_aquatic_visitor_support_20261124.mjs
import {execFileSync} from 'node:child_process';
import {gunzipSync} from 'node:zlib';
const ref='b780370';
const base='results/all_species_forest/wet_network/';
function read(name) {
  const path=ref+':'+base+name;
  const b=execFileSync('git',['show',path],{maxBuffer:20*1024*1024});
  return (name.endsWith('.gz')?gunzipSync(b):b).toString('utf8');
}
function csv(s){
  let rows=[],r=[],f='',q=false;
  for(let i=0;i<s.length;i++){
    let c=s[i];
    if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}
    else if(c==='"'){q=true}else if(c===','){r.push(f);f=''}
    else if(c==='\n'){r.push(f.replace(/\r$/,''));rows.push(r);r=[];f=''}else f+=c;
  }
  if(f||r.length){r.push(f.replace(/\r$/,''));rows.push(r)}
  const h=rows.shift();
  if(rows.some(r=>r.length!==h.length))throw Error('CSV shape');
  return rows.map(r=>Object.fromEntries(h.map((k,i)=>[k,r[i]])));
}
const F=csv(read('flower_sources.csv.gz')), W=csv(read('species.csv'));
const cnt=new Map();
for(const r of W)cnt.set(r.accepted,(cnt.get(r.accepted)||0)+1);
const U=W.filter(r=>r.accepted&&cnt.get(r.accepted)===1);
const by=new Map(U.map(r=>[r.accepted,r]));
const aquatic=new Set(F.filter(r=>r.aquatic_hoverfly==='True').map(r=>r.visitor));
const yes=x=>x.candidate==='True';
if(F.length!==3864||W.length!==363||U.length!==359||U.filter(yes).length!==58||aquatic.size!==29)throw Error('Frozen input changed');
if(F.some(r=>aquatic.has(r.visitor)&&r.aquatic_hoverfly!=='True'))throw Error('Inconsistent aquatic flags');
const group=(a,key)=>{let m=new Map();for(let x of a){let k=key(x);if(!m.has(k))m.set(k,[]);m.get(k).push(x)}return m};
const mean=x=>x.length?x.reduce((a,b)=>a+b,0)/x.length:null;
const count=(a,f)=>a.filter(f).length;
const panelKey=x=>JSON.stringify([x.archive,x.study,x.habitat]);
function agg(rows){
  const d=[];
  for(const [p,a] of group(rows,r=>r.plant)){
    const w=by.get(p);if(!w)continue;
    const v=[...new Set(a.map(x=>x.visitor))],A=v.filter(x=>aquatic.has(x)).length;
    const panels=group(a,panelKey),same=[...panels.values()].filter(g=>g.some(x=>aquatic.has(x.visitor)));
    const mixed=same.filter(g=>g.some(x=>!aquatic.has(x.visitor))).length;
    d.push({p,sp:w.sp,candidate:yes(w),region:w.site_region,island:w.island_band,sites:w.site_band,indicator:+w.indicator,n:v.length,A,share:A/v.length,nstudies:panels.size,nA:same.length,nMixed:mixed,alt:A>0?v.length-A:null});
  }
  return d;
}
function summary(d){
 return [true,false].map(c=>{
  const z=d.filter(x=>x.candidate===c),p=z.filter(x=>x.A>0);
  return {case:c,n:z.length,aquatic_positive:p.length,mean_recorded_aquatic_share:mean(z.map(x=>x.share)),positive_same_panel_other:p.filter(x=>x.nMixed>0).length,positive_no_same_panel_other:p.filter(x=>x.nMixed===0).length,positive_other_any_panel:p.filter(x=>x.alt>0).length,positive_species:p.map(x=>({p:x.p,a:x.A,n:x.n,nA:x.nA,nMixed:x.nMixed}))};
 });
}
function restricted(rows){
  return [...group(rows.filter(x=>by.has(x.plant)),panelKey).values()].filter(g=>g.some(x=>aquatic.has(x.visitor))&&g.some(x=>yes(by.get(x.plant)))&&g.some(x=>!yes(by.get(x.plant)))).flat();
}
const sets={all:F.filter(r=>by.has(r.plant)),wet:F.filter(r=>r.wet_site==='True'&&by.has(r.plant))};
sets.all_informative=restricted(sets.all);sets.wet_informative=restricted(sets.wet);
const degree=n=>n===1?0:n<=4?1:n<=15?2:3;
function stratified(d,field,kind){
 const ss=[...group(d,x=>(kind==='region_degree'?[x.region,degree(x.n)]:[x.region,x.island,x.sites,degree(x.n)]).join('|')).values()].filter(g=>g.some(x=>x.candidate)&&g.some(x=>!x.candidate));
 if(!ss.length)return null;
 const total=ss.reduce((s,g)=>s+count(g,x=>x.candidate),0);
 const obs=ss.reduce((s,g)=>{let k=g.filter(x=>x.candidate),j=g.filter(x=>!x.candidate);return s+k.length*(mean(k.map(field))-mean(j.map(field)))},0)/total;
 const refmean=ss.reduce((s,g)=>s+count(g,x=>x.candidate)*mean(g.filter(x=>!x.candidate).map(field)),0)/total;
 let state=7125;const rand=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296};let extreme=0;
 for(let t=0;t<9999;t++){
   let v=0;
   for(const g of ss){const a=g.map(field),k=count(g,x=>x.candidate);
     for(let i=a.length-1;i>0;i--){const j=Math.floor(rand()*(i+1));[a[i],a[j]]=[a[j],a[i]]}
     v+=k*(mean(a.slice(0,k))-mean(a.slice(k)));
   }
   if(Math.abs(v/total)>=Math.abs(obs)-1e-12)extreme++;
 }
 return {strata:ss.length,case:total,ref:ss.reduce((s,g)=>s+count(g,x=>!x.candidate),0),case_mean:refmean+obs,ref_standardized:refmean,difference:obs,p:(extreme+1)/10000};
}
const result={frozen_ref:ref,source_blobs:{flowers:'8d75a38af58c9fd63a1d1c99b3e69cf6e39f830e',species:'388a950643cda5769a41fbcecc8b4d4194448f7c'},total_rows:F.length,unique_unambiguous_cohort:U.length,case:count(U,yes),aquatic_taxa:[...aquatic].sort(),analyses:{}};
for(const [k,rows] of Object.entries(sets)){
 const d=agg(rows);
 result.analyses[k]={panels:group(rows,panelKey).size,records:rows.length,summary:summary(d),strict_share:stratified(d,x=>x.share,'strict'),relaxed_share:stratified(d,x=>x.share,'region_degree'),strict_positive:stratified(d,x=>+(x.A>0),'strict'),relaxed_positive:stratified(d,x=>+(x.A>0),'region_degree')};
}
console.log(JSON.stringify(result,null,2));
