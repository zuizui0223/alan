// Node.js replay from ALAN root. Both source blobs are fixed.
import {execFileSync} from 'node:child_process';
function csv(s){let rows=[],r=[],f="",q=false;for(let i=0;i<s.length;i++){let c=s[i];if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}else if(c==='"'){q=true}else if(c===','){r.push(f);f=""}else if(c==='\n'){r.push(f.replace(/\r$/,""));rows.push(r);r=[];f=""}else f+=c}if(f||r.length){r.push(f.replace(/\r$/,""));rows.push(r)}let h=rows.shift();if(rows.some(r=>r.length!==h.length))throw Error("CSV shape");return rows.map(r=>Object.fromEntries(h.map((k,i)=>[k,r[i]])))}
function analyze(raw,macro){
 const a=csv(raw),w=csv(macro),wm=new Map(w.map(x=>[x.sp,x])),dups=new Map();for(const x of w)dups.set(x.accepted,(dups.get(x.accepted)||0)+1);
 if(a.length!==363||w.length!==363||w.filter(x=>x.candidate==="True").length!==58)throw Error("frozen input size");
 const rows=a.filter(x=>x.status==="verified"&&x.accepted&&dups.get(x.accepted)===1&&wm.has(x.sp)).map(x=>{
 const z=wm.get(x.sp);if(x.accepted!==z.accepted||x.candidate!==z.candidate)throw Error("Mismatch");
 const p=x["Pollination syndrome"],b=x["Generative reproduction type"],r=x["Reproduction type"];
 const insect=p?p.includes("insect-pollination")&&!p.includes("wind-pollination")&&!p.includes("water-pollination"):null;
 const strict=b?["allogamy","allogamy self-incompatibility"].includes(b):null;
 const SI=b?b==="allogamy self-incompatibility":null;
 const seed=r?r==="only by seed/spores":null;
 const veg=r?r.includes("vegetativ"):null;
 const both=insect===null||strict===null?null:insect&&strict;
 const joint=both===null||seed===null?null:both&&seed;
 const siconstraint=insect===null||SI===null||seed===null?null:insect&&SI&&seed;
 const group=both===null||seed===null?null:both?(seed?"insect_strict_seed_only":"insect_strict_with_vegetative"):(seed?"other_seed_only":"other_with_vegetative");
 return {sp:x.sp,name:x.accepted,c:z.candidate==="True",reg:z.site_region,site:z.site_band,island:z.island_band,m:+z.indicator,p,b,r,insect,strict,SI,seed,veg,both,joint,siconstraint,group};
 });
 const mean=x=>x.length?x.reduce((a,b)=>a+b,0)/x.length:null;
 const grp=(a,key)=>{let m=new Map;for(let x of a){let k=key(x);if(!m.has(k))m.set(k,[]);m.get(k).push(x)}return [...m.values()]};
 function test(d,k,mode){const bins=grp(d,x=>mode==="strict"?[x.reg,x.site,x.island].join("|"):mode==="region"?x.reg:"all").filter(v=>v.some(x=>x.c)&&v.some(x=>!x.c));const n=bins.reduce((s,x)=>s+x.filter(y=>y.c).length,0),nR=bins.reduce((s,x)=>s+x.filter(y=>!y.c).length,0);if(!n)return {mode,n,nR,bins:0,diff:null,p:null};let observed=bins.reduce((s,x)=>{let c=x.filter(y=>y.c),r=x.filter(y=>!y.c);return s+c.length*(mean(c.map(y=>+y[k]))-mean(r.map(y=>+y[k])))},0)/n;let state=7125;const rand=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state/4294967296};let extreme=0;for(let j=0;j<9999;j++){let t=0;for(let x of bins){let a=x.map(y=>+y[k]),nc=x.filter(y=>y.c).length;for(let i=a.length-1;i>0;i--){let j=Math.floor(rand()*(i+1));[a[i],a[j]]=[a[j],a[i]]}t+=nc*(mean(a.slice(0,nc))-mean(a.slice(nc)))}if(Math.abs(t/n)>=Math.abs(observed)-1e-12)extreme++}return {mode,n,nR,bins:bins.length,diff:observed,p:(extreme+1)/10000};}
 let result={inputs:{cohort:w.length,verified_unique:rows.length,trait_blob:"c603e16c96f5ee90d06df05205836dba9d206198",cohort_blob:"388a950643cda5769a41fbcecc8b4d4194448f7c"},traits:{},completeProfiles:[],rows};
 for(const k of ["insect","strict","SI","seed","both","joint","siconstraint"]){const d=rows.filter(x=>x[k]!==null);result.traits[k]={counts:[true,false].map(c=>({case:c,n:d.filter(x=>x.c===c).length,yes:d.filter(x=>x.c===c&&x[k]).length})),tests:["all","region","strict"].map(m=>test(d,k,m)),positives:d.filter(x=>x[k]).map(x=>({name:x.name,case:x.c,p:x.p,b:x.b,r:x.r}))}}
 const prof=rows.filter(x=>x.group!==null);result.completeProfiles=[true,false].map(c=>({case:c,n:prof.filter(x=>x.c===c).length,groups:Object.fromEntries(["insect_strict_seed_only","insect_strict_with_vegetative","other_seed_only","other_with_vegetative"].map(g=>[g,prof.filter(x=>x.c===c&&x.group===g).length]))}));return result;
}
const blob=sha=>execFileSync('git',['show',sha],{encoding:'utf8',maxBuffer:10*1024*1024});
console.log(JSON.stringify(analyze(blob('c603e16c96f5ee90d06df05205836dba9d206198'),blob('388a950643cda5769a41fbcecc8b4d4194448f7c')),null,2));
