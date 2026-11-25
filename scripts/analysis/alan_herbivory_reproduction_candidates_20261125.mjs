// Replay: node scripts/analysis/alan_herbivory_reproduction_candidates_20261125.mjs /path/to/author/seedset.csv
// Requires fixed PR101 objects and exact Li source CSV. No new original Grognuz download: fixed derived site means.
import {execFileSync} from 'node:child_process';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
function csv(s){let rows=[],r=[],f="",q=false;for(let i=0;i<s.length;i++){let c=s[i];if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}else if(c==='"'){q=true}else if(c===','){r.push(f);f=""}else if(c==='\n'){r.push(f.replace(/\r$/,""));rows.push(r);r=[];f=""}else f+=c}if(f||r.length){r.push(f.replace(/\r$/,""));rows.push(r)}let h=rows.shift();if(rows.some(r=>r.length!==h.length))throw Error("CSV shape");return rows.map(r=>Object.fromEntries(h.map((k,i)=>[k,r[i]])))}
function analyze(s,m,l){
let raw=csv(s),old=csv(m),li=csv(l);const mean=x=>x.reduce((a,b)=>a+b,0)/x.length;const varS=x=>x.reduce((a,b)=>a+(b-mean(x))**2,0)/(x.length-1);const group=(rows,f)=>{let a=new Map();for(let x of rows){let k=f(x);if(!a.has(k))a.set(k,[]);a.get(k).push(x)}return [...a.values()]};
const t2=4.302652729696142,t3=3.182446305284263;
const pairs=group(raw.filter(x=>x.stable_pair),x=>[x.stable_pair,x.P_PLANT_SP,x.date].join("|")).filter(z=>z.length===2&&new Set(z.map(x=>x.arm)).size===2).map(z=>({pair:z[0].stable_pair,species:z[0].P_PLANT_SP,date:z[0].date,lit:z.find(x=>x.arm==="LIT"),dark:z.find(x=>x.arm==="DARK")}));
if(pairs.length!==7||new Set(pairs.map(x=>x.pair)).size!==3)throw Error("support changed");
function slugFit(ep,latest=false){let use=latest?group(pairs,x=>x.pair+"|"+x.species).map(z=>z.sort((a,b)=>b.date.localeCompare(a.date))[0]):pairs;
if(use.some(x=>x.lit[ep]===""||x.dark[ep]===""||!(+x.lit[ep]>0)||!(+x.dark[ep]>0)))throw Error("unexpected nonpositive matched cell; no pseudocount");
const z=use.map(x=>({pair:x.pair,species:x.species,date:x.date,lit:+x.lit[ep],dark:+x.dark[ep],y:Math.log(+x.lit[ep]/+x.dark[ep])}));
const by=group(z,x=>x.pair).map(a=>{let sp=group(a,x=>x.species).map(g=>({species:g[0].species,dates:g.length,y:mean(g.map(x=>x.y))}));return {pair:a[0].pair,species:sp,y:mean(sp.map(x=>x.y))}});
let n=by.length,y=mean(by.map(x=>x.y)),v=varS(by.map(x=>x.y))/n;
return {ep,latest,n,y,v,rr:Math.exp(y),ci:[Math.exp(y-t2*Math.sqrt(v)),Math.exp(y+t2*Math.sqrt(v))],pair:by,cells:z,leaveOne:by.map((_,i)=>Math.exp(mean(by.filter((_,j)=>i!==j).map(x=>x.y))))};}
let leaf=slugFit("leaf_mean_damage_all_observed"),flower=slugFit("flower_mean_damage_all_observed"),leaf_latest=slugFit("leaf_mean_damage_all_observed",true),flower_latest=slugFit("flower_mean_damage_all_observed",true);
function pool(p){let k=p.length;function objective(t){let w=p.map(x=>1/(x.v+t)),sw=w.reduce((a,b)=>a+b,0),b=p.reduce((s,x,i)=>s+w[i]*x.y,0)/sw;return p.reduce((s,x,i)=>s+Math.log(x.v+t)+w[i]*(x.y-b)**2,0)+Math.log(sw)}
let a=0,b=10;for(let i=0;i<150;i++){let c=a+(b-a)/3,d=b-(b-a)/3;if(objective(c)<objective(d))b=d;else a=c}let tau=(a+b)/2;if(objective(0)<=objective(tau))tau=0;let w=p.map(x=>1/(x.v+tau)),sw=w.reduce((a,b)=>a+b,0),mu=p.reduce((a,x,i)=>a+w[i]*x.y,0)/sw,Q=p.reduce((a,x,i)=>a+w[i]*(x.y-mu)**2,0),se=Math.sqrt(Math.max(1,Q/(k-1))/sw),t=k===3?t2:t3;return {k,tau,rr:Math.exp(mu),ci:[Math.exp(mu-t*se),Math.exp(mu+t*se)],se,weights:w.map(x=>x/sw)}}
let p=old.map(x=>({name:x.programme,y:+x.log_ratio,v:+x.sampling_variance}));if(p.length!==3)throw Error("prior study count");
let meta3=pool(p),meta4=pool([...p,{name:"Grognuz2025_Seeland2023",y:leaf.y,v:leaf.v}]),meta4_latest=pool([...p,{name:"Grognuz2025_Seeland2023",y:leaf_latest.y,v:leaf_latest.v}]);
if(Math.abs(meta3.rr-.9287555116)>1e-7)throw Error("old meta replay mismatch");
function difference(a,b){const matched=a.pair.map(z=>{let w=b.pair.find(x=>x.pair===z.pair);if(!w)throw Error("pair mismatch");return {pair:z.pair,d:w.y-z.y}});const y=mean(matched.map(x=>x.d)),v=varS(matched.map(x=>x.d))/3;let reference=Array.from({length:8},(_,k)=>mean(matched.map((x,i)=>(k>>i)&1?x.d:-x.d))).filter(x=>Math.abs(x)>=Math.abs(y)-1e-12).length/8;return {y,v,rr:Math.exp(y),ci:[Math.exp(y-t2*Math.sqrt(v)),Math.exp(y+t2*Math.sqrt(v))],matched,signflip_reference_p:reference}}
if(li.length!==79)throw Error("Li size");
const lr=li.filter(x=>x.Seeds!=="NA"&&x.Fruit!=="NA"&&x.pollinated!=="NA");
if(lr.length!==78||lr.some(x=>+x.pollinated!==+(+x.Seeds>0)||+x.Starting<=0||+x.Seeds<0||+x.Fruit<0))throw Error("Li outcome source mismatch");
const fields=group(lr,x=>x.replicate).map(a=>{let arms={};for(let t of ["C","P","L"]){let z=a.filter(x=>x.Treatment===t);if(!z.length)throw Error("missing arm");arms[t]={n:z.length,seed:mean(z.map(x=>+x.Seeds)),fruit:mean(z.map(x=>+x.Fruit)),presence:mean(z.map(x=>+x.pollinated)),seed_per_starting:z.reduce((s,x)=>s+ +x.Seeds,0)/z.reduce((s,x)=>s+ +x.Starting,0),fruit_per_starting:z.reduce((s,x)=>s+ +x.Fruit,0)/z.reduce((s,x)=>s+ +x.Starting,0)}}return {field:a[0].replicate,arms}});
if(fields.length!==3)throw Error("field count");
function fieldRR(ep,t,ctrl="C"){let a=fields.map(x=>x.arms[t][ep]),b=fields.map(x=>x.arms[ctrl][ep]);let ma=mean(a),mb=mean(b);if(ma<=0||mb<=0)throw Error("nonpositive mean");let y=Math.log(ma/mb),scores=a.map((x,i)=>x/ma-b[i]/mb),v=varS(scores)/3;return {ep,t,ctrl,ma,mb,y,v,rr:Math.exp(y),ci:[Math.exp(y-t2*Math.sqrt(v)),Math.exp(y+t2*Math.sqrt(v))],scores,field_differences:a.map((x,i)=>x-b[i]),leaveOne:fields.map((_,i)=>mean(a.filter((_,j)=>i!==j))/mean(b.filter((_,j)=>i!==j)))}}
let liEffects=["seed","fruit","presence","seed_per_starting","fruit_per_starting"].flatMap(ep=>[fieldRR(ep,"L"),fieldRR(ep,"P"),fieldRR(ep,"L","P")]);
let liCov=liEffects.map(a=>liEffects.map(b=>a.scores.reduce((s,x,i)=>s+(x-mean(a.scores))*(b.scores[i]-mean(b.scores)),0)/6));
return {slug:{matchedCells:pairs.length,leaf,flower,leaf_latest,flower_latest,flower_vs_leaf:difference(leaf,flower),flower_vs_leaf_latest:difference(leaf_latest,flower_latest),prior_programmes:p,meta3,meta4,meta4_latest},li:{raw_rows:li.length,observed:lr.length,fields,effects:liEffects,covariance:liCov,covariance_order:liEffects.map(x=>[x.ep,x.t,x.ctrl].join("|"))}};
}
const source={"slug":{"commit":"437011edd2de22904535d648d47e56d94ff48ad2","path":"database/alan_experiments/slug_sources/wild_site_date_means.csv","blob":"d9f87775f99427eb089b3988a70ba7c830dbabe8"},"prior_leaf":{"commit":"437011edd2de22904535d648d47e56d94ff48ad2","path":"database/alan_experiments/leaf_damage_synthesis/programme_estimates.csv","blob":"d3e35d28325000fd6bb7e9a7ad9d019f2729a4f8"},"Li":{"repository":"Dongboli/experimental-data","path":"Artificial lights at night/seedset.csv","blob":"b616aca5beeba1457549841477d61ca4f796a147","doi":"10.64898/2026.07.10.737511"}};
if(!process.argv[2])throw Error('provide author seedset.csv');
const li=readFileSync(process.argv[2]),hash=createHash('sha1').update('blob '+li.length+'\0').update(li).digest('hex');
if(hash!==source.Li.blob)throw Error('Li blob mismatch '+hash);
const a=[source.slug,source.prior_leaf].map(x=>execFileSync('git',['show',x.blob],{encoding:'utf8',maxBuffer:10*1024*1024}));
const result=analyze(...a,li.toString('utf8'));result.source=source;result.status='EXPLORATORY_NOT_CANONICAL';result.new_formally_adopted_programmes=0;
console.log(JSON.stringify(result,null,2));
