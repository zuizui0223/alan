// Node replay of executed analysis. Fixed repository blobs; no coordinates printed.
import {execFileSync} from 'node:child_process';
import {gunzipSync} from 'node:zlib';
function csv(s){let a=s.trim().split(/\r?\n/).map(x=>x.split(",")),h=a.shift();if(a.some(x=>x.length!==h.length))throw Error("CSV columns");return a.map(r=>Object.fromEntries(h.map((k,i)=>[k,r[i]])))}
function analyze(raw,dist,...geo){
 const A=csv(raw),D=csv(dist),dm=new Map(D.filter(x=>x.distance_m!==""&&Number.isFinite(+x.distance_m)).map(x=>[x.unit_id,+x.distance_m]));
 const g=A.filter(x=>dm.has(x.unit_id)).map(x=>({...Object.fromEntries(Object.entries(x).map(([k,v])=>[k,k==="state"?v:v===""?NaN:+v])),distance_m:dm.get(x.unit_id)}));
 if(g.length!==1706||new Set(g.map(x=>x.unit_id)).size!==1706||g.some(x=>x.area0<0))throw Error("source");
 const env=["light0","temp0","rain0","temp_trend","rain_trend","built0","built","effort","days","spring","summer","autumn","start","span","years","island","nc","nr","area0","freq0","area_trend","freq_trend"];
 const sum=x=>x.reduce((s,v)=>s+v,0),mean=x=>sum(x)/x.length,dot=(a,b)=>sum(a.map((x,i)=>x*b[i]));
 const byBlock=(rows)=>{let b=new Map();for(let x of rows){if(!b.has(x.block))b.set(x.block,[]);b.get(x.block).push(x)}return [...b.values()]};
 for(let x of g){x.L=x.light/10;x.N=1-x.distance_m/1000;x.W=Math.log2(1+x.area0);x.LN=x.L*x.N;x.LW=x.L*x.W;x.NW=x.N*x.W;x.LNW=x.L*x.N*x.W;}
 function fit(frame,three,target){
  const blocks=byBlock(frame).filter(z=>z.length>1&&Math.max(...z.map(x=>x.light))-Math.min(...z.map(x=>x.light))>1e-8),n=blocks.flat().length,kb=blocks.length;
  const center=f=>blocks.flatMap(b=>{let z=b.map(f),m=mean(z);return z.map(v=>v-m)});
  const pred=["L","N","W","LN","LW","NW",...(three?["LNW"]:[]),...env],cols=pred.map(k=>center(x=>x[k])),x=cols[pred.indexOf(target)],control=cols.filter((_,i)=>pred[i]!==target);
  if(cols.some(x=>x.some(v=>!Number.isFinite(v))))throw Error("nonfinite");
  const Q=[];
  function residual(x){let v=x.slice();for(let rep=0;rep<2;rep++)for(let q of Q){let c=dot(v,q);for(let i=0;i<v.length;i++)v[i]-=c*q[i]}return v}
  for(let c of control){let norm=Math.sqrt(dot(c,c));if(norm<1e-10)continue;let r=residual(c.map(v=>v/norm)),nr=Math.sqrt(dot(r,r));if(nr>1e-9)Q.push(r.map(v=>v/nr))}
  let xr=residual(x),xx=dot(xr,xr);if(xx<1e-12)return {target,three,status:"no support"};
  let v=xr.map(z=>z/xx),df=n-kb-Q.length-1;if(df<1)return {target,three,status:"no residual df"};
  let shares=[],off=0;for(let b of blocks){shares.push(sum(v.slice(off,off+b.length).map(x=>x*x)));off+=b.length}
  let sm=sum(shares);shares=shares.map(x=>x/sm);let k=shares.filter(x=>x>1e-9).length,eff=1/sum(shares.map(x=>x*x));
  const t0=1.95996398454005,tdf=k-1,t=t0+(t0**3+t0)/(4*tdf)+(5*t0**5+16*t0**3+3*t0)/(96*tdf**2)+(3*t0**7+19*t0**5+17*t0**3-15*t0)/(384*tdf**3);
  let estimates={};
  for(let outcome of ["candidate_rate","reference_rate","rate_difference"]){let yr=residual(center(z=>z[outcome])),b=dot(xr,yr)/xx,resid=yr.map((y,i)=>y-b*xr[i]);let scores=[],off=0;for(let a of blocks){let sc=0;for(let i=off;i<off+a.length;i++)sc+=v[i]*resid[i];scores.push(sc);off+=a.length}let se=Math.sqrt(sum(scores.map(x=>x*x))*k/(k-1)*(n-1)/df);estimates[outcome]={b,se,ci:[b-t*se,b+t*se]}}
  return {target,three,available:frame.length,n,blocks:kb,effective_blocks:eff,df,rank:Q.length+1,estimates};
 }
 const q=(a,p)=>{let x=a.slice().sort((a,b)=>a-b),j=(x.length-1)*p,k=Math.floor(j);return x[k]+(x[Math.min(k+1,x.length-1)]-x[k])*(j-k)};
 const corr=(x,y)=>{let a=x.map(v=>v-mean(x)),b=y.map(v=>v-mean(y));return dot(a,b)/Math.sqrt(dot(a,a)*dot(b,b))};
 const geom=new Map();for(const s of geo)for(let x of csv(s)){if(x.latitude===""||x.longitude===""||!Number.isFinite(+x.latitude)||!Number.isFinite(+x.longitude))continue;const id=+x.unit_id,z={lat:+x.latitude,lon:+x.longitude};if(geom.has(id)&&(Math.abs(geom.get(id).lat-z.lat)>1e-7||Math.abs(geom.get(id).lon-z.lon)>1e-7))throw Error("geo conflict");geom.set(id,z)}
 const observedGeo=g.filter(x=>geom.has(x.unit_id));
 let out={source_gate:{real_coordinate_ids:geom.size,matched_panel:observedGeo.length,matched_brightened:observedGeo.filter(x=>x.light>1e-8).length,matched_blocks:new Set(observedGeo.map(x=>x.block)).size,matched_near100:observedGeo.filter(x=>x.distance_m<=100).length,raw_raster_retrieved:false,actual_patch_metrics_estimated:false},scopes:{}};
 const scopes={all:g,water_detected:g.filter(x=>x.distance_m<1000),water_rich:g.filter(x=>x.area0>=1&&x.distance_m<1000),bright_water_rich:g.filter(x=>x.area0>=1&&x.distance_m<1000&&x.light>1e-8)};
 for(let [name,frame] of Object.entries(scopes)){
  let support=["near","far"].map(k=>{let z=frame.filter(x=>k==="near"?x.distance_m<=100:x.distance_m>100);return {k,n:z.length,bright:z.filter(x=>x.light>1e-8).length,water_quantiles:[.1,.5,.9].map(p=>q(z.map(x=>x.area0),p))}});
  out.scopes[name]={n:frame.length,bright:frame.filter(x=>x.light>1e-8).length,support,models:[fit(frame,false,"LW"),fit(frame,false,"LN"),fit(frame,true,"LNW")]};
 }
 out.common_support_100m_1pct=["near_low","near_high","far_low","far_high"].map(k=>{let z=g.filter(x=>(x.distance_m<=100)==k.startsWith("near")&&(x.area0>=1)==k.endsWith("high"));return {cell:k,n:z.length,bright:z.filter(x=>x.light>1e-8).length,blocks:new Set(z.map(x=>x.block)).size}});
 out.correlations={all:corr(g.map(x=>x.N),g.map(x=>x.W)),detected:corr(scopes.water_detected.map(x=>x.N),scopes.water_detected.map(x=>x.W))};
 return out;
}
const blobs={"results/all_species_forest/wet_water/coordinate_trends.csv.gz":"1c025bcf6c28a6a3cd125ad05628169e126aef0d","results/all_species_forest/wet_shore/distance.csv":"c2bd22ab52c4674fedfb39a12f67cf847a5e2c3d","results/all_species_forest/wet_annual/habitats.csv":"6c464c2c2d2a9ddc71a95cab912175c263a90ff0","results/all_species_forest/wet_shore/nearest_wells.csv":"142bafb45550f73124e6342d694d2871ad3b5f15","results/all_species_forest/seven_moth/nearest_wells.csv":"f36def0fefdff9ea3c4486183cc0804a84e1116a"};
const a=Object.values(blobs).map(s=>execFileSync('git',['show',s],{maxBuffer:10*1024*1024}));
const result=analyze(gunzipSync(a[0]).toString('utf8'),...a.slice(1).map(x=>x.toString('utf8')));
result.source_blobs=blobs;
console.log(JSON.stringify(result,null,2));
