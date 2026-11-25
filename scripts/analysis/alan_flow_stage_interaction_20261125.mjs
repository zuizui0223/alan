// Run: node scripts/analysis/alan_flow_stage_interaction_20261125.mjs [path/to/source.json]
// With no argument, read exact source commit/path (fetch PR101 objects first if absent).
import {execFileSync} from 'node:child_process';
import {readFileSync} from 'node:fs';
function analyze(j){
 const list=["drift","benthos"].flatMap(stage=>["abund","rich"].flatMap(metric=>j.effect_summaries[stage][metric].map((e,i)=>({stage,metric,cond:i===0?"normal":"ponded",...e}))));
 if(list.length!==8||j.n_independent_header_blocks!==6||list.some(e=>e.block_treatment_values.length!==6||e.block_control_values.length!==6||e.block_treatment_values.concat(e.block_control_values).some(x=>!Number.isFinite(x)||x<=0)))throw Error("source size or zero");
 const C=j.joint_shared_HT.covariance_logRR;
 const mu=list.map(x=>x.lnRR), n=6, crit=2.570581835636314;
 const mean=x=>x.reduce((a,b)=>a+b,0)/x.length, varS=x=>x.reduce((s,v)=>s+(v-mean(x))**2,0)/(x.length-1);
 const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
 const order=list.map(e=>`${e.stage}|${e.metric}|${e.cond==="normal"?"L_vs_CT":"LF_vs_F"}`);
 if(JSON.stringify(order)!==JSON.stringify(j.joint_shared_HT.source_outcomes))throw Error("order mismatch");
 for(let e of list)if(Math.abs(mean(e.block_treatment_values)-e.mean_treatment)>1e-12||Math.abs(mean(e.block_control_values)-e.mean_control)>1e-12||Math.abs(Math.log(e.mean_treatment/e.mean_control)-e.lnRR)>1e-12)throw Error("mean mismatch");
 const scores=list.map(e=>Array.from({length:n},(_,i)=>e.block_treatment_values[i]/e.mean_treatment-e.block_control_values[i]/e.mean_control));
 const cc=scores.map(a=>scores.map(b=>a.reduce((s,v,i)=>s+(v-mean(a))*(b[i]-mean(b)),0)/(n*(n-1))));
 const maxdiff=Math.max(...cc.flatMap((a,i)=>a.map((v,k)=>Math.abs(v-C[i][k]))));
 if(maxdiff>1e-12)throw Error("cov mismatch");
 const result={};
 function pT(x){x=Math.abs(x);let steps=20000,h=x/steps,fn=t=>8/(3*Math.PI*Math.sqrt(5))*Math.pow(1+t*t/5,-3);let s=fn(0)+fn(x);for(let i=1;i<steps;i++)s+=(i%2?4:2)*fn(i*h);return Math.max(0,1-2*s*h/3);}
 for(let [metric,ix] of [["abundance",[0,1,4,5]],["richness",[2,3,6,7]]]){
  const w=Array(8).fill(0);ix.forEach((k,i)=>w[k]=[-1,1,1,-1][i]);
  const b=dot(w,mu),variance=dot(w,C.map(a=>dot(w,a))),se=Math.sqrt(variance);
  const perBlock=Array.from({length:6},(_,i)=>dot(w,list.map(e=>Math.log(e.block_treatment_values[i]/e.block_control_values[i]))));
  const loo=Array.from({length:6},(_,omit)=>Math.exp(dot(w,list.map(e=>Math.log(mean(e.block_treatment_values.filter((_,i)=>i!==omit))/mean(e.block_control_values.filter((_,i)=>i!==omit)))))));
  const bm=mean(perBlock),bs=Math.sqrt(varS(perBlock)/6);
  let pSign=Array.from({length:64},(_,mask)=>mean(perBlock.map((v,i)=>(mask>>i)&1?v:-v))).filter(x=>Math.abs(x)>=Math.abs(bm)-1e-12).length/64;
  result[metric]={normalStageRatio:Math.exp(mu[ix[0]]-mu[ix[2]]),pondedStageRatio:Math.exp(mu[ix[1]]-mu[ix[3]]),interaction_logRR:b,se,ratio:Math.exp(b),ci:[Math.exp(b-crit*se),Math.exp(b+crit*se)],p_t5:pT(b/se),leaveOneBlockOut:loo,blockLogRatios:perBlock,blockWeighted:{ratio:Math.exp(bm),ci:[Math.exp(bm-crit*bs),Math.exp(bm+crit*bs)],p_sign_reference:pSign}};
 }
 const p=Object.entries(result).sort((a,b)=>a[1].p_t5-b[1].p_t5);p[0][1].holm=Math.min(1,2*p[0][1].p_t5);p[1][1].holm=Math.min(1,Math.max(p[0][1].holm,p[1][1].p_t5));
 return {source:j.source_sha256,source_programme:j.physical_programme,source_order:order,cov_max_abs_diff:maxdiff,result};
}
const source={"commit":"437011edd2de22904535d648d47e56d94ff48ad2","path":"database/alan_experiments/arthropod_expansion_20261122/eeres2021_source_reanalysis/original_reanalysis_20261123.json","blob":"068b7879feaca1d884e50b97bf21aa3b059ebe6f"};
const j=JSON.parse(process.argv[2]?readFileSync(process.argv[2],'utf8'):execFileSync('git',['show',source.commit+':'+source.path],{encoding:'utf8',maxBuffer:1024*1024}));
const result=analyze(j);result.source_git=source;
console.log(JSON.stringify(result,null,2));
