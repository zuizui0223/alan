// node scripts/analysis/alan_stressor_arrangement_bounds_20261125.mjs [source JSON]
import {execFileSync} from 'node:child_process';
import {readFileSync} from 'node:fs';
function analyze(j){
const s=j.effect_summaries.drift.abund;if(j.n_independent_header_blocks!==6||s[0].treatment!=="L"||s[0].control!=="CT"||s[1].treatment!=="L+F"||s[1].control!=="F")throw Error("source treatment mismatch");
const L=s[0].block_treatment_values,C=s[0].block_control_values,F=s[1].block_control_values,LF=s[1].block_treatment_values;
const mean=x=>x.reduce((a,b)=>a+b,0)/x.length;const sd=x=>Math.sqrt(x.reduce((a,b)=>a+(b-mean(x))**2,0)/(x.length-1));const t=2.570581835636314;
if([L,C,F,LF].some(z=>z.length!==6||z.some(v=>!Number.isFinite(v)||v<=0)))throw Error("source values");
let b=C.map((c,i)=>({block:i+1,c,L:L[i],F:F[i],LF:LF[i],SS:L[i]+F[i],otherTotal:3*c+LF[i],thresholdShare:(L[i]+F[i]-LF[i])/(3*c),thresholdVsPooled:(L[i]+F[i]-LF[i])/c,CS_lower:LF[i],CS_upper:LF[i]+3*c}));
let S=mean(b.map(z=>z.SS)),O=mean(b.map(z=>z.otherTotal)),M=mean(C),K=mean(LF),crit=(S-K)/(3*M);
const total=4*b.reduce((s,z)=>s+z.SS+z.otherTotal,0);if(Math.abs(total-5564)>1e-8)throw Error("published total mismatch");
function pairedRatio(a,b){let ma=mean(a),mb=mean(b),l=Math.log(ma/mb),se=sd(a.map((x,i)=>x/ma-b[i]/mb))/Math.sqrt(6);return {rr:Math.exp(l),se,ci:[Math.exp(l-t*se),Math.exp(l+t*se)]}}
let SS=b.map(x=>x.SS),other=b.map(x=>x.otherTotal/2),dif=b.map(z=>z.SS-z.otherTotal/2);
let sens=[0,.25,1/3,.4,crit,.5,1].map(q=>{let cs=b.map(z=>z.LF+3*z.c*q),ct=b.map(z=>3*z.c*(1-q));return {q,CS_mean:mean(cs),CT_mean:mean(ct),CS_SS:mean(cs)/S,SS_CT:mean(ct)>0?S/mean(ct):null}});
return {source_sha256:j.source_sha256,independent_headers:6,raw_count_reconstruction:total,SS_raw_total:4*6*S,CT_CS_raw_total:4*6*O,block_bounds:b,means:{C:M,L:mean(L),F:mean(F),LF:K,SS:S,CT_CS_total:O,CT_CS_mean:O/2,identified_SS_vs_other_mean:pairedRatio(SS,other),difference:mean(dif),ci_difference:[mean(dif)-t*sd(dif)/Math.sqrt(6),mean(dif)+t*sd(dif)/Math.sqrt(6)],CS_over_SS_identification_interval:[K/S,O/S],threshold_share:crit,threshold_relative_to_pooled:3*crit,threshold_share_CT_gt_SS:1-S/(3*M)},conditional_scenarios:sens,source_limit:"CT pooled over two CT-network branches plus one CS-network branch. Unknown allocation cannot be recovered. Bounds are identification bounds, not confidence intervals; q scenarios are not observed outcomes."};
}
const source={"commit":"437011edd2de22904535d648d47e56d94ff48ad2","path":"database/alan_experiments/arthropod_expansion_20261122/eeres2021_source_reanalysis/original_reanalysis_20261123.json","blob":"068b7879feaca1d884e50b97bf21aa3b059ebe6f"};
const r=analyze(JSON.parse(process.argv[2]?readFileSync(process.argv[2],'utf8'):execFileSync('git',['show',source.commit+':'+source.path],{encoding:'utf8'})));
r.source_git=source;r.original_workbook_retrieved_this_turn=false;
console.log(JSON.stringify(r,null,2));
