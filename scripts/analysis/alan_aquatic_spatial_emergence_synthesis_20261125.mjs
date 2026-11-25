// Recalculate after fetching pinned PR101 blobs. Source tables are rounded published summaries.
import {execFileSync} from 'node:child_process';
function analyze(source,sullivan,ogbeide,origin){
const sum=a=>a.reduce((s,v)=>s+v,0),mean=a=>sum(a)/a.length,variance=a=>sum(a.map(v=>(v-mean(a))**2))/(a.length-1),gm=a=>Math.exp(mean(a.map(Math.log)));
const crit={1:12.706204736432095,2:4.302652729696142,3:3.182446305284263,5:2.570581835636314,6:2.4469118511449692};function interval(y,v,df){return [Math.exp(y-crit[df]*Math.sqrt(v)),Math.exp(y+crit[df]*Math.sqrt(v))]}
function pool(z){let k=z.length;function f(t){let w=z.map(x=>1/(x.v+t)),sw=sum(w),mu=sum(z.map((x,i)=>x.y*w[i]))/sw;return sum(z.map((x,i)=>Math.log(x.v+t)+w[i]*(x.y-mu)**2))+Math.log(sw)}let a=0,b=100;for(let i=0;i<200;i++){let c=a+(b-a)/3,d=b-(b-a)/3;if(f(c)<f(d))b=d;else a=c}let tau=(a+b)/2;if(tau<1e-10||f(0)<=f(tau))tau=0;let w=z.map(x=>1/(x.v+tau)),sw=sum(w),y=sum(z.map((x,i)=>x.y*w[i]))/sw,Q=sum(z.map((x,i)=>w[i]*(x.y-y)**2)),v=Math.max(1,Q/(k-1))/sw;return {k,tau,y,v,RR:Math.exp(y),CI:interval(y,v,k-1),weights:w.map(x=>x/sw)}}
const month=["June","July","August","September"],dist=[0,3,40,75];
const table=source.spree_table2.map(([Al,Ad,Tl,Td],i)=>({month:month[Math.floor(i/4)],distance_m:dist[i%4],Al,Ad,Tl,Td,aquatic_RR:Al/Ad,terrestrial_RR:Tl/Td,RRR:(Al/Ad)/(Tl/Td)}));
for(let i=0;i<4;i++)if(JSON.stringify([...source.spree_table2[i*4+1],...source.spree_table2[i*4+3]])!==JSON.stringify(source.czarnecka_spree_rows[i]))throw Error("source extraction mismatch");
const distance=dist.map(d=>{let a=table.filter(x=>x.distance_m===d);let w=a.map(x=>x.month==="August"?2:3);let wm=k=>sum(a.map((x,i)=>x[k]*w[i]))/sum(w);return {distance_m:d,aquatic_gm_RR:gm(a.map(x=>x.aquatic_RR)),terrestrial_gm_RR:gm(a.map(x=>x.terrestrial_RR)),equal_month_RRR:gm(a.map(x=>x.RRR)),night_pair_weighted_RRR:(wm("Al")/wm("Ad"))/(wm("Tl")/wm("Td")),all_month_RRR:a.map(x=>x.RRR)}});
const same_shore=month.map((m,i)=>({month:m,RRR_0_over_40:table[4*i].RRR/table[4*i+2].RRR}));
const stage=source.spree_table4.map(([el,ed,bl,bd],i)=>({month:month[i],capture_RR:table[4*i].aquatic_RR,emergence_RR_direct:el/ed,emergence_RR_between:bl/bd,stage_ratio:table[4*i].aquatic_RR/(el/ed)}));
const s=sullivan.source_files.find(x=>x.provisional_block_contrasts).provisional_block_contrasts.results.AquaticEmergent;
const cs=Object.values(s.estimands);if(cs.length!==4||cs.some(x=>x.raw_block_lnRR.length!==3))throw Error("Sullivan support");
const blocks=Array.from({length:3},(_,i)=>mean(cs.map(x=>x.raw_block_lnRR[i])));const sy=mean(blocks),sv=variance(blocks)/3;
const S={programme:"Ohio2013_Sullivan",n:3,blocks,y:sy,v:sv,RR:Math.exp(sy),CI:interval(sy,sv,2),endpoint:"Emergence-trap count, mean of four colour log ratios within each block"};
const rows=ogbeide.trim().split(/\r?\n/).slice(1).map(x=>{let [week,treat,flume,count,rate,traps]=x.split(",");return {week,treat,flume,rate:+rate,traps:+traps}});if(rows.length!==32||rows.some(x=>x.traps!==6))throw Error("Ogbeide source support");
function og(weeks){let z=rows.filter(x=>["alan","control"].includes(x.treat)&&weeks.includes(x.week)),g={};for(let r of z){let k=r.treat+"|"+r.flume;(g[k]??=[]).push(r.rate)}let a=Object.entries(g);if(a.length!==8||a.some(([k,x])=>x.length!==weeks.length))throw Error("Ogbeide balance");let L=a.filter(([k,x])=>k.startsWith("alan|")).map(([k,x])=>mean(x)),D=a.filter(([k,x])=>k.startsWith("control|")).map(([k,x])=>mean(x)),l=mean(L),d=mean(D),y=Math.log(l/d),v=variance(L)/(4*l*l)+variance(D)/(4*d*d);return {programme:"Landau2022_Ogbeide",weeks,n_per_arm:4,L,D,y,v,RR:Math.exp(y),CI_conservative_t3:interval(y,v,3),endpoint:"Emergence trap rates, no crayfish, equal week mean per flume"}}
const O=og(["Week 1","Week 4"]),O4=og(["Week 4"]);
function can(rho,kind){const x=kind==="thesis"?source.canada_thesis_emergence:source.canada_compiled_emergence;const a=(x.se_lit??x.sd_lit/Math.sqrt(x.n))/x.lit,b=(x.se_dark??x.sd_dark/Math.sqrt(x.n))/x.dark,y=Math.log(x.lit/x.dark),v=a*a+b*b-2*rho*a*b;return {programme:"Canada2011_Perkin",kind,rho,y,v,RR:Math.exp(y),CI_conditional_t3:interval(y,v,3),status:"hypothetical covariance and marginal SE compatibility, not verified cluster variance"}}
const emergence={S,O,O4,k2:pool([S,O]),k2_week4:pool([S,O4]),k3_sensitivity:["thesis","compiled"].flatMap(kind=>[-1,0,.5,.9,1].map(rho=>{const p=can(rho,kind);return {kind,rho,perkin:p,pooled:pool([S,O,p]),week4_sensitivity:pool([S,O4,p])}}))};
const old=origin.trim().split(/\r?\n/).slice(1).map(x=>x.split(",")).filter(x=>x[0]==="source_defined_menu_mean_counts").map(x=>({programme:x[1],y:+x[2],v:+x[3],RR:+x[4]}));if(old.length!==2)throw Error("origin source");
const broad=distance.map(d=>({distance_m:d,programme_RRR:d.equal_month_RRR,three_programme_equal_geometric_mean:gm([...old.map(x=>x.RR),d.equal_month_RRR]),scope:"point summary only, no new uncertainty or independent meta k"}));
return {source,spree:{table,distance,same_shore,same_shore_equal_month_gm:gm(same_shore.map(x=>x.RRR_0_over_40)),stage,stage_equal_month_gm:gm(stage.map(x=>x.stage_ratio)),source_097_reconciled:"3m and 75m placements, aquatic/terrestrial origin in both, not aquatic/terrestrial traps; 0/40 omitted; one site, 11 night pairs, dark first, no independent site covariance"},emergence,origin_descriptive:{old,broad,existing_k2:pool(old)},formal_new_k:0,status:"EXPLORATORY_NOT_CANONICAL"};
}
const source={
  "thesis": "https://refubium.fu-berlin.de/bitstream/handle/fub188/10784/EK_Perkin_art_light_stream_ecosys.pdf?isAllowed=y&sequence=1",
  "perkin_spree_doi": "10.1111/fwb.12270",
  "perkin_canada_doi": "10.1111/fwb.12426",
  "spree_table2": [
    [
      34.55,
      0.44,
      1.19,
      0.25
    ],
    [
      3.05,
      0.47,
      0.37,
      0.4
    ],
    [
      0.58,
      0.32,
      0.26,
      0.34
    ],
    [
      0.26,
      0.27,
      0.28,
      0.21
    ],
    [
      17.71,
      0.61,
      1.62,
      0.35
    ],
    [
      4.26,
      0.43,
      0.92,
      0.38
    ],
    [
      0.62,
      0.38,
      0.33,
      0.38
    ],
    [
      0.21,
      0.09,
      0.2,
      0.23
    ],
    [
      2.52,
      0.23,
      0.27,
      0.11
    ],
    [
      0.65,
      0.24,
      0.13,
      0.12
    ],
    [
      0.17,
      0.11,
      0.1,
      0.13
    ],
    [
      0.05,
      0.13,
      0.08,
      0.08
    ],
    [
      0.69,
      0.16,
      0.14,
      0.1
    ],
    [
      0.16,
      0.15,
      0.06,
      0.09
    ],
    [
      0.12,
      0.16,
      0.08,
      0.08
    ],
    [
      0.1,
      0.12,
      0.04,
      0.06
    ]
  ],
  "spree_table4": [
    [
      31.83,
      14.85,
      4.88,
      4.12
    ],
    [
      36.4,
      14.49,
      20.79,
      16.57
    ],
    [
      10.37,
      4.08,
      13.44,
      13.16
    ]
  ],
  "canada_thesis_emergence": {
    "lit": 1.22,
    "dark": 1.25,
    "se_lit": 0.21,
    "se_dark": 0.17,
    "source_model_additive_coef": -0.02,
    "source_model_se": 0.233,
    "n_stream_pairs": 4
  },
  "canada_compiled_emergence": {
    "lit": 1.22,
    "dark": 1.3,
    "sd_lit": 0.4,
    "sd_dark": 0.4,
    "n": 4
  },
  "czarnecka_spree_rows": [
    [
      3.05,
      0.47,
      0.37,
      0.4,
      0.26,
      0.27,
      0.28,
      0.21
    ],
    [
      4.26,
      0.43,
      0.92,
      0.38,
      0.21,
      0.09,
      0.2,
      0.23
    ],
    [
      0.65,
      0.24,
      0.13,
      0.12,
      0.05,
      0.13,
      0.08,
      0.08
    ],
    [
      0.16,
      0.15,
      0.06,
      0.09,
      0.1,
      0.12,
      0.04,
      0.06
    ]
  ],
  "primary_page_indices": {
    "spree_table2": 84,
    "spree_table4": 87,
    "canada_results": 114,
    "canada_model": 129
  },
  "source_blobs": {
    "database/alan_experiments/arthropod_expansion_20261122/sullivan2019_source_reanalysis/dryad_provenance_and_frozen_block_results.json": "349ecebb168f1bfdaedf533c7066a764d69c4f89",
    "database/alan_experiments/ogbeide_sources/flume_emergence_rates.csv": "a774508cdbe30641e97de61a7f55d437e15b9aec",
    "database/alan_experiments/capture_origin_synthesis/programme_estimates.csv": "45271835c14bfbf5bf510bd078d547b379068e1e"
  }
};
const [s,o,c]=Object.values(source.source_blobs).map(b=>execFileSync('git',['show',b],{encoding:'utf8',maxBuffer:10*1024*1024}));
console.log(JSON.stringify(analyze(source,JSON.parse(s),o,c),null,2));
