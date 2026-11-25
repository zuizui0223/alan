// Replay fixed, audited source-summary values. Optional: supply exact local extracted field CSV to cross-check input.
import {execFileSync} from 'node:child_process';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
function csv(s){let rows=[],r=[],f="",q=false;for(let i=0;i<s.length;i++){let c=s[i];if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}else if(c==='"'){q=true}else if(c===','){r.push(f);f=""}else if(c==='\n'){r.push(f.replace(/\r$/,""));rows.push(r);r=[];f=""}else f+=c}if(f||r.length){r.push(f.replace(/\r$/,""));rows.push(r)}let a=rows.shift();if(rows.some(x=>x.length!==a.length))throw Error("shape");return rows.map(r=>Object.fromEntries(a.map((k,i)=>[k,r[i]])))}
function analyze(s,es,cs){
const sum=a=>a.reduce((x,y)=>x+y,0),mean=a=>sum(a)/a.length,odds=x=>x/(1-x),crit={1:12.706204736432095,4:2.7764451051977987,5:2.570581835636314};
const R=Object.fromEntries(s.source_extracted_rows.map(([id,d,l])=>[id,{d,l,rr:l/d}]));
const rr=id=>R[id].rr;
function sumrr(ids){let dark=sum(ids.map(i=>R[i].d)),lit=sum(ids.map(i=>R[i].l));return {dark,lit,rr:lit/dark}}
const a=sumrr(["175.13","175.14","175.15","175.16"]),b=sumrr(["175.17","175.18"]),all=sumrr(["175.13","175.14","175.15","175.16","175.17","175.18"]);
const fish={cladocera:a,copepoda:b,six_taxa:all,ratio_of_response_ratios:a.rr/b.rr,dark_cladocera_share:a.dark/all.dark,lit_cladocera_share:a.lit/all.lit,allocation_hold:"single platform lamp; lit and dark on opposite sides, dark cover; no independently allocated lamp replication"};
const seasons=["spring","summer","autumn"],man=seasons.map((season,i)=>{let a=R["253."+(10+i)],b=R["253."+(13+i)],[cd,cl]=s.manfrin_primary_A3_mean_catch_proportions[i];return {season,aquatic_rr:a.rr,terrestrial_rr:b.rr,relative_origin_rr:a.rr/b.rr,lit_share_ratio_of_means:a.l/(a.l+b.l),dark_share_ratio_of_means:a.d/(a.d+b.d),reported_mean_share_dark:cd,reported_mean_share_lit:cl,reported_mean_share_OR:odds(cl)/odds(cd)}});
const diet=s.manfrin_primary_A2_aquatic_carbon.flatMap(([taxon,vs])=>vs.map(([dm,lm,do_,lo],i)=>{let m=man[i],b=odds(lm)/odds(dm);return {taxon,season:seasons[i],diet_dark:dm,diet_lit:lm,diet_mode_dark:do_,diet_mode_lit:lo,diet_mean_OR:b,diet_mode_OR:odds(lo)/odds(do_),diet_vs_mean_catch_share_OR_ratio:b/m.reported_mean_share_OR,diet_vs_CPUE_OR_ratio:b/m.relative_origin_rr}}));
let ds={n_cells:9,posterior_mean_diet_decline:diet.filter(x=>x.diet_lit<x.diet_dark).length,posterior_mode_diet_decline:diet.filter(x=>x.diet_mode_lit<x.diet_mode_dark).length,diet_mean_OR_below_mean_catch_share_OR:diet.filter(x=>x.diet_vs_mean_catch_share_OR_ratio<1).length,diet_mean_OR_below_CPUE_OR:diet.filter(x=>x.diet_vs_CPUE_OR_ratio<1).length,independent_sites_per_arm:1};
const E=es.map(x=>+x.lnRR),V=es.map(x=>es.map(y=>+cs.find(z=>z.effect_id_1===x.effect_id&&z.effect_id_2===y.effect_id).sampling_covariance));
if(es.length!==6||es.some((x,i)=>Math.abs(V[i][i]-+x.sampling_variance)>1e-12||+x.n_lit!==6||+x.n_dark!==6))throw Error("Penryn variance mismatch");
function interval(y,v,df){return [Math.exp(y-crit[df]*Math.sqrt(v)),Math.exp(y+crit[df]*Math.sqrt(v))]}
function con(w){const y=sum(w.map((a,i)=>a*E[i])),v=sum(w.flatMap((a,i)=>w.map((b,j)=>a*b*V[i][j])));return {y,v,rr:Math.exp(y),ci:interval(y,v,5)}}
const penryn={herbivore_vs_plants:con([-.5,-.5,.5,.5,0,0]),herbivore_predator_vs_plants:con([-.5,-.5,0,0,.5,.5]),predator_vs_herbivore:con([0,0,-.5,-.5,.5,.5])};
function arm([d,ds,l,ls]){return {y:Math.log(l/d),v:(ds/d)**2+(ls/l)**2,d,l}}
const calibration=[0,2,4,1,3,5].map((j,i)=>{let x=es[j],r=arm(s.penryn_review_rows[i]);return {id:"191."+(i+1),raw_source_id:x.effect_id,raw_lnRR:+x.lnRR,compiled_lnRR:r.y,raw_v:+x.sampling_variance,compiled_v:r.v,variance_ratio:r.v/+x.sampling_variance}});
let arms=s.liu_arm_rows.map(arm),g=arms.slice(0,4).map((a,i)=>({y:arms[i+4].y-a.y,v:arms[i+4].v+a.v})),y=mean(g.map(x=>x.y));
let liu=[0,.5,1].map(rho=>{let v=sum(g.flatMap((a,i)=>g.map((b,j)=>(i===j?1:rho)*Math.sqrt(a.v*b.v))))/16;return {rho,y,v,rr:Math.exp(y),ci:interval(y,v,4)}});
function pool([a,b]){let tau=Math.max(0,((a.y-b.y)**2-a.v-b.v)/2),w=[1/(a.v+tau),1/(b.v+tau)],sw=sum(w),y=(w[0]*a.y+w[1]*b.y)/sw,Q=w[0]*(a.y-y)**2+w[1]*(b.y-y)**2,v=Math.max(1,Q)/sw;return {k:2,tau,rr:Math.exp(y),ci:interval(y,v,1),weights:w.map(x=>x/sw)}}
const out={source:s,fish,manfrin2018:{catch:man,diet,summary:ds,scope:"one dark/lit site pair reused across seasons and taxa; isotope mixing-model means vs catch counts; not a measured conversion efficiency, no SE or pooled effect"},yuen2017:{summary_point_ratios:["024.1","024.2","024.3","024.4","024.5"].map(id=>({id,rr:rr(id)})),primary_approaches:5375,primary_web_interceptions:300,primary_consumed:10,consumed_per_intercept:10/300,primary_model_coef:-2.47,primary_model_se:.88,complete_webs:21,sites:6,scope:"large/small groups in primary figure refer to spiders, not prey size; separate approaches, interceptions and actual consumption; coefficient kept on source scale; supplement not retrieved"},willmott2019:{summary_point_ratios:["049.7","049.8","049.9"].map(id=>({id,rr:rr(id)})),primary_weighted_capture_rr:1.58/.23,primary_reported_SE_lit:.31,primary_reported_SE_dark:.08,primary_25web_implied_SD_lit:.31*Math.sqrt(25),compiled_SD_lit_049_9:7.75,sites:10,nights:3,webs:48,moth_captured:0,spiders_removed_before_field:true,scope:"not actual consumption, size-weighted prey count index not mass; sites switched light condition, unknown site/night covariance"},penryn_validation:{rows:calibration,max_abs_logRR_diff:Math.max(...calibration.map(x=>Math.abs(x.raw_lnRR-x.compiled_lnRR))),variance_ratio_range:[Math.min(...calibration.map(x=>x.variance_ratio)),Math.max(...calibration.map(x=>x.variance_ratio))]},penryn,liu,plant_consumer_access_sensitivity:liu.map(x=>({rho:x.rho,herbivore_package:pool([penryn.herbivore_vs_plants,x]),herbivore_plus_predator_package:pool([penryn.herbivore_predator_vs_plants,x])})),formal_new_k:0,status:"EXPLORATORY_NOT_CANONICAL"};
return out
}
const s={
  "workbook_sha256": "8fae8a0244047ab057808ea321962bf848cff42b6fd521f34b37a360efa154d6",
  "field_csv_sha256": "b0856e58503b7f3bdb74d75f93d15c110b8794f780529fa4513b0a7b13ee9b5a",
  "source_extracted_rows": [
    [
      "024.1",
      30.605,
      31.601
    ],
    [
      "024.2",
      11.957,
      29.324
    ],
    [
      "024.3",
      3.333,
      2.266
    ],
    [
      "024.4",
      0.142,
      0.074
    ],
    [
      "024.5",
      0.174,
      0.264
    ],
    [
      "049.7",
      0.15,
      0.82
    ],
    [
      "049.8",
      0.3,
      1.29
    ],
    [
      "049.9",
      0.23,
      1.58
    ],
    [
      "175.13",
      77.65,
      140.96
    ],
    [
      "175.14",
      85.72,
      122.34
    ],
    [
      "175.15",
      5.85,
      10.06
    ],
    [
      "175.16",
      0.08,
      0.81
    ],
    [
      "175.17",
      52.21,
      41.65
    ],
    [
      "175.18",
      5.9,
      2.82
    ],
    [
      "253.10",
      0.44,
      10.3
    ],
    [
      "253.11",
      0.87,
      281.69
    ],
    [
      "253.12",
      0.09,
      1.24
    ],
    [
      "253.13",
      0.12,
      1.73
    ],
    [
      "253.14",
      0.75,
      35.26
    ],
    [
      "253.15",
      0.17,
      0.73
    ]
  ],
  "liu_arm_rows": [
    [
      7.07,
      0.85,
      7.62,
      0.83
    ],
    [
      5.25,
      0.44,
      5.5,
      0.4
    ],
    [
      8.5,
      1.36,
      9.01,
      1.36
    ],
    [
      6.03,
      0.47,
      5.37,
      0.47
    ],
    [
      6.88,
      0.89,
      7.5,
      0.91
    ],
    [
      5.59,
      0.46,
      5.65,
      0.49
    ],
    [
      8.77,
      1.4,
      9.24,
      1.36
    ],
    [
      6.05,
      0.49,
      5.8,
      0.49
    ]
  ],
  "penryn_review_rows": [
    [
      341.93,
      22.59,
      452.69,
      64.51
    ],
    [
      326.88,
      20.97,
      332.26,
      26.88
    ],
    [
      315.59,
      19.36,
      339.25,
      38.17
    ],
    [
      341.93,
      22.59,
      387.1,
      32.25
    ],
    [
      326.88,
      20.97,
      333.87,
      47.31
    ],
    [
      315.59,
      19.36,
      346.24,
      38.17
    ]
  ],
  "manfrin_primary_A3_mean_catch_proportions": [
    [
      0.72,
      0.85
    ],
    [
      0.49,
      0.64
    ],
    [
      0.26,
      0.6
    ]
  ],
  "manfrin_primary_A2_aquatic_carbon": [
    [
      "Pachygnatha clercki",
      [
        [
          0.29,
          0.26,
          0.3,
          0.26
        ],
        [
          0.21,
          0.35,
          0.21,
          0.34
        ],
        [
          0.34,
          0.33,
          0.33,
          0.33
        ]
      ]
    ],
    [
      "Pardosa prativaga",
      [
        [
          0.31,
          0.23,
          0.32,
          0.23
        ],
        [
          0.25,
          0.24,
          0.24,
          0.23
        ],
        [
          0.32,
          0.27,
          0.33,
          0.28
        ]
      ]
    ],
    [
      "Opiliones",
      [
        [
          0.32,
          0.28,
          0.32,
          0.28
        ],
        [
          0.21,
          0.2,
          0.2,
          0.2
        ],
        [
          0.31,
          0.26,
          0.32,
          0.25
        ]
      ]
    ]
  ],
  "primary_sources": {
    "yuen": "https://doi.org/10.7717/peerj.4070",
    "willmott": "https://doi.org/10.1111/eth.12852",
    "talanda": "https://doi.org/10.1093/plankt/fbac046",
    "manfrin": "https://doi.org/10.1111/oik.04696",
    "manfrin_accepted_manuscript": "https://qmro.qmul.ac.uk/xmlui/bitstream/handle/123456789/34664/Manfrin%20Dietary%20changes%20in%20predators%20and%20scavengers%20in%20a%20nocturnally%20illuminated%20riparian%20ecosystem%202018%20Accepted.pdf?sequence=1",
    "manfrin_data": "https://doi.org/10.5061/dryad.3hf58",
    "liu": "https://doi.org/10.1111/gcb.16126",
    "anic": "https://doi.org/10.1002/ece3.9157"
  },
  "penryn_source_blobs": {
    "database/alan_experiments/plant_sources/anic2022_effects.csv": "8a5bcd69c10cf682e0a0504251d640747260dad2",
    "database/alan_experiments/plant_sources/anic2022_covariance.csv": "a81712616d04909ae0b2842e12bff792fdcbd9c2"
  }
};
if(process.argv[2]){
 const bytes=readFileSync(process.argv[2]);
 if(createHash('sha256').update(bytes).digest('hex')!==s.field_csv_sha256)throw Error('field CSV hash mismatch');
 const r=Object.fromEntries(csv(bytes.toString('utf8')).map(x=>[x.ID,x]));
 for(const [id,d,l] of s.source_extracted_rows)if(+r[id].Control_Mean!==d||+r[id].Experimental_Mean!==l)throw Error('summary input mismatch '+id);
 for(const [prefix,values] of [['188',s.liu_arm_rows],['191',s.penryn_review_rows]])for(let i=0;i<values.length;i++){let x=r[prefix+'.'+(i+1)];let a=[+x.Control_Mean,+x.Control_SE,+x.Experimental_Mean,+x.Experimental_SE];if(JSON.stringify(a)!==JSON.stringify(values[i]))throw Error('summary arm mismatch');}
}
const [e,c]=Object.values(s.penryn_source_blobs).map(b=>csv(execFileSync('git',['show',b],{encoding:'utf8'})));
console.log(JSON.stringify(analyze(s,e,c),null,2));
