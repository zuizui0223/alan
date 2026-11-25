// Replay on locally extracted source field study-label inventory and frozen 732-key routing table.
import {readFileSync} from 'node:fs';
function csv(s){let rows=[],r=[],f="",q=false;for(let i=0;i<s.length;i++){let c=s[i];if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}else if(c==='"'){q=true}else if(c===','){r.push(f);f=""}else if(c==='\n'){r.push(f.replace(/\r$/,""));rows.push(r);r=[];f=""}else f+=c}if(f||r.length){r.push(f.replace(/\r$/,""));rows.push(r)}let hh=rows.shift();if(rows.some(r=>r.length!==hh.length))throw Error("shape");return rows.map(r=>Object.fromEntries(hh.map((k,i)=>[k,r[i]])))}
function crosswalk(raw,route){
const src=raw.trim().split("\n").map(x=>{let [id,n,d]=x.split(" ");return {id,n:+n,dois:d.split("|")}});if(src.length!==128||new Set(src.map(x=>x.id)).size!==128||src.reduce((s,z)=>s+z.n,0)!==1676||route.length!==732)throw Error("source totals");
const lookup=new Map(route.map(x=>[x.report_key.toLowerCase(),x])),owners=new Map();for(const x of src)for(let d of x.dois){if(!owners.has(d))owners.set(d,[]);owners.get(d).push(x.id)}
const rows=src.map(x=>{let matched=x.dois.filter(d=>lookup.has("doi:"+d)),status,canonical_hint="";
if(x.dois.length>1)status="multiple_DOIs_hold";
else if(owners.get(x.dois[0]).length>1)status="shared_DOI_across_study_IDs_hold";
else if(x.id==="101")status="malformed_DOI_hold";
else status=matched.length?"exact_DOI_candidate_match":"unmatched_DOI_candidate";
if(x.id==="011")canonical_hint="10.1111/1365-2664.12927";
if(x.id==="188")canonical_hint="10.1111/gcb.16126";
if(x.id==="192")canonical_hint="10.1101/2022.12.11.519667";
const old=matched.length===1?lookup.get("doi:"+matched[0]):null;
return {study_id:x.id,field_rows:x.n,source_doi_strings:x.dois.join("|"),status,matched_report_keys:matched.map(d=>"doi:"+d).join("|"),possible_canonical_doi_manual_hint:canonical_hint,matched_priority:old?.priority||"",matched_programme_labels:old?.programme_labels||"",matched_title:old?.title||"",final_eligibility:"not_adjudicated",new_meta_k:0}});
const tally=Object.fromEntries([...new Set(rows.map(x=>x.status))].sort().map(k=>[k,{study_ids:rows.filter(x=>x.status===k).length,field_rows:rows.filter(x=>x.status===k).reduce((s,x)=>s+x.field_rows,0)}]));
return {rows,audit:{field:tally,field_with_any_exact_doi_match:rows.filter(x=>x.matched_report_keys).length,source_hash:"8fae8a0244047ab057808ea321962bf848cff42b6fd521f34b37a360efa154d6",source_study_labels:314,source_effect_rows:4905,field_study_labels:128,field_rows:1676,old_route_keys:732,source_doi_label_conflicts:[...owners].filter(([d,v])=>v.length>1),final_eligibility:false,warning:"DOI overlap candidates only; primary title and physical allocation not adjudicated; unmatched labels are not new independent studies"}};
}
const dir='database/alan_experiments/unified_meta_intake_20261125/czarnecka2025/';
const result=crosswalk(readFileSync(dir+'field_study_doi_inventory.txt','utf8'),csv(readFileSync('database/alan_experiments/unified_meta_intake_20261125/report_routing.csv','utf8')));
console.log(JSON.stringify(result.audit,null,2));
