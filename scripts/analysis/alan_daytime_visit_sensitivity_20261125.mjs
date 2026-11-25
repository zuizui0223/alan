// Separate diagnostic; never changes Lu source no_pool status.
import {execFileSync} from 'node:child_process';
function csv(s){let rows=[],r=[],f="",q=false;for(let i=0;i<s.length;i++){let c=s[i];if(q){if(c==='"'){if(s[i+1]==='"'){f+='"';i++}else q=false}else f+=c}else if(c==='"'){q=true}else if(c===','){r.push(f);f=""}else if(c==='\n'){r.push(f.replace(/\r$/,""));rows.push(r);r=[];f=""}else f+=c}if(f||r.length){r.push(f.replace(/\r$/,""));rows.push(r)}let h=rows.shift();if(rows.some(r=>r.length!==h.length))throw Error("CSV shape");return rows.map(r=>Object.fromEntries(h.map((k,i)=>[k,r[i]])))}
const [g,l,gc,lc]=Object.values({"giavi_effects":"2089e1a59ef4802b94d10d9fe813e053f221093c","lu_effects":"747c00e0ce26ddd78350d33c6da83fef25ba70ae","giavi_cov":"62c02ecf2672a18309336947a641da99904f76db","lu_cov":"2ee8bfeb0e58fa7b32e2cda1064f0559c4497f33"}).map(s=>csv(execFileSync('git',['show',s],{encoding:'utf8'})));
const G=Object.fromEntries(g.map(x=>[x.effect_id,x]));
const L=Object.fromEntries(l.filter(x=>x.diel==='Diurnal'&&x.metric==='source_block_total').map(x=>[x.endpoint,x]));
if(!L.three_orders_DHC.status.includes('no_pool'))throw Error('gate changed');
const vG=(a,b)=>+gc.find(x=>x.effect_id_1===a&&x.effect_id_2===b).sampling_covariance;
const vL=(a,b)=>+lc.find(x=>x.endpoint===a)[b];
function pool(a,b){const tau2=Math.max(0,((a.y-b.y)**2-a.v-b.v)/2),w=[1/(a.v+tau2),1/(b.v+tau2)],m=(w[0]*a.y+w[1]*b.y)/(w[0]+w[1]),Q=w[0]*(a.y-m)**2+w[1]*(b.y-m)**2,se=Math.sqrt(Math.max(1,Q)/(w[0]+w[1])),t=12.706204736432095;return {tau2,rr:Math.exp(m),ci:[Math.exp(m-t*se),Math.exp(m+t*se)]}}
const a=G.Giavi2021_all_three_orders,b=L.three_orders_DHC,gd='Giavi2021_Diptera',gh='Giavi2021_Hymenoptera',ld='D_abundance',lh='H_abundance';
console.log(JSON.stringify({status:'arithmetic_sensitivity_only',total:pool({y:+a.lnRR,v:+a.sampling_variance},{y:+b.lnRR,v:+b.paired_delta_variance}),guild:pool({y:+G[gd].lnRR-+G[gh].lnRR,v:vG(gd,gd)+vG(gh,gh)-2*vG(gd,gh)},{y:+L[ld].lnRR-+L[lh].lnRR,v:vL(ld,ld)+vL(lh,lh)-2*vL(ld,lh)})},null,2))
