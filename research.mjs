// Shared, independently testable cohort and expression calculations.
export const palette = ['#237968','#c77148','#527ca6','#9a6eab','#b3953e','#52776e','#bc6984','#75899b'];
export function humanText(value) {
  let text=String(value??'');
  if(/[ÃÂ]/.test(text)&&[...text].every(c=>c.charCodeAt(0)<256)) try {text=new TextDecoder('utf-8',{fatal:true}).decode(Uint8Array.from(text,c=>c.charCodeAt(0)));} catch {}
  return text;
}
export const esc=value=>humanText(value).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
export const fmt=n=>Number(n).toLocaleString('en-US');
export const cellType=name=>humanText(name).replace(/^C\d+\s*[·:—-]?\s*/,'').trim()||'Unassigned';
export const species=s=>humanText(s).replace('Homo sapiens','Human').replace('Mus musculus','Mouse').replace('Rattus norvegicus','Rat').replace('Macaca fascicularis','Macaque');
export function studyCategory(row) {
  if(['GSE175499','HRA000150','EGAS00001004561'].includes(row.id)) return 'reference';
  if(['GSE150703','GSE152928','GSE199792','PRJNA864092','GSE319256'].includes(row.id)) return 'related';
  if(/diabet|akimba/i.test(row.title||'')||['PRJCA006081','PRJNA653629'].includes(row.id)) return 'dr';
  if(/retin|ocular|macular|microglia/i.test(row.title||'')) return 'related';
  return 'reference';
}
export function sampleInfo(meta) {
  return meta.samples.map((id,index)=>{
    const m=meta.sampleMetadata?.[id]||{},fields={};
    for(const item of m.characteristics_ch1||[]) {const p=item.indexOf(':');if(p>0)fields[item.slice(0,p).trim().toLowerCase()]=humanText(item.slice(p+1).trim());}
    const title=humanText(m.title?.[0]||id), text=[title,...Object.values(fields)].join(' ').toLowerCase();
    let condition='Unassigned',evidence='No unambiguous condition found. Assign this sample in Custom groups.';
    const assign=(label,why)=>{condition=label;evidence=why;};
    const source=Object.entries(fields).filter(([k])=>/disease|condition|treatment|genotype/.test(k)).map(([k,v])=>`${k}: ${v}`).join('; ')||title;
    // Deposited archive filenames explicitly identify NDR1–3 and DR1–3.
    const pbmc={'GSM7910930':'NDR1','GSM7910931':'NDR2','GSM7910932':'NDR3','GSM7910933':'DR1','GSM7910934':'DR2','GSM7910935':'DR3'};
    if(meta.id==='GSE248284'&&pbmc[id]) assign(pbmc[id].startsWith('NDR')?'Diabetes without DR':'Diabetic retinopathy',`Deposited matrix filename: ${id}_${pbmc[id]}matrix.mtx.gz. Both cohorts have type 1 diabetes.`);
    else if(meta.id==='GSE245561') assign('PDR',`PDR fibrovascular membrane study; ${source}. No healthy control deposited.`);
    else if(/proliferative vitreoretinopathy|\bpvr\b/.test(text)) assign('PVR',source);
    else if(/proliferative diabetic|\bpdr\b/.test(text)) assign('PDR',source);
    else if(/oxygen.induced|\boir\b|(?:^|[ _])oir_/i.test(text)&&!(/\bnoir\b/.test(text))) assign('OIR',source);
    else if(/vldlr.?(?:ko|knockout)/.test(text)) assign('Vldlr knockout',source);
    else if(/sfts|severe fever/.test(text)) assign('SFTS',source);
    else if(/non.diabet|normoxia|normoxic|\bnoir\b|norm_|control|\bctrl\b|sham|wild.type|wide.type|\bwt\b|wt\/wt|leprwt|bkswt/.test(text)) assign(meta.id==='GSE150703'?'Normoxia':'Control',source);
    else if(/diabet|akimba|\bstz\b|db\/db|leptin receptor.deficient|leprdb|bksdb|\bdb-retina\b|\bdb-kidney\b/.test(text)) assign('Diabetes',source);
    let tissue=fields.tissue||humanText(m.source_name_ch1?.[0]||'');
    if(!tissue&&/retina|retinal/i.test(meta.title)) tissue='Retina';
    if(meta.id==='GSE204880') tissue=/kidney/i.test(title)?'Kidney':'Retina';
    return {index,id,title,condition,evidence,tissue:humanText(tissue),fields};
  });
}
export function typeMap(meta) {return new Map(meta.groups.map(g=>[g.id,cellType(g.name)]));}
export function filterCells(meta,cells,samples,filters,selectionA=new Set(),selectionB=new Set(),subset=null) {
  const types=typeMap(meta), result=[];
  for(let i=0;i<meta.cells;i++) {
    const sample=samples[cells.sample[i]];
    if(!sample||!filters.samples.has(sample.index)||!filters.types.has(types.get(cells.cluster[i]))||!filters.clusters.has(cells.cluster[i]))continue;
    if(filters.tissue!=='all'&&sample.tissue.toLowerCase()!==filters.tissue.toLowerCase())continue;
    if(subset&&!subset.has(i))continue;
    if(filters.scope==='a'&&!selectionA.has(i)||filters.scope==='b'&&!selectionB.has(i)||filters.scope==='selected'&&!selectionA.has(i)&&!selectionB.has(i))continue;
    result.push(i);
  }
  return result;
}
export function groupCells(meta,cells,samples,ids,mode,custom={},selectionA=new Set(),selectionB=new Set(),subset=null) {
  const types=typeMap(meta),groups=new Map(),names=new Map(meta.groups.map(g=>[String(g.id),humanText(g.name)]));
  if(mode==='selection') {groups.set('a',{key:'a',name:'Selected cells A',ids:[]});groups.set('b',{key:'b',name:'Selected cells B',ids:[]});}
  if(mode==='custom') {groups.set('a',{key:'a',name:custom.nameA||'Group A',ids:[]});groups.set('b',{key:'b',name:custom.nameB||'Group B',ids:[]});}
  for(const i of ids) {
    const s=samples[cells.sample[i]],t=types.get(cells.cluster[i]);let key,name;
    if(mode==='selection'){key=selectionA.has(i)?'a':selectionB.has(i)?'b':null;}
    else if(mode==='custom'){key=custom.assignment?.[s.index];if(!['a','b'].includes(key))key=null;}
    else if(mode==='subset'){const at=subset?.lookup.get(i);if(at===undefined)continue;key=String(subset.cluster[at]);name='Reclustered C'+key;}
    else if(mode==='sample'){key=String(s.index);name=s.title===s.id?s.id:`${s.title} · ${s.id}`;}
    else if(mode==='celltype'){key=t;name=t;}
    else if(mode==='cluster'){key=String(cells.cluster[i]);name=names.get(key);}
    else if(mode==='type-condition'){key=`${t} / ${s.condition}`;name=key;}
    else {key=s.condition;name=s.condition;}
    if(key==null)continue;
    if(!groups.has(key))groups.set(key,{key,name,ids:[]});groups.get(key).ids.push(i);
  }
  return [...groups.values()].map(g=>({...g,samples:[...new Set(g.ids.map(i=>cells.sample[i]))]})).sort((a,b)=>{
    const order=['Control','Normoxia','Diabetes without DR','Diabetes','Diabetic retinopathy','PDR','PVR','OIR','Unassigned'];
    if(mode==='condition')return (order.indexOf(a.name)<0?50:order.indexOf(a.name))-(order.indexOf(b.name)<0?50:order.indexOf(b.name));
    return 0;
  });
}
export function quantile(sorted,p) {if(!sorted.length)return null;const index=(sorted.length-1)*p,lo=Math.floor(index),hi=Math.ceil(index);return sorted[lo]+(sorted[hi]-sorted[lo])*(index-lo);}
export function summarize(ids,values) {
  const sorted=ids.map(i=>values[i]).sort((a,b)=>a-b),n=sorted.length;
  if(!n)return {n:0,mean:null,pct:null,median:null,q1:null,q3:null,low:null,high:null,min:null,max:null,linearMean:null};
  let sum=0,linearSum=0,detected=0;
  for(const v of sorted){sum+=v;linearSum+=Math.expm1(v);if(v>0)detected++;}
  const q1=quantile(sorted,.25),q3=quantile(sorted,.75),iqr=q3-q1;
  let l=0,h=n-1;while(sorted[l]<q1-1.5*iqr)l++;while(sorted[h]>q3+1.5*iqr)h--;
  return {n,mean:sum/n,pct:detected/n*100,median:quantile(sorted,.5),q1,q3,low:sorted[l],high:sorted[h],min:sorted[0],max:sorted[n-1],linearMean:linearSum/n};
}
function erfc(x) {
  const z=Math.abs(x),t=1/(1+z/2),p=t*Math.exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*(1.48851587+t*(-.82215223+t*.17087277)))))))));
  return x>=0?p:2-p;
}
export function mannWhitney(a,b,values) {
  if(a.length<3||b.length<3)return null;
  const bset=new Set(b);if(a.some(i=>bset.has(i)))throw Error('Comparison populations overlap.');
  const sorted=[...a.map(i=>[values[i],1]),...b.map(i=>[values[i],0])].sort((x,y)=>x[0]-y[0]);
  let rankA=0,ties=0;
  for(let i=0;i<sorted.length;){let j=i+1,countA=sorted[i][1];while(j<sorted.length&&sorted[j][0]===sorted[i][0])countA+=sorted[j++][1];const t=j-i;rankA+=countA*(i+1+j)/2;ties+=t*t*t-t;i=j;}
  const n=a.length,m=b.length,N=n+m,u=rankA-n*(n+1)/2,variance=n*m/12*(N+1-ties/(N*(N-1)));
  if(variance<=0)return 1;
  const z=Math.max(0,Math.abs(u-n*m/2)-.5)/Math.sqrt(variance);
  return Math.min(1,Math.max(0,erfc(z/Math.SQRT2)));
}
export function adjustBH(ps) {const valid=ps.map((p,i)=>({p,i})).filter(x=>x.p!==null).sort((a,b)=>a.p-b.p),out=ps.map(()=>null);let last=1;for(let k=valid.length-1;k>=0;k--){last=Math.min(last,valid[k].p*valid.length/(k+1));out[valid[k].i]=last;}return out;}
export function compareGenes(geneNames,vectors,a,b) {
  const rows=vectors.map((v,j)=>{const A=summarize(a,v),B=summarize(b,v);return {gene:geneNames[j],nA:A.n,nB:B.n,meanA:A.mean,meanB:B.mean,pctA:A.pct,pctB:B.pct,log2FC:A.n&&B.n?Math.log2((A.linearMean+1)/(B.linearMean+1)):null,deltaPct:A.n&&B.n?A.pct-B.pct:null,pvalue:mannWhitney(a,b,v)};});
  const qs=adjustBH(rows.map(r=>r.pvalue));return rows.map((r,i)=>({...r,qvalue:qs[i]}));
}
export function insidePolygon(x,y,poly) {let inside=false;for(let i=0,j=poly.length-1;i<poly.length;j=i++){const[xi,yi]=poly[i],[xj,yj]=poly[j];if((yi>y)!==(yj>y)&&x<(xj-xi)*(y-yi)/(yj-yi)+xi)inside=!inside;}return inside;}
