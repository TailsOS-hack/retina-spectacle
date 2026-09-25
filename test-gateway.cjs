const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict'),zlib=require('node:zlib');
const source=fs.readFileSync('api/gateway.js','utf8');
function setup({failDispatch=false}={}){
 const files=new Map(),calls=[];let seq=0;
 const c={studies:[{id:'GSE245561',status:'ready',cells:8960,activeVersion:'legacy',publishedVersions:['legacy','new']},{id:'GSE175499',status:'ready',cells:150504,activeVersion:'legacy',publishedVersions:['legacy']}],freshness:{lastCheckedAt:0}};
 files.set('state/catalog.json',{value:c,sha:'first'});
 const fetch=async(url,options={})=>{
  calls.push([url,options.method||'GET']);
  if(url.includes('/contents/')){
   const path=url.split('/contents/')[1],f=files.get(path),method=options.method||'GET';
   if(method==='GET')return f?new Response(JSON.stringify({content:Buffer.from(JSON.stringify(f.value)).toString('base64'),sha:f.sha}),{status:200}):new Response('',{status:404});
   const body=JSON.parse(options.body);
   if((f&&body.sha!==f.sha)||(!f&&body.sha))return new Response('',{status:409});
   const sha=String(++seq);files.set(path,{sha,value:JSON.parse(Buffer.from(body.content,'base64'))});return new Response(JSON.stringify({content:{sha}}),{status:201});
  }
  if(url.includes('/dispatches'))return new Response(null,{status:failDispatch?503:204});
  if(url.includes('/releases/download/')){
   const name=url.split('/').pop(),base='../free-releases/'+(url.includes('gse175499')?'GSE175499':'GSE245561')+'/';
   let b=fs.readFileSync(base+name);const range=options.headers?.Range;
   if(range){const[a,end]=range.slice(6).split('-').map(Number);b=b.subarray(a,end+1)}
   return new Response(b,{status:range?206:200});
  }
  throw Error('Unexpected service URL');
 };
 const sandbox={require,Buffer,URL,console,fetch,process:{env:{GH_TOKEN:'test-repository-key',ADMIN_KEY:'test-access-code'}},module:{exports:{}},setTimeout};vm.runInNewContext(source,sandbox);
 const request=async(method,url,body,owner=false)=>{let statusCode=200,data,headers={};const res={set statusCode(n){statusCode=n},setHeader(k,v){headers[k]=v},end(b){data=b}};await sandbox.module.exports({method,url,body,headers:owner?{'x-admin-key':'test-access-code'}:{}},res);return {status:statusCode,headers,value:JSON.parse(zlib.gunzipSync(data))}};
 return {request,files,calls};
}
(async()=>{
 const a=setup(),checks=[];
 const ok=(name,test)=>{assert(test,name);checks.push(name)};
 let r=await a.request('GET','/api/health');ok('Health response',r.status===200&&r.value.version==='3.0-free');
 r=await a.request('GET','/api/data/GSE245561/meta.json?revision=unknown');ok('Unpublished version rejected',r.status===404);
 r=await a.request('GET','/api/data/not-a-study/meta.json');ok('Invalid accession rejected',r.status===400);
 r=await a.request('GET','/api/data/GSE245561/gene/-1?revision=legacy');ok('Negative gene rejected',r.status===404);
 r=await a.request('GET','/api/data/GSE245561/gene/100?revision=legacy');ok('Real ranged vector decoded',r.status===200&&r.value.indices.length===r.value.values.length&&r.value.cells===8960);ok('Immutable version caching',r.headers['Cache-Control'].includes('31536000'));
 r=await a.request('GET','/api/data/GSE245561/gene/100');ok('Unpinned vectors expire promptly',r.headers['Cache-Control'].includes('max-age=30'));
 r=await a.request('GET','/api/data/GSE175499/gene/100?revision=legacy');ok('Compressed large-study vector decoded',r.status===200&&r.value.cells===150504&&r.value.indices.length===r.value.values.length);
 r=await a.request('GET','/api/data/GSE175499/cells.json?revision=legacy');ok('Large study uses bounded shards',r.status===200&&r.value.parts.length>1&&r.value.rows===150504);
 r=await a.request('POST','/api/analysis/GSE245561/differential',{a:[1,2,3],b:[4,5,6],minPct:.1,logfc:.25});ok('Analysis access protected',r.status===401);
 r=await a.request('POST','/api/analysis/GSE245561/differential',{a:[1,2,3],b:[3,4,5],minPct:.1,logfc:.25},true);ok('Overlapping groups rejected',r.status===400);
 r=await a.request('POST','/api/analysis/GSE245561/differential',{a:[1,2,3],b:[4,5,99999],minPct:.1,logfc:.25},true);ok('Out-of-range cell rejected',r.status===400);
 r=await a.request('POST','/api/analysis/GSE245561/recluster',{cells:[1,2,3]},true);ok('Tiny reclustering rejected',r.status===400);
 r=await a.request('POST','/api/analysis/GSE245561/differential',{a:[1,2,3],b:[4,5,6],minPct:2,logfc:.25},true);ok('Invalid thresholds rejected',r.status===400);
 const payload={a:[1,2,3],b:[4,5,6],minPct:.1,logfc:.25,revision:'legacy'};
 const jobs=await Promise.all(Array.from({length:8},()=>a.request('POST','/api/analysis/GSE245561/differential',{...payload},true)));ok('Concurrent identical analyses share one job',new Set(jobs.map(x=>x.value.id)).size===1&&a.calls.filter(x=>x[0].includes('/dispatches')).length===1);ok('Public job excludes payload and key',jobs.every(x=>!x.value.payload)&&!JSON.stringify(jobs).includes('test-access-code'));
 r=await a.request('GET','/api/jobs/'+jobs[0].value.id);ok('Queued job is recoverable by id',r.value.status==='queued');
 const b=setup(),refreshes=await Promise.all(Array.from({length:12},()=>b.request('POST','/api/refresh')));ok('Concurrent refreshes dispatch once',b.calls.filter(x=>x[0].includes('/dispatches')).length===1);ok('All refresh callers receive success',refreshes.every(x=>x.status===200));
 const f=setup({failDispatch:true});r=await f.request('POST','/api/refresh');ok('Failed dispatch is reported',r.status===503);ok('Failed dispatch releases the retry gate',f.files.get('control/refresh-lock.json').value.status==='failed');
 fs.writeFileSync('../free-gateway-test-report.json',JSON.stringify({checks,passed:checks.length},null,2));console.log(checks.length+' gateway checks passed');
})().catch(e=>{console.error(e);process.exitCode=1});
