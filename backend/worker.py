import os,sys,json,time,hashlib,gzip,base64,traceback,tempfile,shutil,re
from pathlib import Path
import numpy as np
from release_data import read_state,write_state,pack,publish,tag,retrieve,api
from discovery import review_accessions,latest_geo_candidates,candidate_is_single_cell,geo_meta,arrayexpress_meta,source_fingerprint,bioproject_geo,download,expand,read_matrices,Unavailable
from analysis import process,differential,embedding
ROOT=Path(tempfile.gettempdir())/'retina-spectacle';ROOT.mkdir(exist_ok=True)
def refresh():
 c=read_state('state/catalog.json',{'review':{},'studies':[],'catalog':{},'freshness':{}})
 c['freshness']={**c.get('freshness',{}),'status':'checking','startedAt':time.time()};write_state('state/catalog.json',c)
 errors=[];added=[];changed=[];checked=[];retry=os.environ.get('RETRY_SOURCE','')
 try:
  ids,sha=review_accessions();c['review']={'id':'PMC11214886','accessions':ids,'sha256':sha,'discoveredAt':time.time()}
  known=c.get('catalog',{}).get('accessions',[]);screened=c.get('catalog',{}).get('screened',{});mc={}
  for acc in latest_geo_candidates():
   if acc in ids or acc in known:continue
   try:
    m=geo_meta(acc);mc[acc]=m
    if candidate_is_single_cell(m):known.append(acc);added.append(acc);screened[acc]='Single-cell retinal diabetes record'
    else:screened[acc]='Excluded: not established as retinal-diabetes single-cell data'
   except Exception as e:errors.append({'study':acc,'error':str(e)[:200]})
   time.sleep(.5)
  c['catalog']={'accessions':known,'screened':screened,'queryCheckedAt':time.time()}
  records=json.loads(Path(__file__).with_name('access_records.json').read_text())
  targets=list(dict.fromkeys(ids+known));completed=0
  for acc in targets:
   row=next((x for x in c['studies'] if x['id']==acc),None)
   if row is None:row={'id':acc,'status':'discovering','discoveryOrigin':'review' if acc in ids else 'live GEO search'};c['studies'].append(row)
   if acc!=retry and row.get('status') in ['controlled_access','raw_reads','unavailable','failed'] and time.time()-row.get('lastAttemptAt',0)<86400:continue
   if acc in records:row.update(records[acc],lastAttemptAt=time.time());continue
   if acc.startswith('PRJNA'):
    linked=bioproject_geo(acc);row.update(status='linked' if linked else 'raw_reads',linked=linked,message='Linked GEO sources: '+', '.join(linked) if linked else 'Public raw reads require alignment; no processed matrix discovered.',lastAttemptAt=time.time())
    for child in linked:
     if child not in targets:targets.append(child)
    continue
   if not re.fullmatch(r'GSE\d+|E-MTAB-\d+',acc):row.update(status='access_review',message='Repository access review required');continue
   previous=dict(row)
   try:
    m=mc.get(acc) or (geo_meta(acc) if acc.startswith('GSE') else arrayexpress_meta(acc));fp,headers=source_fingerprint(m);checked.append(acc);row.update(lastCheckedAt=time.time(),sourceHeaders=headers)
    if row.get('status')=='ready' and row.get('sourceFingerprint')==fp:continue
    if completed>=3:row.update(refreshStatus='queued',refreshMessage='Queued for the next source check');continue
    row.update(title=m['title'],species=m['species'],source=m['source'],lastAttemptAt=time.time())
    if previous.get('status')=='ready':row.update(refreshStatus='running',refreshMessage='Importing changed repository data')
    else:row.update(status='running',message='Downloading public processed matrices')
    write_state('state/catalog.json',c)
    if not m['files']:raise Unavailable('No public processed matrix is listed')
    rev=hashlib.sha256((fp+str(time.time())).encode()).hexdigest()[:16];d=ROOT/acc/rev;raw=d/'raw';raw.mkdir(parents=True,exist_ok=True);files=[];m['liveDiscovery']=row.get('discoveryOrigin')=='live GEO search'
    for url in m['files']:
     if not re.search(r'\.(tar|gz|zip|h5|h5ad|csv|tsv|txt|mtx|mtx_rows|mtx_cols)$',url,re.I):continue
     name=url.rsplit('/',1)[-1]
     if not re.fullmatch(r'[A-Za-z0-9_.%+\-]+',name):raise Unavailable('Unsupported source filename')
     p=raw/name;m['provenance'].append(download(url,p,print));files.extend(expand(p,raw/(name+'.expanded')))
    a,norm=read_matrices(files,m,print);process(a,m,d,print,source_normalized=norm);del a
    stage=ROOT/('packed-'+acc);shutil.rmtree(stage,ignore_errors=True);stage.mkdir();meta=pack(d,acc,rev,stage);publish(stage,tag(acc,rev),acc+' · '+rev)
    row.update(status='ready',message='Ready — partial import' if meta.get('partialImport') else 'Ready for analysis',cells=meta['cells'],genes=meta['genes'],activeVersion=rev,publishedVersions=list(dict.fromkeys(previous.get('publishedVersions',['legacy'] if previous.get('status')=='ready' else [])+[rev])),sourceFingerprint=fp,partialImport=meta.get('partialImport',False),refreshStatus='idle',refreshMessage='',completedAt=time.time())
    changed.append(acc);completed+=1;write_state('state/catalog.json',c);shutil.rmtree(d);shutil.rmtree(stage)
   except Exception as e:
    traceback.print_exc();errors.append({'study':acc,'error':str(e)[:200]})
    if previous.get('status')=='ready':row.update(previous);row.update(refreshStatus='failed',refreshMessage=str(e)[:400],lastAttemptAt=time.time())
    else:row.update(status='unavailable' if isinstance(e,Unavailable) else 'failed',message=str(e)[:400],lastAttemptAt=time.time())
   time.sleep(.5)
  c['freshness'].update(status='ready',message='Source check complete',updatedAt=time.time(),lastCheckedAt=time.time(),newStudies=added,updatedStudies=changed,checkedStudies=checked,errors=errors)
 except Exception as e:c['freshness'].update(status='failed',message=str(e)[:400],errors=errors);traceback.print_exc()
 write_state('state/catalog.json',c)

def write_workbook(result,path):
 from openpyxl import Workbook
 wb=Workbook();ws=wb.active;ws.title='Exploratory markers';rows=result['rows']
 fields=list(rows[0]) if rows else ['gene','log2FC','deltaPct','pctA','pctB','meanA','meanB','pvalue','qvalue']
 ws.append(fields)
 for row in rows:ws.append([row[key] for key in fields])
 methods=wb.create_sheet('Methods');methods.append(['Method',result['method']]);methods.append(['Caution',result['warning']]);wb.save(path)

def run_analysis(job_id):
 import scanpy as sc
 if not re.fullmatch('[a-f0-9]{24}',job_id):raise ValueError('Invalid job')
 previous=read_state('state/jobs/'+job_id+'.json');
 if previous and previous.get('status')=='ready':return
 record=read_state('control/jobs/'+job_id+'.json');payload=json.loads(gzip.decompress(base64.b64decode(record.pop('payload'))));d=ROOT/job_id;d.mkdir(exist_ok=True)
 def progress(message):record.update(status='running',message=message,updatedAt=time.time());write_state('state/jobs/'+job_id+'.json',record)
 try:
  # A runner may stop after publishing the result but before updating its status.
  # Reuse those immutable results, including the original workbook ZIP bytes.
  import requests
  try:existing=api('GET','releases/tags/job-'+job_id)
  except requests.HTTPError as e:
   if e.response.status_code!=404:raise
   existing=None
  assets={item['name']:item for item in (existing or {}).get('assets',[]) if item['state']=='uploaded'}
  if 'result.json.gz' in assets:
   response=requests.get(assets['result.json.gz']['browser_download_url'],timeout=120);response.raise_for_status();result=json.loads(gzip.decompress(response.content))
   if record['kind']=='differential' and 'result.xlsx' not in assets:
    out=d/'result';out.mkdir(exist_ok=True);write_workbook(result,out/'result.xlsx');publish(out,'job-'+job_id,'Analysis '+job_id)
   record.update(status='ready',message='Analysis complete',updatedAt=time.time());write_state('state/jobs/'+job_id+'.json',record);return
  progress('Retrieving the versioned expression matrix');path=retrieve(record['study'],record['revision'],d);a=sc.read_h5ad(path)
  if record['kind']=='differential':result=differential(a,payload['a'],payload['b'],payload.get('minPct',.1),payload.get('logfc',.25))
  else:
   ix=np.unique(np.asarray(payload['cells'],dtype=int))
   if len(ix)<12 or len(ix)>50000 or ix.min()<0 or ix.max()>=a.n_obs:raise ValueError('Select 12–50,000 valid cells')
   b=a[ix].copy();del a
   if 'counts' in b.layers:b.X=b.layers['counts'].copy();sc.pp.normalize_total(b,target_sum=10000);sc.pp.log1p(b)
   sc.pp.filter_genes(b,min_cells=3);b=embedding(b,progress,tsne=True)
   result={'indices':ix.tolist(),'x':np.round(b.obsm['X_umap'][:,0],4).tolist(),'y':np.round(b.obsm['X_umap'][:,1],4).tolist(),'tx':np.round(b.obsm['X_tsne'][:,0],4).tolist(),'ty':np.round(b.obsm['X_tsne'][:,1],4).tolist(),'cluster':b.obs.cluster.astype(str).tolist(),'method':'Subset variable genes, scaled PCA, neighbors, UMAP, t-SNE and Leiden; seed 42. Expression overlays retain the parent matrix.'}
  out=d/'result';out.mkdir(exist_ok=True)
  with gzip.open(out/'result.json.gz','wt') as f:json.dump(result,f,separators=(',',':'))
  if record['kind']=='differential':
   write_workbook(result,out/'result.xlsx')
  publish(out,'job-'+job_id,'Analysis '+job_id);record.update(status='ready',message='Analysis complete',updatedAt=time.time())
 except Exception as e:traceback.print_exc();record.update(status='failed',message=str(e)[:400],updatedAt=time.time())
 write_state('state/jobs/'+job_id+'.json',record)

if __name__=='__main__':
 if os.environ.get('JOB_KIND','refresh')=='analysis':run_analysis(os.environ['JOB_ID'])
 else:refresh()
