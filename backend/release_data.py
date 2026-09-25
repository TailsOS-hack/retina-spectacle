"""Immutable public dataset releases and compact ranged expression vectors."""
import os,json,gzip,hashlib,time,shutil
from pathlib import Path
import requests,numpy as np,h5py
try:from anndata.io import read_elem
except ImportError:from anndata._io.specs import read_elem
REPO=os.environ.get('GITHUB_REPOSITORY','TailsOS-hack/retina-spectacle')
TOKEN=os.environ.get('GH_TOKEN') or os.environ.get('GITHUB_TOKEN','')
def api(method,path,**kw):
 headers={'Authorization':'Bearer '+TOKEN,'Accept':'application/vnd.github+json','X-GitHub-Api-Version':'2022-11-28'};headers.update(kw.pop('headers',{}))
 for attempt in range(4):
  r=requests.request(method,'https://api.github.com/repos/'+REPO+'/'+path,headers=headers,timeout=120,**kw)
  if r.status_code not in [429,500,502,503] or attempt==3:break
  time.sleep(2**attempt)
 r.raise_for_status();return r.json() if r.content else None

def read_state(path,default=None):
 import base64
 try:return json.loads(base64.b64decode(api('GET','contents/'+path)['content']))
 except requests.HTTPError as e:
  if e.response.status_code==404:return default
  raise

def write_state(path,obj):
 import base64
 for _ in range(5):
  try:old=api('GET','contents/'+path);sha=old['sha']
  except requests.HTTPError as e:
   if e.response.status_code!=404:raise
   sha=None
  body={'message':'Update '+path,'content':base64.b64encode(json.dumps(obj,separators=(',',':')).encode()).decode()}
  if sha:body['sha']=sha
  try:return api('PUT','contents/'+path,json=body)
  except requests.HTTPError as e:
   if e.response.status_code not in [409,422]:raise
  time.sleep(1)
 raise RuntimeError('Concurrent state update failed')

def tag(acc,rev):return 'dataset-'+acc.lower()+'-'+rev

def browser_json(out,name,obj):
 """Shard large column-oriented data so every Vercel response fits its limit."""
 out=Path(out)
 def encoded(value):return gzip.compress(json.dumps(value,separators=(',',':')).encode(),compresslevel=6,mtime=0)
 data=encoded(obj)
 if len(data)<=3_000_000:
  (out/(name+'.json.gz')).write_bytes(data);return
 if name not in ('cells','genes'):raise ValueError('Metadata exceeds response limit')
 length=len(next(iter(obj.values())))
 if not all(isinstance(v,list) and len(v)==length for v in obj.values()):raise ValueError('Invalid column lengths')
 parts=[];start=0;chunk=20000 if name=='cells' else 5000
 while start<length:
  end=min(length,start+chunk);data=encoded({k:v[start:end] for k,v in obj.items()})
  if len(data)>3_000_000:
   if chunk==1:raise ValueError('One data row exceeds response limit')
   chunk=max(1,chunk//2);continue
  filename=f'{name}-{len(parts):03}.json';(out/(filename+'.gz')).write_bytes(data);parts.append(filename);start=end
 (out/(name+'.json.gz')).write_bytes(encoded({'parts':parts,'rows':length}))

def pack(directory,acc,rev,out):
 directory=Path(directory);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 meta=json.loads((directory/'meta.json').read_text());meta['revision']=rev;cells=json.loads((directory/'cells.json').read_text())
 with h5py.File(directory/'study.h5ad','r') as h:X=read_elem(h['X']).tocsc();names=read_elem(h['var']).index.astype(str).tolist()
 X.sort_indices();X.eliminate_zeros()
 try:genes=json.loads((directory/'genes.json').read_text());assert len(genes['genes'])==X.shape[1]
 except Exception:
  sums=[np.asarray(X.mean(axis=0)).ravel(),np.asarray((X>0).mean(axis=0)).ravel()*100]
  for i in range(len(meta['samples'])):sums.append(np.asarray(X[np.asarray(cells['sample'])==i].mean(axis=0)).ravel())
  genes={'genes':names,'summary':np.column_stack(sums).round(5).tolist()}
 for name,obj in [('meta',meta),('cells',cells),('genes',genes)]:browser_json(out,name,obj)
 index=[];part=0;size=0;f=open(out/'expression-000.bin','wb')
 try:
  for g in range(X.shape[1]):
   start,end=X.indptr[g:g+2];n=int(end-start)
   records=np.empty(n,dtype=[('i','<u4'),('v','<f4')]);records['i']=X.indices[start:end];records['v']=np.round(X.data[start:end],5)
   encoded=gzip.compress(records.tobytes(),compresslevel=3,mtime=0) if n else b''
   if size+len(encoded)>128*1024**2 and size:f.close();part+=1;size=0;f=open(out/f'expression-{part:03}.bin','wb')
   index.append([part,size,n,len(encoded)]);f.write(encoded);size+=len(encoded)
 finally:f.close()
 parts=[]
 with open(directory/'study.h5ad','rb') as src:
  k=0
  while True:
   data=src.read(512*1024**2)
   if not data:break
   name=f'study-{k:03}.h5part';(out/name).write_bytes(data);parts.append({'file':name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()});k+=1
 obj={'encoding':'gzip-per-gene','cells':X.shape[0],'genes':names,'index':index,'matrixParts':parts,'revision':rev}
 with gzip.open(out/'matrix-index.json.gz','wt',compresslevel=6) as f:json.dump(obj,f,separators=(',',':'))
 print(acc,'packed',X.shape,'nonzeros',X.nnz,flush=True)
 return meta

def publish(directory,release_tag,title):
 try:rel=api('GET','releases/tags/'+release_tag)
 except requests.HTTPError as e:
  if e.response.status_code!=404:raise
  rel=api('POST','releases',json={'tag_name':release_tag,'name':title,'body':'Public research data generated by Retina Spectacle. See the included study metadata for sources, methods and provenance.','make_latest':'false'})
 existing={a['name']:a for a in rel['assets']}
 for p in sorted(Path(directory).iterdir()):
  if not p.is_file():continue
  if p.name in existing and existing[p.name]['size']==p.stat().st_size:continue
  if p.name in existing:raise ValueError('Immutable release asset differs: '+p.name)
  url=rel['upload_url'].split('{')[0]
  for attempt in range(4):
   with open(p,'rb') as data:r=requests.post(url,params={'name':p.name},headers={'Authorization':'Bearer '+TOKEN,'Content-Type':'application/octet-stream'},data=data,timeout=600)
   if r.ok:break
   if attempt==3:r.raise_for_status()
   time.sleep(2**attempt)
 print(release_tag,'published',len(list(Path(directory).iterdir())),'assets',flush=True)
 return rel

def retrieve(acc,rev,directory):
 directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);base='https://github.com/'+REPO+'/releases/download/'+tag(acc,rev)+'/'
 r=requests.get(base+'matrix-index.json.gz',timeout=120);r.raise_for_status();index=json.loads(gzip.decompress(r.content))
 p=directory/'study.h5ad'
 with open(p,'wb') as dest:
  for part in index['matrixParts']:
   r=requests.get(base+part['file'],timeout=600);r.raise_for_status()
   if len(r.content)!=part['bytes'] or hashlib.sha256(r.content).hexdigest()!=part['sha256']:raise ValueError('Matrix integrity check failed')
   dest.write(r.content)
 return p
