"""Repository discovery, safe downloads, and processed-matrix import adapters."""
import re,json,hashlib,gzip,csv,tarfile,zipfile,time
from pathlib import Path
from urllib.parse import urlparse,urljoin
import requests
import numpy as np
import scipy.sparse as sp
from scipy.io import mmread

REVIEW='PMC11214886'
ALLOWED={'www.ncbi.nlm.nih.gov','ftp.ncbi.nlm.nih.gov','www.ebi.ac.uk','ftp.ebi.ac.uk','www.ebi.ac.uk','www.ebi.ac.uk','ngdc.cncb.ac.cn','ega-archive.org','www.ena.ebi.ac.uk'}
class Unavailable(Exception):pass

def request(url,**kwargs):
 if url.startswith('ftp://'):url='https://'+url[6:]
 for _ in range(6):
  p=urlparse(url)
  if p.scheme!='https' or p.hostname not in ALLOWED or p.username or p.password:raise ValueError('Repository URL is not permitted')
  r=requests.get(url,timeout=(15,90),allow_redirects=False,**kwargs)
  if r.is_redirect:url=urljoin(url,r.headers['Location']);r.close();continue
  r.raise_for_status();return r
 raise ValueError('Too many redirects')
def get_text(url):return request(url).text

def review_accessions():
 import xml.etree.ElementTree as ET
 xml=get_text(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{REVIEW}/fullTextXML');root=ET.fromstring(xml)
 table=root.find('.//table-wrap');text=' '.join(table.itertext())
 for h in ['‐','‑','–','−']:text=text.replace(h,'-')
 ids=list(dict.fromkeys(re.findall(r'GSE\d+|PRJNA\d+|PRJCA\d+|E-MTAB-\d+|HRA\d+|EGAS\d+',text)))
 return ids,hashlib.sha256(xml.encode()).hexdigest()

def geo_meta(acc):
 if not re.fullmatch(r'GSE\d+',acc):raise ValueError('Invalid GEO accession')
 text=get_text(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full')
 vals=lambda k:[s.strip() for s in re.findall(r'^!Series_'+k+r' = (.*)',text,re.M)]
 if not vals('title'):raise Unavailable('No public GEO record found')
 st=get_text(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=full')
 samples={}
 for part in st.split('^SAMPLE = ')[1:]:
  sid=part.splitlines()[0].strip();sample={}
  for key in ['title','source_name_ch1','characteristics_ch1','data_processing']:
   sample[key]=[s.strip() for s in re.findall(r'^!Sample_'+key+r' = (.*)',part,re.M)]
  samples[sid]=sample
 files=[u.replace('ftp://','https://') for u in vals('supplementary_file') if u.lower()!='none']
 return dict(id=acc,title=vals('title')[0],species=', '.join(dict.fromkeys(vals('sample_organism'))),summary=' '.join(vals('summary')),design=' '.join(vals('overall_design')),samples=samples,files=files,source=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}',provenance=[])

def bioproject_geo(acc):
 s=request('https://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',params={'db':'gds','term':acc,'retmode':'json','retmax':100}).json()
 ids=s.get('esearchresult',{}).get('idlist',[])
 if not ids:return []
 data=request('https://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',params={'db':'gds','id':','.join(ids),'retmode':'json'}).json()
 return list(dict.fromkeys(x.get('accession') for x in data.get('result',{}).values() if isinstance(x,dict) and str(x.get('accession','')).startswith('GSE')))

def arrayexpress_meta(acc):
 obj=request(f'https://www.ebi.ac.uk/biostudies/api/v1/studies/{acc}').json();files=[]
 def walk(x):
  if isinstance(x,dict):
   if 'path' in x and isinstance(x['path'],str):files.append(x['path'])
   for v in x.values():walk(v)
  elif isinstance(x,list):
   for v in x:walk(v)
 walk(obj)
 urls=[f'https://www.ebi.ac.uk/biostudies/files/{acc}/{p}' for p in dict.fromkeys(files) if re.search(r'\.(mtx|tsv|csv|txt|h5|h5ad|tar|zip)(\.gz)?$',p,re.I) and not p.endswith(('idf.txt','sdrf.txt'))]
 attrs={a['name']:a.get('value','') for a in obj.get('section',{}).get('attributes',[])}
 if not urls and 'gxa-sc' in json.dumps(obj):
  base=f'https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/{acc}/'
  listing=get_text(base)
  needed=[acc+'.aggregated_counts.mtx.gz',acc+'.aggregated_counts.decorated.mtx_rows',acc+'.aggregated_counts.mtx_cols',acc+'.cell_metadata.tsv']
  if all('href="'+n+'"' in listing for n in needed):urls=[base+n for n in needed]
 return dict(id=acc,title=attrs.get('Title','ArrayExpress '+acc),species=attrs.get('Organism','See source metadata'),summary=attrs.get('Description',''),design=attrs.get('Description',''),samples={},files=urls,source=f'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{acc}',provenance=[])

def download(url,path,progress):
 if path.exists():return {'url':url,'file':path.name,'bytes':path.stat().st_size,'sha256':hashlib.file_digest(open(path,'rb'),'sha256').hexdigest() if hasattr(hashlib,'file_digest') else hashlib.sha256(path.read_bytes()).hexdigest()}
 r=request(url,stream=True);total=int(r.headers.get('Content-Length',0));limit=3*1024**3
 if total>limit:raise Unavailable('Source file exceeds the current 3 GB per-file import limit')
 h=hashlib.sha256();n=0;t=time.time();temp=path.with_suffix(path.suffix+'.part')
 with open(temp,'wb') as f:
  for chunk in r.iter_content(1024**2):
   if not chunk:continue
   n+=len(chunk)
   if n>limit:raise Unavailable('Source exceeds the 3 GB per-file import limit')
   f.write(chunk);h.update(chunk)
   if time.time()-t>8:progress(f'Downloading {path.name}: {n//1024**2} MB'+(f' / {total//1024**2} MB' if total else ''));t=time.time()
 temp.rename(path);return dict(url=url,file=path.name,bytes=n,sha256=h.hexdigest())

def expand(path,dest):
 dest.mkdir(exist_ok=True,parents=True)
 def target(name):
  p=(dest/name).resolve()
  if not p.is_relative_to(dest.resolve()):raise ValueError('Unsafe archive path')
  return p
 if tarfile.is_tarfile(path):
  with tarfile.open(path) as t:
   total=0
   for item in t:
    if not item.isfile():continue
    total+=item.size
    if total>15*1024**3:raise Unavailable('Expanded archive exceeds 15 GB')
    p=target(item.name);p.parent.mkdir(exist_ok=True,parents=True)
    with t.extractfile(item) as inp,open(p,'wb') as out:
     import shutil;shutil.copyfileobj(inp,out)
 elif zipfile.is_zipfile(path):
  with zipfile.ZipFile(path) as z:
   if sum(i.file_size for i in z.infolist())>15*1024**3:raise Unavailable('Expanded archive exceeds 15 GB')
   for i in z.infolist():
    if i.is_dir():continue
    p=target(i.filename);p.parent.mkdir(exist_ok=True,parents=True);p.write_bytes(z.read(i))
 else:return [path]
 result=[]
 for p in list(dest.rglob('*')):
  if not p.is_file():continue
  if re.search(r'\.(tar|tar.gz|tgz|zip)$',p.name):result.extend(expand(p,p.parent/(p.name+'.expanded')))
  else:result.append(p)
 return list(dict.fromkeys(result))

def open_text(path):return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)
def dense_table(path):
 vals=[];inds=[];ptr=[0];genes=[];bad=False
 with open_text(path) as f:
  line=f.readline()
  while line.startswith('#'):line=f.readline()
  sep='\t' if '\t' in line else ','
  header=next(csv.reader([line],delimiter=sep));bars=header[1:]
  for line in f:
   if not line.strip():continue
   first,rest=line.split(sep,1);a=np.fromstring(rest,sep=sep,dtype=np.float32)
   # Some R tables omit the empty top-left header cell.
   if len(a)==len(header) and not genes:bars=header
   if len(a)!=len(bars):raise Unavailable('Ambiguous expression-table orientation or nonnumeric matrix')
   if not np.isfinite(a).all() or a.min()<0:raise Unavailable('Matrix contains invalid or negative expression values')
   if np.any(np.abs(a-np.round(a))>1e-4):bad=True
   ix=np.flatnonzero(a);genes.append(first.strip('"'));vals.extend(a[ix]);inds.extend(ix);ptr.append(len(vals))
 if not genes:raise Unavailable('Empty matrix')
 X=sp.csr_matrix((np.asarray(vals,dtype=np.float32),np.asarray(inds),np.asarray(ptr)),shape=(len(genes),len(bars))).T.tocsr()
 return X,genes,bars,bad

def read_matrices(files,meta,progress):
 import anndata as ad,scanpy as sc,pandas as pd
 files=list(dict.fromkeys(files));mats=[];used=set()
 atlas=next((p for p in files if p.name.endswith('.aggregated_counts.mtx.gz')),None)
 if atlas:
  root=atlas.parent;prefix=meta['id'];progress('Reading Expression Atlas raw count matrix')
  with gzip.open(atlas,'rb') as f:X=mmread(f).T.tocsr().astype(np.float32)
  rows=pd.read_csv(root/(prefix+'.aggregated_counts.decorated.mtx_rows'),sep='\t',header=None,keep_default_na=False)
  bars=(root/(prefix+'.aggregated_counts.mtx_cols')).read_text().splitlines()
  cm=pd.read_csv(root/(prefix+'.cell_metadata.tsv'),sep='\t').set_index('id')
  samplemap={}
  for barcode,row in cm.iterrows():samplemap[barcode.split('-')[0]]=str(row.get('disease',''))+' · '+str(row.get('genotype',''))
  samples=[b.split('-')[0]+' · '+samplemap.get(b.split('-')[0],'unmapped') for b in bars]
  g=[r[1] or r[0] for r in rows.values.tolist()]
  a=ad.AnnData(X,obs=pd.DataFrame({'sample':samples,'source_barcode':bars},index=bars),var=pd.DataFrame({'gene_id':rows[0].tolist()},index=g));a.var_names_make_unique()
  meta['samples']={s:{'title':[s],'characteristics_ch1':['Metadata from Expression Atlas cell_metadata.tsv']} for s in set(samples)}
  return a,False
 for p in sorted(files):
  if not re.search(r'matrix\.mtx(\.gz)?$',p.name,re.I):continue
  prefix=re.split(r'matrix\.mtx',p.name,flags=re.I)[0]
  options=[x for x in files if x.parent==p.parent and x.name.startswith(prefix)]
  bf=next((x for x in options if re.search(r'barcodes.*tsv',x.name,re.I)),None)
  gf=next((x for x in options if re.search(r'(features|genes).*tsv',x.name,re.I)),None)
  if bf is None or gf is None:raise Unavailable('MatrixMarket file lacks matching barcode or gene files')
  progress('Reading '+p.name)
  try:
   with (gzip.open(p,'rb') if p.name.endswith('.gz') else open(p,'rb')) as f:X=mmread(f).T.tocsr().astype(np.float32)
  except (EOFError,gzip.BadGzipFile) as e:
   meta.setdefault('excludedMatrices',[]).append({'file':p.name,'sample':(re.search(r'GSM\d+',p.name).group() if re.search(r'GSM\d+',p.name) else p.name),'reason':'Deposited compressed matrix is truncated or corrupt: '+str(e)})
   meta['partialImport']=True;used.update([p,bf,gf]);progress('Omitting corrupt source matrix: '+p.name);continue
  with open_text(gf) as f:rows=list(csv.reader(f,delimiter='\t'))
  genes=[r[1] if len(r)>1 else r[0] for r in rows]
  with open_text(bf) as f:bars=[s.strip() for s in f]
  gene_ids=[r[0] for r in rows]
  if len(rows[0])>=3:
   sel=np.array([r[2]=='Gene Expression' for r in rows]);X=X[:,sel];genes=list(np.array(genes)[sel]);gene_ids=list(np.array(gene_ids)[sel])
  mats.append((p,X,genes,bars,False,gene_ids));used.update([p,bf,gf])
 for p in sorted(files):
  if p in used:continue
  if p.name.endswith(('.h5','.h5ad')):
   progress('Reading '+p.name)
   a=sc.read_h5ad(p) if p.name.endswith('.h5ad') else sc.read_10x_h5(p)
   X=sp.csr_matrix(a.layers['counts'] if 'counts' in a.layers else a.X,dtype=np.float32)
   normalized=bool(np.any(np.abs(X.data-np.round(X.data))>1e-4))
   mats.append((p,X,a.var_names.tolist(),a.obs_names.tolist(),normalized,a.var.get('gene_ids',a.var_names).tolist()))
  elif re.search(r'\.(csv|tsv|txt)(\.gz)?$',p.name,re.I) and re.search(r'matrix|count|dge',p.name,re.I):
   progress('Reading '+p.name);X,g,b,norm=dense_table(p);mats.append((p,X,g,b,norm,g))
 if not mats:raise Unavailable('No supported processed expression matrices found. Raw sequencing reads require a separate alignment pipeline.')
 ads=[];normalization=[]
 for p,X,g,b,norm,gids in mats:
  sid=re.search(r'GSM\d+',str(p));sid=sid.group() if sid else p.stem
  if meta.get('liveDiscovery'):
   sample=meta.get('samples',{}).get(sid,{})
   text=' '.join(' '.join(v) if isinstance(v,list) else str(v) for v in sample.values())
   if sample and re.search(r'bulk.?RNA|bulk sequencing',text,re.I) and not re.search(r'single.cell|single.nucle|cell.?ranger|10x',text,re.I):continue
   if re.search(r'\.(csv|tsv|txt)(\.gz)?$',p.name,re.I):
    if len(b)<200 or sum(bool(re.search(r'[ACGT]{12,}',str(x))) for x in b)<len(b)*.8:
     raise Unavailable('Newly discovered dense matrix needs a verified cell-barcode adapter before automatic import; bulk samples will not be treated as cells')
  if X.shape!=(len(b),len(g)):raise ValueError('Matrix dimensions disagree with labels')
  a=ad.AnnData(X,obs=pd.DataFrame(index=pd.Index(b)),var=pd.DataFrame({'gene_id':gids},index=pd.Index(g)));a.var_names_make_unique();a.obs_names_make_unique()
  # Preserve library identities. Only use documented, unambiguous barcode suffixes.
  if meta['id']=='GSE178121':a.obs['sample']=['Control pool' if x.endswith('_con') else 'STZ diabetes pool' if x.endswith('_stz') else 'Unmapped' for x in b]
  elif meta['id']=='GSE150703':a.obs['sample']=['Barcode prefix: '+x.rsplit('_',1)[0] for x in b]
  else:a.obs['sample']=sid
  a.obs['source_barcode']=b;a.obs_names=[sid+':'+x for x in a.obs_names];ads.append(a);normalization.append(norm)
 if any(normalization) and not all(normalization):raise Unavailable('Files mix raw counts and noninteger values; scale review is required before merging')
 if not ads:raise Unavailable('No verified single-cell matrices remain after excluding bulk sample files')
 merged=ad.concat(ads,join='outer',merge='first',fill_value=0);merged.obs_names_make_unique()
 return merged,any(normalization)

LATEST_QUERY='(diabetic retinopathy[All Fields] OR diabetic retina[All Fields] OR diabetic macular edema[All Fields]) AND ("single cell"[All Fields] OR "single-cell"[All Fields] OR "single nucleus"[All Fields]) AND gse[Entry Type]'
def latest_geo_candidates():
 obj=request('https://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',params={'db':'gds','term':LATEST_QUERY,'retmode':'json','retmax':1000,'sort':'pdat'}).json()
 ids=obj.get('esearchresult',{}).get('idlist',[])
 if not ids:return []
 time.sleep(.4)
 data=request('https://www.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',params={'db':'gds','id':','.join(ids),'retmode':'json'}).json().get('result',{})
 return [v['accession'] for v in data.values() if isinstance(v,dict) and re.fullmatch(r'GSE\d+',v.get('accession',''))]

def candidate_is_single_cell(meta):
 sc=r'single[\s-]*(cell|nucle[ui])|scRNA|snRNA'
 primary=meta.get('title','')+' '+meta.get('design','')
 topic=meta.get('title','')+' '+meta.get('summary','')
 processing=' '.join(' '.join(m.get('data_processing',[])) for m in meta.get('samples',{}).values())
 relevant=bool(re.search(r'diabet.*(retin|macular)|(retin|macular).*diabet',topic,re.I|re.S))
 technical=bool(re.search(r'cell.?ranger|10x|drop.?seq|cell barcode|single.cell',processing,re.I))
 return relevant and bool(re.search(sc,primary,re.I) or (re.search(sc,meta.get('summary',''),re.I) and technical))

def source_fingerprint(meta):
 headers=[]
 for url in meta.get('files',[]):
  # Repository timestamps/ETags avoid repeatedly downloading unchanged matrices.
  r=request(url,stream=True)
  headers.append({'url':url,'bytes':r.headers.get('Content-Length'),'etag':r.headers.get('ETag'),'modified':r.headers.get('Last-Modified')});r.close()
 semantic={k:meta.get(k) for k in ['id','title','species','summary','design','samples','files']}
 value=json.dumps({'metadata':semantic,'fileHeaders':headers},sort_keys=True,separators=(',',':'))
 return hashlib.sha256(value.encode()).hexdigest(),headers
