import json,time,gzip,hashlib
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
import scanpy as sc
import pandas as pd
from sklearn.manifold import TSNE

PANELS={'Rod':['RHO','GNAT1','NRL','PDE6A'],'Cone':['ARR3','GNAT2','OPN1SW','PDE6C'],'Bipolar':['VSX2','TRPM1','GRM6','CABP5'],'Müller glia':['RLBP1','GLUL','SLC1A3'],'Amacrine':['GAD1','GAD2','TFAP2A'],'Ganglion':['RBPMS','SNCG','POU4F2'],'Endothelial':['PECAM1','VWF','KDR','CLDN5'],'Pericyte / smooth muscle':['RGS5','PDGFRB','CSPG4','ACTA2'],'Fibroblast':['COL1A1','COL1A2','DCN','LUM'],'Myeloid':['LST1','TYROBP','AIF1','C1QA'],'T cell':['CD3D','CD3E','TRAC'],'B cell':['MS4A1','CD79A','CD79B'],'RPE':['RPE65','BEST1','LRAT']}

def write_json(path,obj):
 temp=Path(str(path)+'.tmp');temp.write_text(json.dumps(obj,separators=(',',':'),allow_nan=False));temp.replace(path)

def save_h5ad_atomic(a,path):
 path=Path(path);pending=path.with_name(path.stem+'.pending.h5ad')
 a.write_h5ad(pending,compression='gzip');pending.replace(path)

def embedding(a,progress,tsne=True,checkpoint=None,checkpoint_commit=None):
 if a.n_obs<12 or a.n_vars<10:raise ValueError('At least 12 cells and 10 genes are needed for reclustering')
 if checkpoint and Path(checkpoint).exists():
  progress('Resuming saved PCA, UMAP and clusters')
  a=sc.read_h5ad(checkpoint)
 else:
  progress('Selecting variable genes and computing scaled PCA')
  sc.pp.highly_variable_genes(a,n_top_genes=min(2000,a.n_vars),flavor='seurat')
  b=a[:,a.var.highly_variable].copy();sc.pp.scale(b,max_value=10)
  sc.tl.pca(b,n_comps=min(30,b.n_vars-1,b.n_obs-1),svd_solver='arpack',random_state=42)
  a.obsm['X_pca']=b.obsm['X_pca'];del b
  progress('Computing neighbors, UMAP and Leiden clusters')
  sc.pp.neighbors(a,n_neighbors=min(15,a.n_obs-1),use_rep='X_pca',random_state=42)
  sc.tl.umap(a,min_dist=.3,random_state=42)
  sc.tl.leiden(a,resolution=.6,key_added='cluster',random_state=42,flavor='igraph',n_iterations=2,directed=False)
  if checkpoint:
   save_h5ad_atomic(a,checkpoint)
   if checkpoint_commit:checkpoint_commit()
 if tsne:
  progress('Computing t-SNE projection')
  if a.n_obs>50000:
   from openTSNE import TSNE as FastTSNE
   a.obsm['X_tsne']=np.asarray(FastTSNE(perplexity=30,random_state=42,initialization='pca',negative_gradient_method='fft',early_exaggeration_iter=250,n_iter=250,n_jobs=2).fit(a.obsm['X_pca']))
   a.uns['tsne_method']='openTSNE FFT approximation; 250 early-exaggeration plus 250 optimization iterations.'
  else:
   a.obsm['X_tsne']=TSNE(n_components=2,perplexity=min(30,max(3,(a.n_obs-1)/3)),random_state=42,init='pca',learning_rate='auto',max_iter=500,n_jobs=2).fit_transform(a.obsm['X_pca'])
   a.uns['tsne_method']='scikit-learn Barnes-Hut t-SNE; 500 total iterations.'
 return a

def process(a,meta,out,progress,source_normalized=False,checkpoint_commit=None):
 out.mkdir(exist_ok=True,parents=True);input_cells=a.n_obs;a.uns['pipeline_version']='2.1.2'
 if source_normalized:
  if meta['id']!='GSE150703' or 'log(1+10,000*' not in json.dumps(meta.get('samples',{})):
   raise ValueError('Noninteger expression data need scale verification before normalization. Source files retained.')
  if not np.isfinite(a.X.data).all() or (a.X.data<0).any() or (a.X.data>25).any():raise ValueError('Invalid deposited log-normalized values')
  # GEO explicitly documents log(1 + 10000 * per-cell normalized DGE).
  # Do not reconstruct or claim raw counts; retain the deposited scale.
  a.obs['n_genes_by_counts']=np.asarray((a.X>0).sum(1)).ravel()
  a.obs['total_counts']=0;a.obs['pct_counts_mt']=0
  a=a[a.obs.n_genes_by_counts>=200].copy()
  if a.n_obs>250000:raise ValueError('Study exceeds the current 250,000-cell processing limit')
  sc.pp.filter_genes(a,min_cells=3)
  a.uns['source_normalized']=True
  meta['normalization']='Deposited ln(1 + CP10K), retained unchanged as documented by GEO. Raw counts unavailable.'
  meta['qc']='At least 200 detected genes per cell; genes detected in at least 3 cells. Raw-count totals and mitochondrial count QC are unavailable and are not inferred.'
  meta['rawCountsAvailable']=False
  a.uns['input_cells']=input_cells;a=embedding(a,progress,tsne=True,checkpoint=out/'embedding-checkpoint.h5ad',checkpoint_commit=checkpoint_commit);save_h5ad_atomic(a,out/'study.h5ad');export(a,meta,out,input_cells);return a
 if not np.isfinite(a.X.data).all() or (a.X.data<0).any():raise ValueError('Invalid counts')
 progress('Filtering cells and normalizing deposited counts')
 a.var['mt']=a.var_names.str.upper().str.startswith('MT-')
 sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
 a=a[(a.obs.n_genes_by_counts>=200)&(a.obs.pct_counts_mt<=20)].copy()
 if a.n_obs<12:raise ValueError('Too few cells passed the documented quality filters')
 if a.n_obs>250000:raise ValueError('Study exceeds the current 250,000-cell processing limit; files retained for a larger worker')
 # Preserve full-feature library totals before filtering genes.
 totals=np.asarray(a.X.sum(axis=1)).ravel()
 sc.pp.filter_genes(a,min_cells=3);a.layers['counts']=a.X.copy()
 a.X=sp.diags(10000/np.maximum(totals,1)).dot(a.X).tocsr();a.X.data=np.log1p(a.X.data)
 a.uns['input_cells']=input_cells
 a=embedding(a,progress,tsne=True,checkpoint=out/'embedding-checkpoint.h5ad',checkpoint_commit=checkpoint_commit)
 progress('Saving expression matrices and analysis provenance')
 a.uns['input_cells']=input_cells;save_h5ad_atomic(a,out/'study.h5ad')
 export(a,meta,out,input_cells)
 return a

def export(a,meta,out,input_cells=None):
 X=sp.csr_matrix(a.X);g=a.var_names.tolist();clusters=sorted(a.obs.cluster.unique(),key=int);samples=list(dict.fromkeys(a.obs['sample'].astype(str)));sample_codes=[samples.index(s) for s in a.obs['sample'].astype(str)];cluster_codes=[clusters.index(s) for s in a.obs.cluster.astype(str)];lookup={s.upper():i for i,s in enumerate(g)}
 overall=np.asarray(X.mean(0)).ravel();groups=[]
 for i,k in enumerate(clusters):
  mask=a.obs.cluster.astype(str)==k;means=np.asarray(X[mask.to_numpy()].mean(0)).ravel()
  scores={p:float(np.mean([means[lookup[v]] for v in gs if v in lookup])) for p,gs in PANELS.items() if any(v in lookup for v in gs)}
  top=max(scores,key=scores.get) if scores else 'Unassigned';groups.append(dict(id=i,name=f'C{i} · {top}',n=int(mask.sum()),markers=[g[j] for j in np.argsort(means-overall)[-5:][::-1]],provisional=True))
 cell=dict(x=np.round(a.obsm['X_umap'][:,0],4).tolist(),y=np.round(a.obsm['X_umap'][:,1],4).tolist(),tx=np.round(a.obsm['X_tsne'][:,0],4).tolist(),ty=np.round(a.obsm['X_tsne'][:,1],4).tolist(),cluster=cluster_codes,sample=sample_codes,barcode=a.obs_names.tolist(),counts=a.obs.total_counts.astype(int).tolist(),detected=a.obs.n_genes_by_counts.astype(int).tolist(),mito=np.round(a.obs.pct_counts_mt,2).tolist())
 write_json(out/'cells.json',cell)
 det=np.asarray((X>0).mean(0)).ravel()*100;sample_means=np.vstack([np.asarray(X[(a.obs['sample'].astype(str)==s).to_numpy()].mean(0)).ravel() for s in samples]).T
 summaries=np.column_stack([overall,det,sample_means]);write_json(out/'genes.json',dict(genes=g,summary=np.round(summaries,5).tolist()))
 m={**meta,'inputCells':input_cells or int(a.uns.get('input_cells',a.n_obs)),'cells':a.n_obs,'genes':a.n_vars,'samples':samples,'sampleMetadata':meta.get('samples',{}),'sampleCounts':[sample_codes.count(i) for i in range(len(samples))],'groups':groups,'removedCells':(input_cells or int(a.uns.get('input_cells',a.n_obs)))-a.n_obs,'rawCountsAvailable':meta.get('rawCountsAvailable',True),'normalization':meta.get('normalization','ln(1 + 10,000 × raw count / total full-feature cell counts)'),'qc':meta.get('qc','Cells: at least 200 detected genes, at most 20% mitochondrial counts. Genes: detected in at least 3 retained cells.'),'embedding':'2,000 highly variable genes (Seurat dispersion method), scaled and clipped at 10; 30-component PCA; 15-neighbor graph; UMAP and t-SNE; seed 42.','clustering':'Leiden resolution 0.6. Marker-panel suggestions are provisional, not author annotations.','processedAt':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'pipelineVersion':str(a.uns.get('pipeline_version','2.0.0')),'tsneMethod':str(a.uns.get('tsne_method','scikit-learn Barnes-Hut t-SNE; 500 total iterations.')),'conditions':samples}
 write_json(out/'meta.json',m)

def differential(a,one,two,min_pct=.1,logfc=.25):
 one=np.unique(np.asarray(one,dtype=int));two=np.unique(np.asarray(two,dtype=int))
 if min(len(one),len(two))<3:raise ValueError('Each group needs at least 3 cells')
 if np.intersect1d(one,two).size:raise ValueError('Comparison groups must not overlap')
 if min(one.min(),two.min())<0 or max(one.max(),two.max())>=a.n_obs:raise ValueError('Cell index is out of range')
 X=sp.csc_matrix(a.X);A=X[one];B=X[two];pa=np.asarray((A>0).mean(0)).ravel();pb=np.asarray((B>0).mean(0)).ravel()
 # Fold change uses arithmetic means on the linear CP10K scale, not mean log values.
 LA=A.copy();LB=B.copy();LA.data=np.expm1(LA.data);LB.data=np.expm1(LB.data)
 ma=np.asarray(LA.mean(0)).ravel();mb=np.asarray(LB.mean(0)).ravel();fc=np.log2((ma+1)/(mb+1))
 eligible=np.flatnonzero((np.maximum(pa,pb)>=min_pct)&(np.abs(fc)>=logfc));rows=[]
 for j in eligible:
  av=A[:,j].toarray().ravel();bv=B[:,j].toarray().ravel();pv=float(mannwhitneyu(av,bv,alternative='two-sided',method='asymptotic').pvalue)
  rows.append(dict(gene=str(a.var_names[j]),log2FC=float(fc[j]),deltaPct=float(100*(pa[j]-pb[j])),pctA=float(pa[j]*100),pctB=float(pb[j]*100),meanA=float(ma[j]),meanB=float(mb[j]),pvalue=pv))
 # BH across tested genes; disclose filtering and family size.
 if rows:
  order=np.argsort([r['pvalue'] for r in rows]);pvals=np.array([rows[i]['pvalue'] for i in order]);adj=np.minimum.accumulate((pvals*len(rows)/np.arange(1,len(rows)+1))[::-1])[::-1]
  for idx,q in zip(order,adj):rows[idx]['qvalue']=float(min(q,1))
 rows.sort(key=lambda r:(r['qvalue'],-abs(r['log2FC'])))
 return dict(rows=rows,nA=len(one),nB=len(two),testedGenes=len(rows),minPct=min_pct,log2FCThreshold=logfc,method='Two-sided Mann–Whitney U (Wilcoxon rank-sum), asymptotic tie correction; Benjamini–Hochberg across tested genes. log2FC uses mean CP10K with pseudocount 1.',warning='Exploratory cell-level statistics. Cells are not independent biological replicates. Results do not establish a replicated disease effect.')
