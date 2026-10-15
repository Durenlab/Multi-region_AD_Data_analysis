"""Small EC comparison from existing results and RNA observation metadata only."""
from pathlib import Path
import json
import numpy as np,pandas as pd,h5py
from scipy.stats import spearmanr,pearsonr,rankdata,t as student_t
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
R=Path(__file__).resolve().parents[2];B=R.parent/'AD_REGA_region_wise';D=R/'experiments/EXP004_results';O=D/'ec_comparison';T=O/'tables';T.mkdir(parents=True,exist_ok=True)
CTS=['Exc','Inh','Ast','Oli','OPC','Mic_Immune'];MOD=['EC_M14','EC_M17'];checks=[];sources=[]
def save(d,name):d.to_csv(T/(name+'.tsv'),sep='\t',index=False)
def check(name,ok,detail=''):checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
def read(p,**kw):
 p=Path(p);d=pd.read_csv(p,**kw);sources.append(dict(path=str(p),rows=len(d),columns=';'.join(map(str,d.columns))));return d
def genes(x):return set(str(x).split(';')) if pd.notna(x) and str(x) else set()
def js(x):return ';'.join(sorted(x))
config=dict(modules=MOD,benchmarks=['PFC_M14','MTC_M13'],pathway_rule='Complete GSEA rows; compare same term/context when positive source FDR < 0.05 in either EC module; threshold differences are descriptive.',activity_numerical_zero_relative=1e-8,activity_sensitivity_relative=[0,1e-12,1e-10,1e-8,1e-6],composition_levels=['RNA.Subclass','RNA.Subtype'],composition_denominator='Nuclei in same donor and RNA.Class; available RNA object before CellChat sampling; Vasc excluded.',composition_eligibility='Parent >=20 nuclei in >=40 matched donors; subtype/subclass >=10 nuclei in at least half those donors; variable fractions.',composition_outcomes=['M14','M17','rank_contrast_M17_minus_M14'],composition_tests='Exploratory Spearman and partial rank correlation adjusting categorical pathology, age capped at90, sex, PMI; BH across all eligible features ×3 outcomes separately for raw and adjusted families.',resampling='2000 donor bootstrap samples for primary correlation; leave-one-donor-out sensitivity',seed=400417)
(O/'analysis_configuration.json').write_text(json.dumps(config,indent=2))
# Same-cell pathway and leading-edge comparisons, independent of hard assignments.
g=read(D/'regional/tables/full_GSEA.tsv.gz',sep='\t');g=g[(g.kind=='cross')&g.module_id.isin(MOD)].copy();g['key']=g.term.str.strip().str.upper()
assign=read(D/'tables/cross_gene_assignments.tsv',sep='\t');assign=assign[assign.candidate.isin(MOD)]
focus_keys={'RIBOSOME','OXIDATIVE PHOSPHORYLATION','ALZHEIMER DISEASE','ANTIGEN PROCESSING AND PRESENTATION','COMPLEMENT AND COAGULATION CASCADES','CYTOKINE-CYTOKINE RECEPTOR INTERACTION','ECM-RECEPTOR INTERACTION','NEUROACTIVE LIGAND-RECEPTOR INTERACTION'}
focus_ids=['GO:0003341','GO:0060285','GO:0002399','GO:0002483','GO:1903053','GO:0010043']
focus_keys.update(g.loc[g.key.str.contains('|'.join(focus_ids)), 'key'])
comparisons=[];detail=[];weight_tables=[]
for ct in CTS:
 ctx=ct+'_EC';p=B/f'part3_process_REGA_output/09extracted_matrices_region_shareH/EC/{ctx}_WBV_genes_by_module.tsv';w=read(p,sep='\t',index_col=0);check('Unique WBV genes '+ctx,w.index.is_unique)
 wn=w.div(w.abs().sum(axis=0),axis=1);wn=wn.div(wn.abs().sum(axis=1).replace(0,np.nan),axis=0).fillna(0);rp=wn[['M14','M17']].rank(pct=True)
 ws=pd.DataFrame({'gene':w.index,'context':ctx,'normalized_WBV_M14':wn.M14,'normalized_WBV_M17':wn.M17,'rank_percentile_M14':rp.M14,'rank_percentile_M17':rp.M17,'WBV_source':str(p),'WBV_source_row':np.arange(2,len(w)+2)}).reset_index(drop=True)
 amap=assign[assign.celltype_region==ctx].set_index('gene').candidate.to_dict();ws['hard_assignment']=ws.gene.map(amap).fillna('other_or_unassigned');weight_tables.append(ws);wi=ws.set_index('gene')
 a=g[(g.context==ctx)&(g.module_id=='EC_M14')].set_index(['library','key']);b=g[(g.context==ctx)&(g.module_id=='EC_M17')].set_index(['library','key'])
 keys=set(a[a.positive_fdr05].index)|set(b[b.positive_fdr05].index)
 for lib,key in sorted(keys):
  x=a.loc[(lib,key)] if (lib,key) in a.index else None;y=b.loc[(lib,key)] if (lib,key) in b.index else None
  gx=genes(x.leading_genes) if x is not None else set();gy=genes(y.leading_genes) if y is not None else set();both=gx&gy
  rec=dict(context=ctx,library=lib,term=key,focus_panel=key in focus_keys,positive_M14=bool(x.positive_fdr05) if x is not None else False,positive_M17=bool(y.positive_fdr05) if y is not None else False,n_leading_M14=len(gx),n_leading_M17=len(gy),n_shared=len(both),leading_jaccard=len(both)/len(gx|gy) if gx|gy else np.nan,shared_genes=js(both),M14_only_leading=js(gx-gy),M17_only_leading=js(gy-gx))
  for label,z in [('M14',x),('M17',y)]:
   for col in ['nes','fdr','pval','source_path','source_row']:rec[label+'_'+col]=z[col] if z is not None else np.nan
  comparisons.append(rec)
  if key not in focus_keys:continue
  for gene in sorted(gx|gy):
   ww=wi.loc[gene].to_dict() if gene in wi.index else {}
   detail.append(dict(context=ctx,library=lib,term=key,gene=gene,leading_M14=gene in gx,leading_M17=gene in gy,positive_M14=rec['positive_M14'],positive_M17=rec['positive_M17'],M14_GSEA_source=rec['M14_source_path'],M14_GSEA_row=rec['M14_source_row'],M17_GSEA_source=rec['M17_source_path'],M17_GSEA_row=rec['M17_source_row'],**{k:v for k,v in ww.items() if k!='context'}))
pc=pd.DataFrame(comparisons);gd=pd.DataFrame(detail);weights=pd.concat(weight_tables);save(pc,'pathway_leading_edge_comparison');save(gd,'focus_pathway_gene_details');save(weights,'EC_normalized_gene_weights')
# Shared H and matching metadata. Preserve source values; numerical zero handling only in new comparisons.
hp=B/'part2_run_REGA/07REGA_region_shareH_output/AD_region_shareH_EC_K20/H_with_samples.csv';h=read(hp,index_col=0).T;h.index.name='donor_id';h.columns=h.columns.str.replace('module','M',regex=False)
metap=B/'part1_prepare_REGA/01region_consistent_metadata/Ast_EC_region_consistent_metadata.tsv';meta=read(metap,sep='\t');check('Unique metadata donors',meta.ROSMAP_IndividualID.is_unique)
meta['metadata_source_row']=np.arange(2,len(meta)+2);meta=meta.set_index('ROSMAP_IndividualID').reindex(h.index);check('All H donors have metadata',meta.Pathology.notna().all());check('46 EC donors',len(h)==46)
don=meta.copy();don['H_source']=str(hp);don['H_source_column']=np.arange(2,len(h)+2);don['metadata_source']=str(metap);don['H_M14_raw']=h.M14;don['H_M17_raw']=h.M17;don['age_capped90']=pd.to_numeric(don.age_death).clip(upper=90);don['pmi']=pd.to_numeric(don.pmi,errors='coerce')
for mod in ['M14','M17']:
 don[mod]=h[mod].mask(h[mod]<h[mod].max()*1e-8,0);don[mod+'_rank_pct']=don[mod].rank(pct=True);don[mod+'_z']=(don[mod]-don[mod].mean())/don[mod].std(ddof=1)
don['rank_contrast_M17_minus_M14']=don.M17_rank_pct-don.M14_rank_pct
X=pd.DataFrame({'intercept':1.,'earlyAD':(don.Pathology=='earlyAD').astype(float),'lateAD':(don.Pathology=='lateAD').astype(float),'age_capped90':don.age_capped90,'female':(don.Sex=='Female').astype(float),'pmi':don.pmi},index=don.index)
check('Complete available covariates',X.notna().all().all())
def partial(x,y,design):
 x=np.asarray(x,float);y=np.asarray(y,float);z=np.asarray(design,float);rx=rankdata(x);ry=rankdata(y);ex=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0];ey=ry-z@np.linalg.lstsq(z,ry,rcond=None)[0];r=float(np.corrcoef(ex,ey)[0,1]);df=len(x)-np.linalg.matrix_rank(z)-1;pv=2*student_t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df);return r,float(pv),int(df)
cor=[]
for eps in config['activity_sensitivity_relative']:
 z=h[['M14','M17']].mask(h[['M14','M17']]<h[['M14','M17']].max()*eps,0)
 for method,fn in [('Spearman',spearmanr),('Pearson',pearsonr)]:
  r,p=fn(z.M14,z.M17);cor.append(dict(analysis=method,relative_zero_threshold=eps,n=len(z),effect=r,p=p,n_zero_M14=int((z.M14==0).sum()),n_zero_M17=int((z.M17==0).sum()),covariates='none'))
for label,design in [('partial_rank_pathology',X.iloc[:,:3]),('partial_rank_pathology_age_sex_pmi',X)]:
 r,p,df=partial(don.M14,don.M17,design);cor.append(dict(analysis=label,relative_zero_threshold=1e-8,n=len(don),effect=r,p=p,covariates=';'.join(design.columns),df=df))
for group,z in don.groupby('Pathology'):
 r,p=spearmanr(z.M14,z.M17);cor.append(dict(analysis='within_'+group,relative_zero_threshold=1e-8,n=len(z),effect=r,p=p,covariates='within pathology group'))
co=pd.DataFrame(cor);save(co,'donor_activity_correlations')
rng=np.random.default_rng(config['seed']);boot=[]
for _ in range(2000):
 idx=rng.integers(0,len(don),len(don));z=don.iloc[idx];boot.append(spearmanr(z.M14,z.M17).statistic)
ci=np.nanquantile(boot,[.025,.975]);loo=[]
for i in range(len(don)):
 z=don.drop(don.index[i]);loo.append(spearmanr(z.M14,z.M17).statistic)
save(pd.DataFrame([dict(analysis='floored Spearman',bootstrap_replicates=2000,CI025=ci[0],CI975=ci[1],leave_one_out_min=min(loo),leave_one_out_max=max(loo))]),'donor_correlation_sensitivity')
# Reproduce original linear slopes/p to check donor alignment; original selection thresholds unchanged.
mp=read(D/'tables/module_pathology.tsv',sep='\t');stage=don.Pathology.map({'nonAD':0,'earlyAD':1,'lateAD':2})
for mod in ['M14','M17']:
 fit=sm.OLS(h[mod],sm.add_constant(stage)).fit();z=mp[(mp.candidate=='EC_'+mod)&(mp.test=='Linear_Regression')].iloc[0];check('Original slope/p reproduced '+mod,np.isclose(fit.params.iloc[1],z.slope,rtol=1e-6) and np.isclose(fit.pvalues.iloc[1],z.pval,rtol=1e-6),f'{fit.params.iloc[1]}, {fit.pvalues.iloc[1]}')
# Read observation metadata only; never load RNA X or expression layers.
rnap=B/'part8_cellchat/AD_scRNA_used.h5ad'
with h5py.File(rnap,'r') as f:
 obs=f['obs'];available=list(obs.keys())
 def col(key):
  obj=obs[key]
  if isinstance(obj,h5py.Group):
   cat=obj['categories'].asstr()[:] if obj['categories'].dtype.kind in 'OS' else obj['categories'][:];codes=obj['codes'][:];return np.where(codes>=0,cat[np.maximum(codes,0)],None)
  return obj.asstr()[:] if obj.dtype.kind in 'OS' else obj[:]
 mask=col('BrainRegion')=='EC';cells=pd.DataFrame({key:col(key)[mask] for key in ['ROSMAP_IndividualID','RNA.Class','RNA.Subclass','RNA.Subtype','Pathology','Study']})
sources.append(dict(path=str(rnap),rows=len(cells),columns='obs only: '+ ';'.join(available)));save(pd.DataFrame({'available_obs_field':available}),'RNA_metadata_fields')
cells=cells[cells.ROSMAP_IndividualID.isin(don.index)];check('All matched donors in RNA metadata',set(cells.ROSMAP_IndividualID)==set(don.index));check('RNA pathology matches H metadata',all(cells.groupby('ROSMAP_IndividualID').Pathology.nunique()==1) and (cells.groupby('ROSMAP_IndividualID').Pathology.first().reindex(don.index)==don.Pathology).all())
classcounts=cells.groupby(['ROSMAP_IndividualID','RNA.Class'],observed=True).size().unstack(fill_value=0).reindex(don.index);save(classcounts.reset_index(),'donor_class_nucleus_counts')
counts=[];fractions=[];elig=[];tests=[]
for level in config['composition_levels']:
 grouped=cells.groupby(['RNA.Class',level,'ROSMAP_IndividualID'],observed=True).size()
 for (ct,label),part in grouped.groupby(level=[0,1]):
  if ct not in CTS:continue
  vals=part.droplevel([0,1]).reindex(don.index,fill_value=0);parent=classcounts[ct];valid=parent>=20;prop=vals/parent.replace(0,np.nan);enough=int(((vals>=10)&valid).sum());eligible=valid.sum()>=40 and enough>=np.ceil(valid.sum()/2) and prop[valid].nunique()>1
  fid=level+'|'+ct+'|'+str(label);elig.append(dict(feature=fid,level=level,celltype=ct,label=label,n_parent20=int(valid.sum()),n_subpopulation10=enough,eligible=eligible))
  for donor in don.index:counts.append(dict(donor_id=donor,feature=fid,level=level,celltype=ct,label=label,n_subpopulation=int(vals[donor]),n_parent=int(parent[donor]),fraction=prop[donor],eligible_feature=eligible,parent_coverage=bool(valid[donor])))
  if not eligible:continue
  for outcome in config['composition_outcomes']:
   r,p=spearmanr(don.loc[valid,outcome],prop[valid]);pr,pp,df=partial(don.loc[valid,outcome],prop[valid],X.loc[valid]);tests.append(dict(feature=fid,level=level,celltype=ct,label=label,outcome=outcome,n=int(valid.sum()),rho=r,p=p,partial_r=pr,partial_p=pp,partial_df=df))
ctab=pd.DataFrame(tests)
for pcol,qcol in [('p','q_BH'),('partial_p','partial_q_BH')]:
 finite=ctab[pcol].notna();ctab[qcol]=np.nan;ctab.loc[finite,qcol]=multipletests(ctab.loc[finite,pcol],method='fdr_bh')[1]
countframe=pd.DataFrame(counts);aliases=[]
wide=countframe.pivot(index='donor_id',columns='feature',values='fraction')
for fid in ctab.feature.unique():
 same=[other for other in ctab.feature.unique() if np.allclose(wide[fid],wide[other],equal_nan=True,rtol=0,atol=1e-14)]
 aliases.append(dict(feature=fid,identical_fraction_labels=';'.join(sorted(same)),n_identical_labels=len(same)))
save(pd.DataFrame(aliases),'subpopulation_fraction_aliases')
save(pd.DataFrame(counts),'donor_subpopulation_counts');save(pd.DataFrame(elig),'subpopulation_eligibility');save(ctab,'subpopulation_activity_associations');save(don.reset_index(),'EC_donor_activities')
# All positively associated assigned TFs: retain motif aliases and edge specificity limits.
tf=read(D/'tables/assigned_TF_activity.tsv',sep='\t');tf=tf[tf.candidate.isin(MOD)&(tf.positive_LM_FDR05|tf.positive_MW_FDR05)].copy()
ed=read(D/'regional/tables/TF_edges_to_pathway_leading_genes.tsv',sep='\t');ed=ed[ed.module_id.isin(MOD)].copy();ed['key']=ed.term.str.strip().str.upper();ed['focus_panel']=ed.key.isin(focus_keys)
motifrows=[]
for ctx,tt in tf.groupby('celltype_region'):
 root=B/f'part4_TF/03TF_activity/{ctx}/TF_activity_pcc_pos_only';pccp=root/f'{ctx}_TF_motif_pcc_M2.csv';bp=root/f'{ctx}_TF_motif_binary_M1.csv';pcc=read(pccp,index_col=0);binary=read(bp,index_col=0).reindex(index=pcc.index,columns=pcc.columns).fillna(0)
 for z in tt.itertuples():
  if z.tf not in pcc.index:continue
  active=(pcc.loc[z.tf]>0)&(binary.loc[z.tf]>0)
  for motif,value in pcc.loc[z.tf,active].items():motifrows.append(dict(module_id=z.candidate,context=ctx,TF=z.tf,motif=motif,pcc_weight=value,PCC_source=str(pccp),binary_source=str(bp)))
mr=pd.DataFrame(motifrows);save(mr,'EC_TF_retained_motifs')
summary=[]
for z in tf.itertuples():
 ee=ed[(ed.module_id==z.candidate)&(ed.context==z.celltype_region)&(ed.TF==z.tf)];focused=ee[ee.focus_panel];mot=mr[(mr.context==z.celltype_region)&(mr.TF==z.tf)]
 summary.append(dict(module_id=z.candidate,context=z.celltype_region,TF=z.tf,LM_beta=z.lm_beta_pathology,LM_p=z.lm_pval_pathology,LM_q=z.lm_fdr_pathology,late_minus_non=z.lateAD_minus_nonAD,MW_p=z.mannwhitney_lateAD_vs_nonAD_pval,MW_q=z.mannwhitney_lateAD_vs_nonAD_fdr,n_pathway_targets=ee.TG.nunique(),n_focus_targets=focused.TG.nunique(),n_focus_targets_topdecile=focused[focused.target_percentile>=.9].TG.nunique(),max_focus_target_percentile=focused.target_percentile.max(),retained_motifs=js(mot.motif),n_retained_motifs=len(mot),activity_source=z.source_path_activity,activity_row=z.source_row_activity,assignment_source=z.source_path_assignment,assignment_row=z.source_row_assignment))
sumtf=pd.DataFrame(summary);sumtf['same_context_motif_aliases']=''
for (ctx,mot),ii in sumtf.groupby(['context','retained_motifs']).groups.items():
 if mot:sumtf.loc[ii,'same_context_motif_aliases']=';'.join(sorted(sumtf.loc[ii,'TF']))
save(sumtf,'EC_TF_pathway_comparison');save(ed,'EC_TF_pathway_edges')
activity_pairs=[]
for ctx,tt in tf.groupby('celltype_region'):
 ap=B/f'part4_TF/03TF_activity/{ctx}/TF_activity_pcc_pos_only/{ctx}_TF_activity_pcc_cutoff.csv';act=read(ap,index_col=0)
 names=[n for n in tt.tf if n in act.index]
 for i,left in enumerate(names):
  for right in names[i+1:]:
   x=act.loc[left,h.index].to_numpy(float);y=act.loc[right,h.index].to_numpy(float);r=pearsonr(x,y).statistic
   activity_pairs.append(dict(context=ctx,TF1=left,TF2=right,n=len(x),pearson_r=r,nearly_identical_shape=r>1-1e-10,source_path=str(ap),TF1_source_row=int(act.index.get_loc(left))+2,TF2_source_row=int(act.index.get_loc(right))+2))
save(pd.DataFrame(activity_pairs),'EC_TF_activity_profile_comparison')
# Candidate communication must intersect significant focus-pathway leading genes in assigned cell context.
cc=read(D/'tables/cellchat_module_links.tsv',sep='\t');cc=cc[cc.candidate.isin(MOD)].copy();gl={}
for mod in ['M14','M17']:
 for ct in CTS:
  s=gd[(gd.context==ct+'_EC')&gd['leading_'+mod]&gd['positive_'+mod]];gl[('EC_'+mod,ct)]=set(s.gene)
cc['pathway_ligand_genes']=[js(genes(z.ligand_module_genes)&gl.get((z.candidate,z.source),set())) for z in cc.itertuples()];cc['pathway_receptor_genes']=[js(genes(z.receptor_module_genes)&gl.get((z.candidate,z.target),set())) for z in cc.itertuples()]
cc=cc[(cc.pathway_ligand_genes!='')|(cc.pathway_receptor_genes!='')];save(cc,'EC_pathway_compatible_CellChat')
bench=read(D/'tables/ldsc_all.tsv',sep='\t');save(bench[(bench.kind=='cross')&bench.candidate.isin(['EC_M14','EC_M17','PFC_M14','MTC_M13'])],'genetic_benchmarks')
check('All six contexts in pathway comparison',set(pc.context)==set(c+'_EC' for c in CTS));check('All eligible EC TF candidates checked against available source',len(sumtf)==len(tf));check('No infinite fraction-test statistics',not np.isinf(ctab[['rho','p','partial_r','partial_p']].to_numpy()).any())
save(pd.DataFrame(sources).drop_duplicates('path'),'source_manifest');save(pd.DataFrame(checks),'analysis_checks')
print('Pathway comparisons',len(pc),'focus gene rows',len(gd),'donors',len(don),'eligible subpopulation features',sum(x['eligible'] for x in elig),'TF candidates',len(tf),'compatible CellChat rows',len(cc));print(co.to_string(index=False));print('Subpopulation top adjusted associations');print(ctab.sort_values('partial_q_BH').head(12).to_string(index=False));print('Checks',all(x['passed'] for x in checks))
if not all(x['passed'] for x in checks):raise RuntimeError('Source consistency check failed')
