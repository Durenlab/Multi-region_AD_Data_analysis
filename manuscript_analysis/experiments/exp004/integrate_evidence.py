#!/usr/bin/env python
"""Extract existing EXP004 evidence. No raw data changes or new model fitting."""
from pathlib import Path
import json, re
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[2]; BASE=ROOT.parent
OUT=ROOT/'experiments/EXP004_results'; TABLE=OUT/'tables'
TABLE.mkdir(parents=True,exist_ok=True)
CAND=['AG_M13','EC_M14','EC_M17','MTC_M6','MTC_M13','PFC_M14','PFC_M18']
SPEC='Exc_AG_M8 Inh_AG_M8 Ast_AG_M2 Oli_AG_M6 Vasc_AG_M11 Exc_EC_M13 Ast_EC_M17 Oli_EC_M14 Oli_EC_M16 Mic_immune_EC_M16 Vasc_EC_M11 Ast_HC_M13 OPC_HC_M13 Exc_MTC_M8 Exc_MTC_M19 Ast_MTC_M8 Oli_MTC_M8 Oli_MTC_M13 OPC_MTC_M8 Mic_immune_MTC_M8 Exc_PFC_M1 OPC_PFC_M14 Vasc_PFC_M20'.split()
manifest=[]; checks=[]
def read(path,unit,**kwargs):
 p=Path(path);p=p if p.is_absolute() else BASE/p
 if not p.exists():
  manifest.append(dict(path=str(p),available=False,unit=unit));return pd.DataFrame()
 d=pd.read_csv(p,sep=kwargs.pop('sep','\t'),**kwargs)
 manifest.append(dict(path=str(p),available=True,unit=unit,rows=len(d),columns=';'.join(d.columns)))
 d['source_path']=str(p);start=1 if kwargs.get('header', 'infer') is None else 2;d['source_row']=np.arange(start,len(d)+start)
 return d

def save(d,name):d.to_csv(TABLE/(name+'.tsv'),sep='\t',index=False,na_rep='NA')
def bhcheck(d,groups,p,q,label):
 errors=[]
 for _,a in d.groupby(groups,dropna=False):
  ok=np.isfinite(a[p]); vals=multipletests(a.loc[ok,p],method='fdr_bh')[1] if ok.any() else []
  errors.extend(np.abs(vals-a.loc[ok,q]).tolist())
 err=max(errors,default=0);checks.append(dict(check=label,passed=bool(err<1e-9),max_abs_error=err))
 return err
cross=BASE/'AD_REGA_region_wise';single=BASE/'AD_REGA'
mp=read(cross/'part2_run_REGA/08region_sharedH_pathology_significant_test/ALL_regions_sharedH_All_four_tests_H_pathology_results.tsv','region x module x test/contrast')
bhcheck(mp,['Region','test','Comparison'],'pval','pval_adj','module BH across 20 modules per region/test/contrast')
bhcheck(mp,['Region','test'],'pval','pval_adj_global','module global column BH within region/test')
mp['candidate']=mp.Region+'_'+mp.Module
mp['late_minus_non_rank_biserial']=np.where(mp.Comparison=='nonAD vs lateAD',1-2*mp.statistic/(mp.n1*mp.n2),np.nan)
primary=mp[(mp.test=='Linear_Regression')|(mp.Comparison=='nonAD vs lateAD')]
selected=set(primary[(primary.pval_adj<.05)&(((primary.test=='Linear_Regression')&(primary.slope>0))|((primary.Comparison=='nonAD vs lateAD')&(primary.late_minus_non_rank_biserial>0)))].candidate)
checks.append(dict(check='positive union equals seven specified candidates',passed=selected==set(CAND),max_abs_error=0))
save(mp[mp.candidate.isin(CAND)],'module_pathology');print('UNION',selected,flush=True)
sp=read(single/'part2_run_REGA/H_pathology_significant_test/module_significancy/selected_module_significant_test_details.tsv','celltype x region x module x test/contrast');save(sp,'single_module_pathology')
lds=[]; acts=[]
for kind,base,ld,tf in [('cross',cross,'part5_ldsc/h2_result_phase3_single_module_extend2kb_FC/AD_Bellenguez_single_module_extend2kb_LDSC_results.tsv','part4_TF/04DriverTF_lm/all_celltype_region_TF_pathology_association.tsv'),('single',single,'part6_ldsc/h2_result_phase3_single_module_extend2kb/AD_Bellenguez_single_module_extend2kb_LDSC_results.tsv','part7_TF/02TF_pathology_association_lm/all_celltype_region_TF_pathology_association.tsv')]:
 d=read(base/ld,'celltype x region x module LDSC annotation');d['kind']=kind
 bhcheck(d,['kind'],'Enrichment_p','FDR_BH_Enrichment_p',kind+' LDSC BH across full source table '+str(len(d)))
 d['candidate']=d.celltype_region.str.rsplit('_',n=1).str[-1]+'_'+d.module;lds.append(d)
 t=read(base/tf,'celltype x region x TF activity association');t['kind']=kind
 for p,q in [('lm_pval_pathology','lm_fdr_pathology'),('mannwhitney_lateAD_vs_nonAD_pval','mannwhitney_lateAD_vs_nonAD_fdr')]:bhcheck(t,['celltype_region'],p,q,kind+' TF '+q+' within context')
 acts.append(t)
ldsc=pd.concat(lds);save(ldsc,'ldsc_all');activity=pd.concat(acts)
assignments=[];single_sets={};cross_sets={};overlaps=[];counts=[]
for kind in ['gene','RE','tf']:
 allrows=[]
 for region in ['AG','EC','MTC','PFC']:
  sub={'gene':f'part3_process_REGA_output/10-1gene_module_assignment_region_shareH_FC/{region}/{region}_ALL_cellgene_FC_assignment_with_unassigned.tsv','RE':f'part3_process_REGA_output/10-2RE_module_assignment_region_shareH_FC/{region}/{region}_ALL_cellRE_FC_assignment_with_unassigned.tsv','tf':f'part4_TF/02tf_assignment_FC/{region}/{region}_ALL_cell_TF_FC_assignment_with_unassigned.tsv'}[kind]
  d=read(cross/sub,'celltype x region x '+kind+' hard assignment')
  d['candidate']=region+'_'+d.assigned_module.fillna('unassigned')
  d=d[d.candidate.isin(CAND)].copy();d['entity']=d[kind]
  allrows.append(d)
  for (cand,pair),a in d.groupby(['candidate','celltype_region']):cross_sets[(kind,cand,pair)]=set(a.entity)
  print('read',region,kind,len(d),flush=True)
 allrows=pd.concat(allrows);save(allrows,'cross_'+kind+'_assignments')
 assignments.append(allrows.assign(entity_kind=kind))
 for original in SPEC:
  spec=original.replace('Mic_immune','Mic_Immune');pair,m=spec.rsplit('_',1)
  sub=f'part7_TF/04tf_assignment/{pair}_TF_module_assignments.tsv' if kind=='tf' else f'part3_process_REGA_output/09module_assignments/{pair}/{"genes" if kind=="gene" else "RE"}_module_assignments.tsv'
  d=read(single/sub,'celltype x region x '+kind+' hard assignment');d=d[d.assigned_module==m].copy();single_sets[(kind,spec)]=set(d.entity)
  for (k,cand,cpair),cset in cross_sets.items():
   if k!=kind or cpair!=pair:continue
   sset=set(d.entity);shared=cset&sset
   overlaps.append(dict(entity_kind=kind,candidate=cand,single_candidate=spec,request_label=original,celltype_region=pair,n_cross=len(cset),n_single=len(sset),n_shared=len(shared),jaccard=len(shared)/len(cset|sset) if cset|sset else np.nan,shared_entities=';'.join(sorted(shared)),single_source_path=str(BASE/sub) if False else str(single/sub),cross_source_path=allrows.loc[(allrows.candidate==cand)&(allrows.celltype_region==pair),'source_path'].iloc[0]))
 for cand in CAND:
  r=cand.rsplit('_',1)[0]
  for ct in ['Exc','Inh','Ast','Oli','OPC','Mic_Immune']:
   pair=ct+'_'+r;counts.append(dict(candidate=cand,celltype_region=pair,entity_kind=kind,n_assigned=len(cross_sets.get((kind,cand,pair),set()))))
save(pd.DataFrame(overlaps),'single_cross_overlaps');save(pd.DataFrame(counts),'assignment_counts')
ta=assignments[-1].merge(activity[activity.kind=='cross'],left_on=['celltype_region','entity'],right_on=['celltype_region','TF'],how='left',suffixes=('_assignment','_activity'))
ta['positive_MW_FDR05']=(ta.mannwhitney_lateAD_vs_nonAD_fdr<.05)&(ta.lateAD_minus_nonAD>0)
ta['positive_LM_FDR05']=(ta.lm_fdr_pathology<.05)&(ta.lm_beta_pathology>0)
save(ta,'assigned_TF_activity')
# Preserve both supplied overlap tables; audit against current assignment/activity joins.
for kind,base,rel in [('cross',cross,'part4_TF/06driverTF_moduleTF_overlap_FC/mannwhitney_fdr_0.05/positive_AD_TF_module_overlap_details_fdr_lt_0.05.tsv'),('single',single,'part7_TF/05DriverTF_moduleTF_overlap_FC/mannwhitney_fdr_0.05/positive_AD_TF_module_overlap_details_fdr_lt_0.05.tsv')]:
 d=read(base/rel,'celltype x region x module x positive MW-associated TF');d['candidate']=d.region+'_'+d.module;save(d,'supplied_'+kind+'_TF_overlap')
 if kind=='cross':
  a=set(zip(ta.loc[ta.positive_MW_FDR05,'candidate'],ta.loc[ta.positive_MW_FDR05,'celltype_region'],ta.loc[ta.positive_MW_FDR05,'tf']))
  dd=d[d.candidate.isin(CAND)];b=set(zip(dd.candidate,dd.celltype_region,dd.tf));checks.append(dict(check='supplied overlap vs current cross assignments and MW activity',passed=a==b,max_abs_error=len(a^b)))
# Source GSEA summaries contain the original nominal p-value and GSEA FDR.
gsea=[]
for kind,base,stage in [('cross',cross,'part6_gsea'),('single',single,'part4_gsea_analysis')]:
 for lib in ['GOBP','KEGG']:
  folder=f'GSEA_from_WBV_{lib}2026'+('_all_celltype_region' if kind=='cross' else '')
  fn=f'region_wise_{lib}_positive_p005_hits.tsv' if kind=='cross' else ('GO' if lib=='GOBP' else lib)+'_positive_p005_hits.tsv'
  h=read(base/stage/folder/fn,'celltype x region x module x gene set');h['candidate']=h.celltype_region.str.rsplit('_',n=1).str[-1]+'_'+h.module
  if kind=='cross':h=h[h.candidate.isin(CAND)]
  else:h=h[(h.celltype_region+'_'+h.module).isin([s.replace('Mic_immune','Mic_Immune') for s in SPEC])]
  h['library']=lib;h['kind']=kind
  gsea.append(h)
save(pd.concat(gsea),'gsea_hits')
# Donor CellChat links: preserve actual sender/receiver membership, including complexes.
ccroot=ROOT/'experiments/EXP002_results/tables'
tr=read(ccroot/'donor_pathology_progression_tests.tsv','region x sender x receiver x ligand-receptor donor trend')
pw=read(ccroot/'donor_pathology_pairwise_tests.tsv','region x sender x receiver x ligand-receptor x donor contrast')
bhcheck(tr,['region'],'p_value','q_value','CellChat donor Spearman BH within region')
bhcheck(pw,['region','comparison'],'p_value','q_value','CellChat donor MW BH within region/contrast')
save(tr.groupby('region').agg(n_rows=('p_value','size'),n_finite=('p_value','count'),min_q=('q_value','min'),n_fdr05=('q_value',lambda v:(v<.05).sum())).reset_index(),'cellchat_family_summary')
links=[]
for r in ['AG','EC','MTC','PFC']:
 meta=read(ccroot/(r+'_candidate_interactions.tsv'),'ligand-receptor candidate metadata')
 d=tr[tr.region==r].merge(meta,on='interaction_id',suffixes=('_trend','_metadata'))
 late=pw[(pw.region==r)&(pw.comparison=='lateAD_minus_nonAD')]
 d=d.merge(late,on=['region','source','target','interaction_id'],suffixes=('','_pairwise'))
 for cand in [c for c in CAND if c.startswith(r+'_')]:
  for _,row in d.iterrows():
   ls=set(re.findall(r'[A-Za-z0-9.-]+',str(row['ligand.symbol']))) if pd.notna(row['ligand.symbol']) else {row.ligand}
   rs=set(re.findall(r'[A-Za-z0-9.-]+',str(row['receptor.symbol']))) if pd.notna(row['receptor.symbol']) else {row.receptor}
   lg=ls&cross_sets.get(('gene',cand,row.source+'_'+r),set());rg=rs&cross_sets.get(('gene',cand,row.target+'_'+r),set())
   if not lg and not rg:continue
   rec=row.to_dict();rec.update(candidate=cand,ligand_module_genes=';'.join(sorted(lg)),receptor_module_genes=';'.join(sorted(rg)),all_ligand_components_assigned=bool(ls and ls<=lg),all_receptor_components_assigned=bool(rs and rs<=rg))
   links.append(rec)
 links_df=pd.DataFrame(links);print('CellChat region',r,'linked rows cumulative',len(links),flush=True)
links_df['both_sides_assigned']=links_df.all_ligand_components_assigned&links_df.all_receptor_components_assigned
save(links_df,'cellchat_module_links')
# Small request file for targeted RDS lookup. All assigned TFs passing either positive test.
queries=[]
for cand in CAND:
 d=links_df[(links_df.candidate==cand)&(links_df.source!=links_df.target)&(links_df.status=='tested')].sort_values(['q_value','p_value'])
 # top ten donor trends per module; exploratory ranking explicitly preserved.
 for _,a in d.head(10).iterrows():
  for pair,genes in [(a.source+'_'+a.region,a.ligand_module_genes),(a.target+'_'+a.region,a.receptor_module_genes)]:
   ts=ta[(ta.candidate==cand)&(ta.celltype_region==pair)&(ta.positive_MW_FDR05|ta.positive_LM_FDR05)]
   for gene in str(genes).split(';'):
    if not gene:continue
    for tf in ts.tf:queries.append(dict(candidate=cand,celltype_region=pair,TF=tf,TG=gene,interaction_id=a.interaction_id,sender=a.source,receiver=a.target))
 # Also inspect selected module gene connections even without CellChat support.
 for pair,ts in ta[(ta.candidate==cand)&(ta.positive_MW_FDR05|ta.positive_LM_FDR05)].groupby('celltype_region'):
  genes=sorted(cross_sets.get(('gene',cand,pair),set()))
  for tf in ts.tf:
   for gene in genes:queries.append(dict(candidate=cand,celltype_region=pair,TF=tf,TG=gene,interaction_id='',sender='',receiver=''))
save(pd.DataFrame(queries).drop_duplicates(),'grn_queries')
# MAGMA audit: notebook V9 is P in the prefiltered 9-column file, original merged file has 10 columns.
sig=read(single/'part7_TF/sig_Genes_pval.tsv','prefiltered MAGMA gene',header=None,names=['GENE','CHR','START','STOP','NSNPS','NPARAM','N','ZSTAT','P'])
mg=read(single/'part5_magma/03merged_output/AD_Bellenguez_hg19_merged.genes.out','MAGMA gene',sep=r'\s+')
with (single/'part7_TF/motif_by_TF_matrix.tsv').open() as f:tfs=f.readline().strip().split('\t')[1:]
predefined='FOSB KLF16 ZFP3 SPI1 MAF ETV4 ZKSCAN4 CCDC6 ZSCAN26 ZNF691 FOXA3 ZKSCAN3 SIX5 SPDEF'.split()
sigtf=sig[sig.GENE.isin(tfs)].copy();sigtf['passes_notebook_threshold']=sigtf.P<.000515;sigtf['predefined_set']=sigtf.GENE.isin(predefined);sigtf['bonferroni_97_for_audit']=np.minimum(sigtf.P*97,1);save(sigtf,'magma_TF_prefilter_audit');save(mg[mg.GENE.isin(predefined)],'magma_predefined_TFs_original')
print('MAGMA',len(sig),len(sigtf),'threshold hits',sigtf[sigtf.passes_notebook_threshold].GENE.tolist(),flush=True)
save(pd.DataFrame(manifest).drop_duplicates('path'),'source_manifest');save(pd.DataFrame(checks),'validation_checks')
print(pd.DataFrame(checks).to_string(index=False),flush=True)
