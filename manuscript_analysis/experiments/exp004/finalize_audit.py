from pathlib import Path
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];T=R/'experiments/EXP004_results/tables';B=R.parent
# Mirror the exact notebook's syntactic-name intersection, confirmed with R.
d=pd.read_csv(T/'magma_TF_prefilter_audit.tsv',sep='\t');# source_row correction handled below only for the initial extraction
d['in_notebook_97_family']=~d.GENE.isin(['NKX2-1','NKX2-8']);d['bonferroni_97_for_audit']=np.where(d.in_notebook_97_family,np.minimum(d.P*97,1),np.nan)
d.to_csv(T/'magma_TF_prefilter_audit.tsv',sep='\t',index=False,na_rep='NA')
# Read original single-context activity evidence for explicit cross-model comparison.
a=pd.read_csv(B/'AD_REGA/part7_TF/02TF_pathology_association_lm/all_celltype_region_TF_pathology_association.tsv',sep='\t');a['source_row']=np.arange(2,len(a)+2);a['source_path']=str(B/'AD_REGA/part7_TF/02TF_pathology_association_lm/all_celltype_region_TF_pathology_association.tsv')
keep=((a.celltype_region=='Oli_EC')&a.TF.isin(['SPI1','CEBPD']))|((a.celltype_region=='OPC_MTC')&(a.TF=='MAF'))
a[keep].to_csv(T/'single_TF_comparators.tsv',sep='\t',index=False)
# Exact pooled LR rows, interpreted descriptively and not as pathology tests.
p=B/'AD_REGA_region_wise/part8_cellchat/04_pathology_comparison/EC/tables/EC_all_conditions_ligand_receptor_communications.tsv';d=pd.read_csv(p,sep='\t');d['source_row']=np.arange(2,len(d)+2);d['source_path']=str(p)
x=d[(d.source=='Ast')&(d.target=='OPC')&d.interaction_name.isin(['TNC_ITGA8_ITGB1','TNC_ITGA9_ITGB1'])];x.to_csv(T/'pooled_TNC_context.tsv',sep='\t',index=False)
# Reconstruct the selected donor display with the published EXP002 eligibility and zero rules.
base=R/'experiments/EXP002_results';cp=base/'tables/EC_donor_celltype_coverage.tsv';sp=base/'intermediate/EC_donor_nonzero_probabilities.tsv.gz'
c=pd.read_csv(cp,sep='\t',dtype={'donor_id':str});s=pd.read_csv(sp,sep='\t',dtype={'donor_id':str});s['source_row']=np.arange(2,len(s)+2)
valid=c[(c.donor_network_status=='completed')&(c.n_cells>20)]
ids=set(valid[valid.celltype=='Ast'].donor_id)&set(valid[valid.celltype=='OPC'].donor_id)
donors=valid[valid.donor_id.isin(ids)][['donor_id','pathology']].drop_duplicates()
rows=[]
for interaction in ['TNC_ITGA8_ITGB1','TNC_ITGA9_ITGB1']:
 q=s[(s.source=='Ast')&(s.target=='OPC')&(s.interaction_id==interaction)]
 d=donors.merge(q[['donor_id','communication_probability','source_row']],on='donor_id',how='left',validate='one_to_one');d['score_present_in_sparse_file']=d.communication_probability.notna();d['communication_probability']=d.communication_probability.fillna(0);d['interaction_id']=interaction;d['source_path']=str(sp);d['coverage_source']=str(cp);rows.append(d)
pd.concat(rows).to_csv(T/'TNC_donor_scores.tsv',sep='\t',index=False,na_rep='NA')
print('TNC donor groups:',donors.groupby('pathology').size().to_dict());print('Pooled TNC rows',len(x));print(x[['condition','interaction_name','prob','pval']].to_string(index=False))
# Preserve all assignment-context counts and representative pathway rows for seven-module audit.
h=pd.read_csv(T/'gsea_hits_audited.tsv',sep='\t');h=h[h.kind=='cross'].sort_values(['fdr','nes'],ascending=[True,False]);h.groupby(['candidate','celltype_region','library']).head(3).to_csv(T/'leading_pathways_by_context.tsv',sep='\t',index=False)
# Compact family counts / source availability, including RDS metadata and original reports.
m=pd.read_csv(T/'source_manifest.tsv',sep='\t');extra=[]
for p in [cp,sp,B/'AD_REGA/part7_TF/TF_with_SNP.ipynb',B/'AD_Driver_TF/cell_type/GRN/Oli_EC_TG_TF.rds']:
 extra.append(dict(path=str(p),available=p.exists(),unit='supporting source; see report'))
for p in h.raw_report_path.dropna().unique():extra.append(dict(path=p,available=Path(p).exists(),unit='GSEA gene set within one module/context/library'))
md=pd.read_csv(T/'grn_metadata.tsv',sep='\t')
for col in ['matrix_source','metadata_source']:
 for p in md[col]:extra.append(dict(path=p,available=Path(p).exists(),unit='GRN target x TF weights / pair metadata'))
pd.concat([m,pd.DataFrame(extra)]).drop_duplicates('path').to_csv(T/'source_manifest.tsv',sep='\t',index=False)
# Compare the notebook input with the requested original MAGMA output, retaining missing genes.
a=pd.read_csv(T/'magma_TF_prefilter_audit.tsv',sep='\t');b=pd.read_csv(T/'magma_predefined_TFs_original.tsv',sep='\t')
a[a.predefined_set].merge(b,on='GENE',how='left',suffixes=('_notebook','_original')).to_csv(T/'magma_source_comparison.tsv',sep='\t',index=False,na_rep='NA')
from statsmodels.stats.multitest import multipletests
checks=[]
for fn in ['LinearRegression_all_celltype_region.tsv','Pairwise_MannWhitneyU_all_celltype_region.tsv']:
 p=B/'AD_REGA/part2_run_REGA/H_pathology_significant_test/module_significancy'/fn;d=pd.read_csv(p,sep='\t');errors=[]
 for _,a in d.groupby(['celltype_region','Comparison']):
  ok=np.isfinite(a.pval);q=multipletests(a.loc[ok,'pval'],method='fdr_bh')[1];errors.extend(abs(q-a.loc[ok,'pval_adj']))
 checks.append(dict(check='single module BH within celltype_region/contrast: '+fn,passed=max(errors)<1e-9,max_abs_error=max(errors)))
# Verify the exact finite counts and no missing source files among inventoried inputs.
m=pd.read_csv(T/'source_manifest.tsv',sep='\t');checks.append(dict(check='all inventoried source files available',passed=bool(m.available.all()),max_abs_error=int((~m.available).sum())))
h=pd.read_csv(T/'gsea_hits_audited.tsv',sep='\t');err=(h.fdr-h.raw_GSEA_FDR).abs().max();checks.append(dict(check='GSEA summary FDR matches raw report within 1e-7',passed=bool(err<1e-7 and h.raw_report_available.all()),max_abs_error=err))
x=pd.read_csv(T/'grn_edges.tsv',sep='\t');checks.append(dict(check='queried GRN identifiers found and finite weights',passed=bool(x.edge_value.notna().all()),max_abs_error=int(x.edge_value.isna().sum())))
c=pd.read_csv(T/'validation_checks.tsv',sep='\t');pd.concat([c,pd.DataFrame(checks)]).drop_duplicates('check',keep='last').to_csv(T/'validation_checks.tsv',sep='\t',index=False)
# Reproduce the motif-specificity evidence supporting the recommended TF discussion.
rows=[]
for pair,tfs in [('Ast_EC',['AR','PGR','BRCA1']),('Oli_EC',['SPI1','CEBPD'])]:
 base=B/'AD_REGA_region_wise/part4_TF/03TF_activity'/pair/'TF_activity_pcc_pos_only'
 pp=base/(pair+'_TF_motif_pcc_M2.csv');mm=base/(pair+'_TF_motif_binary_M1.csv')
 pcc=pd.read_csv(pp,index_col=0);m1=pd.read_csv(mm,index_col=0)
 for tf in tfs:
  weights=pcc.loc[tf]*m1.loc[tf]
  for motif,value in weights[weights!=0].items():rows.append(dict(celltype_region=pair,TF=tf,motif=motif,PCC_times_binary=value,PCC_source=str(pp),M1_source=str(mm)))
pd.DataFrame(rows).to_csv(T/'recommended_TF_motif_support.tsv',sep='\t',index=False)
