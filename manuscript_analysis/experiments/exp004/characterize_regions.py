"""Characterize existing full GSEA reports and aging tables; no new GSEA fits."""
from pathlib import Path
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[2];B=R.parent;T=R/'experiments/EXP004_results/regional/tables';T.mkdir(parents=True,exist_ok=True)
CAND=['AG_M13','EC_M14','EC_M17','MTC_M6','MTC_M13','PFC_M14','PFC_M18'];CT=['Exc','Inh','Ast','Oli','OPC','Mic_Immune']
SPEC='Exc_AG_M8 Inh_AG_M8 Ast_AG_M2 Oli_AG_M6 Vasc_AG_M11 Exc_EC_M13 Ast_EC_M17 Oli_EC_M14 Oli_EC_M16 Mic_Immune_EC_M16 Vasc_EC_M11 Exc_MTC_M8 Exc_MTC_M19 Ast_MTC_M8 Oli_MTC_M8 Oli_MTC_M13 OPC_MTC_M8 Mic_Immune_MTC_M8 Ast_HC_M13 OPC_HC_M13 Exc_PFC_M1 OPC_PFC_M14 Vasc_PFC_M20'.split()
manifest=[];rows=[]
for kind,ids,base in [('cross',CAND,B/'AD_REGA_region_wise/part6_gsea'),('single',SPEC,B/'AD_REGA/part4_gsea_analysis')]:
 for cid in ids:
  if kind=='cross':region,module=cid.rsplit('_',1);pairs=[c+'_'+region for c in CT]
  else:pair,module=cid.rsplit('_',1);pairs=[pair]
  for pair in pairs:
   for lib in ['GOBP','KEGG']:
    folder=f'GSEA_from_WBV_{lib}2026'+('_all_celltype_region' if kind=='cross' else '')
    p=base/folder/pair/module/'gseapy.gene_set.prerank.report.csv'
    manifest.append(dict(kind=kind,module_id=cid,context=pair,library=lib,path=str(p),exists=p.exists()))
    if not p.exists():continue
    d=pd.read_csv(p);d['source_row']=np.arange(2,len(d)+2);d['source_path']=str(p);d['kind']=kind;d['module_id']=cid;d['context']=pair;d['celltype']=pair.rsplit('_',1)[0];d['region']=pair.rsplit('_',1)[1];d['module']=module;d['library']=lib
    d=d.rename(columns={'Term':'term','NES':'nes','NOM p-val':'pval','FDR q-val':'fdr','Lead_genes':'leading_genes'})
    d['positive_fdr05']=(d.nes>0)&(d.fdr<.05);rows.append(d)
g=pd.concat(rows,ignore_index=True);g.to_csv(T/'full_GSEA.tsv.gz',sep='\t',index=False)
sig=g[g.positive_fdr05].sort_values(['fdr','nes'],ascending=[True,False]);sig.groupby(['kind','module_id','context','library']).head(12).to_csv(T/'leading_GSEA_by_context.tsv',sep='\t',index=False)
terms=['RIBOSOME','OXIDATIVE PHOSPHORYLATION','ALZHEIMER DISEASE','PARKINSON DISEASE','PATHWAYS OF NEURODEGENERATION - MULTIPLE DISEASES','ANTIGEN PROCESSING AND PRESENTATION','COMPLEMENT AND COAGULATION CASCADES','CYTOKINE-CYTOKINE RECEPTOR INTERACTION','PHAGOSOME','ECM-RECEPTOR INTERACTION','CELL ADHESION MOLECULES','NEUROACTIVE LIGAND-RECEPTOR INTERACTION','GLYCOLYSIS / GLUCONEOGENESIS']
g['term_key']=g.term.str.strip().str.upper();panel=g[(g.library=='KEGG')&g.term_key.isin(terms)].copy();panel.to_csv(T/'pathway_panel_exact_terms.tsv',sep='\t',index=False)
# Cross/single same-celltype pathway matches; biological similarity is separate from hard-set overlap.
co=pd.read_csv(R/'experiments/EXP004_results/tables/single_cross_overlaps.tsv',sep='\t');co=co[(co.candidate.isin(CAND))&(co.entity_kind=='gene')]
shared=[]
for _,o in co.iterrows():
 a=sig[(sig.kind=='cross')&(sig.module_id==o.candidate)&(sig.context==o.celltype_region)]
 b=sig[(sig.kind=='single')&(sig.module_id==o.single_candidate)&(sig.context==o.celltype_region)]
 d=a.merge(b,on=['context','library','term'],suffixes=('_cross','_single'))
 if d.empty:continue
 d['n_hard_shared_genes']=o.n_shared;d['hard_gene_jaccard']=o.jaccard
 d['shared_leading_genes']=[ ';'.join(sorted(set(str(x).split(';'))&set(str(y).split(';')))) for x,y in zip(d.leading_genes_cross,d.leading_genes_single)]
 shared.append(d)
pd.concat(shared).to_csv(T/'cross_single_shared_pathways.tsv',sep='\t',index=False)
# Aging: preserve negative NES, zero permutation estimates, and source rows.
ages=[];checks=[]
for p in sorted((B/'AD_REGA_region_wise/part10_Age_related/GSEA_all_celltype_region').glob('*/*_aging_GSEA.tsv')):
 d=pd.read_csv(p,sep='\t');d['source_path']=str(p);d['source_row']=np.arange(2,len(d)+2);d['module_id']=d.context.str.rsplit('_',n=1).str[-1]+'_'+d.module
 ok=np.isfinite(d.pval);q=multipletests(d.loc[ok,'pval'],method='fdr_bh')[1];err=np.max(np.abs(q-d.loc[ok,'FDR']))
 checks.append(dict(check='aging BH over 20 modules: '+d.context.iloc[0],passed=len(d)==20 and err<1e-10,max_error=err));ages.append(d)
age=pd.concat(ages);age.to_csv(T/'aging_all_contexts.tsv',sep='\t',index=False);age[age.module_id.isin(CAND)].to_csv(T/'aging_requested_modules.tsv',sep='\t',index=False)
pd.DataFrame(manifest).to_csv(T/'GSEA_source_manifest.tsv',sep='\t',index=False)
checks.append(dict(check='all full GSEA sources available',passed=all(x['exists'] for x in manifest),max_error=sum(not x['exists'] for x in manifest)))
pd.DataFrame(checks).to_csv(T/'regional_checks.tsv',sep='\t',index=False)
print('Full GSEA rows',len(g),'sources',len(manifest),'aging contexts',len(ages),flush=True)
for cid in CAND:
 print('\nMODULE',cid)
 a=panel[(panel.kind=='cross')&(panel.module_id==cid)];print(a.pivot(index='term_key',columns='celltype',values='fdr').reindex(columns=CT).round(4).to_string())
 for cell in CT:
  a=sig[(sig.kind=='cross')&(sig.module_id==cid)&(sig.celltype==cell)]
  print(cell, '; '.join(f'{z.term} [NES {z.nes:.2f}; q {z.fdr:.3g}]' for z in a[a.library=='GOBP'].head(4).itertuples()))
# Connect existing positive GRN weights to genes in significant pathway leading edges.
# This reuses the previously extracted matrices; it does not infer new edges.
edges=pd.read_csv(R/'experiments/EXP004_results/tables/grn_edges.tsv',sep='\t')
linked=[]
for (cid,ctx),e in edges[edges.candidate.isin(CAND)].groupby(['candidate','celltype_region']):
 e=e.drop_duplicates(['TF','TG']);e=e[e.edge_value>0]
 for z in sig[(sig.kind=='cross')&(sig.module_id==cid)&(sig.context==ctx)].itertuples():
  for x in e[e.TG.isin(str(z.leading_genes).split(';'))].itertuples():
   linked.append(dict(module_id=cid,context=ctx,TF=x.TF,TG=x.TG,edge_value=x.edge_value,target_percentile=x.within_TF_target_percentile,term=z.term,library=z.library,nes=z.nes,fdr=z.fdr,GSEA_source=z.source_path,GSEA_row=z.source_row,GRN_source=x.matrix_source,GRN_TG_row=x.TG_row_index,GRN_TF_column=x.TF_column_index))
pd.DataFrame(linked).to_csv(T/'TF_edges_to_pathway_leading_genes.tsv',sep='\t',index=False)
