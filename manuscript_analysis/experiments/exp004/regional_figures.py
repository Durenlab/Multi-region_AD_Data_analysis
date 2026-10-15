"""Draw descriptive panels from existing GSEA, aging and LDSC results."""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[2];D=R/'experiments/EXP004_results';T=D/'regional/tables';F=D/'regional/figures';F.mkdir(exist_ok=True)
g=pd.read_csv(T/'full_GSEA.tsv.gz',sep='\t');g['key']=g.term.str.strip().str.upper()
a=pd.read_csv(T/'aging_requested_modules.tsv',sep='\t');l=pd.read_csv(D/'tables/ldsc_all.tsv',sep='\t')
C=['AG_M13','EC_M14','EC_M17','MTC_M6','MTC_M13','PFC_M14','PFC_M18'];CT=['Exc','Inh','Ast','Oli','OPC','Mic_Immune'];short=['Exc','Inh','Ast','Oli','OPC','Mic/Imm']
plt.get_cmap('RdBu_r').set_bad('#dddddd')
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
terms=[('RIBOSOME','Ribosome [K]'),('OXIDATIVE PHOSPHORYLATION','Oxidative phosphorylation [K]'),('ALZHEIMER DISEASE','Alzheimer disease [K]'),('PATHWAYS OF NEURODEGENERATION - MULTIPLE DISEASES','Neurodegeneration, multiple diseases [K]'),('GLYCOLYTIC PROCESS (GO:0006096)','Glycolytic process [GO]'),('TRIGLYCERIDE CATABOLIC PROCESS (GO:0019433)','Triglyceride catabolic process [GO]'),('ANTIGEN PROCESSING AND PRESENTATION','Antigen processing / presentation [K]'),('MHC CLASS II PROTEIN COMPLEX ASSEMBLY (GO:0002399)','MHC class II complex assembly [GO]'),('COMPLEMENT AND COAGULATION CASCADES','Complement / coagulation [K]'),('CYTOKINE-CYTOKINE RECEPTOR INTERACTION','Cytokine–receptor interaction [K]'),('ECM-RECEPTOR INTERACTION','ECM–receptor interaction [K]'),('EXTRACELLULAR MATRIX ASSEMBLY (GO:0085029)','Extracellular matrix assembly [GO]'),('NEUROACTIVE LIGAND-RECEPTOR INTERACTION','Neuroactive ligand–receptor interaction [K]')]
figure_rows=[]
def row(cid,ctx,key,kind='cross'):
 z=g[(g.kind==kind)&(g.module_id==cid)&(g.context==ctx)&(g.key==key)]
 if len(z)>1:raise ValueError((cid,ctx,key,len(z)))
 if z.empty:return pd.Series(dict(kind=kind,module_id=cid,context=ctx,term=key,nes=np.nan,fdr=np.nan,source_path='',source_row=np.nan,availability='not scored in full source report'))
 return z.iloc[0]
def save(fig,name):
 for ext in ['pdf','png','svg']:fig.savefig(F/f'{name}.{ext}',dpi=180,bbox_inches='tight',facecolor='white')
 shutil.copy2(F/f'{name}.pdf',R/'manuscript'/f'EXP004_{name}.pdf');plt.close(fig)
def draw_matrix(ax,mods,subset,title):
 vals=[];qs=[]
 for key,label in subset:
  zs=[row(m,c+'_'+m.rsplit('_',1)[0],key) for m in mods for c in CT]
  vals.append([z.nes for z in zs]);qs.append([z.fdr for z in zs])
  for z in zs:figure_rows.append(dict(z,figure_panel=title,display_label=label))
 vals=np.asarray(vals);qs=np.asarray(qs)
 cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
 im=ax.imshow(vals,cmap=cmap,vmin=-3.5,vmax=3.5,aspect='auto')
 for i,j in np.argwhere(~np.isfinite(vals)):ax.text(j,i,'—',ha='center',va='center',fontsize=8)
 for i,j in np.argwhere(qs<.05):ax.text(j,i,'*',ha='center',va='center',fontsize=12,color='black' if abs(vals[i,j])<2.7 else 'white')
 ax.set_yticks(range(len(subset)),[x[1] for x in subset]);ax.set_xticks(range(len(mods)*6),short*len(mods),rotation=55,ha='right');ax.tick_params(length=0)
 for j in range(6,len(mods)*6,6):ax.axvline(j-.5,color='white',lw=3)
 for i,m in enumerate(mods):ax.text(i*6+2.5,-1.0,m,ha='center',fontweight='bold',fontsize=11)
 ax.set_title(title,loc='left',pad=40,fontweight='bold')
 return im
fig=plt.figure(figsize=(24,10.8));gs=fig.add_gridspec(3,1,height_ratios=[13,1.5,1.5],hspace=.7)
ax=fig.add_subplot(gs[0]);im=draw_matrix(ax,C,terms,'a  Pathway architecture across all seven pathology-selected modules')
fig.colorbar(im,ax=ax,fraction=.016,pad=.015,label='Signed pathway NES (color capped at ±3.5)',extend='both')
for n,typ in enumerate(['aging','ldsc'],1):
 ax=fig.add_subplot(gs[n]);zs=[]
 for m in C:
  for c in CT:
   ctx=c+'_'+m.rsplit('_',1)[0]
   z=(a[(a.module_id==m)&(a.context==ctx)].iloc[0] if typ=='aging' else l[(l.kind=='cross')&(l.candidate==m)&(l.celltype_region==ctx)].iloc[0]);zs.append(z)
 vals=[z.NES if typ=='aging' else z.Enrichment for z in zs];qs=[z.FDR if typ=='aging' else z.FDR_BH_Enrichment_p for z in zs]
 im=ax.imshow([vals],aspect='auto',cmap='RdBu_r',vmin=-3 if typ=='aging' else -15,vmax=3 if typ=='aging' else 15)
 for j,(v,q,z) in enumerate(zip(vals,qs,zs)):
  mark='*' if q<.05 else ('†' if typ=='ldsc' and z.Enrichment_p<.05 else '')
  ax.text(j,0,f'{v:.1f}{mark}',ha='center',va='center',fontsize=8,color='white' if abs(v)>(2.1 if typ=='aging' else 11) else 'black')
 ax.set_yticks([0],['Aging NES' if typ=='aging' else 'LDSC enrichment']);ax.set_xticks(range(len(C)*6),short*len(C),rotation=55,ha='right');ax.tick_params(length=0)
 for j in np.arange(5.5,len(C)*6-1,6):ax.axvline(j,color='white',lw=3)
 ax.set_title('b  Age-up gene sets: signed NES' if typ=='aging' else 'c  AD heritability: enrichment estimate',loc='left',fontweight='bold',pad=10)
 fig.colorbar(im,ax=ax,fraction=.016,pad=.015,label='NES' if typ=='aging' else 'Enrichment')
fig.text(.12,-.005,'* FDR < 0.05 within the source analysis family; † LDSC nominal P < 0.05 only. K = KEGG; GO = GO Biological Process.\nGray / dash = term not scored; otherwise all values are shown, including nonsignificant and negative NES. Rows were selected for regional interpretation; this is not an exhaustive pathway screen.',fontsize=9)
save(fig,'01_regional_patterns')
# EC detail: identical term order and color scales, with exact q values printed.
subset=[terms[i] for i in [0,1,2,6,7,8,9,10,12]]
fig,axes=plt.subplots(1,2,figsize=(14,6.7),sharey=True)
for ax,m in zip(axes,['EC_M14','EC_M17']):
 zs=[[row(m,c+'_EC',k) for c in CT] for k,_ in subset];vals=np.array([[z.nes for z in rr] for rr in zs]);im=ax.imshow(vals,cmap='RdBu_r',vmin=-3.5,vmax=3.5,aspect='auto')
 for i,rr in enumerate(zs):
  for j,z in enumerate(rr):
   label=(f'{z.fdr:.3f}'+('*' if z.fdr<.05 else '')) if np.isfinite(z.fdr) else '—'
   ax.text(j,i,label,ha='center',va='center',fontsize=8,color='white' if abs(z.nes)>2.7 else 'black')
 ax.set_xticks(range(6),short);ax.set_yticks(range(len(subset)),[x[1] for x in subset]);ax.set_title(m,fontweight='bold',pad=12);ax.tick_params(length=0)
fig.subplots_adjust(wspace=.1,right=.88,bottom=.16);cax=fig.add_axes([.90,.19,.014,.66]);fig.colorbar(im,cax=cax,label='Signed NES',extend='both')
fig.text(.125,.05,'Numbers = GSEA FDR; * FDR < 0.05. Colors = signed NES; dash = term not scored.\nEC_M14: Oli/OPC bioenergetic components within a broader immune / interaction program.\nEC_M17: glial antigen presentation with weaker genetic evidence; see regional panel for opposite Inh aging direction.',fontsize=10)
save(fig,'02_EC_M14_vs_M17')
# Same-context comparisons deliberately include selected nonsignificant counterparts.
comparisons=[('AG_M13','Exc_AG_M8','RIBOSOME'),('AG_M13','Exc_AG_M8','OXIDATIVE PHOSPHORYLATION'),('AG_M13','Inh_AG_M8','OXIDATIVE PHOSPHORYLATION'),('AG_M13','Ast_AG_M2','EXTRACELLULAR MATRIX ASSEMBLY (GO:0085029)'),('AG_M13','Oli_AG_M6','OXIDATIVE PHOSPHORYLATION'),('EC_M14','Ast_EC_M17','COMPLEMENT AND COAGULATION CASCADES'),('EC_M14','Oli_EC_M16','ANTIGEN PROCESSING AND PRESENTATION'),('EC_M14','Mic_Immune_EC_M16','ECM-RECEPTOR INTERACTION'),('EC_M17','Oli_EC_M14','ANTIGEN PROCESSING AND PRESENTATION'),('EC_M17','Ast_EC_M17','COMPLEMENT AND COAGULATION CASCADES'),('MTC_M13','Exc_MTC_M19','GLYCOLYTIC PROCESS (GO:0006096)'),('MTC_M13','Exc_MTC_M19','GLYCOLYSIS / GLUCONEOGENESIS'),('MTC_M13','Oli_MTC_M8','NEUROACTIVE LIGAND-RECEPTOR INTERACTION'),('MTC_M13','OPC_MTC_M8','OXIDATIVE PHOSPHORYLATION'),('MTC_M6','Oli_MTC_M8','NEUROACTIVE LIGAND-RECEPTOR INTERACTION'),('PFC_M14','Exc_PFC_M1','CYTOPLASMIC TRANSLATION (GO:0002181)'),('PFC_M18','Exc_PFC_M1','CYTOPLASMIC TRANSLATION (GO:0002181)')]
ov=pd.read_csv(D/'tables/single_cross_overlaps.tsv',sep='\t');out=[];labels=[]
comparisons.sort(key=lambda x:C.index(x[0]))
fig,ax=plt.subplots(figsize=(12,11))
for i,(m,s,key) in enumerate(comparisons):
 ctx=s.rsplit('_',1)[0];x=row(m,ctx,key);y=row(s,ctx,key,'single');o=ov[(ov.candidate==m)&(ov.single_candidate==s)&(ov.entity_kind=='gene')].iloc[0]
 ax.plot([x.nes,y.nes],[i,i],color='#acacac',lw=1.5)
 for z,col,marker in [(x,'#197c80','o'),(y,'#ba5932','s')]:ax.scatter(z.nes,i,s=58,marker=marker,edgecolor=col,facecolor=col if z.fdr<.05 else 'white',zorder=3)
 term=dict(terms).get(key,{'GLYCOLYSIS / GLUCONEOGENESIS':'Glycolysis / gluconeogenesis [K]','CYTOPLASMIC TRANSLATION (GO:0002181)':'Cytoplasmic translation [GO]'}.get(key,key));labels.append(f'{m} ↔ {s}\n{term}')
 ax.text(1.02,i,f'{int(o.n_shared):3d}   /   {o.jaccard:.3f}',transform=ax.get_yaxis_transform(),va='center',fontfamily='monospace',fontsize=9)
 out.append(dict(cross_module=m,single_module=s,context=ctx,term=x.term,cross_NES=x.nes,cross_FDR=x.fdr,single_NES=y.nes,single_FDR=y.fdr,shared_hard_genes=o.n_shared,gene_jaccard=o.jaccard,cross_source=x.source_path,cross_row=x.source_row,single_source=y.source_path,single_row=y.source_row))
ax.set_yticks(range(len(labels)),labels,fontsize=9);ax.invert_yaxis();ax.set_xlabel('Signed pathway NES');ax.axvline(0,color='#bbb',lw=.8);ax.grid(axis='x',alpha=.2);ax.set_title('Cross-cell programs and independently fitted cell-specific modules',loc='left',fontweight='bold',pad=65)
ax.text(1.02,1.005,'Shared genes / Jaccard',transform=ax.transAxes,fontsize=9)
ax.legend(handles=[Line2D([],[],marker='o',color='#197c80',ls='',label='Cross-cell module'),Line2D([],[],marker='s',color='#ba5932',ls='',label='Cell-specific module'),Line2D([],[],marker='o',mfc='white',color='#555',ls='',label='Hollow: FDR ≥ 0.05')],loc='lower left',bbox_to_anchor=(-.02,1.035),ncol=3,frameon=False)
fig.text(.125,.015,'Gene overlap refers to hard assignments across the two fits. Shared pathway enrichment does not imply module equivalence.\nDisplayed comparisons were chosen to show both concordance and divergence; significance is evaluated separately in each source GSEA.',fontsize=9)
fig.subplots_adjust(left=.49,right=.85,bottom=.10,top=.90)
save(fig,'03_cross_vs_cell_specific')
pd.DataFrame(out).to_csv(T/'figure_cross_single_comparisons.tsv',sep='\t',index=False)
pd.DataFrame(figure_rows).to_csv(T/'figure_pathway_values.tsv',sep='\t',index=False)
print('Saved three regional figures in PDF, PNG and SVG; manuscript PDF copies; source tables.')
