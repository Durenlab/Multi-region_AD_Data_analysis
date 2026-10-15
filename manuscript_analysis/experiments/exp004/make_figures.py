"""Publication-size summaries of existing EXP004 statistics; no new tests."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
R=Path(__file__).resolve().parents[2];T=R/'experiments/EXP004_results/tables';F=R/'experiments/EXP004_results/figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
order=['EC_M14','EC_M17','MTC_M13','MTC_M6','PFC_M14','AG_M13','PFC_M18'];cts=['Exc','Inh','Ast','Oli','OPC','Mic_Immune']
mp=pd.read_csv(T/'module_pathology.tsv',sep='\t');ld=pd.read_csv(T/'ldsc_all.tsv',sep='\t');ct=pd.read_csv(T/'assignment_counts.tsv',sep='\t');tf=pd.read_csv(T/'assigned_TF_activity.tsv',sep='\t')
def save(fig,name):
 for ext in ['pdf','png','svg']:fig.savefig(F/(name+'.'+ext),bbox_inches='tight')
 plt.close(fig)
def heat(ax,values,labels,rows,title,cmap='Blues',fmt=None,vmax=None):
 a=np.array(values,float);cm=plt.get_cmap(cmap).copy();cm.set_bad('#dddddd')
 im=ax.imshow(a,aspect='auto',cmap=cm,vmin=0,vmax=vmax)
 ax.set_xticks(range(len(labels)),labels,rotation=45,ha='right');ax.set_yticks(range(len(rows)),rows);ax.set_title(title,loc='left',fontweight='bold',pad=12)
 for i in range(a.shape[0]):
  for j in range(a.shape[1]):
   if np.isfinite(a[i,j]):
    text=fmt(i,j,a[i,j]) if fmt else f'{a[i,j]:.2g}';color='white' if a[i,j]>.62*(vmax or np.nanmax(a)) else '#222222';ax.text(j,i,text,ha='center',va='center',color=color,fontsize=9)
   else:ax.text(j,i,'NA',ha='center',va='center',fontsize=8,color='#555555')
 return im
fig,axs=plt.subplots(1,4,figsize=(16,5.7),gridspec_kw={'width_ratios':[1.2,2.4,2.4,1.1]})
q=[]
for c in order:
 a=mp[mp.candidate==c];q.append([a[a.test=='Linear_Regression'].pval_adj.iloc[0],a[a.Comparison=='nonAD vs lateAD'].pval_adj.iloc[0]])
q=np.array(q);heat(axs[0],-np.log10(q),['Linear','Late vs nonAD'],order,'A  Pathology',fmt=lambda i,j,v:f'{q[i,j]:.3g}',vmax=4)
l=[]
for c in order:
 r,m=c.split('_');a=ld[(ld.kind=='cross')&(ld.candidate==c)].set_index('celltype_region');l.append([a.loc[x+'_'+r,'FDR_BH_Enrichment_p'] if x+'_'+r in a.index else np.nan for x in cts])
l=np.array(l);heat(axs[1],-np.log10(l),cts,order,'B  AD LDSC',fmt=lambda i,j,v:f'{l[i,j]:.3g}',vmax=4)
g=ct[ct.entity_kind=='gene'].copy();g['celltype']=g.celltype_region.str.rsplit('_',n=1).str[0];g=g.pivot(index='candidate',columns='celltype',values='n_assigned').reindex(index=order,columns=cts)
heat(axs[2],np.log10(g+1),cts,order,'C  Assigned genes',fmt=lambda i,j,v:str(int(g.iloc[i,j])),vmax=3)
b=tf.groupby('candidate')[['positive_LM_FDR05','positive_MW_FDR05']].sum().reindex(order,fill_value=0)
heat(axs[3],b.values,['Linear','Late vs nonAD'],order,'D  Positive TFs',fmt=lambda i,j,v:str(int(v)),vmax=26)
for ax in axs[1:]:ax.set_yticklabels([])
fig.suptitle('Seven candidate modules: evidence supports one focused EC_M14 story',fontsize=15,y=1.02)
fig.text(.04,-.03,'A–B: cell labels are supplied FDR values; darker = smaller FDR. C: cell labels are counts; color = log10(count + 1).\nD: assigned TFs with positive effect and FDR < 0.05, separately by test. PFC_M14 has no assigned TFs. Order is a qualitative integration rank.',fontsize=10)
fig.tight_layout();save(fig,'01_seven_module_evidence')
# Focused panel: overlap and regulators, with exact cell types.
fig,axs=plt.subplots(1,3,figsize=(16,5.5),gridspec_kw={'width_ratios':[1.15,1.3,1.5]})
a=ld[(ld.kind=='cross')&(ld.candidate=='EC_M14')].copy();a['ct']=a.celltype_region.str.replace('_EC','',regex=False);a=a.set_index('ct').reindex(cts)
y=np.arange(6);axs[0].errorbar(a.Enrichment,y,xerr=a.Enrichment_std_error,fmt='o',color='#3969AC',capsize=3);axs[0].axvline(1,color='gray',ls=':');axs[0].set_yticks(y,cts);axs[0].invert_yaxis();axs[0].set_xlabel('AD LDSC enrichment ± source SE');axs[0].set_title('A  Genetic support is cell-specific',loc='left',fontweight='bold')
for i,row in enumerate(a.itertuples()):axs[0].annotate(f'q={row.FDR_BH_Enrichment_p:.3g}',(row.Enrichment,i),xytext=(7,5),textcoords='offset points',fontsize=8)
o=pd.read_csv(T/'single_cross_overlaps.tsv',sep='\t');o=o[(o.candidate=='EC_M14')&o.single_candidate.isin(['Ast_EC_M17','Oli_EC_M16','Oli_EC_M14','Mic_Immune_EC_M16'])]
names=['Ast_EC_M17','Oli_EC_M16','Oli_EC_M14','Mic_Immune_EC_M16'];v=o.pivot(index='single_candidate',columns='entity_kind',values='jaccard').reindex(index=names,columns=['gene','RE','tf']);n=o.pivot(index='single_candidate',columns='entity_kind',values='n_shared').reindex(index=names,columns=['gene','RE','tf'])
heat(axs[1],v,['Genes','REs','TFs'],names,'B  Actual overlap',fmt=lambda i,j,v:f'{int(n.iloc[i,j])}\nJ={v:.3f}',vmax=.25)
select=[('Oli_EC','SPI1'),('Oli_EC','CEBPD'),('Ast_EC','AR'),('Ast_EC','PGR'),('Ast_EC','BRCA1'),('Exc_EC','IRF5'),('Mic_Immune_EC','HMGN3')]
rows=pd.DataFrame([tf[(tf.candidate=='EC_M14')&(tf.celltype_region==p)&(tf.tf==t)].iloc[0] for p,t in select]);y=np.arange(len(rows))
axs[2].scatter(-np.log10(rows.lm_fdr_pathology),y-.12,label='Linear FDR',s=45,color='#3969AC');axs[2].scatter(-np.log10(rows.mannwhitney_lateAD_vs_nonAD_fdr),y+.12,label='Late/nonAD FDR',s=45,color='#E69F00');axs[2].axvline(-np.log10(.05),ls=':',color='gray');axs[2].set_yticks(y,[p.replace('_EC','')+' / '+t for p,t in select]);axs[2].invert_yaxis();axs[2].set_xlabel('−log10(FDR)');axs[2].set_title('C  Assigned TF activity',loc='left',fontweight='bold');axs[2].legend(loc='lower right',fontsize=8)
fig.suptitle('EC_M14: convergent evidence, with regulator identity limits',fontsize=15,y=1.04)
fig.text(.03,-.04,'B: labels show intersection size and Jaccard; overlaps are descriptive, without a new enrichment test.\nC: all displayed activity effects are positive. AR/PGR/BRCA1 share one retained ZBTB33-labelled motif; these are not three independent signals.',fontsize=10)
fig.tight_layout(w_pad=3);save(fig,'02_EC_M14_support')
# Existing donor scores, plus explicit inferred network and evidence gaps.
d=pd.read_csv(T/'TNC_donor_scores.tsv',sep='\t');cc=pd.read_csv(T/'cellchat_module_links.tsv',sep='\t');fig=plt.figure(figsize=(13,8));gs=fig.add_gridspec(2,2,height_ratios=[1.15,1]);rng=np.random.default_rng(4004)
for j,inter in enumerate(['TNC_ITGA8_ITGB1','TNC_ITGA9_ITGB1']):
 ax=fig.add_subplot(gs[0,j]);a=d[d.interaction_id==inter];st=cc[(cc.candidate=='EC_M14')&(cc.source=='Ast')&(cc.target=='OPC')&(cc.interaction_id==inter)].iloc[0]
 for i,(group,color) in enumerate(zip(['nonAD','earlyAD','lateAD'],['#66AABB','#EEAA55','#CC6677'])):
  vals=a[a.pathology==group].communication_probability.to_numpy();ax.scatter(i+rng.uniform(-.12,.12,len(vals)),vals,s=25,color=color,alpha=.75,edgecolors='none');ax.plot([i-.19,i+.19],[np.median(vals)]*2,color='black',lw=2)
 ax.set_xticks([0,1,2],['nonAD\nn=25','earlyAD\nn=14','lateAD\nn=7']);ax.set_ylabel('Donor CellChat score');ax.set_title(f'{chr(65+j)}  Ast TNC → OPC '+inter.replace('TNC_','').replace('_',' + '),loc='left',fontweight='bold')
 ax.text(.02,.97,f'Trend ρ={st.spearman_rho:.3f}, p={st.p_value:.3g}, q={st.q_value:.2g}\nLate/nonAD p={st.p_value_pairwise:.3g}, q={st.q_value_pairwise:.3g}',transform=ax.transAxes,va='top',fontsize=10)
 ax.margins(y=.25)
ax=fig.add_subplot(gs[1,:]);ax.axis('off');ax.set_xlim(0,12);ax.set_ylim(0,3.2)
def box(x,y,w,h,text,color):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.10',facecolor=color,edgecolor='#444444'));ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=10)
box(.2,1.1,3,.95,'Ast_EC / EC_M14\nshared motif activity\nAR · PGR · BRCA1','#DCEEE5');box(4.35,1.1,2,.95,'Ast_EC / EC_M14\nTNC gene assignment','#DCEEE5');box(8.1,1.1,3.3,.95,'OPC_EC\nITGA8/ITGA9 + ITGB1\nnot assigned to EC_M14','#F3E5EF')
for x1,x2 in [(3.3,4.2),(6.5,7.95)]:ax.annotate('',xy=(x2,1.6),xytext=(x1,1.6),arrowprops={'arrowstyle':'->','ls':'--','color':'#555555','lw':1.5})
ax.text(3.8,.75,'GRN weights\nno binding test',ha='center',fontsize=9);ax.text(7.3,.75,'CellChat inference\nFDR not significant',ha='center',fontsize=9)
ax.text(.1,2.7,'C  A testable communication branch; the upstream regulator remains unresolved',fontsize=12,fontweight='bold')
ax.text(.1,.05,'SPI1/CEBPD → ATP6V0D2 is a separate Oli_EC GRN branch; no measured edge joins it to astrocyte TNC.',fontsize=10)
fig.suptitle('EC_M14 communication hypothesis: display the uncertainty',fontsize=15,y=.99);fig.tight_layout(rect=[0,0,1,.96]);save(fig,'03_TNC_communication_hypothesis')
print('Created 3 figures in PDF, PNG and SVG.')
