"""Figures for the focused EC comparison; all values read from saved tables."""
from pathlib import Path
import shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[2];O=R/'experiments/EXP004_results/ec_comparison';T=O/'tables';F=O/'figures'
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
read=lambda n:pd.read_csv(T/(n+'.tsv'),sep='\t')
P=read('pathway_leading_edge_comparison');W=read('EC_normalized_gene_weights');D=read('EC_donor_activities');C=read('donor_activity_correlations');S=read('subpopulation_activity_associations');N=read('donor_subpopulation_counts')
colors={'nonAD':'#6d9ab5','earlyAD':'#dea549','lateAD':'#b85465'}
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(F/f'{name}.{ext}',bbox_inches='tight',dpi=180,facecolor='white')
 shutil.copy2(F/f'{name}.pdf',R/'manuscript'/f'EXP004_EC_{name}.pdf');plt.close(fig)
pairs=[('Oli_EC','ANTIGEN PROCESSING AND PRESENTATION'),('Ast_EC','COMPLEMENT AND COAGULATION CASCADES'),('Ast_EC','CILIUM MOVEMENT (GO:0003341)'),('Mic_Immune_EC','ANTIGEN PROCESSING AND PRESENTATION'),('Exc_EC','CYTOKINE-CYTOKINE RECEPTOR INTERACTION'),('Oli_EC','OXIDATIVE PHOSPHORYLATION'),('OPC_EC','OXIDATIVE PHOSPHORYLATION')]
labels=['Oli / antigen presentation','Ast / complement','Ast / cilium movement','Mic/Imm / antigen presentation','Exc / cytokine–receptor','Oli / oxidative phosphorylation','OPC / oxidative phosphorylation']
fig,(ax,bx)=plt.subplots(1,2,figsize=(14,8.8),gridspec_kw={'width_ratios':[1.7,1]})
selected=[]
for i,((ctx,term),label) in enumerate(zip(pairs,labels)):
 z=P[(P.context==ctx)&(P.term==term)].iloc[0];selected.append(z);sizes=[z.n_leading_M14-z.n_shared,z.n_shared,z.n_leading_M17-z.n_shared];left=0
 for n,col in zip(sizes,['#228b8d','#b5b5b5','#ba6848']):
  ax.barh(i,n,left=left,height=.48,color=col)
  if n>=4:ax.text(left+n/2,i,str(n),ha='center',va='center',fontsize=9,color='white' if col!='#b5b5b5' else 'black')
  left+=n
 ax.text(0,i+.35,f'q: M14 {z.M14_fdr:.3g} / M17 {z.M17_fdr:.3g}; J = {z.leading_jaccard:.2f}',fontsize=8,color='#444')
ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis();ax.set_xlabel('Genes in the union of source leading edges');ax.set_xlim(0,175);ax.set_title('a  Shared pathways can contain shared and distinct leading genes',loc='left',pad=40,fontweight='bold',fontsize=10)
ax.legend(handles=[Line2D([],[],color=c,lw=7,label=l) for c,l in zip(['#228b8d','#b5b5b5','#ba6848'],['M14 only','Shared','M17 only'])],loc='lower left',bbox_to_anchor=(0,1.01),ncol=3,frameon=False,fontsize=8)
pg=[('Oli_EC',g) for g in ['CD74','CIITA','CTSS','HLA-B','HLA-F','HLA-DMA','HLA-DOB','ATP6V0D2']]+[('Mic_Immune_EC',g) for g in ['CD74','CIITA','HLA-DMB','CTSS','NPY1R','MT1E','MT1G']]+[('Ast_EC',g) for g in ['CFI','C3','DNAH5','DNAH11']]
plotgenes=pd.DataFrame([W[(W.context==c)&(W.gene==g)].iloc[0] for c,g in pg]);val=plotgenes[['rank_percentile_M14','rank_percentile_M17']].to_numpy();im=bx.imshow(val,aspect='auto',cmap='YlGnBu',vmin=0,vmax=1)
for i in range(len(pg)):
 for j in range(2):bx.text(j,i,f'{val[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if val[i,j]>.65 else 'black')
bx.set_yticks(range(len(pg)),[c.replace('_EC','').replace('Mic_Immune','Mic/Imm')+' / '+g for c,g in pg]);bx.set_xticks([0,1],['EC_M14','EC_M17']);bx.set_title('b  Normalized gene-weight percentiles',loc='left',pad=15,fontweight='bold',fontsize=10);fig.colorbar(im,ax=bx,fraction=.04,pad=.04,label='Within-module gene percentile')
fig.subplots_adjust(wspace=.85,bottom=.13)
fig.text(.08,.025,'Leading-edge membership is descriptive and is not differential expression or exclusive module membership.\nCounterpart leading edges are retained even when their GSEA FDR is nonsignificant; q values are shown.\nSelected genes illustrate pathway overlap and candidate regulatory targets. Gene ranks are computed before source tie-breaking jitter.',fontsize=9)
save(fig,'01_leading_genes');pd.DataFrame(selected).to_csv(T/'figure_leading_edge_values.tsv',sep='\t',index=False);plotgenes.to_csv(T/'figure_gene_rank_values.tsv',sep='\t',index=False)
# Donor relationships: display tied ranks, numerical sensitivities and donor profiles.
fig=plt.figure(figsize=(13,8));gs=fig.add_gridspec(2,2,height_ratios=[3,1],hspace=.6,wspace=.45);ax=fig.add_subplot(gs[0,0]);rng=np.random.default_rng(417)
for group,col in colors.items():
 z=D[D.Pathology==group];ax.scatter(z.M14_rank_pct+rng.uniform(-.005,.005,len(z)),z.M17_rank_pct+rng.uniform(-.005,.005,len(z)),color=col,s=35,alpha=.8,label=f'{group} (n={len(z)})')
ax.set(xlim=(0,1.04),ylim=(0,1.04),xlabel='EC_M14 donor rank percentile',ylabel='EC_M17 donor rank percentile');ax.legend(frameon=False,fontsize=8,loc='upper left');ax.set_title('a  Same donors, partly distinct activity patterns',loc='left',fontweight='bold',fontsize=10)
ci=read('donor_correlation_sensitivity').iloc[0];z=C[(C.analysis=='Spearman')&(C.relative_zero_threshold==1e-8)].iloc[0];ax.text(.03,.02,f'Spearman ρ = {z.effect:.3f}, P = {z.p:.3f}\nBootstrap 95% CI [{ci.CI025:.3f}, {ci.CI975:.3f}]',transform=ax.transAxes,fontsize=9)
bx=fig.add_subplot(gs[0,1]);rows=[C[(C.analysis=='Spearman')&(C.relative_zero_threshold==0)].iloc[0],z,C[C.analysis=='partial_rank_pathology'].iloc[0],C[C.analysis=='partial_rank_pathology_age_sex_pmi'].iloc[0]]
for i,r in enumerate(rows):bx.scatter(r.effect,i,s=50,color='#228b8d');bx.text(.6,i,f'{r.effect:.3f}; P={r.p:.3g}',va='center',fontsize=8)
bx.axvline(0,color='#aaa',lw=1);bx.set_yticks(range(4),['Raw ranks','Numerical zeros tied','Adjusted for pathology','Pathology + age/sex/PMI']);bx.invert_yaxis();bx.set(xlim=(-.3,1.02),xlabel='Correlation / partial rank correlation');bx.set_title('b  Dependence on numerical and pathology handling',loc='left',fontweight='bold',fontsize=10)
cx=fig.add_subplot(gs[1,:]);z=D.copy();z['order']=pd.Categorical(z.Pathology,['nonAD','earlyAD','lateAD'],ordered=True);z=z.sort_values(['order','rank_contrast_M17_minus_M14']);im=cx.imshow(z[['M14_rank_pct','M17_rank_pct']].to_numpy().T,aspect='auto',vmin=0,vmax=1,cmap='YlGnBu');cx.set_yticks([0,1],['EC_M14','EC_M17']);cx.set_xticks([]);cx.set_title('c  Donor ranks, ordered within pathology by M17 − M14 rank contrast',loc='left',fontweight='bold',fontsize=10,pad=20)
start=0
for group in ['nonAD','earlyAD','lateAD']:
 n=int((z.Pathology==group).sum());cx.text(start+(n-1)/2,-.8,f'{group}: {n}',ha='center',fontsize=8)
 if start:cx.axvline(start-.5,color='white',lw=3)
 start+=n
fig.colorbar(im,ax=cx,fraction=.02,pad=.015,label='Rank percentile');fig.text(.1,.015,'46 matched donors. For new comparisons only, H < 10⁻⁸ × module maximum is tied at zero; raw source H and original selection statistics are preserved.\nScatter jitter (±0.005) separates tied points for display only. Secondary P values are descriptive sensitivity results, not independent discoveries.',fontsize=8)
fig.subplots_adjust(bottom=.15);save(fig,'02_donor_patterns')
# Show both corrected M14 associations and the relevant nonrobust immune-mixture association.
features=[('RNA.Subtype|Inh|Inh FBN2 EPB41L4A','M14','Inh FBN2 EPB41L4A'),('RNA.Subtype|Inh|Inh L1-6 LAMP5 CA13','M14','Inh L1-6 LAMP5 CA13'),('RNA.Subclass|Mic_Immune|T','M17','T within Mic/Imm')]
fig,axs=plt.subplots(1,3,figsize=(14,4.8))
for ax,(feature,outcome,label) in zip(axs,features):
 st=S[(S.feature==feature)&(S.outcome==outcome)].iloc[0];z=N[(N.feature==feature)&N.parent_coverage].merge(D,on='donor_id')
 for group,col in colors.items():
  zz=z[z.Pathology==group];ax.scatter(zz[outcome+'_rank_pct'],zz.fraction*100,color=col,s=28,alpha=.7,label=group)
 ax.set_xlabel(f'EC_{outcome} donor rank percentile');ax.set_ylabel('Nuclei within parent cell class (%)');ax.set_title(label,fontweight='bold',fontsize=10)
 ax.text(.02,.98,f'n={int(st.n)}\nRaw ρ={st.rho:.3f}; q={st.q_BH:.3g}\nAdjusted r={st.partial_r:.3f}; q={st.partial_q_BH:.3g}',transform=ax.transAxes,va='top',fontsize=9);ax.margins(y=.3)
axs[-1].legend(frameon=False,fontsize=8,loc='lower right');fig.subplots_adjust(wspace=.35,bottom=.27)
fig.text(.08,.045,'Exploratory selection for display: two covariate-adjusted FDR hits and the T-cell fraction relevant to Mic/Imm gene annotations.\nPartial rank correlations adjust categorical pathology, capped age, sex and PMI; BH across 48 eligible features × 3 outcomes (144 tests) per family.\nFractions measure sampled RNA nuclei, not absolute cell abundance. Matching subclass/subtype labels and complementary fractions are not independent findings.\nNo association with the M17 − M14 rank contrast passes adjusted FDR < 0.05.',fontsize=8)
save(fig,'03_composition_context')
print('Saved three focused EC figures with manuscript PDF copies.')
