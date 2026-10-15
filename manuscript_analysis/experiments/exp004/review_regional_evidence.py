"""Audit existing source rows and regional deliverables without fitting models."""
from pathlib import Path
import ast,re,hashlib
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[2];D=R/'experiments/EXP004_results';T=D/'regional/tables';checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
for p in (R/'experiments/exp004').glob('*.py'):
 ast.parse(p.read_text());check('Syntax: '+p.name,True)
g=pd.read_csv(T/'full_GSEA.tsv.gz',sep='\t');g['key']=g.term.str.strip().str.upper();manifest=pd.read_csv(T/'GSEA_source_manifest.tsv',sep='\t')
check('130 complete reports; 229446 GSEA rows',len(manifest)==130 and len(g)==229446 and all(Path(p).exists() for p in manifest.path))
check('GSEA source key unique',not g.duplicated(['kind','module_id','context','library','key']).any())
old=pd.read_csv(T/'regional_checks.tsv',sep='\t');check('All 36 aging BH families and source availability',old.passed.all() and len(old)==37)
# Independently reread every full report and compare all extracted values by original row.
err=0;termsok=True
for path,part in g.groupby('source_path'):
 raw=pd.read_csv(path);sel=raw.iloc[part.source_row.astype(int).to_numpy()-2]
 for ext,col in [('nes','NES'),('fdr','FDR q-val'),('pval','NOM p-val')]:err=max(err,float(np.nanmax(np.abs(part[ext].to_numpy()-sel[col].to_numpy()))))
 termsok &= np.array_equal(part.term.str.strip().to_numpy(),sel.Term.str.strip().to_numpy())
check('All GSEA numeric fields and term identities reproduce original rows',err<1e-12 and termsok,err)
ld=pd.read_csv(D/'tables/ldsc_all.tsv',sep='\t')
for kind,n in [('cross',60),('single',27)]:
 z=ld[ld.kind==kind];e=np.max(np.abs(multipletests(z.Enrichment_p,method='fdr_bh')[1]-z.FDR_BH_Enrichment_p));check('LDSC BH family: '+kind,len(z)==n and e<1e-12,e)
a=pd.read_csv(T/'aging_requested_modules.tsv',sep='\t');check('All seven modules × six aging contexts',len(a)==42 and len(a.groupby(['module_id','context']))==42)
x=a[(a.module_id=='EC_M17')&(a.context=='Inh_EC')].iloc[0];check('EC_M17 Inh significant negative aging retained',x.NES<0 and x.FDR<.05)
# Figure numbers are directly aligned to exact term keys.
f=pd.read_csv(T/'figure_pathway_values.tsv',sep='\t');merged=f.merge(g,on=['kind','module_id','context','term'],suffixes=('_fig','_raw'),how='left');valid=merged.nes_fig.notna();err=float(np.max(np.abs(merged.loc[valid,'nes_fig']-merged.loc[valid,'nes_raw'])));qerr=float(np.max(np.abs(merged.loc[valid,'fdr_fig']-merged.loc[valid,'fdr_raw'])));check('Regional figure values match full reports',valid.sum()==f.nes.notna().sum() and merged.loc[valid,'nes_raw'].notna().all() and err<1e-12 and qerr<1e-12,f'{valid.sum()} scored, {(~valid).sum()} unscored')
# Every plotted unscored cell must really be absent from the full report.
absent=True
for z in f[f.nes.isna()].itertuples():absent &= g[(g.kind==z.kind)&(g.module_id==z.module_id)&(g.context==z.context)&(g.key==z.term.strip().upper())].empty
check('Unscored figure cells are absent, not failed significance',absent)
e=pd.read_csv(T/'selected_TF_pathway_evidence.tsv',sep='\t');assign=pd.read_csv(D/'tables/cross_gene_assignments.tsv',sep='\t');grn=pd.read_csv(D/'tables/grn_edges.tsv',sep='\t');tf=pd.read_csv(D/'tables/assigned_TF_activity.tsv',sep='\t');ok=True
for z in e.itertuples():
 hit=g[(g.source_path==z.GSEA_source)&(g.source_row==z.GSEA_row)].iloc[0]
 assigned=assign[(assign.candidate==z.module_id)&(assign.celltype_region==z.context)&(assign.gene==z.TG)]
 edge=grn[(grn.candidate==z.module_id)&(grn.celltype_region==z.context)&(grn.TF==z.TF)&(grn.TG==z.TG)]
 activity=tf[(tf.candidate==z.module_id)&(tf.celltype_region==z.context)&(tf.tf==z.TF)]
 ok &= hit.nes>0 and hit.fdr<.05 and z.TG in hit.leading_genes.split(';') and not assigned.empty and not edge.empty and np.isclose(edge.iloc[0].edge_value,z.edge_value) and not activity.empty and (activity.iloc[0].positive_MW_FDR05 or activity.iloc[0].positive_LM_FDR05)
check('Selected TF targets are assigned genes and significant pathway leading genes with source GRN/activity support',ok,len(e))
# Sanity checks for main distinctions, without inventing a test of module differences.
keys=['OXIDATIVE PHOSPHORYLATION','ALZHEIMER DISEASE'];x=g[(g.kind=='cross')&(g.module_id=='EC_M14')&g.celltype.isin(['Exc','Inh'])&g.key.isin(keys)];check('EC_M14 neuronal OxPhos/AD do not pass positive GSEA FDR',len(x)==4 and not x.positive_fdr05.any())
x=g[(g.kind=='cross')&(g.module_id=='EC_M14')&g.celltype.isin(['Oli','OPC'])&(g.key=='OXIDATIVE PHOSPHORYLATION')];check('EC_M14 Oli/OPC OxPhos pass positive GSEA FDR',len(x)==2 and x.positive_fdr05.all())
check('Three figures × three formats',all((D/'regional/figures'/f'{n}.{ext}').stat().st_size>1000 for n in ['01_regional_patterns','02_EC_M14_vs_M17','03_cross_vs_cell_specific'] for ext in ['pdf','png','svg']))
check('Manuscript figure copies identical',all(hashlib.sha256((D/'regional/figures'/f'{n}.pdf').read_bytes()).digest()==hashlib.sha256((R/'manuscript'/f'EXP004_{n}.pdf').read_bytes()).digest() for n in ['01_regional_patterns','02_EC_M14_vs_M17','03_cross_vs_cell_specific']))
missing=[]
for p in [R/'experiments/EXP004_integrated_results.md',R/'manuscript/EXP004_figure_legends.md']:
 for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if not target.startswith(('http:','https:')) and not (p.parent/target).exists():missing.append(str(p)+': '+target)
check('Report and legend relative links resolve',not missing,missing)
s=(R/'experiments/EXP004_integrated_results.md').read_text();check('Final module selection deferred in current report','Final module selection is deferred' in s and 'Select EC_M14 alone' not in s)
check('Manuscript backup and diff present',all((R/p).exists() for p in ['experiments/EXP004_results/archive/manuscript_before_regional_revision.txt','manuscript/EXP004_manuscript_changes.diff']))
expected={'AG_M13','EC_M14','EC_M17','MTC_M6','MTC_M13','PFC_M14','PFC_M18'}
check('All seven modules represented in main pathway figure',set(f.module_id)==expected)
comp=pd.read_csv(T/'figure_cross_single_comparisons.tsv',sep='\t')
check('All seven modules represented in cross/single comparison figure',set(comp.cross_module)==expected)
manuscript=(R/'manuscript/draft.txt').read_text()
check('Inclusive manuscript synthesis follows broad shared results',manuscript.index('Cross-cell-type modules reveal complementary regional programs')<manuscript.index('Integrated functional and genetic evidence distinguishes multicellular response patterns')<manuscript.index('Regulatory-element-derived transcription factor activity'))
check('Premature focused regulatory section removed', 'Pathway-linked regulatory candidates and exploratory communication' not in manuscript)
z=pd.DataFrame(checks);z.to_csv(T/'final_review_checks.tsv',sep='\t',index=False);print(z.to_string(index=False));print('PASS',int(z.passed.sum()),'/',len(z))
if not z.passed.all():raise SystemExit('An evidence check failed; see final_review_checks.tsv')
