"""Small output and provenance checks; does not refit the biological analyses."""
from pathlib import Path
import ast,re
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];T=R/'experiments/EXP004_results/tables';checks=[]
for p in (R/'experiments/exp004').glob('*.py'):ast.parse(p.read_text());checks.append(dict(check='Python syntax: '+p.name,passed=True,max_abs_error=0))
p=R/'experiments/EXP004_integrated_results.md';missing=[]
for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
 if not target.startswith(('http:','https:')) and not (p.parent/target).exists():missing.append(target)
checks.append(dict(check='report relative links resolve',passed=not missing,max_abs_error=len(missing)))
c=pd.read_csv(T/'cross_gene_assignments.tsv',sep='\t');a=pd.read_csv(T/'assignment_counts.tsv',sep='\t');
checks.append(dict(check='gene counts match assignment rows',passed=len(c)==a[a.entity_kind=='gene'].n_assigned.sum(),max_abs_error=abs(len(c)-a[a.entity_kind=='gene'].n_assigned.sum())))
scores=pd.read_csv(T/'TNC_donor_scores.tsv',sep='\t');links=pd.read_csv(T/'cellchat_module_links.tsv',sep='\t');maxerr=0
for inter,d in scores.groupby('interaction_id'):
 row=links[(links.candidate=='EC_M14')&(links.source=='Ast')&(links.target=='OPC')&(links.interaction_id==inter)].iloc[0]
 diff=d[d.pathology=='lateAD'].communication_probability.median()-d[d.pathology=='nonAD'].communication_probability.median();maxerr=max(maxerr,abs(diff-row.median_difference))
 assert d.groupby('pathology').size().to_dict()=={'earlyAD':14,'lateAD':7,'nonAD':25}
checks.append(dict(check='TNC plotted donor counts and median contrast match source',passed=maxerr<1e-12,max_abs_error=maxerr))
figs=[p for p in (R/'experiments/EXP004_results/figures').glob('*') if p.is_file() and p.suffix in {'.pdf','.png','.svg'}];checks.append(dict(check='three figures in PDF PNG SVG',passed=len(figs)==9 and all(p.stat().st_size>1000 for p in figs),max_abs_error=abs(len(figs)-9)))
z=pd.read_csv(T/'validation_checks.tsv',sep='\t');z=pd.concat([z,pd.DataFrame(checks)]).drop_duplicates('check',keep='last');z.to_csv(T/'validation_checks.tsv',sep='\t',index=False)
print(z.to_string(index=False));assert z.passed.all(), 'A validation failed'
