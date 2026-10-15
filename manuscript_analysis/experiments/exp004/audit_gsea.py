from pathlib import Path
import pandas as pd, numpy as np
root=Path(__file__).resolve().parents[2];tab=root/'experiments/EXP004_results/tables'
h=pd.read_csv(tab/'gsea_hits.tsv',sep='\t'); allrows=[]
for (source,m),a in h.groupby(['source_file','module']):
 p=Path(source).parent/m/'gseapy.gene_set.prerank.report.csv'
 if not p.exists():
  a['raw_report_available']=False;allrows.append(a);continue
 raw=pd.read_csv(p);raw['raw_report_path']=str(p);raw['raw_report_row']=np.arange(2,len(raw)+2)
 raw=raw.rename(columns={'Term':'term','NOM p-val':'nominal_p','FDR q-val':'raw_GSEA_FDR','Lead_genes':'raw_lead_genes'})
 a=a.merge(raw[['term','nominal_p','raw_GSEA_FDR','raw_lead_genes','raw_report_path','raw_report_row']],on='term',how='left',validate='many_to_one');a['raw_report_available']=True
 allrows.append(a)
d=pd.concat(allrows);d.to_csv(tab/'gsea_hits_audited.tsv',sep='\t',index=False,na_rep='NA')
print('GSEA rows',len(d),'raw found',d.raw_report_available.sum(),'max FDR disagreement',(d.fdr-d.raw_GSEA_FDR).abs().max())
for cand,a in d[d.kind=='cross'].groupby('candidate'):
 print('\n',cand)
 for lib in ['GOBP','KEGG']:
  b=a[a.library==lib].sort_values(['fdr','nes'],ascending=[True,False]);print(b[['celltype_region','term','nes','nominal_p','fdr']].head(9).to_string(index=False))
