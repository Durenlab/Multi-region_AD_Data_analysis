"""Synthesize all selected modules before question-driven mechanistic narrowing."""
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2];D=R/'experiments/EXP004_results';T=D/'tables';U=D/'regional/tables'
read=lambda n:pd.read_csv(T/(n+'.tsv'),sep='\t')
mp=read('module_pathology');ld=read('ldsc_all');age=pd.read_csv(U/'aging_requested_modules.tsv',sep='\t');manifest=pd.read_csv(U/'GSEA_source_manifest.tsv',sep='\t')
order=['AG_M13','EC_M14','EC_M17','MTC_M6','MTC_M13','PFC_M14','PFC_M18'];fmt=lambda x:f'{x:.4g}'
def table(headers,rows):
 clean=lambda x:str(x).replace('|','/').replace('\n',' ')
 return '\n| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(clean(x) for x in row)+' |' for row in rows)+'\n\n'
parts=['''# EXP004: Integrated synthesis across the complete candidate landscape

Revised for the clarified narrative: **cell-specific LDSC/GSEA → cross-cell LDSC/GSEA → inclusive synthesis → justified narrowing → focused TF/CellChat mechanism**. All seven pathology-selected shared modules are considered on the same basis. **This section characterizes all seven modules before narrowing.** AG_M13, EC_M14/EC_M17 and MTC_M13 are not a preselected subset.

The integration reads 130 full GSEA reports covering all seven cross-cell modules × six cell types × two libraries and all 23 cell-specific candidates × two libraries (229,446 rows). No new model fitting, large analysis or Slurm job was run. [Revised manuscript](../manuscript/draft.txt); [figure legends](../manuscript/EXP004_figure_legends.md).

The earlier single-module recommendation and four-module revision are historical, superseded drafts. Their underlying source audits remain available, but their editorial selection is not the premise of this synthesis.

The current [whole-module profile](EXP004_module_profiles.md) characterizes EC_M14 across genes/GSEA, genetics, aging, within-module cell-type intersections, TFs and CellChat after summarizing all seven candidates. The [earlier EC comparison](EXP004_EC_comparison.md) remains a supporting resource. The primary organizing question is what EC_M14 represents as a whole.
''']
parts.append((R/'experiments/exp004/regional_narrative.md').read_text())
parts.append((R/'experiments/exp004/module_profile_narrative.md').read_text())
parts.append('## Quantitative support: identical coverage for all seven modules\n\n### Pathology and genetic context\n\nThe pathology screen is a positive-direction union of separate tests, not a combined test. Linear and pairwise BH families contain 20 modules per source region/test/contrast. No rank is assigned below.\n')
rows=[];ledger=[]
for c in order:
 x=mp[(mp.candidate==c)&(mp.test=='Linear_Regression')].iloc[0]
 y=mp[(mp.candidate==c)&(mp.test=='Pairwise_MannWhitneyU')&mp.Comparison.str.contains('nonAD')&mp.Comparison.str.contains('lateAD')].iloc[0]
 z=ld[(ld.kind=='cross')&(ld.candidate==c)].sort_values('FDR_BH_Enrichment_p').iloc[0]
 qualify=[]
 if x.pval_adj<.05 and x.slope>0:qualify.append('Linear')
 if y.pval_adj<.05 and y.late_minus_non_rank_biserial>0:qualify.append('MW')
 rows.append([c,fmt(x.slope),fmt(x.pval)+' / '+fmt(x.pval_adj),fmt(y.late_minus_non_rank_biserial),fmt(y.pval)+' / '+fmt(y.pval_adj),'+'.join(qualify),z.celltype_region+': '+fmt(z.Enrichment)+'; p '+fmt(z.Enrichment_p)+'; q '+fmt(z.FDR_BH_Enrichment_p)])
 ledger.extend([dict(layer='pathology',module_id=c,context=c,test=t.test,source_path=t.source_path,source_row=t.source_row) for t in [x,y]])
parts.append(table(['Module','Linear slope','Linear p / q','Late−non rank effect','MW p / q','Qualifying test','Lowest LDSC q: context and statistics'],rows));pd.DataFrame(ledger).to_csv(U/'pathology_source_ledger.tsv',sep='\t',index=False)
parts.append('### Full LDSC coverage\n\nBH correction uses all 60 cross-cell annotations; the independent cell-specific source uses 27 annotations. Negative estimates reflect regression uncertainty, not negative biological heritability.\n')
x=ld[(ld.kind=='cross')&ld.candidate.isin(order)].copy();x['order']=pd.Categorical(x.candidate,order,ordered=True);x=x.sort_values(['order','celltype_region']);x.drop(columns='order').to_csv(U/'regional_LDSC.tsv',sep='\t',index=False)
parts.append(table(['Module','Cell context','Enrichment','p','BH q'],[[z.candidate,z.celltype_region,fmt(z.Enrichment),fmt(z.Enrichment_p),fmt(z.FDR_BH_Enrichment_p)] for z in x.itertuples()]))
parts.append('### Full signed aging coverage\n\nBH correction is across 20 modules within each cell type–region context. Source zero p/q values are finite-permutation estimates from the supplied 1,000-permutation workflow. Negative enrichment is distinct from positive aging-up enrichment.\n')
a=age.copy();a['order']=pd.Categorical(a.module_id,order,ordered=True);a=a.sort_values(['order','context'])
parts.append(table(['Module','Context','NES','p','BH q','Interpretation'],[[z.module_id,z.context,fmt(z.NES),fmt(z.pval),fmt(z.FDR),'Positive aging-up' if z.NES>0 and z.FDR<.05 else 'Negative aging-up enrichment' if z.NES<0 and z.FDR<.05 else 'Not FDR significant'] for z in a.itertuples()]))
parts.append('''## Source provenance and retained evidence

Each extracted pathway/aging row retains its absolute source path and original row number (header = row 1). Evidence remains available at both modeling levels; the inclusive synthesis does not discard the earlier audits.

| Analysis unit | Evidence table | Definition |
| --- | --- | --- |
| Module × pathology test | [Pathology source ledger](EXP004_results/regional/tables/pathology_source_ledger.tsv), [statistics](EXP004_results/tables/module_pathology.tsv) | Linear ordinal-stage association and late/non-AD MW; separate source BH families |
| Module × cell context × pathway | [Full reports](EXP004_results/regional/tables/full_GSEA.tsv.gz), [130-source manifest](EXP004_results/regional/tables/GSEA_source_manifest.tsv), [leading positive terms](EXP004_results/regional/tables/leading_GSEA_by_context.tsv) | NES, nominal p, source GSEA FDR; missing/unscored terms differ from nonsignificant results |
| Same-context cross/single module pair | [Shared significant terms](EXP004_results/regional/tables/cross_single_shared_pathways.tsv), [gene/RE/TF overlaps](EXP004_results/tables/single_cross_overlaps.tsv) | Exact terms and leading-gene intersection; hard-set overlap is separate from ranked enrichment |
| Cell context × module RE annotation | [42 candidate annotations](EXP004_results/regional/tables/regional_LDSC.tsv), [all source annotations](EXP004_results/tables/ldsc_all.tsv) | Enrichment, nominal p, BH q; original families preserved |
| Cell context × module × age-up gene set | [All 42 candidate contexts](EXP004_results/regional/tables/aging_requested_modules.tsv), [all 36 contexts × 20 modules](EXP004_results/regional/tables/aging_all_contexts.tsv) | Signed NES and BH over 20 modules per context |
| Cell-tagged entity × module | [Assignment counts](EXP004_results/tables/assignment_counts.tsv) | Hard-assigned genes/REs/TFs; not the GSEA-ranked universe |

The aging source scripts are `AD_REGA_region_wise/part10_Age_related/01making_region_union_aging_gene_list.R`, `02making_up_gene_gmt.R` and `03aging_gsea_from_wbv.py` in the sibling project. The supplied GMT pools matched age-up genes across source regions and excludes discordant directions. Aging GSEA used column-then-row L1-normalized WBV, min/max gene-set sizes 5/5,000, 1,000 permutations and seed 123, with small tie-breaking jitter enabled by the launcher. Standalone GMT reconstruction needs clarification because the second R script references an inherited object name; it was not rerun here.

### TF and communication evidence reserved for the focused stage

The [assigned-TF activity](EXP004_results/tables/assigned_TF_activity.tsv), [pathway-linked existing GRN edges](EXP004_results/regional/tables/TF_edges_to_pathway_leading_genes.tsv), [CellChat module-gene links](EXP004_results/tables/cellchat_module_links.tsv), and [MAGMA source audit](EXP004_results/tables/magma_source_comparison.tsv) remain available for follow-up after the biological question is selected. Their availability does not determine which modules enter the synthesis. The previous five-TF shortlist is a historical exploratory subset, not the outcome of this inclusive comparison.

Existing limitations also remain: discussed TNC donor contrasts are not FDR-significant; receptor genes are outside the paired OPC module; shared motif mapping does not uniquely resolve the proposed Ast TFs. MAGMA supplies gene-level evidence only, and the notebook input differs from the original merged table. None of these gaps justifies omitting a module with informative GSEA/LDSC/aging results from the overall summary.

## Figure and manuscript sequence

1. **Existing Figure 2:** establish the cell-specific module landscape and its LDSC/GSEA patterns, including HC and vascular context.
2. **Existing Figure 3 overview:** present all seven pathology-associated shared modules and their composition, LDSC and GSEA.
3. **Integrated synthesis at the end of the shared-module Results:** explain recurring themes, distinct cellular locations of genetic/functional evidence, and differences within regions. Use the [expanded all-seven pathway/aging/LDSC panel](EXP004_results/regional/figures/01_regional_patterns.pdf) and [cross/single comparisons spanning all seven](EXP004_results/regional/figures/03_cross_vs_cell_specific.pdf). These figures support synthesis without an arbitrary four-module subset.
4. **Subsequent narrowing:** choose a question and justify a module or contrast using the inclusive evidence. The [EC comparison](EXP004_results/regional/figures/02_EC_M14_vs_M17.pdf) is a retained optional detail for an explicitly justified within-region immune question, not the premise of the broad summary.
5. **TF activity and focused TF/CellChat interpretation:** retain the existing Figure 4 activity-versus-expression analysis as a distinct analysis layer; connect it to a specific mechanistic story after focus is justified. The manuscript follows the inclusive synthesis with an integrated EC_M14 profile; comparison plots remain supporting material.

PDF, PNG and SVG outputs have exact [pathway](EXP004_results/regional/tables/figure_pathway_values.tsv) and [cross/single comparison](EXP004_results/regional/tables/figure_cross_single_comparisons.tsv) source tables. PDF copies are in `manuscript/`; `Main_Figures.pdf` is preserved. [Legends](../manuscript/EXP004_figure_legends.md), [revision notes](../manuscript/EXP004_revision_notes.md), and the [manuscript diff](../manuscript/EXP004_manuscript_changes.diff) document the revised structure.

## Reproduction and checks

Use the existing `graph_clean` environment:

```bash
ml anaconda3
source activate graph_clean
python experiments/exp004/characterize_regions.py
python experiments/exp004/regional_figures.py
python experiments/exp004/write_report.py
python experiments/exp004/review_regional_evidence.py
```

These commands summarize source results and produce descriptive figures; they do not rerun REGA, GSEA, LDSC or CellChat. [Extraction checks](EXP004_results/regional/tables/regional_checks.tsv) and [source/artifact checks](EXP004_results/regional/tables/final_review_checks.tsv) record the audit. Source-version reconciliation for historical Figure 4 values, LDSC bookkeeping and upstream assignment scripts remains necessary before publication. The existing TF classification calculations and cohort-wide counts were not independently revalidated in this synthesis revision.
''')
parts.append('## Supporting molecular comparisons and literature context\n\nThe following detail is retained as supporting evidence for the whole-module synthesis.\n\n'+(R/'experiments/exp004/biological_synthesis.md').read_text())
(R/'experiments/EXP004_integrated_results.md').write_text('\n\n'.join(parts));print('Wrote inclusive synthesis and evidence tables for all seven modules.')

(R/'experiments/EXP004_module_profiles.md').write_text('# EXP004: biological landscape and integrated EC_M14 profile\n\n'+(R/'experiments/exp004/module_profile_narrative.md').read_text())
