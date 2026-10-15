"""Plot existing CellChat evidence beside EC hypotheses; no inferential tests/refits.

Keep the prior hard-assignment intersection and the full GSEA leading-edge
intersection distinct. Pooled rows are descriptive; donor P/q values are read
from EXP002 without recomputing correction on the selected examples.
"""
from pathlib import Path
import shutil
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

R = Path(__file__).resolve().parents[2]
O = R / 'experiments/EXP004_results/ec_comparison'
T, F = O / 'tables', O / 'figures'
B = R / 'experiments/EXP002_results'
sources = []
def read(path):
    path = Path(path)
    x = pd.read_csv(path, sep='\t')
    x['input_source_path'] = str(path)
    x['input_source_row'] = np.arange(2, len(x) + 2)
    sources.append({'path': str(path), 'rows': len(x)})
    return x

pooled = read(R.parent / 'AD_REGA_region_wise/part8_cellchat/04_pathology_comparison/EC/tables/EC_all_conditions_ligand_receptor_communications.tsv')
meta = read(B / 'tables/EC_candidate_interactions.tsv')
genes = read(T / 'focus_pathway_gene_details.tsv')
assigned = read(R / 'experiments/EXP004_results/tables/cross_gene_assignments.tsv')
strict = read(T / 'EC_pathway_compatible_CellChat.tsv')
trend = read(B / 'tables/donor_pathology_progression_tests.tsv')
pairwise = read(B / 'tables/donor_pathology_pairwise_tests.tsv')
scores = read(B / 'intermediate/EC_donor_nonzero_probabilities.tsv.gz')
coverage = read(B / 'tables/EC_donor_celltype_coverage.tsv')
edges = read(T / 'EC_TF_pathway_edges.tsv')
trend = trend[trend.region.eq('EC')]
pairwise = pairwise[pairwise.region.eq('EC') & pairwise.comparison.eq('lateAD_minus_nonAD')]
keys = ['source', 'target', 'interaction_id']
order = ['nonAD', 'earlyAD', 'lateAD']
colors = ['#66AABB', '#EEAA55', '#CC6677']
cts = ['Exc', 'Inh', 'Ast', 'Oli', 'OPC', 'Mic_Immune', 'Vasc']
def symbols(v):
    return set(str(v).replace(' ', '').split(',')) if pd.notna(v) else set()
def joined(v):
    return ';'.join(sorted(v))

pool = pooled.rename(columns={'interaction_name': 'interaction_id', 'input_source_path': 'pooled_source_path', 'input_source_row': 'pooled_source_row'})
pool = pool.merge(meta[['interaction_id', 'ligand.symbol', 'receptor.symbol', 'input_source_path', 'input_source_row']], on='interaction_id', validate='many_to_one')
pool = pool.rename(columns={'input_source_path': 'metadata_source_path', 'input_source_row': 'metadata_source_row'})
mapped = []
for module in ['M14', 'M17']:
    lead = {c: set(z.gene) for c, z in genes[genes['leading_' + module] & genes['positive_' + module]].groupby('context')}
    hard = {c: set(z.gene) for c, z in assigned[assigned.candidate.eq('EC_' + module)].groupby('celltype_region')}
    z = pool[pool.source.ne(pool.target)].copy()
    z['module_id'] = 'EC_' + module
    z['leading_ligand_genes'] = [joined(symbols(g) & lead.get(c + '_EC', set())) for g, c in zip(z['ligand.symbol'], z.source)]
    z['leading_receptor_genes'] = [joined(symbols(g) & lead.get(c + '_EC', set())) for g, c in zip(z['receptor.symbol'], z.target)]
    z['assigned_leading_ligand_genes'] = [joined(set(g.split(';')) & hard.get(c + '_EC', set())) for g, c in zip(z.leading_ligand_genes, z.source)]
    z['assigned_leading_receptor_genes'] = [joined(set(g.split(';')) & hard.get(c + '_EC', set())) for g, c in zip(z.leading_receptor_genes, z.target)]
    mapped.append(z[z.leading_ligand_genes.ne('') | z.leading_receptor_genes.ne('')])
mapped = pd.concat(mapped, ignore_index=True)
tr = trend[keys + ['n_nonAD', 'n_earlyAD', 'n_lateAD', 'spearman_rho', 'p_value', 'q_value', 'status', 'input_source_path', 'input_source_row']].rename(columns={c: 'trend_' + c for c in ['p_value', 'q_value', 'status', 'input_source_path', 'input_source_row']})
pw = pairwise[keys + ['n_condition_1', 'n_condition_2', 'median_difference', 'rank_biserial_delta', 'p_value', 'q_value', 'status', 'input_source_path', 'input_source_row']].rename(columns={c: 'pairwise_' + c for c in ['p_value', 'q_value', 'status', 'input_source_path', 'input_source_row']})
mapped = mapped.merge(tr, on=keys, how='left', validate='many_to_one').merge(pw, on=keys, how='left', validate='many_to_one')
mapped.to_csv(T / 'cellchat_ranked_pathway_pooled_links.tsv', sep='\t', index=False)
family = mapped.groupby(['module_id', 'pathway_name']).agg(pooled_condition_rows=('prob', 'size')).reset_index()
unique = mapped.drop_duplicates(['module_id'] + keys)
unique.groupby(['module_id', 'pathway_name']).size().rename('unique_directed_interactions').reset_index().merge(family).to_csv(T / 'cellchat_ranked_pathway_family_counts.tsv', sep='\t', index=False)

plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False})
def save(fig, name):
    for ext in ['pdf', 'png', 'svg']:
        fig.savefig(F / f'{name}.{ext}', dpi=170, bbox_inches='tight', facecolor='white')
    shutil.copy2(F / f'{name}.pdf', R / 'manuscript' / f'EXP004_EC_{name}.pdf')
    plt.close(fig)

# Broad pathway coverage, followed by matched receptor hypotheses in two contexts.
fig = plt.figure(figsize=(15, 7.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=.6)
ax = fig.add_subplot(gs[0, 0])
counts = unique.groupby(['pathway_name', 'module_id']).size().unstack(fill_value=0).sort_index()
y = np.arange(len(counts))
for j, (m, color) in enumerate([('EC_M14', '#228b8d'), ('EC_M17', '#ba6848')]):
    ax.barh(y + (j - .5) * .36, counts[m], height=.34, color=color, label=m)
ax.set_yticks(y, counts.index)
ax.invert_yaxis()
ax.set_xlabel('Unique sender–receiver–interaction links')
ax.set_title('a  Pooled links to ranked pathway genes', loc='left', fontsize=11, pad=24)
ax.legend(frameon=False, fontsize=9)
for panel, module, target, interaction, title in [(1, 'EC_M17', 'Mic_Immune', 'APP_CD74', 'b  EC_M17: APP → CD74 in Mic/Imm'), (2, 'EC_M14', 'Exc', 'APP_TNFRSF21', 'c  EC_M14: APP → TNFRSF21 in Exc')]:
    ax = fig.add_subplot(gs[0, panel])
    z = mapped[mapped.module_id.eq(module) & mapped.target.eq(target) & mapped.interaction_id.eq(interaction)]
    sender = [c for c in cts if c != target and c in set(z.source)]
    a = z.pivot(index='source', columns='condition', values='prob').reindex(index=sender, columns=order)
    im = ax.imshow(a, vmin=0, vmax=.16, cmap='YlGnBu', aspect='auto')
    for i in range(len(sender)):
        for j in range(3):
            v = a.iloc[i, j]
            ax.text(j, i, f'{v:.3f}' if pd.notna(v) else 'NR', ha='center', va='center', fontsize=10, color='white' if v > .085 else 'black')
    ax.set_yticks(range(len(sender)), [s.replace('Mic_Immune', 'Mic/Imm') for s in sender])
    ax.set_xticks(range(3), order, rotation=30)
    ax.set_ylabel('Sender cell class')
    ax.set_title(title, loc='left', fontsize=11, pad=24)
    fig.colorbar(im, ax=ax, fraction=.046, pad=.04, label='Pooled communication probability')
fig.subplots_adjust(bottom=.25, top=.87)
fig.text(.07, .045, 'GSEA extension: a ligand or receptor is in a significant focused-pathway leading edge; hard module assignment is not required.\n'
         'Panels b–c include every reported cross-cell sender for these receptors. CD74 (Mic/Imm) and TNFRSF21 (Exc) are not hard-assigned\n'
         'to the displayed module. Pooled condition networks are descriptive; their probabilities and within-network P values do not test\n'
         'pathology differences between donors. One extended link, OPC → Ast COL6A1–CD44, passes donor BH correction (Figures 5–7).', fontsize=10)
save(fig, '04_cellchat_pathway_context')

# Display examples chosen for biological context, not smallest P values.
examples = [
    ('EC_M17 antigen pathway', 'Exc', 'Mic_Immune', 'APP_CD74', 'Ranked-pathway extension'),
    ('EC_M14 interaction pathway', 'Oli', 'Exc', 'APP_TNFRSF21', 'Ranked-pathway extension'),
    ('Immune receptor context', 'Inh', 'Mic_Immune', 'CX3CL1_CX3CR1', 'Context; see pathway mapping table'),
    ('EC_M14 Oli CD74 hypothesis', 'Mic_Immune', 'Oli', 'APP_CD74', 'Original assignment + leading-edge filter'),
    ('EC_M14 Oli CD74 hypothesis', 'Mic_Immune', 'Oli', 'MIF_CD74_CD44', 'Original assignment + leading-edge filter'),
    ('EC_M14 Ast ECM pathway', 'OPC', 'Ast', 'COL6A1_CD44', 'Only corrected link in ranked-pathway pooled extension'),
]
selected_stats, donor_rows = [], []
for title, source, target, interaction, reason in examples:
    a = coverage[coverage.celltype.eq(source)].copy()
    b = coverage[coverage.celltype.eq(target)].copy()
    v = a.merge(b, on=['donor_id', 'pathology'], suffixes=('_sender', '_receiver'), validate='one_to_one')
    v['eligible'] = v.donor_network_status_sender.eq('completed') & v.donor_network_status_receiver.eq('completed') & (v.n_cells_sender >= v.min_cells_required_sender) & (v.n_cells_receiver >= v.min_cells_required_receiver)
    raw = scores[scores.source.eq(source) & scores.target.eq(target) & scores.interaction_id.eq(interaction)]
    v = v.merge(raw[['donor_id', 'communication_probability', 'input_source_row']], on='donor_id', how='left', validate='one_to_one')
    v['communication_probability'] = v.communication_probability.where(~v.eligible, v.communication_probability.fillna(0))
    v.loc[~v.eligible, 'communication_probability'] = np.nan
    v['zero_filled_for_eligible_pair'] = v.eligible & v.input_source_row.isna()
    v['source'] = source; v['target'] = target; v['interaction_id'] = interaction
    v['score_source_path'] = str(B / 'intermediate/EC_donor_nonzero_probabilities.tsv.gz')
    donor_rows.append(v)
    st = tr[tr.source.eq(source) & tr.target.eq(target) & tr.interaction_id.eq(interaction)].merge(pw, on=keys).iloc[0].to_dict()
    st.update(display_title=title, selection_reason=reason)
    selected_stats.append(st)
display = pd.concat(donor_rows, ignore_index=True)
st = pd.DataFrame(selected_stats)
display.to_csv(T / 'cellchat_figure_donor_scores.tsv', sep='\t', index=False)
st.to_csv(T / 'cellchat_figure_selected_statistics.tsv', sep='\t', index=False)
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
rng = np.random.default_rng(417)
for ax, spec, stat in zip(axs.flat, examples, selected_stats):
    title, source, target, interaction, reason = spec
    d = display[display.source.eq(source) & display.target.eq(target) & display.interaction_id.eq(interaction)]
    labels = []
    for i, (group, color) in enumerate(zip(order, colors)):
        z = d[d.pathology.eq(group)]
        vals = z.loc[z.eligible, 'communication_probability'].to_numpy()
        ax.scatter(i + rng.uniform(-.13, .13, len(vals)), vals, color=color, edgecolor='white', linewidth=.3, s=30, alpha=.85)
        if len(vals): ax.plot([i - .2, i + .2], [np.median(vals)] * 2, color='#333', lw=2)
        labels.append(f'{group}\nn={len(vals)}; missing={int((~z.eligible).sum())}')
    ax.set_xticks(range(3), labels, fontsize=9)
    ax.set_ylabel('Donor communication probability')
    ax.set_title(f'{title}\n{source.replace("Mic_Immune", "Mic/Imm")} → {target.replace("Mic_Immune", "Mic/Imm")}: {interaction}', fontsize=10, loc='left', pad=12)
    pv, qv = stat['pairwise_p_value'], stat['pairwise_q_value']
    caption = f'lateAD vs nonAD: P={pv:.3g}; q={qv:.3g}' if np.isfinite(pv) else 'lateAD vs nonAD: P/q unavailable (constant scores)'
    ax.text(.02, .95, caption + f'\nTrend P={stat["trend_p_value"]:.3g}; q={stat["trend_q_value"]:.3g}', transform=ax.transAxes, va='top', fontsize=8.5)
    maxv = d.communication_probability.max()
    ax.set_ylim(-.03 * maxv if maxv > 0 else -.005, maxv * 1.45 if maxv > 0 else .04)
fig.subplots_adjust(hspace=.6, wspace=.32, bottom=.14)
fig.text(.07, .025, 'Dots are eligible donors; horizontal bars are medians. Cell-class pairs require ≥21 nuclei in both classes and a completed donor network.\n'
         'Absent sparse scores become zero only for eligible pairs; inadequate cell coverage remains missing. No new statistical tests were run.\n'
         'P/q values retain the source EC-wide trend and EC/contrast-wide pairwise BH families. Independent panel scales show the score distributions.', fontsize=10)
save(fig, '05_cellchat_donor_evidence')

# A hypothesis diagram shows regulatory and communication evidence separately.
fig, axs = plt.subplots(2, 1, figsize=(14, 9.5))
def box(ax, x, y, w, h, text, color='#eef3f4', size=10):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=.012', facecolor=color, edgecolor='#52616b', linewidth=1))
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=size)
def arrow(ax, a, b, label='', dashed=False):
    ax.annotate('', xy=b, xytext=a, arrowprops=dict(arrowstyle='->', color='#777' if dashed else '#228b8d', lw=1.8, linestyle='--' if dashed else '-'))
    if label: ax.text((a[0]+b[0])/2, (a[1]+b[1])/2+.035, label, ha='center', fontsize=9)
for ax in axs: ax.set(xlim=(0,1), ylim=(0,1)); ax.axis('off')
ax = axs[0]
ax.set_title('a  EC_M14 Oli: regulatory hypothesis with an unsupported communication extension', loc='left', fontsize=12)
box(ax,.02,.57,.21,.24,'SPI1 / CEBPD\nOli inferred TF activity\nBoth pathology tests q < 0.05')
box(ax,.36,.65,.20,.19,'CD74\nAntigen leading edge\nHard-assigned EC_M14')
box(ax,.36,.23,.20,.19,'ATP6V0D2\nVacuolar V-ATPase\nHard-assigned EC_M14')
arrow(ax,(.23,.73),(.36,.74)); arrow(ax,(.23,.62),(.36,.33))
ax.text(.025,.25,'Existing GRN weights:\nCD74: SPI1 86.70; CEBPD 19.29\nATP6V0D2: SPI1 174.22; CEBPD 41.28\nModel weights do not establish binding.', fontsize=9)
box(ax,.76,.57,.21,.24,'Other cell classes\nAPP or MIF ligand\nCandidate database links', '#f4eeee')
arrow(ax,(.76,.69),(.56,.74),'Hypothesis only',dashed=True)
ax.text(.79,.23,'Incoming APP–CD74 / MIF–CD74:\nno matching pooled Oli rows;\ndonor scores constant.\nMIF receptor complex includes CD74\nwith CD44 / CXCR4 / CXCR2.',ha='center',fontsize=9)
ax.text(.02,.04,'Solid arrows: existing inferred TF→target weights. Dashed arrow: unvalidated ligand–receptor hypothesis.\nNo arrow connects the incoming ligand to TF activation; that downstream mechanism has not been established.',fontsize=10)
ax = axs[1]
ax.set_title('b  EC_M17 Mic/Imm: antigen-pathway receptor context, with unresolved pathology change',loc='left',fontsize=12)
box(ax,.02,.53,.24,.25,'APP-expressing sender classes\nExc, Inh, Ast, Oli, OPC, Vasc\nPooled CellChat rows available')
box(ax,.44,.53,.23,.25,'CD74 in Mic/Imm\nEC_M17 antigen leading edge\nNot hard-assigned EC_M17')
arrow(ax,(.26,.65),(.44,.65),'APP–CD74',dashed=True)
box(ax,.77,.53,.21,.25,'Antigen-associated\nreceiver context\nRegulator unresolved', '#f4eeee')
ax.text(.51,.31,'Observed: pooled APP–CD74 communication probabilities in all three groups.\nDonor example, Exc → Mic/Imm: lateAD/nonAD P=0.218; q=0.793.\nInterpretation: candidate receptor context; no corrected evidence of a pathology difference.',ha='center',fontsize=10)
ax.text(.02,.06,'Mic/Imm and Oli are different receiver contexts. The EC_M14 Oli TF→CD74 weights cannot be transferred to EC_M17 Mic/Imm.\nModule LDSC does not identify a causal ligand, receptor or TF. Each proposed mechanism requires cell-resolved validation.',fontsize=10)
fig.subplots_adjust(hspace=.22)
save(fig, '06_cellchat_hypotheses')
# The corrected extended connection deserves its own evidence panel.
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], wspace=.4, hspace=.48)
ax = fig.add_subplot(gs[0, :]); ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
ax.set_title('EC_M14 pathway-associated hypothesis: OPC → Ast COL6A1–CD44',loc='left',fontsize=14,fontweight='bold')
box(ax,.025,.35,.24,.42,'OPC sender\nCOL6A1 ligand\nNo EC_M14 assignment link')
box(ax,.42,.35,.25,.42,'Ast receiver\nCD44 receptor\nEC_M14 ECM leading edge')
box(ax,.80,.35,.17,.42,'Candidate\nECM response\nTF unresolved', '#f4eeee')
arrow(ax,(.265,.56),(.42,.56),'CellChat',dashed=True)
ax.text(.49,.08,'Hypothesis: increased inferred COL6A1–CD44 communication accompanies an Ast ECM-associated program.\nCD44 is a ranked pathway gene, not an EC_M14 hard-assigned gene; neither a downstream TF nor causation is established.',ha='center',fontsize=10)
ax=fig.add_subplot(gs[1,0]); z=pooled[pooled.source.eq('OPC')&pooled.target.eq('Ast')&pooled.interaction_name.eq('COL6A1_CD44')]
for i,(group,color) in enumerate(zip(order,colors)):
    row=z[z.condition.eq(group)]
    if len(row):
        val=float(row.iloc[0].prob); ax.bar(i,val,color=color,width=.6);ax.text(i,val+.00015,f'{val:.5f}',ha='center',fontsize=10)
    else:ax.text(i,.0003,'Not reported',rotation=45,ha='center',fontsize=10,color='#666')
ax.set_xticks(range(3),order);ax.set(xlim=(-.6,2.6),ylim=(0,.0045),ylabel='Pooled communication probability')
ax.set_title('a  Pooled networks: descriptive',loc='left',fontsize=11)
ax=fig.add_subplot(gs[1,1]); d=display[display.source.eq('OPC')&display.target.eq('Ast')&display.interaction_id.eq('COL6A1_CD44')]
for i,(group,color) in enumerate(zip(order,colors)):
    vals=d.loc[d.pathology.eq(group)&d.eligible,'communication_probability'].to_numpy()
    ax.scatter(i+rng.uniform(-.13,.13,len(vals)),vals,s=32,color=color,alpha=.85)
    ax.plot([i-.2,i+.2],[np.median(vals)]*2,color='#333',lw=2)
ax.set_xticks(range(3),['nonAD\nn=25','earlyAD\nn=14','lateAD\nn=7']);ax.set_ylabel('Donor communication probability')
ax.set_title('b  Donor evidence: source BH correction retained',loc='left',fontsize=11)
stat=st[st.interaction_id.eq('COL6A1_CD44')].iloc[0]
ax.text(.03,.97,f'Trend ρ={stat.spearman_rho:.3f}; P={stat.trend_p_value:.3g}; q={stat.trend_q_value:.3g}\nlateAD/nonAD: P={stat.pairwise_p_value:.3g}; q={stat.pairwise_q_value:.3g}\nRank-biserial effect={stat.rank_biserial_delta:.3f}',transform=ax.transAxes,va='top',fontsize=9)
ax.set_ylim(-.0003,d.communication_probability.max()*1.55)
fig.subplots_adjust(bottom=.16)
fig.text(.07,.025,'This link emerged from an explicitly broader, descriptive leading-edge mapping; the original strict assignment intersection remains unchanged.\nNot reported in a pooled table is not a measured zero. Donor scores and source P/q values are existing EXP002 results, not new tests.\nCellChat infers communication from expression. Spatial localization, composition-aware expression and perturbation are needed to evaluate the mechanism.',fontsize=10)
save(fig,'07_cellchat_COL6A1_CD44')
edges[edges.context.eq('Oli_EC') & edges.TF.isin(['SPI1','CEBPD']) & edges.TG.isin(['CD74','ATP6V0D2'])].to_csv(T/'cellchat_hypothesis_GRN_edges.tsv',sep='\t',index=False)
pd.DataFrame(sources).to_csv(T/'cellchat_figure_sources.tsv',sep='\t',index=False)
summary = {
    'selection': 'Original hard-assignment links preserved; extension includes any significant focused-pathway leading-edge ligand/receptor in the same sender/receiver context and a reported pooled row.',
    'pooled_condition_rows_by_module': mapped.groupby('module_id').size().to_dict(),
    'unique_directed_links_by_module': unique.groupby('module_id').size().to_dict(),
    'ranked_pooled_trend_BH_hits': int((unique.trend_q_value < .05).sum()),
    'ranked_pooled_pairwise_BH_hits': int((unique.pairwise_q_value < .05).sum()),
    'minimum_ranked_pooled_trend_q': float(unique.trend_q_value.min()),
    'minimum_ranked_pooled_pairwise_q': float(unique.pairwise_q_value.min()),
    'strict_pooled_matching_rows': int((mapped.assigned_leading_ligand_genes.ne('') | mapped.assigned_leading_receptor_genes.ne('')).sum()),
    'missingness': 'Sparse missing scores zero-filled only for completed donor networks with >=21 cells in both classes; other pairs remain missing.',
    'tests': 'No statistical tests or BH corrections rerun; original source values retained.',
}
(O/'cellchat_figure_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print(st[['source','target','interaction_id','pairwise_p_value','pairwise_q_value']].to_string(index=False))
print('Saved four CellChat figure sets and manuscript PDF copies.')
