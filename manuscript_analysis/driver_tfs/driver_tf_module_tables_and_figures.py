#!/usr/bin/env python3
"""Summarize disease-module driver TFs and current MAGMA TF-list membership.

Inputs are the current REGA TF results, MAGMA gene results, and the
project's MAGMA FDR-significant TF list. No upstream result is modified.
"""
from pathlib import Path
import argparse
import math
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

ROOT = Path(os.environ.get("AD_MANUSCRIPT_WORK_ROOT", Path(__file__).resolve().parents[3])).expanduser().resolve()
OUT = ROOT / "experiments/driver_TF_disease_modules"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--figures-only", action="store_true", help="Regenerate figures without rewriting tables")
FIGURES_ONLY = parser.parse_args().figures_only

BASE = Path(os.environ.get("REGA_DATA_ROOT", ROOT.parent)).expanduser().resolve()
CELLTYPE_MODULES = (BASE / "REGA_separate/02_run_rega/pathology_association/selected_module.txt").read_text().split()
CROSS_MODULES = (BASE / "REGA_shared/part2_run_REGA/08pathology_tests/selected_module.txt").read_text().split()
EXPLORATORY_CROSS_MODULES = set()

SNP_LIST_PATH = Path(os.environ.get("AD_SNP_TF_LIST", ROOT / "experiments/TF_prioritize/TF_with_SNP/TFs_with_SNP.txt")).expanduser().resolve()
MAGMA_RESULTS = BASE / "REGA_separate/06_magma/results/merged/AD_Bellenguez_hg19_merged.genes.out"
SINGLE_ASSOC_PATH = BASE / "REGA_separate/07_tf/results/activity_pathology/all_celltype_region_TF_pathology_association.tsv"
CROSS_ASSOC_PATH = BASE / "REGA_shared/part7_TF/04DriverTF_lm/all_celltype_region_TF_pathology_association.tsv"
SINGLE_ASSIGN_ROOT = BASE / "REGA_separate/07_tf/results/tf_assignment"
CROSS_ASSIGN_ROOT = BASE / "REGA_shared/part7_TF/02tf_assignment_FC"

magma = pd.read_csv(MAGMA_RESULTS, sep=r"\s+")
magma["source_path"] = str(MAGMA_RESULTS)
magma["source_row"] = range(2, len(magma) + 2)
snp_tfs = {line.strip().upper() for line in SNP_LIST_PATH.read_text().splitlines() if line.strip()}

def read_assoc(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t")
    required = {"celltype_region", "TF", "lm_beta_pathology", "lm_fdr_pathology",
                "lateAD_minus_nonAD", "mannwhitney_lateAD_vs_nonAD_fdr"}
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"{path} lacks columns: {sorted(missing)}")
    d["TF"] = d["TF"].astype(str).str.upper()
    for col in required - {"celltype_region", "TF"}:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["positive_lm_fdr05"] = (d.lm_beta_pathology > 0) & (d.lm_fdr_pathology < 0.05)
    d["positive_mw_fdr05"] = (d.lateAD_minus_nonAD > 0) & (d.mannwhitney_lateAD_vs_nonAD_fdr < 0.05)
    d["is_driver"] = d.positive_lm_fdr05 | d.positive_mw_fdr05
    return d


def assignment_map(path: Path, entity_col: str, pair: str | None = None) -> tuple[dict, dict]:
    d = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if not {entity_col, "assigned_module"}.issubset(d.columns):
        raise ValueError(f"{path} must contain {entity_col} and assigned_module")
    if pair is not None:
        if "celltype_region" not in d.columns:
            raise ValueError(f"{path} must contain celltype_region")
        d = d[d.celltype_region == pair]
    d[entity_col] = d[entity_col].str.strip().str.upper()
    d["assigned_module"] = d["assigned_module"].str.strip()
    d = d[d[entity_col].ne("")]
    module_map, row_map = {}, {}
    for entity, rows in d.groupby(entity_col, sort=False):
        modules = set(rows.assigned_module[rows.assigned_module.ne("")])
        if len(modules) == 1:
            module_map[entity] = next(iter(modules))
            row_map[entity] = int(rows.index[0]) + 2
    return module_map, row_map


def collect_module_drivers(scope: str, module_ids: list[str]) -> pd.DataFrame:
    if scope == "celltype-specific":
        stats_path, assign_root = SINGLE_ASSOC_PATH, SINGLE_ASSIGN_ROOT
    else:
        stats_path, assign_root = CROSS_ASSOC_PATH, CROSS_ASSIGN_ROOT
    stats = read_assoc(stats_path)
    rows = []
    for module_id in module_ids:
        context_key, module = module_id.rsplit("_", 1)
        if scope == "celltype-specific":
            pair = context_key
            assignment_path = assign_root / f"{pair}_TF_module_assignments.tsv"
            module_map, assignment_rows = assignment_map(assignment_path, "entity")
            pair_stats = stats[stats.celltype_region == pair]
            output_region = pair.rsplit("_", 1)[-1]
            output_pair = pair
        else:
            region = context_key
            assignment_path = assign_root / region / f"{region}_ALL_cell_TF_FC_assignment_with_unassigned.tsv"
            pair_stats = stats[stats.celltype_region.str.endswith("_" + region)]
            output_region = region
            module_map, assignment_rows = {}, {}
        for _, stat in pair_stats[pair_stats.is_driver].iterrows():
            pair = str(stat.celltype_region)
            if scope == "cross-celltype":
                # The region-wide assignment table has one row per celltype_region;
                # fetch the provenance for the exact pair rather than another cell.
                pair_map, pair_rows = assignment_map(assignment_path, "tf", pair=pair)
                assigned_module = pair_map.get(stat.TF)
                assignment_source_row = pair_rows.get(stat.TF)
            else:
                assigned_module = module_map.get(stat.TF)
                assignment_source_row = assignment_rows.get(stat.TF)
            if assigned_module != module:
                continue
            evidence = []
            if stat.positive_lm_fdr05:
                evidence.append("positive linear pathology trend; lm_fdr_pathology < 0.05")
            if stat.positive_mw_fdr05:
                evidence.append("positive lateAD vs nonAD Mann–Whitney; FDR < 0.05")
            rows.append({
                "module_scope": scope,
                "module_id": module_id,
                "region": output_region,
                "celltype_region": pair,
                "module": module,
                "TF": stat.TF,
                "module_assignment_and_driver_overlap": True,
                "driver_evidence": " OR ".join(evidence),
                "positive_lm_fdr05": bool(stat.positive_lm_fdr05),
                "positive_mw_fdr05": bool(stat.positive_mw_fdr05),
                "lm_beta_pathology": stat.lm_beta_pathology,
                "lm_pval_pathology": getattr(stat, "lm_pval_pathology", float("nan")),
                "lm_fdr_pathology": stat.lm_fdr_pathology,
                "lateAD_minus_nonAD": stat.lateAD_minus_nonAD,
                "mannwhitney_lateAD_vs_nonAD_fdr": stat.mannwhitney_lateAD_vs_nonAD_fdr,
                "AD_SNP_nearby_TF": "Yes" if stat.TF in snp_tfs else "No",
                "source_path": str(stats_path),
                "source_row": int(stat.name) + 2,
                "assignment_source_path": str(assignment_path),
                "assignment_source_row": assignment_source_row,
            })
    return pd.DataFrame(rows)


single = collect_module_drivers("celltype-specific", CELLTYPE_MODULES)
cross = collect_module_drivers("cross-celltype", CROSS_MODULES)

# Append source MAGMA statistics where the TF occurs in the audit table.
magma["TF"] = magma["GENE"].astype(str).str.upper()
magma["P"] = pd.to_numeric(magma["P"], errors="coerce")
magma = magma.sort_values("P", na_position="last").drop_duplicates("TF", keep="first")
magma_cols = ["TF", "P", "ZSTAT", "source_path", "source_row"]
magma = magma[magma_cols].rename(columns={
    "P": "MAGMA_gene_P", "ZSTAT": "MAGMA_gene_ZSTAT",
    "source_path": "MAGMA_gene_P_source", "source_row": "MAGMA_gene_P_source_row",
})

def finalize(d):
    d = d.merge(magma, on="TF", how="left", validate="many_to_one")
    d["module_assignment_and_driver_overlap"] = True
    # Retain original source provenance, plus the exact MAGMA-list provenance.
    d["MAGMA_TF_SNP_list_source"] = str(SNP_LIST_PATH)
    cols = ["module_scope", "module_id", "region", "celltype_region", "module", "TF",
            "module_assignment_and_driver_overlap", "driver_evidence", "positive_lm_fdr05",
            "positive_mw_fdr05", "lm_beta_pathology", "lm_pval_pathology", "lm_fdr_pathology",
            "lateAD_minus_nonAD", "mannwhitney_lateAD_vs_nonAD_fdr", "AD_SNP_nearby_TF",
            "MAGMA_gene_P", "MAGMA_gene_ZSTAT", "source_path", "source_row",
            "assignment_source_path", "assignment_source_row",
            "MAGMA_gene_P_source", "MAGMA_gene_P_source_row", "MAGMA_TF_SNP_list_source"]
    return d[cols].drop_duplicates().sort_values(["module_id", "celltype_region", "TF"])

single_out = finalize(single)
cross_out = finalize(cross)
if not FIGURES_ONLY:
    single_out.to_csv(TABLES / "celltype_module_driver_TFs.tsv", sep="\t", index=False, na_rep="NA")
    cross_out.to_csv(TABLES / "cross_celltype_module_driver_TFs.tsv", sep="\t", index=False, na_rep="NA")

def module_inventory(d, ids, path):
    counts = d.groupby("module_id").size().to_dict()
    snp_counts = d[d.AD_SNP_nearby_TF.eq("Yes")].groupby("module_id").size().to_dict()
    pd.DataFrame({"module_id": ids,
                  "n_driver_TFs_assigned": [int(counts.get(x, 0)) for x in ids],
                  "n_driver_TFs_in_MAGMA_TF_SNP_list": [int(snp_counts.get(x, 0)) for x in ids],
                  "has_overlapping_driver_TF": ["Yes" if counts.get(x, 0) else "No" for x in ids]}).to_csv(path, sep="\t", index=False)

if not FIGURES_ONLY:
    module_inventory(single_out, CELLTYPE_MODULES, TABLES / "celltype_module_inventory.tsv")
    module_inventory(cross_out, CROSS_MODULES, TABLES / "cross_celltype_module_inventory.tsv")

def wrap_names(names, max_chars=18):
    lines, row, n = [], [], 0
    for name in names:
        add = len(name) + (2 if row else 0)
        if row and n + add > max_chars:
            lines.append(row)
            row, n = [], 0
        row.append(name)
        n += len(name) + (2 if len(row) > 1 else 0)
    if row:
        lines.append(row)
    return lines or [["No overlapping driver TF"]]

def draw_boxes(df, module_ids, scope, output_stem, cross_mode=False):
    grouped = {m: df[df.module_id.eq(m)] for m in module_ids}
    cards = []
    for module in module_ids:
        sub = grouped[module]
        context_lines = []
        if cross_mode:
            contexts = sorted(sub.celltype_region.dropna().unique())
            if not contexts:
                context_lines.append(("No celltype-region overlap", [["No overlapping driver TF"]]))
            for context in contexts:
                names = sorted(sub.loc[sub.celltype_region.eq(context), "TF"].unique())
                context_lines.append((context, wrap_names(names)))
        else:
            context = str(sub.celltype_region.iloc[0]) if not sub.empty else module.rsplit("_", 1)[0]
            names = sorted(sub.TF.unique()) if not sub.empty else []
            context_lines.append((context, wrap_names(names)))
        line_count = sum(1 + len(rows) for _, rows in context_lines)
        cards.append((module, context_lines, max(4.0, 2.35 + line_count * 0.92)))

    ncols = 2
    colcards = [[], []]
    heights = [0.0, 0.0]
    for card in cards:
        col = 0 if heights[0] <= heights[1] else 1
        colcards[col].append(card)
        heights[col] += card[2] + 0.46
    fig_h = max(10, max(heights) * 0.55)
    # Keep the two-column layout and scale its physical width by 1.5 so each
    # module box is half again as wide in the exported figure.
    fig, ax = plt.subplots(figsize=(15, fig_h))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, max(heights))
    ax.axis("off")
    color_no = "#222222"
    color_yes = "#B22222"
    fs = 14.4
    max_y = max(heights)
    for col, collist in enumerate(colcards):
        x0 = 0.025 + col * 0.49
        width = 0.455
        y = max_y - 0.3
        for module, context_lines, card_h in collist:
            ytop = y
            ybottom = ytop - card_h
            ax.add_patch(Rectangle((x0, ybottom), width, card_h,
                                   facecolor="#FAFBFC", edgecolor="#52616B", linewidth=1.6))
            module_color = "#B7791F" if cross_mode and module in EXPLORATORY_CROSS_MODULES else "#16324F"
            ax.text(x0 + 0.012, ytop - 0.52, module, fontsize=20, fontweight="bold", va="top", color=module_color)
            # Increase the vertical gap between the module title and its
            # celltype-region subtitle while keeping the subtitle compact.
            yy = ytop - 1.55
            for context, lines in context_lines:
                ax.text(x0 + 0.014, yy, context, fontsize=16,
                        fontweight="bold", va="top", color="#405261")
                yy -= 0.68
                for names in lines:
                    xx = x0 + 0.025
                    for i, name in enumerate(names):
                        color = color_yes if name.upper() in snp_tfs else color_no
                        t = ax.text(xx, yy, name, fontsize=fs, va="top", color=color,
                                    fontweight="bold" if color == color_yes else "normal")
                        # Approximate text advance in axis-coordinate units; wrap_names
                        # keeps the rendered line within each card.
                        xx += (len(name) + 1.5) * 0.0135
                    yy -= 0.72
                yy -= 0.2
            y = ybottom - 0.46
    title = "Cell-type-specific disease modules" if not cross_mode else "Cross-celltype disease modules"
    fig.suptitle(title + "\nDriver TFs assigned to each module; red labels are in the MAGMA FDR-significant TF list",
                 fontsize=28, fontweight="bold", y=0.995)
    legend_handles = [Patch(facecolor=color_no, label="Not in MAGMA FDR-significant TF list"),
                      Patch(facecolor=color_yes, label="In MAGMA FDR-significant TF list")]
    if cross_mode and EXPLORATORY_CROSS_MODULES:
        legend_handles.append(Patch(facecolor="#B7791F", edgecolor="#B7791F",
                                    label="Amber module name: exploratory; not selected by cross-celltype FDR screen"))
    fig.legend(handles=legend_handles, loc="lower center", ncol=1,
               frameon=False, bbox_to_anchor=(0.5, 0.005), fontsize=20)
    fig.subplots_adjust(top=0.92, bottom=0.12, left=0.02, right=0.98)
    for ext in ["pdf", "png"]:
        fig.savefig(FIGURES / f"{output_stem}.{ext}", dpi=250, bbox_inches="tight")
    plt.close(fig)

draw_boxes(single_out, CELLTYPE_MODULES, "celltype-specific", "celltype_specific_disease_module_driver_TFs")
draw_boxes(cross_out, CROSS_MODULES, "cross-celltype", "cross_celltype_disease_module_driver_TFs", cross_mode=True)

print(f"MAGMA FDR-significant TF list: {len(snp_tfs)} TFs from {SNP_LIST_PATH}")
print(f"Cell-type-specific: {len(CELLTYPE_MODULES)} modules, {len(single_out)} overlapping driver-TF rows, {int((single_out.AD_SNP_nearby_TF == 'Yes').sum())} MAGMA-list overlaps")
print(f"Cross-celltype: {len(CROSS_MODULES)} modules, {len(cross_out)} overlapping driver-TF rows, {int((cross_out.AD_SNP_nearby_TF == 'Yes').sum())} MAGMA-list overlaps")
