#!/usr/bin/env python
"""Widen fixed-feature CellChat plots and add an FDR-significant-only panel.

Uses finalized donor scores/model estimates; it does not rerun any statistics.
"""
from __future__ import annotations

import base64
import hashlib
import html
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "experiments" / "EXP002_results" / "finalized"
FIG = OUT / "figures"
TABLE = OUT / "tables"
REPORT = OUT / "report" / "cellchat_final_report.html"
PATHS = ["nonAD", "earlyAD", "lateAD"]
COLORS = {"nonAD": "#4477AA", "earlyAD": "#EE9944", "lateAD": "#CC6677"}
SEED = 20260929


def plot(df: pd.DataFrame, scores: pd.DataFrame, path: Path, title: str,
         ncols: int, width: float) -> None:
    nrows = max(1, math.ceil(len(df) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(width, max(4.5, 3.55 * nrows)),
                             squeeze=False)
    rng = np.random.default_rng(SEED)
    if df.empty:
        axes.flat[0].text(.5, .5, "No fixed selected features pass pair-specific BH q < 0.05",
                          transform=axes.flat[0].transAxes, ha="center", va="center", wrap=True)
        axes.flat[0].axis("off")
    for ax, row in zip(axes.flat, df.itertuples(index=False)):
        obs = scores[(scores.region == row.region) & (scores.source == row.source) &
                     (scores.target == row.target) &
                     (scores.interaction_id == row.interaction_id)]
        for x, pathology in enumerate(PATHS):
            values = obs.loc[obs.pathology == pathology, "score"].to_numpy(float)
            jitter = rng.normal(0, .055, len(values))
            ax.scatter(x + jitter, values, color=COLORS[pathology], alpha=.8,
                       s=21, edgecolor="white", linewidth=.3, zorder=3)
            ax.text(x, .98, f"n={len(values)}", transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=8.5)
        fitted = [row.prediction_nonAD_mean_covars,
                  row.prediction_earlyAD_mean_covars,
                  row.prediction_lateAD_mean_covars]
        if np.isfinite(fitted).all():
            ax.plot(range(3), fitted, color="#777777", linestyle="--",
                    linewidth=1.8, zorder=2, label="Simple linear pathology trend")
        qtext = "q=NA" if pd.isna(row.q_value_within_family) else f"BH q={row.q_value_within_family:.3g}"
        betatext = "β=NA" if pd.isna(row.beta_pathology_stage) else (
            f"β={row.beta_pathology_stage:.3g} [{row.ci_low:.3g}, {row.ci_high:.3g}]"
        )
        ax.set_title(f"{row.region}: {row.source} → {row.target}\n"
                     f"{row.interaction_id}\n{betatext}; {qtext}", fontsize=9)
        ax.set_xticks(range(3), PATHS, rotation=15, ha="right")
        ax.set_ylabel("CellChat probability")
        ax.grid(axis="y", color="#dddddd", linewidth=.55)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.flat[len(df):]:
        ax.set_visible(False)
    fig.suptitle(title, fontsize=15, y=1.003)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not REPORT.exists():
        raise FileNotFoundError(REPORT)
    selected = pd.read_csv(TABLE / "fixed_selected_features_refreshed.tsv", sep="\t")
    scores = pd.read_csv(TABLE / "fixed_selected_feature_donor_scores.tsv", sep="\t")
    required = {"prediction_nonAD_mean_covars", "prediction_earlyAD_mean_covars",
                "prediction_lateAD_mean_covars", "q_value_within_family"}
    if not required.issubset(selected.columns):
        raise ValueError(f"Missing model columns: {sorted(required - set(selected.columns))}")
    sig = selected[selected.q_value_within_family < .05].copy()
    plot(selected, scores, FIG / "selected_features",
         "Fixed selected CellChat features: donor scores and simple linear pathology trends",
         ncols=3, width=10.25)
    plot(sig, scores, FIG / "selected_features_significant",
         "Fixed selected features with pair-specific BH q < 0.05",
         ncols=2 if len(sig) > 1 else 1, width=8.5 if len(sig) > 1 else 6.0)

    html_report = REPORT.read_text()
    old_title = "Fixed selected communication features"
    pattern = (r"(<h2>" + re.escape(old_title) +
               r"</h2>)<img alt='" + re.escape(old_title) +
               r"' src='data:image/png;base64,[^']+'>")
    image_all = base64.b64encode((FIG / "selected_features.png").read_bytes()).decode("ascii")
    image_sig = base64.b64encode((FIG / "selected_features_significant.png").read_bytes()).decode("ascii")
    block = (f"<h2>{old_title}</h2><img alt='{old_title}' src='data:image/png;base64,{image_all}'>"
             f"<h2>Within-family significant selected features (BH q &lt; 0.05)</h2>"
             f"<img alt='Within-family significant selected features' src='data:image/png;base64,{image_sig}'>")
    updated, n = re.subn(pattern, block, html_report, count=1)
    if n != 1:
        raise ValueError("Could not locate the selected-features image block in the HTML report")
    REPORT.write_text(updated)

    legends = OUT / "methods" / "figure_legends.md"
    text = legends.read_text()
    text = text.replace("the fitted donor-level progression model at the mean covariate profile.",
                        "the fitted donor-level progression model at the mean covariate profile (gray dashed line).")
    if "**Figure 1b." not in text:
        text += ("\n\n**Figure 1b. Within-family significant selected features.** This companion figure shows only fixed-list features with BH q<0.05 within their region × directed sender/receiver × progression family. Points are eligible donors; the gray dashed line is the adjusted linear pathology trend. Donor, cell-count eligibility, HC3 model, covariates, and bootstrap definitions are as in Figure 1.\n")
    legends.write_text(text)

    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    revision = {
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "plot_script": str(Path(__file__).resolve()),
        "plot_script_sha256": script_hash,
        "selected_feature_count": int(len(selected)),
        "within_family_significant_selected_count": int(len(sig)),
        "all_selected_layout_inches": {"width": 10.25, "columns": 3},
        "significant_only_layout_inches": {"width": 8.5 if len(sig) > 1 else 6.0, "columns": 2 if len(sig) > 1 else 1},
        "fitted_slope_line": "gray dashed",
        "statistical_results_recomputed": False,
    }
    (OUT / "config" / "plot_revision.json").write_text(json.dumps(revision, indent=2) + "\n")
    config_path = OUT / "config" / "analysis_config.json"
    config = json.loads(config_path.read_text())
    config["plot_revision"] = revision
    config_path.write_text(json.dumps(config, indent=2, default=str) + "\n")
    print("PLOT_UPDATE_COMPLETE")
    print("fixed selected features:", len(selected))
    print("pair-specific significant selected features:", len(sig))
    print("figures: selected_features.pdf/png and selected_features_significant.pdf/png")
    print("report refreshed:", REPORT)


if __name__ == "__main__":
    main()
