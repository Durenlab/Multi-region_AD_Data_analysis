#!/usr/bin/env python3
"""Create EXP008 pair and module-assignment preflight tables."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ADREGA = Path("/project/zduren/durenlab/naqing_palmetto/AD_REGA")
DRIVER_FILE = ADREGA / "part7_TF/02TF_pathology_association_lm/all_celltype_region_TF_pathology_association.tsv"
GEX_ROOT = ADREGA / "part1_prepare_REGA/04cleaned_gex"
GENE_MODULE_ROOT = ADREGA / "part3_process_REGA_output/09module_assignments"
TF_MODULE_ROOT = ADREGA / "part7_TF/04tf_assignment"
GRN_ROOT = Path("/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/cell_type/GRN")
PATHOLOGY = Path("/project/zduren/durenlab/naqing_palmetto/AD_Driver_TF/individual_pathology.tsv")
OUT = ROOT / "EXP008_results/audit"


def assignment_map(path: Path):
    if not path.is_file():
        return None, "missing_file"
    tab = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if not {"entity", "assigned_module"}.issubset(tab.columns):
        return None, "invalid_schema"
    tab["entity"] = tab.entity.str.strip()
    tab["assigned_module"] = tab.assigned_module.str.strip()
    tab = tab.loc[tab.entity.ne("")]
    nmods = tab.groupby("entity").assigned_module.nunique()
    conflicts = nmods[nmods > 1].index.tolist()
    if conflicts:
        return None, "conflicting_assignment:" + ";".join(conflicts[:5])
    return dict(zip(tab.entity, tab.assigned_module)), "ok"


def usable(module):
    return bool(module) and module.upper() not in {"NA", "NAN", "NULL"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (ROOT / "EXP008_results/logs").mkdir(parents=True, exist_ok=True)
    if not DRIVER_FILE.is_file() or not PATHOLOGY.is_file():
        raise SystemExit("Missing global driver or pathology input")
    drivers = pd.read_csv(DRIVER_FILE, sep="\t", dtype=str)
    drivers["lm_fdr_pathology"] = pd.to_numeric(drivers["lm_fdr_pathology"], errors="coerce")
    drivers = drivers.loc[drivers.lm_fdr_pathology.lt(0.05)].drop_duplicates(["celltype_region", "TF"])
    gex_pairs = sorted(p.stem.removeprefix("cleaned_") for p in GEX_ROOT.glob("cleaned_*.tsv"))
    if not gex_pairs:
        raise SystemExit(f"No cleaned expression matrices under {GEX_ROOT}")

    module_rows = []
    pair_rows = []
    for pair in gex_pairs:
        pair_drivers = drivers.loc[drivers.celltype_region.eq(pair)]
        gm_path = GENE_MODULE_ROOT / pair / "genes_module_assignments.tsv"
        tm_path = TF_MODULE_ROOT / f"{pair}_TF_module_assignments.tsv"
        gm, gm_status = assignment_map(gm_path)
        tm, tm_status = assignment_map(tm_path)
        grn_paths = [GRN_ROOT / f"{pair}_TG_TF.rds"] if len(pair_drivers) else []
        expression_path = GEX_ROOT / f"cleaned_{pair}.tsv"
        missing = [str(p) for p in [expression_path, PATHOLOGY, *grn_paths] if not p.is_file()]
        pair_rows.append({
            "pair": pair, "n_driver_TF": len(pair_drivers),
            "expression_exists": expression_path.is_file(),
            "TG_TF_exists": not grn_paths or grn_paths[0].is_file(),
            "gene_module_table_status": gm_status,
            "TF_module_table_status_QC_only": tm_status,
            "required_missing_inputs": ";".join(missing),
            "status": "ready" if not missing else "missing_required_input",
        })
        for tf in pair_drivers.TF.astype(str):
            tf_gene_module = gm.get(tf, "") if gm is not None else ""
            tf_external_module = tm.get(tf, "") if tm is not None else ""
            good_gm = usable(tf_gene_module)
            good_tm = usable(tf_external_module)
            module_rows.append({
                "pair": pair, "TF": tf,
                "same_gene_table_module": tf_gene_module,
                "same_gene_table_module_usable": good_gm,
                "separate_TF_table_module_QC_only": tf_external_module,
                "separate_TF_table_module_usable": good_tm,
                "module_assignments_differ": bool(good_gm and good_tm and tf_gene_module != tf_external_module),
                "restricted_mode_status": "ready" if good_gm else "unavailable_missing_TF_gene_module",
            })

    pairs = pd.DataFrame(pair_rows)
    module_audit = pd.DataFrame(module_rows)
    pairs.to_csv(OUT / "pairs.tsv", sep="\t", index=False)
    module_audit.to_csv(OUT / "driver_module_audit.tsv", sep="\t", index=False)
    print(f"Observed expression pairs: {len(pairs)}; significant driver TFs: {len(module_audit)}")
    print("Pair input status:", pairs.status.value_counts().to_dict())
    print("Restricted driver module status:", module_audit.restricted_mode_status.value_counts().to_dict())
    print("Separate TF-table module mismatches:", int(module_audit.module_assignments_differ.sum()))
    if (pairs.status != "ready").any():
        raise SystemExit(2)


if __name__ == "__main__":
    main()
