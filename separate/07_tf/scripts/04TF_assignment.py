import os
import numpy as np
import pandas as pd

REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
CELLTYPES = ["Ast", "Exc", "OPC", "Oli", "Inh", "Mic_Immune"]

BASE_DIR = os.environ.get(
    "REGA_TF_BY_MODULE_ROOT",
    "/project/zduren/durenlab/naqing_palmetto/REGA_separate/07_tf/results/tf_by_module",
)
OUT_DIR = os.environ.get(
    "REGA_TF_ASSIGNMENT_ROOT",
    "/project/zduren/durenlab/naqing_palmetto/REGA_separate/07_tf/results/tf_assignment",
)
os.makedirs(OUT_DIR, exist_ok=True)

TOPK_FRAC = 0.10
P_ABS = 0.10

def normalize_l1l1(df: pd.DataFrame) -> pd.DataFrame:
    x = df.copy().astype(float)

    col_sum = x.sum(axis=0)
    col_sum = col_sum.replace(0, np.nan).fillna(1.0)
    x = x.div(col_sum, axis=1)

    row_sum = x.sum(axis=1)
    row_sum = row_sum.replace(0, np.nan).fillna(1.0)
    x = x.div(row_sum, axis=0)

    return x.fillna(0.0)

def assign_topk_sets(M: pd.DataFrame, k: int):
    out = {}
    for mod in M.columns:
        s = M[mod].replace([np.inf, -np.inf], np.nan).dropna()
        out[mod] = list(s.sort_values(ascending=False).head(k).index)
    return out

def harden_sets_by_score(M_score: pd.DataFrame, set_list: dict):
    out = {mod: [] for mod in M_score.columns}

    all_entities = sorted(set(e for v in set_list.values() for e in v))
    all_entities = [e for e in all_entities if e in M_score.index]

    for tf in all_entities:
        cand_mods = [mod for mod, ents in set_list.items() if tf in ents]
        if not cand_mods:
            continue

        scores = M_score.loc[tf, cand_mods].astype(float)
        scores = scores.replace([np.inf, -np.inf], np.nan).fillna(-np.inf)

        if (scores == -np.inf).all():
            continue

        best_mod = scores.idxmax()
        out[best_mod].append(tf)

    for mod in out:
        out[mod] = sorted(set(out[mod]))
    return out

def assign_label_absolute(M: pd.DataFrame, p: float = 0.10):
    M2 = normalize_l1l1(M)
    S = M2.replace([np.inf, -np.inf], np.nan).fillna(-np.inf)

    best_col = S.values.argmax(axis=1)
    best_score = S.values[np.arange(S.shape[0]), best_col]

    thr = p
    labels = np.array([S.columns[i] if score >= thr else None for i, score in zip(best_col, best_score)], dtype=object)

    return labels, M2, thr

def build_union_module_sets(M: pd.DataFrame, topk_frac=0.10, p=0.10):
    M_norm = normalize_l1l1(M)
    topk = max(1, int(np.ceil(M.shape[0] * topk_frac)))

    topk_sets = assign_topk_sets(M_norm, topk)

    labels, _, _ = assign_label_absolute(M, p=p)
    label_sets = {mod: [] for mod in M_norm.columns}
    for tf, lab in zip(M_norm.index, labels):
        if lab is not None:
            label_sets[lab].append(tf)

    union_sets = {}
    for mod in M_norm.columns:
        union_sets[mod] = sorted(set(topk_sets.get(mod, [])).union(label_sets.get(mod, [])))

    return union_sets, M_norm, topk

def assign_tf_to_modules(tf_by_module: pd.DataFrame, topk_frac=0.10, p=0.10):
    union_sets, M_norm, topk = build_union_module_sets(tf_by_module, topk_frac=topk_frac, p=p)
    hardened_sets = harden_sets_by_score(M_norm, union_sets)

    tf2module = pd.Series(index=M_norm.index, dtype="object")
    tf2module[:] = None

    for mod, tfs in hardened_sets.items():
        tf2module.loc[tfs] = mod

    counts = pd.Series({mod: len(tfs) for mod, tfs in hardened_sets.items()}, name="n")
    counts.index.name = "module"

    return tf2module, counts, hardened_sets

for ct in CELLTYPES:
    for region in REGIONS:
        name = f"{ct}_{region}"
        in_file = os.path.join(BASE_DIR, f"{name}_TF_by_module.tsv")

        if not os.path.exists(in_file):
            print(f"[SKIP] missing {in_file}")
            continue

        print(f"[READ] {in_file}")
        M = pd.read_csv(in_file, sep="\t", index_col=0)
        M = M.apply(pd.to_numeric, errors="coerce").fillna(0.0)

        tf2module, counts, hardened_sets = assign_tf_to_modules(
            M,
            topk_frac=TOPK_FRAC,
            p=P_ABS,
        )

        out_assign = os.path.join(OUT_DIR, f"{name}_TF_module_assignments.tsv")
        out_counts = os.path.join(OUT_DIR, f"{name}_TF_module_counts.tsv")
        out_sets = os.path.join(OUT_DIR, f"{name}_TF_sets.txt")

        pd.DataFrame({
            "entity": tf2module.index,
            "assigned_module": tf2module.values
        }).to_csv(out_assign, sep="\t", index=False)

        counts.reset_index().to_csv(out_counts, sep="\t", index=False)

        with open(out_sets, "w") as f:
            for mod in M.columns:
                tfs = hardened_sets.get(mod, [])
                if len(tfs) == 0:
                    continue
                f.write(f"TF_{mod}\t" + "\t".join(tfs) + "\n")

        print(f"[DONE] {name} -> {out_assign}")
