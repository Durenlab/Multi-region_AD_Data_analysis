from pathlib import Path
import os
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ.get("REGA_SEPARATE_ROOT", "/project/zduren/durenlab/naqing_palmetto/REGA_separate"))
import os
import pandas as pd

REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
if os.environ.get("REGA_REGION"):
    if os.environ["REGA_REGION"] not in REGIONS:
        raise ValueError(f"Unknown REGA_REGION: {os.environ['REGA_REGION']}")
    REGIONS = [os.environ["REGA_REGION"]]
CELLTYPES = ["Ast", "Exc", "OPC", "Oli", "Inh", "Mic_Immune"]

base_dir = str(SHARED_ROOT / 'part3_process_REGA_output/09extracted_matrices_region_shareH')
motif_tf_file = str(SHARED_ROOT / "part7_TF/motif_by_TF_matrix.tsv")
out_dir = str(SHARED_ROOT / 'part7_TF/01tf_by_module')
os.makedirs(out_dir, exist_ok=True)

# read motif-by-TF matrix
motif_tf = pd.read_csv(motif_tf_file, sep="\t", index_col=0)

# make sure values are numeric
motif_tf = motif_tf.apply(pd.to_numeric, errors="coerce").fillna(0)

for ct in CELLTYPES:
    for region in REGIONS:
        ct_region = f"{ct}_{region}"
        v_file = os.path.join(base_dir, region, f"{ct_region}_V_TF_by_module.tsv")

        if not os.path.exists(v_file):
            print(f"Missing: {v_file}")
            continue

        # read motif-by-module matrix
        motif_by_module = pd.read_csv(v_file, sep="\t", index_col=0)
        motif_by_module = motif_by_module.apply(pd.to_numeric, errors="coerce").fillna(0)

        # align motifs between the two matrices
        common_motifs = motif_tf.index.intersection(motif_by_module.index)

        if len(common_motifs) == 0:
            print(f"No common motifs for {ct_region}")
            continue

        motif_tf_sub = motif_tf.loc[common_motifs]
        motif_by_module_sub = motif_by_module.loc[common_motifs]

        # TF x motif @ motif x module -> TF x module
        tf_by_module = motif_tf_sub.T @ motif_by_module_sub

        out_file = os.path.join(out_dir, f"{ct_region}_TF_by_module.tsv")
        tf_by_module.to_csv(out_file, sep="\t")

        print(f"Saved {ct_region}: {tf_by_module.shape} -> {out_file}")
