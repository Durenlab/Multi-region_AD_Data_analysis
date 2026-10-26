import os
import pandas as pd

# input file
path = (os.environ["REGA_SEPARATE_ROOT"] + '/07_tf/data/Motif_TF_name_correspondance_uppercase.txt')

# read two columns: motif, TF
df = pd.read_csv(path, sep=r"\s+", header=None, names=["motif", "TF"])

# make binary motif x TF matrix
mat = pd.crosstab(df["motif"], df["TF"])
mat = (mat > 0).astype(int)

# save
motif_tf_file = os.environ.get(
    "REGA_MOTIF_TF_MATRIX",
    (os.environ["REGA_SEPARATE_ROOT"] + '/07_tf/data/motif_by_TF_matrix.tsv'),
)
os.makedirs(os.path.dirname(motif_tf_file), exist_ok=True)
mat.to_csv(motif_tf_file, sep="\t")

print(mat.shape)
print(mat.head())

REGIONS = ["AG", "EC", "HC", "MTC", "PFC", "TH"]
CELLTYPES = ["Ast", "Exc", "OPC", "Oli", "Inh", "Mic_Immune"]

base_dir = os.environ.get(
    "REGA_MATRICES_ROOT",
    (os.environ["REGA_SEPARATE_ROOT"] + '/03_process_rega/matrices'),
)
# motif x TF matrix was written above (or selected through REGA_MOTIF_TF_MATRIX)
out_dir = os.environ.get(
    "REGA_TF_BY_MODULE_ROOT",
    (os.environ["REGA_SEPARATE_ROOT"] + '/07_tf/results/tf_by_module'),
)
os.makedirs(out_dir, exist_ok=True)

# read motif-by-TF matrix
motif_tf = pd.read_csv(motif_tf_file, sep="\t", index_col=0)

# make sure values are numeric
motif_tf = motif_tf.apply(pd.to_numeric, errors="coerce").fillna(0)

for ct in CELLTYPES:
    for region in REGIONS:
        ct_region = f"{ct}_{region}"
        v_file = os.path.join(base_dir, ct_region, f"{ct_region}_V_TF_by_module.tsv")

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
