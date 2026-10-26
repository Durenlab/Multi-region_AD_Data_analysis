import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse

# ============================================================
# Paths
# ============================================================
BASE_DIR = (os.environ["REGA_SEPARATE_ROOT"] + '/')
INPUT_DIR = os.environ.get("REGA_STAGE2_ROOT", os.path.join(BASE_DIR, "02_run_rega/rega_output"))
CONTEXT_SIGNAL_PATH = "selected_context.txt"
TF_binding_DIR = os.environ.get("REGA_INPUTS_ROOT", os.path.join(BASE_DIR, "01_prepare/data/rega_inputs"))
OUT_DIR = os.environ.get("REGA_MATRICES_ROOT", os.path.join(BASE_DIR, "03_process_rega/matrices"))
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# Helpers
# ============================================================

def parse_celltype_region(celltype_region):
    """
    Parse names like:
        Ast_EC
        Oli_MTC
        Mic_Immune_MTC

    Region is the last underscore-separated field.
    Celltype is everything before the last underscore.
    """
    parts = celltype_region.split("_")
    region = parts[-1]
    celltype = "_".join(parts[:-1])
    return celltype, region


def matlab_string_array_to_list(arr):
    """
    Convert MATLAB-loaded string/cell arrays into a clean Python list of strings.

    Handles arrays like:
        shape (1, n)
        shape (n, 1)
        object arrays
        string arrays
    """
    arr = np.asarray(arr).squeeze()

    out = []
    for x in arr:
        if isinstance(x, np.ndarray):
            x = x.squeeze()
            if x.size == 1:
                x = x.item()
        if isinstance(x, bytes):
            x = x.decode("utf-8")
        out.append(str(x).strip())

    return out


def read_name_file(path):
    """
    Read one-column or multi-column text/tsv/csv file and return first column as names.
    """
    df = pd.read_csv(path, sep=None, engine="python", header=None)
    names = df.iloc[:, 0].dropna().astype(str).str.strip().tolist()
    return names


def find_tf_names(group_input_dir, expected_n):
    input_binding_path = os.path.join(group_input_dir, "input_binding.csv")

    if not os.path.exists(input_binding_path):
        print(f"  TF name file not found: {input_binding_path}")
        return [f"TF_{i + 1}" for i in range(expected_n)]

    cols = pd.read_csv(input_binding_path, nrows=0).columns.astype(str).tolist()

    if len(cols) == expected_n:
        tf_names = cols

    elif len(cols) == expected_n + 1:
        # First column is probably RE ID / index column
        tf_names = cols[1:]

    else:
        print(
            f"  Warning: TF name length mismatch from {input_binding_path}. "
            f"Expected {expected_n}, got {len(cols)} columns."
        )
        return [f"TF_{i + 1}" for i in range(expected_n)]

    print(f"  Using TF names from header of: {input_binding_path}")
    return tf_names


def find_motif_names(group_input_dir, expected_n):
    """Use binding-matrix headers for B columns, or numbered motif names."""
    input_binding_path = os.path.join(group_input_dir, "input_binding.csv")
    if not os.path.exists(input_binding_path):
        print(f"  Motif name file not found: {input_binding_path}")
        return [f"Motif_{i + 1}" for i in range(expected_n)]

    columns = pd.read_csv(input_binding_path, nrows=0).columns.astype(str).tolist()
    if len(columns) == expected_n:
        return columns
    if len(columns) == expected_n + 1:
        return columns[1:]
    print(
        f"  Warning: motif name length mismatch in {input_binding_path}. "
        f"Expected {expected_n}, got {len(columns)} columns."
    )
    return [f"Motif_{i + 1}" for i in range(expected_n)]


def check_names(names, expected_n, label):
    """
    Check name length. If mismatched, fall back to generic names.
    """
    if names is None or len(names) != expected_n:
        print(
            f"  Warning: {label} name length mismatch. "
            f"Expected {expected_n}, got {0 if names is None else len(names)}. "
            f"Using generic {label}_1, {label}_2, ..."
        )
        return [f"{label}_{i + 1}" for i in range(expected_n)]
    return names


def write_matrix_tsv(matrix, out_path, row_names, col_names):
    """
    Write dense or sparse matrix to TSV.
    """
    if scipy.sparse.issparse(matrix):
        matrix = matrix.toarray()

    matrix = np.asarray(matrix)

    if matrix.shape[0] != len(row_names):
        raise ValueError(
            f"Row name mismatch for {out_path}: "
            f"matrix has {matrix.shape[0]} rows, names has {len(row_names)}"
        )

    if matrix.shape[1] != len(col_names):
        raise ValueError(
            f"Column name mismatch for {out_path}: "
            f"matrix has {matrix.shape[1]} columns, names has {len(col_names)}"
        )

    df = pd.DataFrame(matrix, index=row_names, columns=col_names)
    df.to_csv(out_path, sep="\t")


def write_binary_sparse_tsv(matrix, out_path, row_names, col_names, chunk_rows=100):
    """Write positive W/B entries as binary row/col/value triplets."""
    if scipy.sparse.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix)
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2D matrix for {out_path}, got {matrix.ndim}D")
    if matrix.shape[0] != len(row_names):
        raise ValueError(
            f"Row name mismatch for {out_path}: "
            f"matrix has {matrix.shape[0]} rows, names has {len(row_names)}"
        )
    if matrix.shape[1] != len(col_names):
        raise ValueError(
            f"Column name mismatch for {out_path}: "
            f"matrix has {matrix.shape[1]} columns, names has {len(col_names)}"
        )

    temporary_path = out_path + ".tmp"
    nonzero_count = 0
    try:
        with open(temporary_path, "w") as handle:
            handle.write("row\tcol\tvalue\n")
            for start in range(0, matrix.shape[0], chunk_rows):
                end = min(start + chunk_rows, matrix.shape[0])
                block = matrix[start:end]
                negative_count = np.count_nonzero(block < 0)
                if negative_count:
                    raise ValueError(
                        f"Found {negative_count} negative values in rows "
                        f"{start}:{end} of {out_path}"
                    )
                local_rows, local_cols = np.nonzero(block > 0)
                if not len(local_rows):
                    continue
                sparse_df = pd.DataFrame({
                    "row": [row_names[start + i] for i in local_rows],
                    "col": [col_names[j] for j in local_cols],
                    "value": np.ones(len(local_rows), dtype=np.uint8),
                })
                sparse_df.to_csv(handle, sep="\t", index=False, header=False)
                nonzero_count += len(local_rows)
        os.replace(temporary_path, out_path)
    except Exception:
        if os.path.exists(temporary_path):
            os.remove(temporary_path)
        raise
    print(f"  Saved {nonzero_count:,} nonzero entries to {out_path}")


# ============================================================
# Discover all available region/cell type REGA inputs.
celltype_regions = []
for region in sorted(os.listdir(INPUT_DIR)):
    region_dir = os.path.join(INPUT_DIR, region)
    if not os.path.isdir(region_dir):
        continue
    for celltype in sorted(os.listdir(region_dir)):
        if celltype == "Vasc":
            continue
        stage2 = os.path.join(region_dir, celltype, "stage2.mat")
        if os.path.isfile(stage2):
            celltype_regions.append(f"{celltype}_{region}")

print(f"Found {len(celltype_regions)} cell type-region inputs.")


# Main loop
# ============================================================

for celltype_region in celltype_regions:

    if celltype_region == "":
        continue

    celltype, region = parse_celltype_region(celltype_region)

    print(f"\nProcessing {celltype_region}")
    print(f"  celltype = {celltype}")
    print(f"  region   = {region}")

    group_input_dir = os.path.join(INPUT_DIR, region, celltype)
    mat_path = os.path.join(group_input_dir, "stage2.mat")

    if not os.path.exists(mat_path):
        print(f"  Skipping: missing file {mat_path}")
        continue

    mat = scipy.io.loadmat(mat_path)

    required_keys = ["W", "B", "V", "H"]
    missing = [k for k in required_keys if k not in mat]

    if len(missing) > 0:
        print(f"  Skipping: missing keys in stage2.mat: {missing}")
        continue

    W = mat["W"]
    B = mat["B"]
    V = mat["V"]
    H = mat["H"]

    print(f"  W shape: {W.shape}")
    print(f"  B shape: {B.shape}")
    print(f"  V shape: {V.shape}")
    print(f"  H shape: {H.shape}")

    # ------------------------------------------------------------
    # Extract real names
    # ------------------------------------------------------------

    gene_names = None
    re_names = None
    sample_names = None

    if "W_rows" in mat:
        gene_names = matlab_string_array_to_list(mat["W_rows"])
    elif "X_rows" in mat:
        gene_names = matlab_string_array_to_list(mat["X_rows"])

    if "W_cols" in mat:
        re_names = matlab_string_array_to_list(mat["W_cols"])

    if "X_cols" in mat:
        sample_names = matlab_string_array_to_list(mat["X_cols"])

    gene_names = check_names(gene_names, W.shape[0], "Gene")
    re_names = check_names(re_names, W.shape[1], "RE")
    sample_names = check_names(sample_names, H.shape[1], "Sample")

    tf_binging_input_dir = os.path.join(TF_binding_DIR, region, celltype)
    tf_names = find_tf_names(tf_binging_input_dir, V.shape[0])
    tf_names = check_names(tf_names, V.shape[0], "TF")
    motif_names = find_motif_names(tf_binging_input_dir, B.shape[1])

    module_names = [f"M{i + 1}" for i in range(V.shape[1])]

    # ------------------------------------------------------------
    # Check dimensions
    # ------------------------------------------------------------

    if W.shape[1] != B.shape[0]:
        print(
            f"  Skipping: dimension mismatch W @ B. "
            f"W shape {W.shape}, B shape {B.shape}"
        )
        continue

    if B.shape[1] != V.shape[0]:
        print(
            f"  Skipping: dimension mismatch B @ V. "
            f"B shape {B.shape}, V shape {V.shape}"
        )
        continue

    if H.shape[0] != V.shape[1]:
        print(
            f"  Warning: H rows != module number from V. "
            f"H shape {H.shape}, V shape {V.shape}"
        )

    # ------------------------------------------------------------
    # Compute BV and WBV
    # ------------------------------------------------------------

    BV = B @ V
    WBV = W @ BV

    print(f"  BV shape:  {BV.shape}")
    print(f"  WBV shape: {WBV.shape}")

    # ------------------------------------------------------------
    # Output folder
    # ------------------------------------------------------------

    group_out_dir = os.path.join(OUT_DIR, celltype_region)
    os.makedirs(group_out_dir, exist_ok=True)

    # ------------------------------------------------------------
    # Write matrices
    # ------------------------------------------------------------

    write_matrix_tsv(
        WBV,
        os.path.join(group_out_dir, f"{celltype_region}_WBV_genes_by_module.tsv"),
        row_names=gene_names,
        col_names=module_names
    )

    write_matrix_tsv(
        BV,
        os.path.join(group_out_dir, f"{celltype_region}_BV_RE_by_module.tsv"),
        row_names=re_names,
        col_names=module_names
    )

    write_matrix_tsv(
        V,
        os.path.join(group_out_dir, f"{celltype_region}_V_TF_by_module.tsv"),
        row_names=tf_names,
        col_names=module_names
    )

    write_matrix_tsv(
        H,
        os.path.join(group_out_dir, f"{celltype_region}_H_module_by_sample.tsv"),
        row_names=module_names,
        col_names=sample_names
    )

    write_binary_sparse_tsv(
        W,
        os.path.join(group_out_dir, f"{celltype_region}_W_binary_sparse.tsv"),
        row_names=gene_names,
        col_names=re_names,
    )
    write_binary_sparse_tsv(
        B,
        os.path.join(group_out_dir, f"{celltype_region}_B_binary_sparse.tsv"),
        row_names=re_names,
        col_names=motif_names,
    )

    print(f"  Saved matrices to: {group_out_dir}")


print("\nDone.")
print(f"All outputs saved under: {OUT_DIR}")
