#!/usr/bin/env python
from pathlib import Path
import os
SHARED_ROOT = Path(os.environ["REGA_SHARED_ROOT"]) if "REGA_SHARED_ROOT" in os.environ else next(
    p for p in Path(__file__).resolve().parents if p.name in {"REGA_shared", "shared"}
)
SEPARATE_ROOT = Path(os.environ["REGA_SEPARATE_ROOT"])

import os
import argparse
import tempfile
import torch
import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

BASE_DIR = str(SHARED_ROOT)

REGA_OUT_ROOT = os.path.join(
    BASE_DIR,
    "part2_run_REGA/07REGA_region_shareH_output"
)

BLOCKS_ROOT = os.path.join(
    BASE_DIR,
    "part1_prepare_REGA/05REGA_blocks_region_consistent"
)

OUT_DIR = os.path.join(
    BASE_DIR,
    "part3_process_REGA_output/09extracted_matrices_region_shareH"

)

os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# Settings
# ============================================================

regions = ["PFC", "TH", "AG", "EC", "HC", "MTC"]
K = 20
outer_step = 3000


# ============================================================
# Helper functions
# ============================================================

def load_checkpoint(pt_path):
    try:
        return torch.load(pt_path, map_location="cpu")
    except Exception:
        return torch.load(pt_path, map_location="cpu", weights_only=False)


def to_str_array(x):
    return np.array([
        v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v)
        for v in x
    ])


def safe_name(x):
    return str(x).replace("/", "_")


def get_block_npz_path(celltype_region, chr_name):
    return os.path.join(
        BLOCKS_ROOT,
        f"{celltype_region}_blocks_aligned",
        f"{chr_name}.npz"
    )


def write_matrix_tsv(matrix, out_path, row_names, col_names):
    matrix = np.asarray(matrix)

    if matrix.shape[0] != len(row_names):
        raise ValueError(
            f"Row name mismatch for {out_path}: "
            f"matrix rows={matrix.shape[0]}, row_names={len(row_names)}"
        )

    if matrix.shape[1] != len(col_names):
        raise ValueError(
            f"Column name mismatch for {out_path}: "
            f"matrix cols={matrix.shape[1]}, col_names={len(col_names)}"
        )

    df = pd.DataFrame(
        matrix,
        index=row_names,
        columns=col_names
    )

    df.to_csv(out_path, sep="\t")


def merge_duplicate_rows(df, method="max"):
    """
    Safety step. Genes/REs should usually be unique after chromosome split,
    but if duplicated names appear, collapse them.
    """
    if not df.index.duplicated().any():
        return df

    if method == "max":
        return df.groupby(df.index).max()
    elif method == "mean":
        return df.groupby(df.index).mean()
    else:
        raise ValueError(f"Unknown duplicate merge method: {method}")


def process_one_region(region):

    print("\n============================================================")
    print(f"Processing region: {region}")
    print("============================================================")

    #run_name = f"AD_region_shareH_{region}_K{K}"
    run_name = f"AD_region_shareH_{region}_K{K}"

    run_dir = os.path.join(REGA_OUT_ROOT, run_name)

    phase2_pt = os.path.join(
        run_dir,
        "phase2",
        f"mu_blocks_outer{outer_step}.pt"
    )

    if not os.path.exists(phase2_pt):
        print(f"WARNING: missing checkpoint, skipping region {region}: {phase2_pt}")
        return None

    region_out_dir = os.path.join(OUT_DIR, region)
    os.makedirs(region_out_dir, exist_ok=True)

    # ------------------------------------------------------------
    # Load checkpoint
    # ------------------------------------------------------------

    ckpt = load_checkpoint(phase2_pt)

    W_blocks = ckpt["W_blocks"]
    B_blocks = ckpt["B_blocks"]
    V_cells = ckpt["V_cells"]
    H = ckpt["H"]

    block_cells = ckpt["block_cells"]
    block_chrs = ckpt["block_chrs"]
    samples = ckpt["samples"]
    motifs = ckpt["motifs"]

    samples = [str(x) for x in samples]
    motifs = [str(x) for x in motifs]

    n_modules = H.shape[0]
    module_names = [f"M{i + 1}" for i in range(n_modules)]

    print("Loaded checkpoint:", phase2_pt)
    print("step:", ckpt.get("step", "NA"))
    print("loss:", ckpt.get("loss", "NA"))
    print("W blocks:", len(W_blocks))
    print("B blocks:", len(B_blocks))
    print("V cells:", list(V_cells.keys()))
    print("H shape:", tuple(H.shape))
    print("samples:", len(samples))
    print("motifs:", len(motifs))

    # ------------------------------------------------------------
    # Save H for this region
    # ------------------------------------------------------------

    H_np = H.detach().cpu().numpy()

    h_out = os.path.join(
        region_out_dir,
        f"{region}_H_module_by_sample.tsv"
    )

    write_matrix_tsv(
        H_np,
        h_out,
        row_names=module_names,
        col_names=samples
    )

    print("Saved H:", h_out)

    # ------------------------------------------------------------
    # Save V for each celltype-region
    # ------------------------------------------------------------

    v_manifest = []

    for cell_key, V in V_cells.items():

        if not str(cell_key).endswith(f"_{region}"):
            raise ValueError(
                f"Region mismatch: processing {region}, but V key is {cell_key}"
            )

        V_np = V.detach().cpu().numpy()

        if V_np.shape[0] != len(motifs):
            raise ValueError(
                f"Motif mismatch for {cell_key}: "
                f"V rows={V_np.shape[0]}, motifs={len(motifs)}"
            )

        v_out = os.path.join(
            region_out_dir,
            f"{cell_key}_V_TF_by_module.tsv"
        )

        write_matrix_tsv(
            V_np,
            v_out,
            row_names=motifs,
            col_names=module_names
        )

        v_manifest.append({
            "region": region,
            "celltype_region": cell_key,
            "V_shape": str(V_np.shape),
            "V_file": v_out
        })

        print("Saved V:", v_out)

    # ------------------------------------------------------------
    # Compute block-level BV / WBV, then merge by celltype-region
    # ------------------------------------------------------------

    bv_parts = {}
    wbv_parts = {}
    w_streams = {}
    b_streams = {}
    block_manifest = []

    with torch.no_grad():

        for i, (W, B, cell_key, chr_name) in enumerate(
            zip(W_blocks, B_blocks, block_cells, block_chrs)
        ):

            print(f"Computing block {i + 1}/{len(W_blocks)}: {cell_key} {chr_name}")

            if not str(cell_key).endswith(f"_{region}"):
                raise ValueError(
                    f"Region mismatch: processing {region}, but block cell_key is {cell_key}"
                )

            V = V_cells[cell_key]

            if B.shape[1] != V.shape[0]:
                raise ValueError(
                    f"B/V mismatch for block {i}: "
                    f"B={tuple(B.shape)}, V={tuple(V.shape)}"
                )

            if W.shape[1] != B.shape[0]:
                raise ValueError(
                    f"W/B mismatch for block {i}: "
                    f"W={tuple(W.shape)}, B={tuple(B.shape)}"
                )

            block_npz_path = get_block_npz_path(cell_key, chr_name)

            if not os.path.exists(block_npz_path):
                raise FileNotFoundError(
                    f"Missing block npz for names: {block_npz_path}"
                )

            z = np.load(block_npz_path, allow_pickle=True)

            genes = to_str_array(z["genes"])
            RE = to_str_array(z["RE"])

            if len(genes) != W.shape[0]:
                raise ValueError(
                    f"Gene count mismatch for {cell_key} {chr_name}: "
                    f"genes={len(genes)}, W rows={W.shape[0]}"
                )

            if len(RE) != W.shape[1]:
                raise ValueError(
                    f"RE count mismatch for {cell_key} {chr_name}: "
                    f"RE={len(RE)}, W cols={W.shape[1]}"
                )

            # Main matrices
            BV = B @ V
            WBV = W @ BV

            BV_np = BV.detach().cpu().numpy()
            WBV_np = WBV.detach().cpu().numpy()

            if cell_key not in w_streams:
                w_streams[cell_key] = tempfile.TemporaryFile(
                    mode="w+t", dir=region_out_dir
                )
                w_streams[cell_key].write("gene\tRE\tvalue\n")
                b_streams[cell_key] = tempfile.TemporaryFile(
                    mode="w+t", dir=region_out_dir
                )
                b_streams[cell_key].write("RE\tmotif\tvalue\n")
            binary_sparse_triplets(
                W.detach().cpu().numpy(), genes, RE, "gene", "RE"
            ).to_csv(w_streams[cell_key], sep="\t", index=False, header=False)
            binary_sparse_triplets(
                B.detach().cpu().numpy(), RE, motifs, "RE", "motif"
            ).to_csv(b_streams[cell_key], sep="\t", index=False, header=False)

            bv_df = pd.DataFrame(
                BV_np,
                index=RE,
                columns=module_names
            )

            wbv_df = pd.DataFrame(
                WBV_np,
                index=genes,
                columns=module_names
            )

            if cell_key not in bv_parts:
                bv_parts[cell_key] = []
            if cell_key not in wbv_parts:
                wbv_parts[cell_key] = []

            bv_parts[cell_key].append(bv_df)
            wbv_parts[cell_key].append(wbv_df)

            block_manifest.append({
                "region": region,
                "block_index": i,
                "celltype_region": cell_key,
                "chromosome": chr_name,
                "W_shape": str(tuple(W.shape)),
                "B_shape": str(tuple(B.shape)),
                "V_shape": str(tuple(V.shape)),
                "BV_shape": str(BV_np.shape),
                "WBV_shape": str(WBV_np.shape),
                "n_genes": len(genes),
                "n_RE": len(RE),
                "source_block_npz": block_npz_path
            })

    # ------------------------------------------------------------
    # Merge and save one BV/WBV matrix per celltype-region
    # ------------------------------------------------------------

    merged_manifest = []

    for cell_key in sorted(wbv_parts.keys()):

        print(f"Merging matrices for {cell_key}")

        wbv_merged = pd.concat(wbv_parts[cell_key], axis=0)
        bv_merged = pd.concat(bv_parts[cell_key], axis=0)

        wbv_merged = merge_duplicate_rows(wbv_merged, method="max")
        bv_merged = merge_duplicate_rows(bv_merged, method="max")

        wbv_out = os.path.join(
            region_out_dir,
            f"{cell_key}_WBV_genes_by_module.tsv"
        )

        bv_out = os.path.join(
            region_out_dir,
            f"{cell_key}_BV_RE_by_module.tsv"
        )

        wbv_merged.to_csv(wbv_out, sep="\t")
        bv_merged.to_csv(bv_out, sep="\t")

        merged_manifest.append({
            "region": region,
            "celltype_region": cell_key,
            "WBV_shape": str(wbv_merged.shape),
            "BV_shape": str(bv_merged.shape),
            "WBV_file": wbv_out,
            "BV_file": bv_out
        })

        print("Saved WBV:", wbv_out)
        print("Saved BV :", bv_out)

    for cell_key in sorted(w_streams):
        prefix = safe_name(cell_key)
        try:
            write_merged_sparse(
                w_streams[cell_key],
                os.path.join(region_out_dir, f"{prefix}_W_binary_sparse.tsv"),
                "gene", "RE"
            )
            write_merged_sparse(
                b_streams[cell_key],
                os.path.join(region_out_dir, f"{prefix}_B_binary_sparse.tsv"),
                "RE", "motif"
            )
        finally:
            w_streams[cell_key].close()
            b_streams[cell_key].close()

    # ------------------------------------------------------------
    # Save manifests
    # ------------------------------------------------------------

    block_manifest_df = pd.DataFrame(block_manifest)
    block_manifest_out = os.path.join(
        region_out_dir,
        f"{region}_block_extraction_manifest.tsv"
    )
    block_manifest_df.to_csv(block_manifest_out, sep="\t", index=False)

    merged_manifest_df = pd.DataFrame(merged_manifest)
    merged_manifest_out = os.path.join(
        region_out_dir,
        f"{region}_merged_matrix_manifest.tsv"
    )
    merged_manifest_df.to_csv(merged_manifest_out, sep="\t", index=False)

    v_manifest_df = pd.DataFrame(v_manifest)
    v_manifest_out = os.path.join(
        region_out_dir,
        f"{region}_V_manifest.tsv"
    )
    v_manifest_df.to_csv(v_manifest_out, sep="\t", index=False)

    print("\nDone region:", region)
    print("Region output:", region_out_dir)
    print("Block manifest:", block_manifest_out)
    print("Merged matrix manifest:", merged_manifest_out)
    print("V manifest:", v_manifest_out)

    return {
        "region": region,
        "checkpoint": phase2_pt,
        "out_dir": region_out_dir,
        "n_blocks": len(block_manifest),
        "n_celltype_regions": len(merged_manifest),
        "n_modules": n_modules,
        "n_samples": len(samples),
        "n_motifs": len(motifs)
    }



def binary_sparse_triplets(matrix, row_names, col_names, row_label, col_label):
    """Return nonzero entries as row name, column name, and binary value."""
    if matrix.ndim != 2:
        raise ValueError(f"Expected a 2D matrix, got shape {matrix.shape}")
    if matrix.shape[0] != len(row_names) or matrix.shape[1] != len(col_names):
        raise ValueError(
            f"Name/shape mismatch: matrix={matrix.shape}, "
            f"rows={len(row_names)}, columns={len(col_names)}"
        )

    row_idx, col_idx = np.nonzero(matrix > 0)
    return pd.DataFrame({
        row_label: np.asarray(row_names)[row_idx],
        col_label: np.asarray(col_names)[col_idx],
        "value": np.ones(len(row_idx), dtype=np.uint8),
    })


def write_merged_sparse(stream, out_path, row_label, col_label):
    stream.seek(0)
    df = pd.read_csv(
        stream, sep="\t",
        dtype={row_label: str, col_label: str, "value": np.uint8},
    )
    # Collapse repeated coordinates if the same feature occurs in multiple blocks.
    df = df.drop_duplicates([row_label, col_label], keep="first")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved {out_path} ({len(df)} nonzero entries)")


# ============================================================
# Main: write H/V/W/B/WBV/BV and binary sparse W/B for one region
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Extract shared REGA matrices")
    parser.add_argument("--region", choices=regions, help="Process just one region")
    args = parser.parse_args()
    selected_regions = [args.region] if args.region else regions
    summary_rows = []
    for region in selected_regions:
        result = process_one_region(region)
        if result is None:
            raise FileNotFoundError(f"Missing REGA checkpoint for {region}")
        summary_rows.append(result)
    summary_name = f"{args.region}_extracted_matrix_summary.tsv" if args.region else "region_extracted_matrix_summary.tsv"
    summary_out = os.path.join(OUT_DIR, summary_name)
    pd.DataFrame(summary_rows).to_csv(summary_out, sep="\t", index=False)
    print("Saved summary:", summary_out)


if __name__ == "__main__":
    main()
