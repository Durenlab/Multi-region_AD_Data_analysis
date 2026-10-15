import os
import glob
import re
import numpy as np
from collections import defaultdict

blocks_root = "05REGA_blocks_region_consistent"

npz_files = sorted(
    path for path in glob.glob(os.path.join(blocks_root, "*_blocks_aligned", "chr*.npz"))
    if not os.path.basename(os.path.dirname(path)).startswith("Vasc_")
)

print("Total npz blocks:", len(npz_files))

if len(npz_files) == 0:
    raise SystemExit("No npz files found.")

# ------------------------------------------------------------
# Helper: extract key and region
# Example:
#   05REGA_blocks_region_consistent/Mic_Immune_AG_blocks_aligned/chr1.npz
# key    = Mic_Immune_AG
# region = AG
# ------------------------------------------------------------

def parse_key_region(npz_path):
    block_dir = os.path.basename(os.path.dirname(npz_path))
    key = re.sub(r"_blocks_aligned$", "", block_dir)

    # region is the last underscore-separated field
    parts = key.split("_")
    region = parts[-1]

    return key, region

# ------------------------------------------------------------
# References
# ------------------------------------------------------------

region_sample_ref = {}
motif_ref = None

bad_samples = []
bad_motifs = []

region_counts = defaultdict(int)
region_sample_counts = {}

# ------------------------------------------------------------
# Check every block
# ------------------------------------------------------------

for f in npz_files:
    z = np.load(f, allow_pickle=True)

    samples = [
        x.decode() if isinstance(x, bytes) else str(x)
        for x in z["samples"]
    ]

    motifs = [
        x.decode() if isinstance(x, bytes) else str(x)
        for x in z["motifs"]
    ]

    key, region = parse_key_region(f)

    region_counts[region] += 1

    # Region-wise sample check
    if region not in region_sample_ref:
        region_sample_ref[region] = samples
        region_sample_counts[region] = len(samples)
    else:
        if samples != region_sample_ref[region]:
            bad_samples.append({
                "file": f,
                "key": key,
                "region": region,
                "sample_count": len(samples),
                "expected_sample_count": len(region_sample_ref[region])
            })

    # Global motif check
    if motif_ref is None:
        motif_ref = motifs
    else:
        if motifs != motif_ref:
            bad_motifs.append({
                "file": f,
                "key": key,
                "region": region,
                "motif_count": len(motifs),
                "expected_motif_count": len(motif_ref)
            })

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

print("\nRegion sample references:")
for region in sorted(region_sample_ref.keys()):
    print(
        f"{region}: "
        f"samples={len(region_sample_ref[region])}, "
        f"blocks={region_counts[region]}"
    )

print("\nMotif count:", len(motif_ref))
print("Blocks with bad region-wise samples:", len(bad_samples))
print("Blocks with bad motifs:", len(bad_motifs))

if bad_samples:
    print("\nExample bad sample blocks:")
    for x in bad_samples[:10]:
        print(
            x["file"],
            "| region:", x["region"],
            "| sample_count:", x["sample_count"],
            "| expected:", x["expected_sample_count"]
        )

if bad_motifs:
    print("\nExample bad motif blocks:")
    for x in bad_motifs[:10]:
        print(
            x["file"],
            "| motif_count:", x["motif_count"],
            "| expected:", x["expected_motif_count"]
        )

# ------------------------------------------------------------
# Final decision
# ------------------------------------------------------------

if len(bad_samples) == 0 and len(bad_motifs) == 0:
    print("\nPASS: samples are identical within each region, and motifs are identical across all blocks.")
else:
    raise SystemExit("\nFAIL: region-wise sample mismatch or motif mismatch.")