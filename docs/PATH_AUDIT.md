# Path audit

The published Python, R, shell, and notebook code no longer embeds the original Palmetto project path. Scripts read model results through `REGA_SEPARATE_ROOT` and `REGA_SHARED_ROOT`; manuscript scripts use `AD_MANUSCRIPT_WORK_ROOT` for external manuscript inputs and tables. The focused SPI1 target-expression scripts take an explicit output directory for each fitted model. The broader shared-model TF–target-gene Slurm wrappers locate their code relative to the wrapper and write to `TF_TARGET_RESULTS_ROOT` or an ignored `results/` directory beside the scripts. Slurm logs go to the directory from which `sbatch` is called unless an output path is supplied at submission.

Set the roots and reference resources described in [EXTERNAL_INPUTS.md](EXTERNAL_INPUTS.md) before running a stage. `SINGLE_GRN_ROOT` must identify the independent-model GRNs used for the reported target-gene analysis; it is deliberately not inferred from an archived folder. MAGMA requires a new `MAGMA_OUTPUT_ROOT` and refuses to reuse an existing directory.

The original numbered paths in `SOURCE_MANIFEST.tsv` are source provenance, not runtime paths. This is a static path audit; the full analysis has not been rerun from the code release.
