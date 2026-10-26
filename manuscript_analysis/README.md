# Analyses used in the manuscript

These scripts assemble the reported tables and figures from the completed independent and shared REGA analyses. The source data and result matrices remain outside this code repository; see [external inputs](../docs/EXTERNAL_INPUTS.md). Set `AD_MANUSCRIPT_WORK_ROOT` to the `Alzheimer_project_explore` directory and `REGA_DATA_ROOT` to its parent when running scripts outside those directories.

| Folder | Manuscript result | Input stage |
| --- | --- | --- |
| `module_tables` | Candidate-module pathology, GO/KEGG, and LDSC supplementary tables (Tables 1–3) | Selected modules and downstream REGA results |
| `tf_comparison` | Significant-TF counts, overall classifier performance, the OPC_AG ROC panel, and Oli_EC TF volcano figures (Figure 4 and supplementary panels) | Separate-model TF pathology and classification results |
| `driver_tfs` | MAGMA TF table, driver-TF Supplementary Table 5, and TF-module figures | MAGMA, TF pathology, and module assignments |
| `spi1_target_expression` | Oli_EC SPI1 target expression (Figure 4i and Supplementary Figure 6d) | Separate and shared GRNs, expression, assignments, and pathology |
| `spp1_grn_overlap` | TF-to-SPP1 GRN edges and driver overlap (Supplementary Table 7) | Separate and shared GRNs and driver TFs |

The three `module_tables` scripts are independent: each reads its own completed source analysis and produces one table family. In `tf_comparison`, `significant_tf_count_comparison.ipynb` plots TF-count results, `overall_classifier_performance.ipynb` plots the all-context classifier, and `opc_ag_roc_comparison.ipynb` makes the selected OPC_AG ROC panel. The context-wise classification script in `separate/07_tf/scripts` already makes the full PR, ROC, and paired-summary figures; the OPC_AG notebook contains only the selected panel. `tf_activity_expression_volcano.R` makes the Oli_EC TF comparison panels.

In `driver_tfs`, run `01_driver_tf_pathology_tables.py` before `02_annotate_driver_tf_tables.py`; the second script reads the first script's Supplementary Table 5 files. `magma_tf_association_table.py` produces Supplementary Table 4, while `driver_tf_module_tables_and_figures.py` creates the reported driver-TF/module tables and figures from the existing TF results. Both annotation and module figures require the MAGMA TF list specified by `AD_SNP_TF_LIST` or their documented default path.

In `spi1_target_expression`, run step 01 for both model types, then step 02 to plot SPI1 from the donor scores. Targets are restricted to SPI1's gene-assignment module separately in each fitted model. [Detailed commands and the saved-target audit](spi1_target_expression/README.md) are in that folder. The single `spp1_grn_overlap` script produces the AG and EC results for Supplementary Table 7.

Numbered names mark a required input-to-output dependency. Folder order follows the Results narrative; independent analyses within a folder can run in any order after their source results exist.
