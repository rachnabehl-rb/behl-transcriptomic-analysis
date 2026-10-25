# Analysis code accompanying Behl et al.

This archive contains the supplied R scripts for the Arabidopsis thaliana WT and gsnor-ko multigenerational study. Scripts are grouped by analysis. Their contents, analytical steps and original input/output paths have been preserved; numbered upload suffixes have been removed from the added script filenames.

## Folder contents and execution order

1. **01_differential_expression/**: six genotype-by-generation DESeq2 scripts for G1, G3 and G5. These generate full differential-expression results and up/down gene lists for downstream analyses.
2. **02_quality_control_and_ordination/**: normalised library distributions, dendrograms, PCA, DAPC and clock-gene analysis.
3. **03_deg_counts_overlap_and_trajectories/**: DEG counts, Venn/UpSet plots, overlap statistics, volcano plots, alluvial states, heatmaps and transcription-factor analysis. Run `Venn Fig3C.R` and then `Fig3C_overlap_statistics.R` in the same R session. The latter requires objects created by the former and generates the three-generation overlap statistics associated with Supplementary File S23. Random DEG sets are sampled independently from each generation's eligible-gene pool.
4. **04_functional_enrichment/**: contrast-specific Mercator enrichment and Figure 4 plotting. Run the updated `Mercator_enrichment_S17_contrast_background.R` with full DESeq2 results and Mercator annotations, then `Go term_logic revived_Fig 4.R` with the generated enrichment results. `Mercator Go analysis.R` analyses WGCNA modules and requires module assignments from folder 06.
5. **05_candidate_genes/**: candidate-gene expression heatmap for Figure 6. The literature-informed regulatory models are not statistically inferred by this plotting script.
6. **06_wgcna/**: preprocessing, network analysis and module ranking. `WGCNA_final_ME40,0.4.R` generates the module output used by downstream scripts. Before running, confirm the intended VST input: `Preprocessing WGCNA.R` writes `vsd_batch.csv`, whereas the network script reads `vsd.csv`.
7. **07_pathogen_assay/**: supplementary CFU analysis using the spreadsheet specified in the script.
8. **08_genotype_treatment_interactions/**: run `DESeq2_genotype_treatment_interactions_G1_G3_G5.R` first, using integer gene counts and sample metadata. It generates generation-specific genotype-by-treatment tests and per-gene WT, mutant and interaction effects. Run `Make_supplementary_interaction_table.R` on `interaction_summary_all_generations.csv` to create the table associated with Supplementary File S25. `Supplementary_interaction_figures_and_normalized_counts.R` is an optional companion for additional plots and example-gene normalised counts; its example genes are user-editable and are not a required part of the manuscript. It expects per-gene interaction result CSVs in `INPUT_DIR`; point this to the generated results directory and separately set the count/metadata paths as needed. Its local figure labels do not define the final manuscript numbering.
9. **09_fold_change_sensitivity/**: `DEG_counts_at_LFC_0_0.2_0.5.R` refits within-genotype contrasts and exports DEG counts at absolute estimated log2 fold-change cutoffs of 0, 0.2 and 0.5 (Supplementary File S26). `Trajectory_LFC_sensitivity_S18_S22.R` reads full within-genotype DESeq2 result files and evaluates the original heatmap persistent/switcher rules at these cutoffs (Supplementary File S28). These are post hoc filters on estimated effects after adjusted P < 0.05, not tests against a nonzero effect-size threshold and not shrunken-fold-change analyses. The trajectory rules include documented two-point exceptions and do not establish uninterrupted expression in unsampled generations.
10. **10_optional_templates/**: `pedigree_sensitivity_template.R` is an optional template requiring pedigree/parent assignments verified against propagation records. It is not evidence that a pedigree sensitivity analysis was completed or reported. Sample-name prefixes alone are insufficient to establish parentage.

## Inputs and configuration

- Extract this archive before use. Adapt local working directories, input paths and output paths in each script to your computer.
- Supply the relevant count matrix, sample metadata, full DESeq2 results, Mercator annotations, phenotype tables, WGCNA inputs and CFU data, as required by the script. Input data are not bundled with this code archive; consult the associated Zenodo data files and supplementary information.
- Count files are referenced by different names in the supplied scripts (`Count_matrix.txt` or `genes_gtf_merged_gene_counts.txt`). Match each path to the deposited gene-count table and confirm its columns and sample identifiers before running.
- The added interaction and DEG-count scripts expect `Metadata_Sample_Area_medians_WITH_GROUPS.csv` with `SampleName`, `Genotype`, `Treatment` and `Generation` columns. Reconcile the deposited metadata filename and column names with these requirements.
- Some scripts consume outputs created by earlier scripts rather than raw counts. Read their configuration and input checks before execution.
- Install the R/Bioconductor packages required by each script. The interaction and supplementary interaction-plot scripts write session-information files when run; no verified record of the original computing environment is bundled here.
- At G5, treatment labels denote treatment history, because G4 and G5 plants were grown under control conditions.

## Scope and verification

The archive contains 37 supplied R scripts, including one optional pedigree template. Script bytes have been preserved and package integrity checked; analyses have not been rerun as part of repackaging. Main phenotyping-analysis scripts and the remeasurement analysis underlying Supplementary File S27 are not included in the supplied files. The optional pedigree template does not replace those analyses. Check the final manuscript figure numbering and the provenance of WGCNA/module-enrichment inputs before reuse.

`CHANGES.md` records additions and replacements in this update. `SHA256SUMS.csv` contains SHA-256 checksums for all packaged files except the checksum table itself.
