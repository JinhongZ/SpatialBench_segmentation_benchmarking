# SpatialBench segmentation benchmarking

This repository contains an R workflow for benchmarking cell segmentation methods in spatial transcriptomics data. It processes MERSCOPE and Xenium samples, builds Seurat objects from different segmentation outputs, performs sample processing and clustering, annotates cell types, computes benchmarking metrics, and generates thesis figures.

## Folder structure

```text
analysis_pipeline
  01_preprocessing
  02_data_processing
  03_downstream_analysis
  04_benchmarking
  helper_functions

thesis_figures
  plot_data
  figures
```

## Workflow

1. Create cell segmented Seurat objects from MERSCOPE, Xenium, or Proseg outputs.
2. Process each sample by adding metadata, filtering cells, and running SCTransform.
3. Merge processed samples for each platform and segmentation method.
4. Run PCA, Harmony correction, clustering, and UMAP.
5. Annotate cell types with ScType marker based annotation.
6. Run pseudo bulk differential expression analysis.
7. Compute segmentation benchmarking metrics.
8. Generate thesis figures from processed summary tables.

## Main scripts

```text
analysis_pipeline/01_preprocessing/create_seurat_segmented_merscope.R
analysis_pipeline/01_preprocessing/create_seurat_segmented_xenium.R
analysis_pipeline/01_preprocessing/create_seurat_segmented_proseg.R
```

These scripts create Seurat objects from platform specific or Proseg segmentation outputs.

```text
analysis_pipeline/02_data_processing/01_sample_processing.R
analysis_pipeline/02_data_processing/02_merge_samples.R
analysis_pipeline/02_data_processing/03_harmony_clustering.R
```

These scripts process samples, merge samples, and run Harmony based clustering.

```text
analysis_pipeline/03_downstream_analysis/01_annotation.R
analysis_pipeline/03_downstream_analysis/02_pseudo_bulk_DE.R
```

These scripts annotate cell types and run pseudo bulk differential expression analysis.

```text
analysis_pipeline/04_benchmarking/create_assignment_info.R
analysis_pipeline/04_benchmarking/01_add_morphological_metrics.R
analysis_pipeline/04_benchmarking/02_create_metadata_df.R
analysis_pipeline/04_benchmarking/run_purity_analysis.R
```

These scripts compute transcript assignment metrics, morphology metrics, sample and cell level metadata, and purity metrics.

## Example commands

Run scripts from the project root so helper functions are found correctly.

```bash
Rscript analysis_pipeline/01_preprocessing/create_seurat_segmented_xenium.R \
  --data_dir path/to/xenium_output \
  --sample_name sample_name \
  --method default \
  --out_dir path/to/preprocessed
```

```bash
Rscript analysis_pipeline/02_data_processing/01_sample_processing.R \
  --input_rds path/to/input.rds \
  --sample_name sample_name \
  --platform xenium \
  --seg default \
  --assay Xenium \
  --qc_min_counts 10 \
  --out_dir path/to/processed
```

```bash
Rscript analysis_pipeline/02_data_processing/02_merge_samples.R \
  --input_dir path/to/processed \
  --method xenium_default \
  --out_file path/to/merged.rds
```

```bash
Rscript analysis_pipeline/02_data_processing/03_harmony_clustering.R \
  --input_rds path/to/merged.rds \
  --method xenium_default \
  --resolution 0.6 \
  --seed 123 \
  --out_file path/to/clustered.rds
```

```bash
Rscript analysis_pipeline/03_downstream_analysis/01_annotation.R \
  --input_rds path/to/clustered.rds \
  --method xenium_default \
  --res 0.6 \
  --out_file path/to/annotated.rds
```

## Inputs

MERSCOPE inputs are expected to include cell by gene counts, detected transcripts, cell metadata, and cell boundary files.

Xenium inputs are expected to include the cell feature matrix, cell metadata, cell boundaries, and transcript files.

Proseg inputs are expected to include expected counts, cell metadata, transcript metadata, and cell polygons.

## Outputs

The pipeline produces Seurat objects, annotated Seurat objects, benchmarking summary tables, pseudo bulk differential expression results, MECR purity outputs, and figure PDFs.

Generated figure files are stored in `thesis_figures/figures`.

Processed plotting data are stored in `thesis_figures/plot_data`.

## R packages

The scripts use packages including `Seurat`, `SeuratWrappers`, `harmony`, `optparse`, `arrow`, `Matrix`, `future`, `tidyverse`, `data.table`, `purrr`, `tidyr`, `dplyr`, `limma`, `edgeR`, `sf`, `sp`, and `SingleCellExperiment`.

## Known consistency notes

1. `analysis_pipeline/02_data_processing/01_sample_processing.R` documents `cellpose2`, but the validation currently accepts `cellpose22`.

2. `analysis_pipeline/04_benchmarking/create_assignment_info.R`, `analysis_pipeline/04_benchmarking/01_add_morphological_metrics.R`, `analysis_pipeline/04_benchmarking/02_create_metadata_df.R`, and `analysis_pipeline/04_benchmarking/run_purity_analysis.R` use `make_option` and `parse_args`, but do not load `optparse`.

3. `analysis_pipeline/04_benchmarking/run_purity_analysis.R` reads `--data_info_dir` into `sample_info`, while the rest of the script uses `data_info`.

4. `analysis_pipeline/04_benchmarking/run_purity_analysis.R` saves files under `out_path`, but the final message refers to `out_file`.

5. `analysis_pipeline/helper_functions/create_sample_info.R` appears to construct `file_path` with `mutate(file_path = segmentation_dir, sub_dir)`. This likely should combine `segmentation_dir` and `sub_dir`.

6. `analysis_pipeline/04_benchmarking/02_create_metadata_df.R` labels the MERSCOPE default segmentation as `Cellpose1`, while other scripts and plotting functions use `Default`.

7. `analysis_pipeline/03_downstream_analysis/02_pseudo_bulk_DE.R` accepts `sample_id_col`, but the helper function currently groups by the fixed column name `sample_id`.

8. `analysis_pipeline/03_downstream_analysis/01_annotation.R` loads the input RDS twice. This is not a naming inconsistency, but it is a small redundant step.

## Notes

This workflow is adapted from SpatialBench and customized for segmentation benchmarking across MERSCOPE and Xenium datasets.
