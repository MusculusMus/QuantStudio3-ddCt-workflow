# QuantStudio 3  ΔΔCt workfloW

This lightweight R workflow imports one or more tab-delimited QuantStudio™ 3 exports, converts them into the format expected by the Bioconductor [`ddCt`](https://bioconductor.org/packages/ddCt/) package, and calculates relative expression with the ΔΔCt method.


## What it does

- Finds date-prefixed `.txt` QuantStudio 3 exports in a folder, avoiding unrelated text files.
- Reads the table after its `[Results]` line.
- Standardizes `Sample Name`, `Target Name`, and `CT` to the column names required by `ddCt`.
- Excludes wells marked `Omit = true` by default.
- Replaces `Undetermined` Ct values with a configurable Ct value (35 by default).
- Lets the user choose exactly one housekeeping gene and one or more calibration samples from menus.
- Returns both the cleaned long Ct data and an analysis-ready, wide relative-expression table.

## Install dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ddCt")
```

## Interactive use

```r
source("R/quantstudio3_ddct_workflow.R")

result <- run_quantstudio3_ddct(
  folder_path = "/path/to/your/QuantStudio3_txt_exports",
  write_expression_csv = TRUE
)

qPCR_df <- result$expression
```

By default, imported filenames must begin with a date such as `2026-07-08_run.txt`, `20260708_run.txt`, or `2026_07_08_run.txt`. To import every `.txt` file in a dedicated export folder, use:

```r
result <- run_quantstudio3_ddct(
  folder_path = "/path/to/your/QuantStudio3_txt_exports",
  file_pattern = "\\.txt$"
)
```

For files beginning with a particular year only, for example 2026, use `file_pattern = "^2026.*\\.txt$"`. Note that `^2026*` is not equivalent: the `*` applies only to the preceding `6`.

The menus first show target names for the housekeeping-gene choice. Then they show sample names repeatedly for calibration selection; choose **Finish selection** when all calibration samples have been selected.

## Reproducible, non-interactive use

For scripts and reproducible analysis, record selections directly instead of using menus:

```r
result <- run_quantstudio3_ddct(
  folder_path = "/path/to/your/QuantStudio3_txt_exports",
  housekeeping_gene = "Cotl1",
  calibration_samples = c("S7", "S28", "S17", "S20"),
  undetermined_ct = 35,
  write_expression_csv = TRUE
)

qPCR_df <- result$expression
```

`qPCR_df` has one row per sample, a `Sample` column, and one expression column per target gene.

## Main outputs

| Object | Contents |
| --- | --- |
| `result$expression` | Wide ΔΔCt relative-expression table, ready to merge with sample metadata. |
| `result$ct_table` | Cleaned long Ct table passed to `ddCt`. |
| `result$import_summary` | Imported files and usable row counts. |
| `result$ddct_result` | The original `ddCtExpression` object. |
| `result$warning_file` | `ddCt` warnings written during calculation. |

## Add experimental metadata before plotting or statistics

QuantStudio exports generally contain sample IDs but not study groups or ages. Join a sample-metadata file before plotting or modelling:

```r
metadata <- read.csv("sample_metadata.csv") # Must contain Sample, Group, Age, etc.
analysis_df <- merge(qPCR_df, metadata, by = "Sample", all.x = TRUE)
```

For example, after sourcing your plotting/statistics helper script (the one that defines `plot_qpcr_interactive()` and `run_qpcr_workflow()`), you can use `IL6` as the response:

```r
plot_qpcr_interactive(
  analysis_df,
  response = "IL6",
  x_factor = "Group",
  facet_factor = "Age"
)

run_qpcr_workflow(
  analysis_df,
  response = "IL6",
  factors = c("Group", "Age")
)
```

## Important analytical choices

- **Undetermined values:** This workflow maps `Undetermined` to Ct 35, matching the original analysis. This is not universally appropriate. Set `undetermined_ct` to the validated detection-limit Ct for your assay, or use a censored-data strategy when scientifically warranted.
- **Omitted wells:** Wells explicitly marked `Omit = true` are excluded by default. Set `exclude_omitted = FALSE` to retain them.
- **Reference gene:** Select a housekeeping gene only after confirming that its expression is stable under the experimental conditions.
- **Calibration samples:** The selected calibration samples form the reference baseline used by `ddCt`.

## Supported QuantStudio 3 export columns

The results table must contain equivalents of:

```text
Sample Name    Target Name    CT
```

The importer also accepts R-sanitized versions (`Sample.Name`, `Target.Name`) and `Ct`. Other QuantStudio columns are ignored.
