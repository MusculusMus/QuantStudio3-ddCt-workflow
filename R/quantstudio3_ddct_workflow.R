# QuantStudio 3 .txt -> ddCt workflow
#
# This script imports one or more QuantStudio 3 tab-delimited exports, finds
# the [Results] table automatically, creates a ddCt InputFrame, and interactively
# selects one housekeeping gene plus one or more calibration samples.
#
# Main interactive use:
#   source("R/quantstudio3_ddct_workflow.R")
#   ddct_result <- run_quantstudio3_ddct("path/to/QuantStudio_exports")
#   qPCR_df <- ddct_result$expression

require_ddct <- function() {
  if (!requireNamespace("ddCt", quietly = TRUE)) {
    stop(
      paste(
        "The Bioconductor package 'ddCt' is required.",
        "Install it with:",
        "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')",
        "BiocManager::install('ddCt')",
        sep = "\n"
      ),
      call. = FALSE
    )
  }
}

locate_quantstudio_results <- function(file) {
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  result_row <- which(tolower(trimws(lines)) == "[results]")
  if (!length(result_row)) {
    stop(
      sprintf("'%s' does not contain a [Results] section.", basename(file)),
      call. = FALSE
    )
  }
  result_row[[1L]]
}

normalise_quantstudio_name <- function(x) {
  x <- trimws(gsub("\ufeff", "", x, fixed = TRUE))
  tolower(gsub("[[:space:]_.]+", "", x))
}

rename_quantstudio_columns <- function(df) {
  compact_names <- normalise_quantstudio_name(names(df))
  aliases <- list(
    Sample = c("sample", "samplename"),
    Detector = c("detector", "target", "targetname"),
    Ct = c("ct", "cт") # The second value contains a Cyrillic small te.
  )
  new_names <- names(df)
  for (destination in names(aliases)) {
    matches <- which(compact_names %in% aliases[[destination]])
    if (length(matches) > 1L) {
      stop(
        sprintf("More than one column in the export matches '%s'.", destination),
        call. = FALSE
      )
    }
    if (length(matches) == 1L) new_names[matches] <- destination
  }
  names(df) <- new_names

  required <- c("Sample", "Detector", "Ct")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop(
      sprintf("Missing required QuantStudio column(s): %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  df
}

clean_quantstudio_table <- function(df, plate_name, undetermined_ct = 35,
                                    exclude_omitted = TRUE) {
  if (!is.numeric(undetermined_ct) || length(undetermined_ct) != 1L ||
      is.na(undetermined_ct) || undetermined_ct <= 0) {
    stop("'undetermined_ct' must be one positive number.", call. = FALSE)
  }
  df <- rename_quantstudio_columns(df)

  if (isTRUE(exclude_omitted)) {
    omit_column <- which(normalise_quantstudio_name(names(df)) == "omit")
    if (length(omit_column) == 1L) {
      omitted <- tolower(trimws(as.character(df[[omit_column]]))) %in% c("true", "yes", "y", "1")
      df <- df[is.na(omitted) | !omitted, , drop = FALSE]
    }
  }

  ct_character <- trimws(as.character(df$Ct))
  undetermined <- !is.na(ct_character) & tolower(ct_character) %in% c("undetermined", "undet")
  ct <- suppressWarnings(as.numeric(ct_character))
  ct[undetermined] <- undetermined_ct

  invalid_ct <- is.na(ct) & !undetermined
  if (any(invalid_ct)) {
    bad_values <- unique(ct_character[invalid_ct])
    bad_values <- bad_values[!is.na(bad_values)]
    stop(
      sprintf(
        "'%s' has missing or non-numeric Ct value(s): %s.",
        plate_name,
        paste(utils::head(bad_values, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  clean <- data.frame(
    Ct = ct,
    Sample = trimws(as.character(df$Sample)),
    Detector = trimws(as.character(df$Detector)),
    Platename = plate_name,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  keep <- !is.na(clean$Sample) & nzchar(clean$Sample) &
    !is.na(clean$Detector) & nzchar(clean$Detector)
  clean <- clean[keep, , drop = FALSE]
  if (!nrow(clean)) {
    stop(sprintf("'%s' has no usable Sample/Target rows.", plate_name), call. = FALSE)
  }
  clean
}

read_quantstudio3_txt <- function(file, undetermined_ct = 35, exclude_omitted = TRUE) {
  if (!file.exists(file)) stop(sprintf("File not found: %s", file), call. = FALSE)
  results_row <- locate_quantstudio_results(file)
  raw <- utils::read.delim(
    file,
    skip = results_row,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE
  )
  clean_quantstudio_table(
    raw,
    plate_name = basename(file),
    undetermined_ct = undetermined_ct,
    exclude_omitted = exclude_omitted
  )
}

# Default pattern accepts date-prefixed names such as 2026-07-08_run.txt,
# 20260708_run.txt, or 2026_07_08_run.txt. Supply file_pattern = "\\.txt$"
# explicitly if every text file in the folder is a QuantStudio export.
find_quantstudio3_files <- function(folder_path,
                                    file_pattern = "^[0-9]{4}[-_]?[0-9]{2}[-_]?[0-9]{2}.*\\.txt$") {
  if (!dir.exists(folder_path)) {
    stop(sprintf("Folder not found: %s", folder_path), call. = FALSE)
  }
  files <- list.files(
    path = folder_path,
    pattern = file_pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(files)) {
    stop(
      sprintf("No files matching '%s' were found in '%s'.", file_pattern, folder_path),
      call. = FALSE
    )
  }
  files
}

read_quantstudio3_folder <- function(folder_path,
                                      file_pattern = "^[0-9]{4}[-_]?[0-9]{2}[-_]?[0-9]{2}.*\\.txt$",
                                      undetermined_ct = 35, exclude_omitted = TRUE) {
  files <- find_quantstudio3_files(folder_path, file_pattern)
  tables <- lapply(
    files,
    read_quantstudio3_txt,
    undetermined_ct = undetermined_ct,
    exclude_omitted = exclude_omitted
  )
  combined <- do.call(rbind, tables)
  rownames(combined) <- NULL
  list(
    ct_table = combined,
    files = files,
    import_summary = data.frame(
      file = basename(files),
      usable_rows = vapply(tables, nrow, integer(1)),
      stringsAsFactors = FALSE
    )
  )
}

choose_one_value <- function(values, title) {
  values <- sort(unique(as.character(values)))
  if (!length(values)) stop("There are no values to select.", call. = FALSE)
  if (!interactive()) {
    stop("Interactive selection requires an interactive R session.", call. = FALSE)
  }
  choice <- utils::menu(values, title = title)
  if (choice == 0L) stop("Selection cancelled.", call. = FALSE)
  values[[choice]]
}

choose_multiple_values <- function(values, title) {
  values <- sort(unique(as.character(values)))
  if (!length(values)) stop("There are no values to select.", call. = FALSE)
  if (!interactive()) {
    stop("Interactive selection requires an interactive R session.", call. = FALSE)
  }

  selected <- character()
  repeat {
    remaining <- setdiff(values, selected)
    if (!length(remaining)) break
    options <- c("Finish selection", remaining)
    choice <- utils::menu(
      options,
      title = paste0(title, "\nSelected: ", if (length(selected)) paste(selected, collapse = ", ") else "none")
    )
    if (choice == 0L || choice == 1L) break
    selected <- c(selected, options[[choice]])
  }
  if (!length(selected)) {
    stop("Select at least one calibration sample.", call. = FALSE)
  }
  selected
}

validate_ddct_choices <- function(ct_table, housekeeping_gene, calibration_samples) {
  detectors <- sort(unique(ct_table$Detector))
  samples <- sort(unique(ct_table$Sample))
  if (!is.character(housekeeping_gene) || length(housekeeping_gene) != 1L ||
      !housekeeping_gene %in% detectors) {
    stop("'housekeeping_gene' must be one detector name found in the imported data.", call. = FALSE)
  }
  if (!is.character(calibration_samples) || !length(calibration_samples) ||
      anyDuplicated(calibration_samples) || !all(calibration_samples %in% samples)) {
    stop("'calibration_samples' must be one or more unique sample names found in the imported data.", call. = FALSE)
  }

  samples_without_hk <- setdiff(samples, unique(ct_table$Sample[ct_table$Detector == housekeeping_gene]))
  if (length(samples_without_hk)) {
    stop(
      sprintf(
        "The housekeeping gene '%s' is missing for: %s.",
        housekeeping_gene,
        paste(samples_without_hk, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ddct_expression_table <- function(ct_table, housekeeping_gene, calibration_samples,
                                  warning_file = "ddct_warnings.txt") {
  require_ddct()
  validate_ddct_choices(ct_table, housekeeping_gene, calibration_samples)
  input <- ddCt::InputFrame(ct_table)
  expression_result <- ddCt::ddCtExpression(
    input,
    warningStream = warning_file,
    algorithm = "ddCt",
    type = "mean",
    efficiencies = NULL,
    calibrationSample = calibration_samples,
    housekeepingGenes = housekeeping_gene
  )
  expression <- as.data.frame(t(expression_result@assayData[["exprs"]]), check.names = FALSE)
  expression$Sample <- rownames(expression)
  rownames(expression) <- NULL
  list(
    input_frame = input,
    ddct_result = expression_result,
    expression = expression
  )
}

print.quantstudio3_ddct <- function(x, ...) {
  cat("\nQuantStudio 3 → ddCt complete\n")
  cat("Files imported: ", length(x$files), "\n", sep = "")
  cat("Usable Ct rows: ", nrow(x$ct_table), "\n", sep = "")
  cat("Housekeeping gene: ", x$housekeeping_gene, "\n", sep = "")
  cat("Calibration samples: ", paste(x$calibration_samples, collapse = ", "), "\n", sep = "")
  cat("Expression table: result$expression\n")
  cat("Long cleaned Ct table: result$ct_table\n")
  cat("ddCt warnings: ", x$warning_file, "\n\n", sep = "")
  invisible(x)
}

run_quantstudio3_ddct <- function(folder_path,
                                  file_pattern = "^[0-9]{4}[-_]?[0-9]{2}[-_]?[0-9]{2}.*\\.txt$",
                                  housekeeping_gene = NULL, calibration_samples = NULL,
                                  undetermined_ct = 35, exclude_omitted = TRUE,
                                  warning_file = NULL, write_expression_csv = FALSE,
                                  output_dir = folder_path, print_summary = TRUE) {
  if (missing(folder_path) || !nzchar(folder_path)) {
    stop("Supply 'folder_path' containing one or more QuantStudio 3 .txt exports.", call. = FALSE)
  }
  imported <- read_quantstudio3_folder(
    folder_path = folder_path,
    file_pattern = file_pattern,
    undetermined_ct = undetermined_ct,
    exclude_omitted = exclude_omitted
  )
  if (is.null(housekeeping_gene)) {
    housekeeping_gene <- choose_one_value(
      imported$ct_table$Detector,
      "Choose ONE housekeeping gene"
    )
  }
  if (is.null(calibration_samples)) {
    calibration_samples <- choose_multiple_values(
      imported$ct_table$Sample,
      "Choose calibration sample(s); select Finish selection when done"
    )
  }
  if (is.null(warning_file)) warning_file <- file.path(folder_path, "ddct_warnings.txt")

  calculated <- ddct_expression_table(
    ct_table = imported$ct_table,
    housekeeping_gene = housekeeping_gene,
    calibration_samples = calibration_samples,
    warning_file = warning_file
  )
  expression_csv <- NULL
  if (isTRUE(write_expression_csv)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    expression_csv <- file.path(output_dir, "ddct_expression.csv")
    utils::write.csv(calculated$expression, expression_csv, row.names = FALSE)
  }

  result <- structure(
    list(
      ct_table = imported$ct_table,
      files = imported$files,
      import_summary = imported$import_summary,
      housekeeping_gene = housekeeping_gene,
      calibration_samples = calibration_samples,
      warning_file = warning_file,
      expression_csv = expression_csv,
      input_frame = calculated$input_frame,
      ddct_result = calculated$ddct_result,
      expression = calculated$expression
    ),
    class = "quantstudio3_ddct"
  )
  if (isTRUE(print_summary)) print(result)
  invisible(result)
}
