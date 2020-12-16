#' Load and check the MSK Prostate REDCap Prostate clinical database
#'
#' @description
#' This package loads, merges, reformats, corrects, labels, and quality controls
#' the MSK Prostate REDCap Prostate clinical database.
#' Run \code{\link{load_prostate_redcap}}, followed by \code{\link{check_prostate_redcap}}.
#'
#' @docType package
#' @name prostateredcap
#' @examples
#' \dontrun{
#' # Load labelled CSV exported from REDCap:
#' pts_smp <- load_prostate_redcap(
#'   labelled_csv = "GUPIMPACTDatabaseFre_DATA_LABELS_2020-01-01_0001.csv")
#'
#' # Perform quality control checks and exclusions:
#' pts_smp_qc <- check_prostate_redcap(pts_smp, recommended_only = TRUE)
#'
#' # Access final patient-level data:
#' pts_smp_qc$pts
#'
#' # Access final sample-level data:
#' pts_smp_qc$smp
#'
#' # View the QC results and exclusion criteria:
#' pts_smp_qc$qc_pts
#' pts_smp_qc$qc_smp
#' }
NULL
#> NULL
