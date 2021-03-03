#' Quality Control Criteria for Patient-Level Data
#'
#' For use in \code{\link{check_prostate_redcap}}.
#'
#' @return Tibble:
#' * \code{label}: Description of the exclusion criterion.
#' * \code{filter_criterion}: Expression that will be used in
#'   \code{\link[dplyr]{filter}}
#'   for restriction to observations that pass the criterion.
#'
#' @export
#' @examples
#' # Show all criteria:
#' qc_criteria_pts()
#'
#' # Show code for criterion #3:
#' qc_criteria_pts()$filter_criterion[[3]]
qc_criteria_pts <- function() {
  tibble::tribble(
    ~label,                                   ~filter_criterion,
    "All patients",                           expr(TRUE),
    "Incomplete record",                      expr(complete_pts == "Complete"),
    "Missing date of birth or diagnosis",     expr(!is.na(age_dx)),
    "Metastatic/CRPC but no associated date", expr(!(((met_date_na == TRUE &
                                                         stage != "M1" &
                                                         is_met == "Yes") |
                                                        (crpc_date_na == TRUE &
                                                           is_crpc == "Yes")) &
                                                       lastvisit_na == TRUE)),
    "No lastvisit+met+CRPC date",             expr(!(met_date_na &
                                                       crpc_date_na &
                                                       lastvisit_na)),
    "Metastases before diagnosis",            expr(!(dx_met_mos <= 0 &
                                                       stage != "M1")),
    "Missing stage",                          expr(!is.na(stage))) %>%
    dplyr::mutate(index = dplyr::row_number())
}

#' Quality Control Criteria for Sample-Level Data
#'
#' For use in \code{\link{check_prostate_redcap}}.
#'
#' @return Tibble:
#' * \code{label}: Description of the exclusion criterion.
#' * \code{filter_criterion}: Expression that will be used in
#'   \code{\link[dplyr]{filter}}
#'   for restriction to observations that pass the criterion.
#' * \code{index}: Level of QC that can be provided to
#'   \code{\link{check_prostate_redcap}} in the \code{qc_level_pts} and
#'   \code{qc_level_smp} arguments.
#'
#' @export
#' @examples
#' # Show all criteria:
#' qc_criteria_smp()
#'
#' # Show code for criterion #4:
#' qc_criteria_smp()$filter_criterion[[4]]
qc_criteria_smp <- function() {
  tibble::tribble(
    ~label,                              ~filter_criterion,
    "All samples",                       expr(TRUE),
    "Samples without qc'd patient data", expr(ptid %in% qc_pts_data$ptid),
    "Incomplete record",                 expr(complete_smp == "Complete"),
    "Missing date of sample",            expr(!is.na(dx_smp_mos)),
    "Missing disease extent",            expr(!is.na(dzextent_smp)),
    "Sample date after last follow-up",  expr(!(smp_os_mos < 0)),
    "Localized/reg. nodes sample; met_date before sample", expr(
      !(primmet_smp == "Primary" & smp_met_mos < 0.5 & is_met_for_qc == "Yes")),
    "Metastatic sample; met_date after sample", expr(!(
      dzextent_smp %in% c("Metastatic castration-resistant",
                          "Metastatic hormone-sensitive",
                          "Metastatic, variant histology") &
        smp_met_mos > 0.5)),
    "Localized sample; stage N1 or M1",  expr(!(dzextent_smp == "Localized" &
                                                  stage_for_qc %in%
                                                  c("N1 M0", "M1"))),
    "Regional nodes sample; stage M1",   expr(!(dzextent_smp == "Regional nodes" &
                                                  stage_for_qc == "M1"))) %>%
    dplyr::mutate(index = dplyr::row_number())
}


#' Quality Control Check for Prostate REDCap database
#'
#' @description Run sequential quality control checks and
#'   thereby restrict the \code{pts} and \code{smp} datasets
#'   to cases that pass certain criteria.
#'
#' @param data List with elements \code{pts} and \code{smp},
#'   returned by \code{\link{load_prostate_redcap}}.
#' @param qc_crit_pts Criteria for checking the \code{pts} tibble.
#'   Defaults to the return of the \code{\link{qc_criteria_pts}} function.
#'   Custom criteria can be supplied instead.
#' @param qc_crit_smp Criteria for checking the \code{smp} tibble.
#'   Defaults to the return of the \code{\link{qc_criteria_smp}} function.
#'   Custom criteria can be supplied instead.
#' @param qc_level_pts Level of QC that the return \code{pts} tibble
#'   will be restricted to. By default, all QC steps in \code{qc_crit_pts}
#'   will be applied. A integer index (row) number of
#'   \code{qc_crit_pts} can be provided instead to perform
#'   less strict exclusions. \code{qc_level_pts = 1} will perform
#'   no exclusions.
#' @param qc_level_smp As \code{qc_level_pts}, for \code{smp}.
#' @param recommended_only Return qc'd \code{pts} and \code{smp}
#'   restricted to variables that are recommended for use in analyses?
#'   Defaults to \code{FALSE} but is recommended for use.
#'
#' @return List:
#'
#' * \code{pts}: Patient-level data after QC.
#' * \code{smp}: Sample-level data after QC.
#' * \code{qc_pts}: Tibble of sequential exclusions for \code{pts}.
#' * \code{qc_smp}: Tibble of sequential exclusions for \code{smp}.
#' The \code{qc_pts} and \code{qc_smp} tibbles can be used to extract
#' information on which records failed which QC step.
#'
#' @import dplyr purrr
#' @export
#' @seealso Overview of analysis-ready data elements:
#'   \url{https://stopsack.github.io/prostateredcap/articles/dataelements.html}
#'
#' @examples
#' \dontrun{
#' # Process output of load_prostate_redcap():
#' pts_smp_qc <- check_prostate_redcap(pts_smp, recommended_only = TRUE)
#' }
check_prostate_redcap <- function(data,
                                  qc_crit_pts = qc_criteria_pts(),
                                  qc_crit_smp = qc_criteria_smp(),
                                  qc_level_pts = NULL,
                                  qc_level_smp = NULL,
                                  recommended_only = FALSE) {
  if(!is.data.frame(data$pts) | !is.data.frame(data$smp))
    stop("Must provide a list with the elements 'pts' and 'smp', both data frames/tibbles.")

  # 'pts' patient-level data: generate tibble with sequential exclusions
  qc_pts <- qc_crit_pts %>%
    mutate(included     = accumulate(.x = filter_criterion,
                                     .f = ~filter(.data = .x, eval(.y)),
                                     .init = data$pts)[-1],
           n            = map_int(.x = included, .f = nrow),
           diff         = lag(n) - n,
           included_lag = lag(included),
           excluded     = map2(.x = included_lag, .y = included, .f = setdiff)) %>%
    select(-filter_criterion, -included_lag)
  # QC'd 'pts' tibble- either the last row of the above or a user-defined level of QC
  qc_pts_data <- qc_pts %>%
    slice(ifelse(is.null(qc_level_pts), nrow(.), qc_level_pts)) %>%
    pull(included) %>%
    pluck(1) %>%
    select(-crpc_date_na, -met_date_na, -lastvisit_na, -complete_pts) %>%
    labelled::copy_labels(from = data$smp)

  # 'smp' sample-level data: generate tibble with sequential exclusions
  qc_smp <- qc_crit_smp %>%
    mutate(included     = accumulate(.x = filter_criterion,
                                     .f = ~filter(.data = .x, eval(.y)),
                                     .init = data$smp)[-1],
           n            = map_int(.x = included, .f = nrow),
           diff         = lag(n) - n,
           included_lag = lag(included),
           excluded     = map2(.x = included_lag, .y = included, .f = setdiff)) %>%
    select(-filter_criterion, -included_lag)
  # QC'd 'smp' tibble- either the last row of the above or a user-defined level of QC
  qc_smp_data <- qc_smp %>%
    slice(ifelse(is.null(qc_level_smp), nrow(.), qc_level_smp)) %>%
    pull(included) %>%
    pluck(1) %>%
    select(-stage_for_qc, -is_met_for_qc, -complete_smp) %>%
    labelled::copy_labels(from = data$smp)

  # Optional: Return only variables in qc'd 'pts' and 'smp' tibbles
  # that are recommended for analyses
  if(recommended_only == TRUE) {
    qc_pts_data <- qc_pts_data %>%
      select(ptid, age_dx, race4, race3, smoking, bx_gl34, psa_dx, psa_dxcat,
             lnpsa_dx, stage, clin_tstage, clin_nstage, mstage, rxprim, rxprim_oth,
             rxprim_rp, rxprim_adt, rxprim_chemo, rxprim_xrt, rxprim_other,
             rp_gl34, path_t, path_n, is_crpc, crpc_event, is_met,
             met_event, is_dead, death_event)
    qc_smp_data <- qc_smp_data %>%
      select(ptid, dmpid, hist_smp, dzextent_smp, dzextent_seq,
             ext_pros, ext_lndis, ext_bone,
             ext_vis, bonevol, cntadt, tissue, smp_pros, smp_tissue, primmet_smp,
             age_smp, age_seq, dx_smp_mos, adt_smp_mos, dx_seq_mos, seq_met_mos,
             seq_crpc_mos, seq_os_mos, dzvol, denovom_smp, denovom_seq)
  }

  list(pts = qc_pts_data, smp = qc_smp_data,
       qc_pts = qc_pts, qc_smp = qc_smp)
}
