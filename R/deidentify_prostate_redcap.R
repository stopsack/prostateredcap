#' Deidentify Prostate REDCap Data
#'
#' @description
#'
#' Deidentifies the dataset for analysis and sharing:
#' * Drops all dates, keeping time intervals only
#' * Rounds age (in years, 1 decimal) and PSA (1 decimal)
#' * Replaces patient ID by a sequential index number that still allows for
#'   merging \code{pts} and \code{smp} datasets.
#'
#' By default, this function is already being called automatically within
#' \code{\link{load_prostate_redcap}}.
#'
#' @param data List with elements \code{pts} and \code{smp}, as generated
#'    by \code{\link{load_prostate_redcap}}.
#'
#' @return List:
#' * \code{pts}: Deidentified patient-level data.
#' * \code{smp}: Deidentified sample-level data.
#' @export
deidentify_prostate_redcap <- function(data) {
  pts <- data$pts %>%
    select(-dob, -dxdate, -met_date, -crpc_date, -lastfu, -adtstart, -bxdate,
           -freezedate, -lastvisit) %>%
    mutate(age_dx   = as.numeric(round(age_dx, digits = 1)),
           psa_dx   = as.numeric(round(psa_dx, digits = 1)),
           lnpsa_dx = as.numeric(round(lnpsa_dx,  digits = 2)),
           ptid2    = row_number())
  smp <- data$smp %>%
    left_join(pts %>% select(ptid, ptid2), by = "ptid") %>%
    select(-ptid, -smpdate) %>%
    mutate(age_smp = as.numeric(round(age_smp, digits = 1)),
           age_seq = as.numeric(round(age_seq, digits = 1))) %>%
    select(ptid = ptid2, everything()) %>%
    labelled::set_variable_labels(ptid    = "Patient ID",
                                  age_smp = "Age at sample (years)",
                                  age_smp = "Age at sequencing (years)")
  pts <- pts %>%
    select(-ptid) %>%
    select(ptid = ptid2, everything()) %>%
    labelled::set_variable_labels(ptid     = "Patient ID",
                                  age_dx   = "Age at diagnosis (years)",
                                  psa_dx   = "PSA at diagnosis (ng/ml)",
                                  lnpsa_dx = "PSA at diagnosis (ng/ml, log)")
  lst(pts, smp)
}
