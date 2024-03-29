#' Deidentify Prostate REDCap Data
#'
#' @description
#'
#' Deidentifies the dataset for analysis and sharing:
#' * Drops all dates, keeping time intervals only
#' * Rounds age (in years, 1 decimal) and PSA (1 decimal)
#' * Replaces patient ID by a sequential index number that still allows for
#'   merging \code{pts}, \code{smp}, and \code{trt} datasets.
#' * Removes free text with potentially identifying information (treatments).
#'
#' By default, this function is already being called automatically within
#' \code{\link{load_prostate_redcap}}.
#'
#' @param data List with elements \code{pts}, \code{smp}, and \code{trt},
#'   as generated by \code{\link{load_prostate_redcap}}.
#' @param trt_freetext Optional. Remove free text treatment names? Defaults to
#'   \code{TRUE}.
#' @param extra Optional. Tibble with additional data set, or list of additional
#'   data sets, in which a \code{ptid} variable of type \code{character} should
#'   be replaced with the same sequential ID as in \code{pts}. Other
#'   deidentification steps of such extra datasets must be done separately if
#'   necessary. Defaults to \code{NULL} (none).
#'
#' @return List:
#' * \code{pts}: Deidentified patient-level data.
#' * \code{smp}: Deidentified sample-level data.
#' * \code{trt}: Deidentified treatment data.
#' * \code{ext}: List of additional data sets (if provided, otherwise
#'   empty, \code{NULL}).
#' @export
deidentify_prostate_redcap <- function(data,
                                       trt_freetext = TRUE,
                                       extra = NULL) {
  pts <- data$pts %>%
    select(-dob, -dxdate, -met_date, -crpc_date, -lastfu, -adtstart, -bxdate,
           -freezedate, -lastvisit) %>%
    mutate(age_dx   = as.numeric(round(age_dx, digits = 1)),
           psa_dx   = as.numeric(round(psa_dx, digits = 1)),
           lnpsa_dx = as.numeric(round(lnpsa_dx,  digits = 2)),
           ptid2    = row_number())

  smp <- data$smp %>%
    left_join(pts %>% select(ptid, ptid2), by = "ptid") %>%
    select(-ptid, -smpdate, -seqdate) %>%
    mutate(age_smp = as.numeric(round(age_smp, digits = 1)),
           age_seq = as.numeric(round(age_seq, digits = 1))) %>%
    select(ptid = ptid2, everything()) %>%
    labelled::set_variable_labels(ptid    = "Patient ID",
                                  age_smp = "Age at sample (years)",
                                  age_smp = "Age at sequencing (years)")
  trt <- data$trt %>%
    left_join(pts %>% select(ptid, ptid2), by = "ptid") %>%
    select(-ptid, -rx_start, -rx_end) %>%
    select(ptid = ptid2, everything()) %>%
    labelled::set_variable_labels(ptid     = "Patient ID")
  if(trt_freetext == TRUE)
    trt <- trt %>% select(-rx_name_other)

  ext <- list()
  if(!is.null(extra)) {
    if(is.data.frame(extra))
      extra <- list(extra)
    ext <- map(.x = extra,
               .f = ~{ .x %>%
                   left_join(pts %>% select(ptid, ptid2), by = "ptid") %>%
                   select(-ptid) %>%
                   select(ptid = ptid2, everything()) %>%
                   labelled::set_variable_labels(ptid     = "Patient ID")
               })
  }

  pts <- pts %>%
    select(-ptid) %>%
    select(ptid = ptid2, everything()) %>%
    labelled::set_variable_labels(ptid     = "Patient ID",
                                  age_dx   = "Age at diagnosis (years)",
                                  psa_dx   = "PSA at diagnosis (ng/ml)",
                                  lnpsa_dx = "PSA at diagnosis (ng/ml, log)")

  lst(pts, smp, trt, ext)
}
