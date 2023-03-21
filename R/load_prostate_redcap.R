#' Guess Dates
#'
#' @description Uses \code{\link[lubridate]{parse_date_time}}
#'   to transform string vectors with different date formats into
#'   an R date vector, printed in ISO format.
#'   If only the year is provided, then
#'   the date will be set as June 30 of that year.
#'
#' @param messydate Date in the format MM/DD/YYYY,
#'   MM/YYYY, or YYYY only.
#'
#' @return Date vector in ISO date notation.
#' @export
#'
#' @examples
#' guessdate(c("03/28/2019", "03/2019", "03/19", "2019"))
guessdate <- function(messydate) {
  # set year-only dates to the mid-point of the year
  messydate <- dplyr::if_else(stringr::str_length(messydate) == 4 &
                                !stringr::str_detect(string = messydate,
                                                     pattern = "/") &
                                !stringr::str_detect(messydate, pattern = "-"),
                              true = paste0("06/30", messydate),
                              false = messydate)
  lubridate::as_date(lubridate::parse_date_time(messydate,
                                                orders = c("mdy", "my", "y"),
                                                quiet = TRUE))
}

#' Recode "N/A", "Unknown", etc. to NA
#'
#' @description Replaces \code{c("Unknown / Not Reported", "N/A", "NA",
#' "Unknown", "X", "x")} with \code{NA}.
#'
#' @param x String or factor vector
#'
#' @return Original vector
#' @noRd
unknowns <- function(x) {
  ifelse(x %in% c("Unknown / Not Reported", "N/A", "NA", "Unknown", "X", "x"),
         yes = NA, no = x)
}

#' Load MSK-IMPACT Prostate REDCap Labeled CSV File
#'
#' @description Loads, merges, reformats, corrects, and labels
#'   the REDCap file used for the MSK-IMPACT Prostate clinical database.
#'   It is recommended that the returned list is next processed by
#'   \code{\link{check_prostate_redcap}}.
#'
#' @param labeled_csv CSV file with labels, exported from REDCap.
#'   Must be the labeled version and must contain dates
#'   in order to derive time intervals.
#' @param deidentify De-identify the returned data set using
#'   \code{\link{deidentify_prostate_redcap}}? Defaults to
#'   \code{TRUE}. Should only be disabled if additional data need
#'   to be merged by identifiers, followed by calling
#'   \code{\link{deidentify_prostate_redcap}} separately.
#' @param keep_also Optional. Additional patient-level variables to keep
#'   without editing. As applicable, they would need to be
#'   deidentified manually.
#'   Provide as list with vectors of variable names for baseline and freeze
#'   forms: \code{list(baseline = c("var1", "var2"), freeze = "varX")}.
#'
#' @details The following edits and assumptions are made:
#' 1. Potentially incomplete date variables are converted to
#'    date format, using \code{\link{guessdate}}.
#' 2. Various missingness indicators in strings and factors,
#'    \code{c("Unknown / Not Reported", "N/A", "NA", "Unknown", "X", "x")},
#'    are converted to \code{NA}.
#' 3. "Undetectable" PSA is set to 0, PSA \code{">x"} is set to \code{x + 1},
#'    PSA \code{"a-b"} (e.g., \code{4.5-4.7}) is set to the mean of the two
#'    values.
#' 4. Clinical T and N stage variables are set to missing if M1.
#' 5. Event dates and follow-up time for metastases (\code{met_date}),
#'    castration resistance (\code{crpc_date}), and death are set:
#'    *  Event date is the last clinic visit (\code{lastvisit})
#'       if a CRPC/metastases event has not occurred.
#'    *  Event date is the last follow up/contact (\code{lastfu})
#'       if last known survival status is alive.
#'    *  If stage is M1 and the recorded metastasis date is no more than
#'       1 month discrepant, \code{met_date} is set to the diagnosis
#'       date (\code{dxdate}).
#'    *  If the sample is a variant histology (e.g., neuroendocrine),
#'       the castration resistance date (\code{crpc_date}) is the date of
#'       diagnosis and the event indicator for survival analyses
#'       (\code{event_crpc}) is \code{NA}.
#'    *  Time intervals for these three survival outcomes are calculated
#'       from the time of sequencing. For late-entry survival models,
#'       time intervals from diagnosis to sequencing and from sample/biopsy to
#'       sequencing are also provided.
#' 6. Disease extent, distinguishing CRPC from castration-sensitive disease,
#'    at sampling is based on the sample date and the date of castration
#'    resistance. If the samples was obtained before the CRPC date, or
#'    CRPC did not occur, the sample is from castration-sensitive disease
#'    by definition.
#'
#' @return List of three labeled tibbles (data frames):
#'
#' * \code{pts}: Patient-level data
#'
#' * \code{smp}: Sample-level data
#'
#' * \code{trt}: Treatment data
#'
#' Access variables labels in RStudio via \code{\link[utils]{View}}
#' or using \code{\link{attr}(., "label")}.
#'
#' The warning message, \code{Duplicated column names deduplicated}, is
#' expected due to the design of the REDCap dataset. Another warning message
#' that a factor does not contain all levels is also possible.
#'
#' @import dplyr stringr purrr forcats
#' @export
#'
#' @seealso Overview of analysis-ready data elements:
#'   \url{https://stopsack.github.io/prostateredcap/articles/dataelements.html}
#'
#' @examples
#' # Get path to toy data provided by the package:
#' example_csv_file <- system.file("extdata",
#'   "SampleGUPIMPACTDatab_DATA_LABELS_2021-05-26.csv",
#'   package = "prostateredcap",
#'   mustWork = TRUE)
#'
#' # Load data:
#' pts_smp <- load_prostate_redcap(labeled_csv = example_csv_file)
#'
#' # Access patient-level data:
#' pts_smp$pts
#'
#' # Access sample-level data:
#' pts_smp$smp
#'
#' # Access treatment data:
#' pts_smp$trt
#'
#' # Pass 'pts_smp' to check_prostate_redcap() next
load_prostate_redcap <- function(labeled_csv,
                                 deidentify = TRUE,
                                 keep_also  = list(baseline = NULL,
                                                   sample   = NULL,
                                                   freeze   = NULL)) {
  # Read REDCap file
  rcclin <- readr::read_csv(
    file = labeled_csv,
    col_types = readr::cols(
      .default = readr::col_character()),
    name_repair = ~vctrs::vec_as_names_legacy(., sep = "_")) %>%
    select(`Record ID`:`Complete?_3`) %>%
    rename(ptid = `Record ID`)

  # Patient-level baseline form
  pts_baseline <- rcclin %>%
    filter(is.na(`Repeat Instrument`)) %>%
    select(ptid,
           complete_pts = `Complete?`,
           dob        = `Birth Date`,
           race       = `Race`,
           ethnicity  = `Ethnicity`,
           smoking    = `Smoking Status at Diagnosis`,
           dxdate     = `Date of Initial Diagnosis`,
           bxdate     = `Date of Initial Prostate Biopsy`,
           bx_gl_sum  = `Sum Gleason at Diagnosis (Biopsy)`,
           bx_gl_maj  = `Primary Gleason Pattern at Diagnosis (Biopsy)`,
           bx_gl_min  = `Secondary Gleason Pattern at Diagnosis (Biopsy)`,
           psa_dx     = `PSA at Diagnosis`,
           clin_t     = `Clinical T Stage at Diagnosis`,
           clin_n     = `Clinical N Stage at Diagnosis (regional lymph node metastases)`,
           clin_m     = `Clinical M Stage at Diagnosis`,
           rxprim     = `Primary Therapy`,
           rxprim_oth = `Other Primary Therapy`,
           rp_gl_sum  = `Sum Gleason at Prostatectomy`,
           rp_gl_maj  = `Primary Gleason Pattern at Prostatectomy`,
           rp_gl_min  = `Secondary Gleason Pattern at Prostatectomy`,
           path_t     = `Pathologic T Stage at Diagnosis`,
           path_n     = `Pathologic N Stage at Diagnosis`,
           any_of(keep_also$baseline))

  # Patient-level "freeze" (outcome) form
  pts_freeze <- rcclin %>%
    filter(`Repeat Instrument` == "Freeze Data") %>%
    select(ptid,
           freezedate = `Freeze Date`,
           adtstart   = `Continuous ADT Start Date`,
           is_crpc    = `Castration Resistant`,
           crpc_date  = `Castration Resistance Date`,
           is_met     = `Metastatic`,
           met_date   = `Date of Metastasis`,
           lastvisit  = `Last MD Visit Date`,
           is_dead    = `Survival Status`,
           lastfu     = `Date of Death/Last Contact`,
           any_of(keep_also$freeze))

  # Combined patient-level data
  pts <- left_join(pts_baseline, pts_freeze, by = "ptid") %>%  # allow for missing freeze form
    mutate_at(.vars = vars(dob, dxdate, bxdate, freezedate, adtstart,
                           crpc_date, met_date, lastvisit, lastfu),
              .funs = guessdate) %>%
    mutate_at(.vars = vars(race, ethnicity, smoking, is_crpc, is_met,
                           starts_with("bx_"), starts_with("rp_"),
                           starts_with("clin_"), starts_with("path_")),
              .funs = unknowns) %>%
    mutate_at(.vars = vars(starts_with("bx_"), starts_with("rp_")),
              .funs = as.numeric) %>%
    mutate_at(.vars = vars(race, ethnicity, smoking, rxprim,
                           is_met, is_crpc, is_dead,
                           starts_with("clin_"), starts_with("path_")),
              .funs = as.factor) %>%
    # PSA
    mutate(
      psa_dx = suppressWarnings(case_when(
        str_sub(psa_dx, start = 1, end = 1) == ">" ~
          as.numeric(str_sub(psa_dx, start = 2, end = 10)) + 1,
        str_detect(string = psa_dx, pattern = "-") ~  # two numbers, e.g. "4.5-5.2"
          mean(as.numeric(str_split(string = psa_dx, pattern = "-", n = 2,
                                    simplify = TRUE)[1]),
               as.numeric(str_split(string = psa_dx, pattern = "-", n = 2,
                                    simplify = TRUE)[2])),
        str_detect(string = psa_dx, pattern = "etect") |
          str_detect(string = psa_dx, pattern = "dec") ~
          0,  # "undetectable", also misspelled "undectable"
        TRUE ~ as.numeric(psa_dx))),
      psa_dxcat = factor(cut(psa_dx, breaks = c(0, 4, 10, 20, Inf),
                            include.lowest = TRUE),
                        labels = c("<4", "4-10", "10-20", ">20")),
      # Gleason grade
      bx_gl34    = factor(case_when(  # grade groups with splitting of score 7
        bx_gl_sum < 7                   ~ "<7",
        bx_gl_maj == 3 & bx_gl_min == 4 ~ "3+4",
        bx_gl_maj == 4 & bx_gl_min == 3 ~ "4+3",
        bx_gl_sum == 8                  ~ "8",
        bx_gl_sum > 8                   ~ "9-10")),
      bx_gl    = factor(case_when(    # Gleason scores, lumping 3+4 and 4+3 as score 7
        bx_gl_sum < 7      ~ "<7",
        bx_gl_sum %in% 7:8 ~ as.character(bx_gl_sum),
        bx_gl_sum == 8     ~ "8",
        bx_gl_sum > 8      ~ "9-10")),
      rp_gl34    = factor(case_when(  # prostatectomy Gleason, treat like bx_gl34
        rp_gl_sum < 7                   ~ "<7",
        rp_gl_maj == 3 & rp_gl_min == 4 ~ "3+4",
        rp_gl_maj == 4 & rp_gl_min == 3 ~ "4+3",
        rp_gl_sum == 8                  ~ "8",
        rp_gl_sum > 8                   ~ "9-10")),
      # Stage
      stage_detailed = factor(case_when(
        clin_t %in% c("1a", "1b", "1c", "2", "2a", "2b", "2c") & clin_n == "0" &
          clin_m == "0" ~ "T1/T2 N0 M0",
        clin_t %in% c("1a", "1b", "1c", "2", "2a", "2b", "2c") & is.na(clin_n) &
          clin_m == "0" ~ "T1/T2 NX M0",
        clin_t %in% c("3",  "3a", "3b")                        & clin_n == "0" &
          clin_m == "0" ~ "T3 N0 M0",
        clin_t %in% c("3",  "3a", "3b")                        & is.na(clin_n) &
          clin_m == "0" ~ "T3 NX M0",
        clin_t == "4"                                                          &
          clin_m == "0" ~ "T4 M0",
        is.na(clin_t)                                          & clin_n == "0" &
          clin_m == "0" ~ "TX N0 M0",
        is.na(clin_t)                                          & is.na(clin_n) &
          clin_m == "0" ~ "TX NX M0",
        clin_n == "1" & clin_m == "0" ~ "N1 M0",
        clin_m %in% c("1", "1a", "1b", "1c") ~ "M1")),
      stage_detailed = fct_relevel(stage_detailed,
                                   "T1/T2 N0 M0",
                                   "T1/T2 NX M0",
                                   "T3 N0 M0",
                                   "T3 NX M0",
                                   "T4 M0",
                                   "TX N0 M0",
                                   "TX NX M0",
                                   "N1 M0",
                                   "M1"),
      stage = fct_other(stage_detailed, keep = c("N1 M0", "M1"),
                        other_level = "N0/NX M0"),
      clin_tstage = factor(case_when(
        clin_t %in% c("1a", "1b", "1c", "2", "2a", "2b", "2c") &
          (clin_n == "0" | is.na(clin_n)) & clin_m == "0" ~ "T1/T2",
        clin_t %in% c("3",  "3a", "3b") &
          (clin_n == "0" | is.na(clin_n)) & clin_m == "0" ~ "T3",
        clin_t == "4" &
          (clin_n == "0" | is.na(clin_n)) & clin_m == "0" ~ "T4"),
        levels = c("T1/T2", "T3", "T4")),
      clin_nstage = factor(case_when(
        clin_n == "0" & clin_m == "0" ~ "N0 M0",
        clin_n == "1" & clin_m == "0" ~ "N1 M0"),
        levels = c("N0 M0", "N1 M0")),
      mstage = factor(case_when(
        clin_m == "0"                        ~ "M0",
        clin_m %in% c("1", "1a", "1b", "1c") ~ "M1"),
        levels = c("M0", "M1")),
      # Primary treatment
      rxprim_rp    = if_else(grepl("RP", rxprim) | grepl("Surgery", rxprim),   TRUE, FALSE),
      rxprim_adt   = if_else(grepl("ADT", rxprim),                             TRUE, FALSE),
      rxprim_chemo = if_else(grepl("Chemo", rxprim),                           TRUE, FALSE),
      rxprim_xrt   = if_else(grepl("RT", rxprim) | grepl("Radiation", rxprim), TRUE, FALSE),
      rxprim_other = if_else(!is.na(rxprim_oth), TRUE, FALSE),
      # If crpc_date, met_date, or last MD visit are (erroneously) after lastfu,
      # then update lastfu
      lastfu = if_else(lastfu > crpc_date | is.na(crpc_date),
                       true = lastfu, false = crpc_date),
      lastfu = if_else(lastfu > met_date | is.na(met_date),
                       true = lastfu, false = met_date),
      lastfu = case_when(lastfu < lastvisit & is_dead == "Alive" ~ lastvisit,
                         TRUE ~ lastfu),
      # Set missing CRPC and met dates to last visit date:
      crpc_date    = case_when(
        # Special case: Histologies that are castration resistant by definition
        # (e.g., neuroendocrine). Set CRPC date as the date of diagnosis.
        str_detect(string = is_crpc, pattern = "istology") |
          str_detect(string = is_crpc, pattern = "euroendo") ~ dxdate,
        is.na(crpc_date) ~ lastvisit,
        TRUE             ~ crpc_date),
      met_date = case_when(
        # Change only if <1 mo discrepancy between met_date and dxdate
        mstage == "M1" &
          abs(lubridate::interval(dxdate, met_date) / months(1)) <= 1 ~
          dxdate,  # Metastatic at diagnosis.
        is.na(met_date) & is_met == "No"   ~ lastvisit,  # censor at last visit
        # this case should not exist:
        is.na(met_date) & is_met == "Yes"  ~ as.Date(NA_real_),
        !is.na(met_date) & is_met == "Yes" ~ met_date),
      is_met   = factor(case_when(
        mstage == "M1" ~ "Yes",
        TRUE           ~ as.character(is_met))),
      # CRPC indicator for survival analyses is missing for variant histologies:
      crpc_event   = case_when(is_crpc == "Yes" ~ 1,
                               is_crpc == "No"  ~ 0),
      met_event    = case_when(is_met == "Yes"  ~ 1,
                               is_met == "No"   ~ 0),
      death_event  = if_else(is_dead == "Dead", 1, 0),
      mfs_event = if_else(is_dead == "Dead" | is_met == "Yes", 1, 0),
      is_mfs = factor(
        mfs_event,
        levels = 0:1,
        labels = c("Event-free", "Metastasis/death")),
      # Time intervals
      age_dx        = lubridate::interval(dob,     dxdate)    / lubridate::years(1),
      dx_bx_mos    = lubridate::interval(dxdate,   bxdate)    / months(1),
      dx_adt_mos   = lubridate::interval(dxdate,   adtstart)  / months(1),
      adt_crpc_mos = lubridate::interval(adtstart, crpc_date) / months(1),
      dx_crpc_mos  = lubridate::interval(dxdate,   crpc_date) / months(1),
      dx_met_mos   = lubridate::interval(dxdate,   met_date)  / months(1),
      dx_os_mos    = lubridate::interval(dxdate,   lastfu)    / months(1),
      adt_os_mos   = lubridate::interval(adtstart, lastfu)    / months(1),
      met_os_mos   = lubridate::interval(met_date, lastfu)    / months(1),
      # Metastasis-free survival: until met or death (if no met). If neither,
      # only until last visit (NOT last contact, where met status is unknown)
      dx_mfs_mos   = case_when(
        is_met == "Yes" ~ dx_met_mos,
        is_met == "No" & is_dead == "Dead" ~ dx_os_mos,
        TRUE ~ dx_met_mos),
      crpc_os_mos  = if_else(is_crpc == "Yes",
                             true = lubridate::interval(crpc_date, lastfu) / months(1),
                             false = NA_real_),
      # Combine race categories and restrict to 3 race categories
      race4        = fct_lump(f = race, n = 3),
      race3        = fct_recode(fct_recode(race4, NULL = "Other"),
                                "Black" = "Black or African American"),
      race3        = fct_relevel(race3, "Asian", "White", "Black"),
      lnpsa_dx        = log(if_else(psa_dx == 0, true = 0.01, false = psa_dx)),
      smoke01      = case_when(smoking == "Current"              ~ 1,
                               smoking %in% c("Former", "Never") ~ 0)) %>%
    # for qc function, allowing to use the deidentified dataset
    mutate_at(.vars = vars(met_date, crpc_date, lastvisit),
              .funs = list(na = is.na))

  # Sample forms
  smp <- rcclin %>%
    select(ptid, `DMP ID`:`Complete?_1`) %>%
    rename(dmpid = `DMP Sample ID`,
           complete_smp = `Complete?_1`) %>%
    filter(!is.na(dmpid)) %>%  # get sample form only
    # fix data entry errors in dmpid, fourth part:
    tidyr::separate(col = "dmpid", into = paste0("dmpid_part", 1:4),
                    sep = "-") %>%
    mutate(dmpid_part4 = if_else(dmpid_part4 %in% c("IM06", "iM6", "Im6"),
                                 true = "IM6", false = dmpid_part4)) %>%
    transmute(
      ptid      = ptid,
      complete_smp = complete_smp,
      dmpid     = paste(dmpid_part1, dmpid_part2, dmpid_part3, dmpid_part4,
                        sep = "-"),
      smpdate   = guessdate(`Date of Collection`),
      seqdate   = guessdate(`DMP Date of Receipt`),
      hist_smp  = factor(unknowns(`Histology for Sample`)),
      hist_cmt  = `Histology (Other)`,
      dzextent_smp = factor(unknowns(`Extent of Disease at Collection`)),
      dzextent2 = dzextent_smp,
      ext_pros  = factor(if_else(
        `Sites of Disease (choice=Prostate/Prostate Bed)` == "Checked",
        TRUE, FALSE)),
      ext_lndis = factor(if_else(
        `Sites of Disease (choice=LN (distant))` ==
          "Checked", TRUE, FALSE)),
      ext_bone  = factor(if_else(
        `Sites of Disease (choice=Bone)` == "Checked",
        TRUE, FALSE)),
      ext_vis   = factor(if_else(
        `Sites of Disease (choice=Liver)` == "Checked" |
          `Sites of Disease (choice=Lung)` == "Checked" |
          `Sites of Disease (choice=Other Soft Tissue)` ==
          "Checked",
        TRUE, FALSE)),
      ext_liver   = factor(if_else(
        `Sites of Disease (choice=Liver)` == "Checked",
        TRUE, FALSE)),
      ext_lung   = factor(if_else(
        `Sites of Disease (choice=Lung)` == "Checked",
        TRUE, FALSE)),
      ext_other   = factor(if_else(
        `Sites of Disease (choice=Other Soft Tissue)` == "Checked",
        TRUE, FALSE)),
      bonevol   = factor(unknowns(`Volume of Bone Metastases at Time of Collection`)),
      cntadt    = factor(unknowns(`Continuous ADT`)),
      tissue    = factor(unknowns(`Sample Type`)),
      smp_pros  = fct_other(tissue, keep = "Prostate",
                            other_level = "Non-prostate"),
      smp_tissue  = fct_other(tissue, keep = c("Prostate", "Bone", "Lymph Node",
                                             "Liver", "Lung"),
                              other_level = "Other soft tissue"),
      smp_tissue  = fct_recode(smp_tissue, `Lymph node` = "Lymph Node"),
      smp_tissue  = fct_collapse(smp_tissue, Visceral = c("Liver", "Lung")),
      smp_tissue  = fct_relevel(smp_tissue,
                                "Prostate", "Lymph node", "Bone",
                                "Visceral", "Other soft tissue"),
      pur_rev   = factor(`Reviewed for Tumor Purity`),
      pur_remov = factor(`Removed for Low Tumor Purity`),
      select(., any_of(keep_also$sample))) %>%
    # Derive variables for which data from "pts" are needed:
    left_join(pts %>% select(ptid, stage, age_dx, dxdate, met_date,
                             is_met_for_qc = is_met, is_dead,
                             crpc_date, is_crpc, lastfu, adtstart),
              by = "ptid") %>%
    mutate(
      # Define disease extent at sequencing based on sample data plus
      # CRPC/mets data between sample and sequencing
      dzextent_seq = factor(case_when(
        dzextent_smp %in% c("Localized", "Regional nodes") &
          is_crpc != "No" & crpc_date <= seqdate &
          (is_met_for_qc != "Yes" | met_date > seqdate) ~
          "Non-metastatic castration-resistant",
        # NB, does not capture progression to Regional Nodes between sample
        # and sequencing:
        dzextent_smp %in% c("Localized", "Regional nodes") &
          (is_crpc == "No" | crpc_date > seqdate) &
          (is_met_for_qc != "Yes" | met_date > seqdate) ~
          as.character(dzextent_smp),
        (dzextent_smp == "Metastatic" |
           (is_met_for_qc == "Yes" & met_date <= seqdate)) &
          grepl("N/A-", is_crpc, fixed = TRUE) ~
          "Metastatic, variant histology",
        (dzextent_smp == "Metastatic" |
           (is_met_for_qc == "Yes" & met_date <= seqdate)) &
          (is_crpc == "No" | (is_crpc == "Yes" & crpc_date > seqdate)) ~
          "Metastatic hormone-sensitive",
        (dzextent_smp == "Metastatic" |
           (is_met_for_qc == "Yes" & met_date <= seqdate)) &
          (is_crpc != "No" & crpc_date <= seqdate) ~
          "Metastatic castration-resistant")),
      dzextent_seq = fct_relevel(dzextent_seq,
                                 "Localized",
                                 "Regional nodes",
                                 "Metastatic hormone-sensitive",
                                 "Non-metastatic castration-resistant",
                                 "Metastatic castration-resistant",
                                 "Metastatic, variant histology"),
      # Define disease extent at sample based on sample data plus CRPC data
      dzextent_smp = factor(case_when(
        dzextent_smp %in% c("Localized", "Regional nodes") &
          is_crpc == "Yes" & crpc_date <= smpdate ~
          "Non-metastatic castration-resistant",
        dzextent_smp %in% c("Localized", "Regional nodes") &
          (is_crpc != "Yes" | crpc_date > smpdate) ~
          as.character(dzextent_smp),
        dzextent_smp == "Metastatic" &
          grepl("N/A-", is_crpc, fixed = TRUE) ~
          "Metastatic, variant histology",
        dzextent_smp == "Metastatic" &
          (is_crpc == "No" | (is_crpc == "Yes" & crpc_date > smpdate)) ~
          "Metastatic hormone-sensitive",
        dzextent_smp == "Metastatic" &
          (is_crpc != "No" & crpc_date <= smpdate) ~
          "Metastatic castration-resistant")),
      dzextent_smp = fct_relevel(dzextent_smp,
                                 "Localized",
                                 "Regional nodes",
                                 "Metastatic hormone-sensitive",
                                 "Non-metastatic castration-resistant",
                                 "Metastatic castration-resistant",
                                 "Metastatic, variant histology"),
      primmet_smp = case_when(
        dzextent_smp %in% c("Localized", "Regional nodes") ~ "Primary",
        dzextent_smp %in% c("Metastatic castration-resistant",
                            "Metastatic hormone-sensitive",
                            "Metastatic, variant histology") ~ "Metastatic"),
      age_smp        = age_dx + lubridate::interval(dxdate, smpdate) /
        lubridate::years(1),
      age_seq        = age_dx + lubridate::interval(dxdate, seqdate) /
        lubridate::years(1),
      dx_smp_mos    = lubridate::interval(dxdate,   smpdate)   / months(1),
      adt_smp_mos   = lubridate::interval(adtstart, smpdate)   / months(1),
      dx_seq_mos    = lubridate::interval(dxdate,   seqdate)   / months(1),
      adt_seq_mos   = lubridate::interval(adtstart, seqdate)   / months(1),
      smp_met_mos   = lubridate::interval(smpdate,  met_date)  / months(1),
      smp_os_mos    = lubridate::interval(smpdate,  lastfu)    / months(1),
      seq_met_mos   = lubridate::interval(seqdate,  met_date)  / months(1),
      seq_crpc_mos  = lubridate::interval(seqdate,  crpc_date) / months(1),
      seq_os_mos    = lubridate::interval(seqdate,  lastfu)    / months(1),
      met_seq_mos   = case_when(
        seqdate > met_date ~
          lubridate::interval(met_date, seqdate)  / months(1)),
      seq_mfs_mos   = case_when(
        # Metastasis-free survival: define only if non-metastatic at sequencing
        !(dzextent_seq %in% c("Localized", "Regional nodes")) ~
            NA_real_,
        # Follow until met or death (if no met). If neither, follow
        # only until last visit (NOT last contact, where met status is unknown)
        is_met_for_qc == "Yes" ~
          lubridate::interval(seqdate, met_date) / months(1),
        is_met_for_qc == "No" & is_dead == "Dead" ~
          lubridate::interval(seqdate, lastfu) / months(1),
        TRUE ~
          lubridate::interval(seqdate, met_date) / months(1)),
      dzvol         = factor(case_when(
        str_detect(string = as.character(dzextent_smp), pattern = "Metastatic") &
          (bonevol != "High-Volume Bone Metastases" | is.na(bonevol)) &
          ext_vis == FALSE ~
          "Low-volume disease",
        str_detect(string = as.character(dzextent_smp), pattern = "Metastatic") &
          (bonevol == "High-Volume Bone Metastases" | ext_vis == TRUE) ~
          "High-volume disease")),
      dzvol = fct_relevel(dzvol, "Low-volume disease", "High-volume disease"),
      denovom_smp   = factor(case_when(
        str_detect(string = as.character(dzextent_smp), pattern = "Metastatic") &
          stage != "M1" ~ "Metastatic recurrence",
        str_detect(string = as.character(dzextent_smp), pattern = "Metastatic") &
          stage == "M1" ~ "De-novo metastatic")),
      denovom_seq   = factor(case_when(
        str_detect(string = as.character(dzextent_seq), pattern = "Metastatic") &
          stage != "M1" ~ "Metastatic recurrence",
        str_detect(string = as.character(dzextent_seq), pattern = "Metastatic") &
          stage == "M1" ~ "De-novo metastatic")),
      stage_for_qc = stage) %>%
    # Left-censor follow-up for times starting at sequencing
    mutate_at(.vars = vars(seq_os_mos, seq_crpc_mos, seq_met_mos, seq_mfs_mos),
              .funs = ~if_else(. <= 0, true = NA_real_, false = .)) %>%
    # remove variables from "pts":
    select(-stage, -age_dx, -dxdate, -met_date, -crpc_date,
           -is_crpc, -is_dead, -lastfu, -adtstart)


  ####
  # Treatment data

  # add empty columns if variables missing in data set
  added_cols <- tibble(
    `Extent of Disease before starting PARPi`           = NA_character_,
    `Sites of Disease (choice=Prostate/Prostate Bed)_1` = NA_character_,
    `Sites of Disease (choice=LN (distant))_1`          = NA_character_,
    `Sites of Disease (choice=Bone)_1`                  = NA_character_,
    `Sites of Disease (choice=Liver)_1`                 = NA_character_,
    `Sites of Disease (choice=Lung)_1`                  = NA_character_,
    `Sites of Disease (choice=Other Soft Tissue)_1`     = NA_character_,
    `Volume of Bone Metastases`                         = NA_character_)

  trt <- rcclin %>%
    filter(`Repeat Instrument` == "Treatment Data") %>%
    select(ptid, `Repeat Instance`, `Treatment Name`:`Complete?_3`)

  trt <- trt %>%
    tibble::add_column(added_cols %>%
                         select(-dplyr::any_of(names(trt)))) %>%
    transmute(
      ptid                 = ptid,
      rx_line              = `Repeat Instance`,
      rx_name              = factor(unknowns(`Treatment Name`)),
      rx_name_other        = `Treatment Name (Other)`,
      rx_name_parpi        = `PARPi Agent Name`,
      rx_start             = guessdate(`Treatment Start Date`),
      rx_end               = guessdate(`Treatment End Date/Last Known Treatment Date`),
      rx_censor            = `Treatment Ongoing`,
      rx_stop_reason       = factor(unknowns(`Reason for Treatment Stop`)),
      rx_stop_reason_other = `Reason for Treatment Stop (Other)`,
      rx_dzextent          = factor(unknowns(`Extent of Disease before starting PARPi`)),
      rx_ext_pros  = factor(if_else(
        `Sites of Disease (choice=Prostate/Prostate Bed)_1` == "Checked",
        TRUE, FALSE)),
      rx_ext_lndis = factor(if_else(
        `Sites of Disease (choice=LN (distant))_1` ==
          "Checked", TRUE, FALSE)),
      rx_ext_bone  = factor(if_else(
        `Sites of Disease (choice=Bone)_1` == "Checked",
        TRUE, FALSE)),
      rx_ext_vis   = factor(if_else(
        `Sites of Disease (choice=Liver)_1` == "Checked" |
          `Sites of Disease (choice=Lung)_1` == "Checked" |
          `Sites of Disease (choice=Other Soft Tissue)_1` ==
          "Checked",
        TRUE, FALSE)),
      rx_ext_liver   = factor(if_else(
        `Sites of Disease (choice=Liver)_1` == "Checked",
        TRUE, FALSE)),
      rx_ext_lung   = factor(if_else(
        `Sites of Disease (choice=Lung)_1` == "Checked",
        TRUE, FALSE)),
      rx_ext_other   = factor(if_else(
        `Sites of Disease (choice=Other Soft Tissue)_1` == "Checked",
        TRUE, FALSE)),
      rx_bonevol   = factor(unknowns(`Volume of Bone Metastases`))) %>%
    # Derive variables for which data from "pts" are needed:
    left_join(pts %>% select(ptid, age_dx, dxdate, met_date,
                             crpc_date, is_crpc, lastfu, adtstart),
              by = "ptid") %>%
    mutate(
      dx_rx_start_mos = lubridate::interval(dxdate,   rx_start) / months(1),
      dx_rx_end_mos   = lubridate::interval(dxdate,   rx_end)   / months(1),
      rx_wks          = lubridate::interval(rx_start, rx_end)   / lubridate::weeks(1),
      rx_dzextent = factor(case_when(
        rx_dzextent %in% c("Localized", "Regional nodes") &
          is_crpc == "Yes" & crpc_date <= rx_start ~
          "Non-metastatic castration-resistant",
        rx_dzextent %in% c("Localized", "Regional nodes") &
          (is_crpc != "Yes" | crpc_date > rx_start) ~
          as.character(rx_dzextent),
        rx_dzextent == "Metastatic" &
          grepl("N/A-", is_crpc, fixed = TRUE) ~
          "Metastatic, variant histology",
        rx_dzextent == "Metastatic" &
          (is_crpc == "No" | (is_crpc == "Yes" & crpc_date > rx_start)) ~
          "Metastatic hormone-sensitive",
        rx_dzextent == "Metastatic" &
          (is_crpc != "No" & crpc_date <= rx_start) ~
          "Metastatic castration-resistant")),
      rx_dzextent = fct_relevel(rx_dzextent,
                                "Localized",
                                "Regional nodes",
                                "Metastatic hormone-sensitive",
                                "Non-metastatic castration-resistant",
                                "Metastatic castration-resistant",
                                "Metastatic, variant histology")) %>%
    # remove variables from "pts":
    select(-age_dx, -dxdate, -met_date, -crpc_date,
           -is_crpc, -lastfu, -adtstart)


  pts <- pts %>% labelled::set_variable_labels(
    ptid           = "Patient ID",
    complete_pts   = "Record completely entered",
    age_dx         = "Age at diagnosis (years)",
    race           = "Self-reported race",
    race4          = "Self-reported race",
    race3          = "Self-reported race",
    ethnicity      = "Ethnicity",
    smoking        = "Smoking status",
    smoke01        = "Ever-smoker",
    bx_gl_sum      = "Gleason score",
    bx_gl          = "Gleason grade",
    bx_gl34        = "Gleason grade",
    bx_gl_maj      = "Gleason grade, major pattern",
    bx_gl_min      = "Gleason grade, minor pattern",
    psa_dx         = "PSA at diagnosis (ng/ml)",
    psa_dxcat      = "PSA at diagnosis (ng/ml)",
    lnpsa_dx       = "PSA at diagnosis (ng/ml, log)",
    clin_t         = "Stage (T)",
    clin_n         = "Stage (N)",
    clin_m         = "Stage (M)",
    stage_detailed = "Stage",
    stage          = "Stage",
    clin_tstage    = "Stage (T)",
    clin_nstage    = "Stage (N)",
    mstage         = "Stage (M)",
    rxprim         = "Primary treatment",
    rxprim_oth     = "Primary treatment (freetext)",
    rxprim_rp      = "Primary therapy: Prostatectomy",
    rxprim_adt     = "Primary therapy: Androgen deprivation",
    rxprim_chemo   = "Primary therapy: Chemotherapy",
    rxprim_xrt     = "Primary therapy: Radiation",
    rxprim_other   = "Primary therapy: Other",
    rp_gl_sum      = "Prostatectomy Gleason score",
    rp_gl34        = "Prostatectomy Gleason grade",
    rp_gl_maj      = "Prostatectomy Gleason grade, major pattern",
    rp_gl_min      = "Prostatectomy Gleason grade, minor pattern",
    path_t         = "Prostatectomy stage (T)",
    path_n         = "Prostatectomy stage (N)",
    is_crpc        = "Castration resistant",
    crpc_event     = "Castration resistant",
    is_met         = "Metastatic",
    met_event      = "Metastatic",
    is_mfs         = "Metastasis or death",
    mfs_event      = "Metastasis or death",
    is_dead        = "Death",
    death_event    = "Death",
    dx_bx_mos      = "Diagnosis to biopsy (months)",
    dx_adt_mos     = "Diagnosis to ADT (months)",
    adt_crpc_mos   = "ADT to castration resistance (months)",
    dx_crpc_mos    = "Diagnosis to castration resistance (months)",
    dx_met_mos     = "Diagnosis to metastasis (months)",
    dx_os_mos      = "Overall survival from diagnosis (months)",
    dx_mfs_mos     = "Metastasis-free survival from diagnosis (months)",
    adt_os_mos     = "Overall survival from ADT initiation (months)",
    crpc_os_mos    = "Overall survival from castration resistance (months)",
    met_os_mos     = "Overall survival from metastasis (months)") %>%
    # reorder variables
    select(ptid, complete_pts, age_dx, race, race4, race3, ethnicity,
           smoking, smoke01, bx_gl_sum, bx_gl, bx_gl34, bx_gl_maj, bx_gl_min,
           psa_dx, psa_dxcat, lnpsa_dx, clin_t, clin_n, clin_m, stage_detailed,
           stage, clin_tstage, clin_nstage, mstage, rxprim, rxprim_oth,
           rxprim_rp, rxprim_adt, rxprim_chemo, rxprim_xrt, rxprim_other,
           rp_gl_sum, rp_gl34, rp_gl_maj, rp_gl_min, path_t, path_n, is_crpc,
           crpc_event, is_met, met_event, is_dead, death_event, is_mfs,
           mfs_event, dx_bx_mos, dx_adt_mos, adt_crpc_mos, dx_crpc_mos,
           dx_met_mos, dx_os_mos, dx_mfs_mos, adt_os_mos, crpc_os_mos,
           met_os_mos,
           everything())

  smp <- smp %>% labelled::set_variable_labels(
    ptid         = "Patient ID",
    complete_smp = "Record completely entered",
    dmpid        = "Sample/molecular pathology ID",
    hist_smp     = "Histology of sample",
    hist_cmt     = "Histology of sample (freetext)",
    dzextent_smp = "Extent of disease at sample",
    dzextent_seq = "Extent of disease at sequencing",
    primmet_smp  = "Extent of disease at sample",
    ext_pros     = "Extent of disease: Prostate",
    ext_lndis    = "Extent of disease: Distant lymph node",
    ext_bone     = "Extent of disease: Bone",
    ext_vis      = "Extent of disease: Visceral/soft tissue",
    ext_liver    = "Extent of disease: Liver",
    ext_lung     = "Extent of disease: Lung",
    ext_other    = "Extent of disease: Other soft tissue",
    bonevol      = "Volume of bone metastases",
    cntadt       = "Sample on continuous ADT",
    tissue       = "Sample tissue",
    smp_tissue   = "Sample tissue",
    smp_pros     = "Prostate sample",
    pur_rev      = "Reviewed for purity",
    pur_remov    = "Removed for purity",
    age_smp      = "Age at sample (years)",
    age_seq      = "Age at sequencing (years)",
    dx_smp_mos   = "Diagnosis to sample (months)",
    adt_smp_mos  = "ADT to sample (months)",
    dx_seq_mos   = "Diagnosis to sequencing (months)",
    adt_seq_mos  = "ADT to sequencing (months)",
    met_seq_mos    = "Metastasis to sequencing (months)",
    smp_met_mos  = "Sample to metastases (months)",  # for QC
    smp_os_mos   = "Overall survival from sample (months)",  # for QC
    seq_met_mos  = "Sequencing to metastases (months)",
    seq_crpc_mos = "Sequencing to castration resistance (months)",
    seq_os_mos   = "Overall survival from sequencing (months)",
    seq_mfs_mos  = "Metastasis-free survival from sequencing (months)",
    dzvol        = "Volume of disease",
    denovom_smp  = "Timing of metastases",
    denovom_seq  = "Timing of metastases")

  trt <- trt %>% labelled::set_variable_labels(
    rx_line              = "Line of therapy",
    rx_name              = "Treatment Name",
    rx_name_other        = "Treatment Name (Other)",
    rx_name_parpi        = "PARPi Agent Name",
    rx_start             = "Treatment Start Date",
    rx_end               = "Treatment End Date/Last Known Treatment Date",
    dx_rx_start_mos      = "Time from diagnosis to treatment start (months)",
    dx_rx_end_mos        = "Time from diagnosis to treatment end (months)",
    rx_wks               = "Duration of treatment (weeks)",
    rx_censor            = "Treatment Ongoing",
    rx_stop_reason       = "Reason for Treatment Stop",
    rx_stop_reason_other = "Reason for Treatment Stop (Other)",
    rx_dzextent          = "Extent of Disease before starting PARPi",
    rx_ext_pros          = "Sites of Disease at PARPi start: Prostate",
    rx_ext_lndis         = "Sites of Disease at PARPi start: Distant lymph node",
    rx_ext_bone          = "Sites of Disease at PARPi start: Bone",
    rx_ext_vis           = "Sites of Disease at PARPi start: Visceral/soft tissue",
    rx_ext_liver         = "Sites of Disease at PARPi start: Liver",
    rx_ext_lung          = "Sites of Disease at PARPi start: Lung",
    rx_ext_other         = "Sites of Disease at PARPi start: Other soft tissue",
    rx_bonevol           = "Volume of Bone Metastases at PARPi start")

  pts_smp <- lst(pts, smp, trt)

  if(deidentify == TRUE)  # Default: deidentify the dataset
    pts_smp <- deidentify_prostate_redcap(pts_smp)

  pts_smp
}
