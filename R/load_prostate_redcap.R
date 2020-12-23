#' Guess Dates
#'
#' @description Uses \code{\link[lubridate]{parse_date_time}}
#'   to transform string vectors with different date formats into
#'   an R date vector, pritedn in ISO format.
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
                                !stringr::str_detect(string = messydate, pattern = "/") &
                                !stringr::str_detect(messydate, pattern = "-"),
                              true = paste0("06/30", messydate),
                              false = messydate)
  lubridate::as_date(lubridate::parse_date_time(messydate, orders = c("mdy", "my", "y")))
}

#' Recode "N/A", "Unknown", etc. to NA
#'
#' @description Replaces \code{c("Unknown / Not Reported", "N/A", "NA", "Unknown", "X", "x")}
#'  with \code{NA}.
#'
#' @param x String or factor vector
#'
#' @return Original vector
#' @noRd
unknowns <- function(x) {
  ifelse(x %in% c("Unknown / Not Reported", "N/A", "NA", "Unknown", "X", "x"),
         NA, x)
}

#' Load MSK-IMPACT Prostate REDCap Labelled CSV File
#'
#' @description Loads, merges, reformats, corrects, and labels
#'   the REDCap file used for the MSK-IMPACT Prostate clinical database.
#'   It is recommended that the returned list is next processed by
#'   \code{\link{check_prostate_redcap}}.
#'
#' @param labelled_csv CSV file with labels, exported from REDCap.
#'   Must be the labelled version and must contain dates
#'   in order to derive time intervals.
#' @param deidentify De-identify the returned data set?
#'   Drops all dates (keeping time intervals), rounds age and PSA,
#'   and replaces patient ID by an index number. Defaults to
#'   \code{TRUE} and should rarely, if ever, be turned off.
#'
#' @details The following edits and assumptions are made:
#' 1. Potentially incomplete date variables are converted to
#'    date vectors, using \code{\link{guessdate}}.
#' 2. Various missingness indicators in strings and factors,
#'    \code{c("Unknown / Not Reported", "N/A", "NA", "Unknown", "X", "x")},
#'    are converted to \code{NA}.
#' 3. "Undetectable" PSA is set to 0, PSA \code{">x"} is set to \code{x + 1}.
#' 4. Clinical T and N stage variables are set to missing if M1.
#' 5. Event dates for metastases (\code{met_date}) and castration resistance
#'    (\code{crpc_date}) are set:
#'    *  To the last clinic visit (\code{lastvisit})
#'       if the event has not occurred.
#'    *  If stage is M1 and the recorded metastasis date is no more than
#'       1 month discrepant, \code{met_date} is set to the diagnosis
#'       date (\code{dxdate}).
#'    *  If the sample is a variant histology (e.g., neuroendocrine),
#'       the castration resistance date (\code{crpc_date}) is the date of
#'       diagnosis and the event indicator for survival analyses
#'       (\code{event_crpc}) is \code{NA}.
#' 6. Disease extent, distinguishing CRPC from castration-sensitive disease,
#'    at sampling is based on the sample date and the date of the sample.
#'    If the samples was obtained before the CRPC date, or CRPC did not occur,
#'    the the sample is from castration-sensitive disease by definition.
#'
#' @return List of two labelled tibbles (data frames):
#'
#' * \code{pts}: Patient-level data
#'
#' * \code{smp}: Sample-level data
#'
#' Access variables labels in RStudio via \code{\link[utils]{View}}
#' or using \code{\link{attr}(., "label")}.
#'
#' @import dplyr stringr purrr forcats
#' @export
#'
#' @examples
#' \dontrun{
#' # Load data:
#' pts_smp <- load_prostate_redcap(
#'   labelled_csv = "GUPIMPACTDatabaseFre_DATA_LABELS_2020-01-01_0001.csv")
#'
#' # Access patient-level data:
#' pts_smp$pts
#'
#' # Access sample-level data:
#' pts_smp$smp
#'
#' # Pass 'pts_smp' to check_prostate_redcap() next
#'
#' }
load_prostate_redcap <- function(labelled_csv, deidentify = TRUE) {
  # Read REDCap file
  rcclin <- readr::read_csv(file = labelled_csv, guess_max = 50000) %>%
    filter(is.na(`Repeat Instrument`) | `Repeat Instrument` != "Treatment Data") %>%
    select(`Record ID`:`Complete?_2`) %>%
    rename(ptid = `Record ID`)

  # Patient-level baseline form
  pts_baseline <- rcclin %>%
    filter(is.na(`Repeat Instrument`)) %>%
    select(ptid       = ptid,
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
           psadx      = `PSA at Diagnosis`,
           clin_t     = `Clinical T Stage at Diagnosis`,
           clin_n     = `Clinical N Stage at Diagnosis (regional lymph node metastases)`,
           clin_m     = `Clinical M Stage at Diagnosis`,
           rxprim     = `Primary Therapy`,
           rxprim_oth = `Other Primary Therapy`,
           rp_gl_sum  = `Sum Gleason at Prostatectomy`,
           rp_gl_maj  = `Primary Gleason Pattern at Prostatectomy`,
           rp_gl_min  = `Secondary Gleason Pattern at Prostatectomy`,
           path_t     = `Pathologic T Stage at Diagnosis`,
           path_n     = `Pathologic N Stage at Diagnosis`)

  # Patient-level "freeze" (outcome) form
  pts_freeze <- rcclin %>%
    filter(`Repeat Instrument` == "Freeze Data") %>%
    select(ptid       = ptid,
           freezedate = `Freeze Date`,
           adtstart   = `Continuous ADT Start Date`,
           is_crpc    = `Castration Resistant`,
           crpc_date  = `Castration Resistance Date`,
           is_met     = `Metastatic`,
           met_date   = `Date of Metastasis`,
           lastvisit  = `Last MD Visit Date`,
           is_dead    = `Survival Status`,
           lastfu     = `Date of Death/Last Contact`)

  # Combined patient-level data
  pts <- left_join(pts_baseline, pts_freeze, by = "ptid") %>%  # allow for missing freeze form
    mutate_at(.vars = vars(dob, dxdate, bxdate, freezedate, adtstart, crpc_date, met_date, lastvisit, lastfu),
              .funs = guessdate) %>%
    mutate_at(.vars = vars(race, ethnicity, smoking, is_crpc, is_met,
                           starts_with("bx_"), starts_with("rp_"), starts_with("clin_"),
                           starts_with("path_")),
              .funs = unknowns) %>%
    mutate_at(.vars = vars(starts_with("bx_"), starts_with("rp_")), as.numeric) %>%
    mutate_at(.vars = vars(race, ethnicity, smoking, rxprim, is_met, is_crpc, is_dead,
                           starts_with("clin_"), starts_with("path_")),
              .funs = as.factor) %>%
    # PSA
    mutate(
      psadx = case_when(
        str_sub(psadx, start = 1, end = 1) == ">" ~
          as.numeric(str_sub(psadx, start = 2, end = 10)) + 1,
        str_detect(string = psadx, pattern = "-") ~  # two numbers, e.g. "4.5-5.2"
          mean(as.numeric(str_split(string = psadx, pattern = "-", n = 2, simplify = TRUE)[1]),
               as.numeric(str_split(string = psadx, pattern = "-", n = 2, simplify = TRUE)[2])),
        str_detect(string = psadx, pattern = "etect") |
          str_detect(string = psadx, pattern = "dec") ~
          0,  # "undetectable", also misspelled "undectable"
        TRUE ~ as.numeric(psadx)),
      psadxcat = factor(cut(psadx, breaks = c(0, 4, 10, 20, Inf), include.lowest = TRUE),
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
        clin_t %in% c("1a", "1b", "1c", "2", "2a", "2b", "2c") & clin_n == "0" & clin_m == "0" ~ "T1/T2 N0 M0",
        clin_t %in% c("1a", "1b", "1c", "2", "2a", "2b", "2c") & is.na(clin_n) & clin_m == "0" ~ "T1/T2 NX M0",
        clin_t %in% c("3",  "3a", "3b")                        & clin_n == "0" & clin_m == "0" ~ "T3 N0 M0",
        clin_t %in% c("3",  "3a", "3b")                        & is.na(clin_n) & clin_m == "0" ~ "T3 NX M0",
        clin_t == "4"                                                          & clin_m == "0" ~ "T4 M0",
        is.na(clin_t)                                          & clin_n == "0" & clin_m == "0" ~ "TX N0 M0",
        is.na(clin_t)                                          & is.na(clin_n) & clin_m == "0" ~ "TX NX M0",
        clin_n == "1" & clin_m == "0"                                                          ~ "N1 M0",
        clin_m %in% c("1", "1a", "1b", "1c")                                                   ~ "M1"),
        levels = c("T1/T2 N0 M0", "T1/T2 NX M0", "T3 N0 M0", "T3 NX M0",
                   "T4 M0", "TX N0 M0", "TX NX M0", "N1 M0", "M1")),
      stage = fct_other(stage_detailed, keep = c("N1 M0", "M1"), other_level = "N0/NX M0"),
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
      # if crpc_date or met_date are (erroneously) after lastfu, then update lastfu
      lastfu = if_else(lastfu > crpc_date | is.na(crpc_date), true = lastfu, false = crpc_date),
      lastfu = if_else(lastfu > met_date | is.na(met_date), true = lastfu, false = met_date),
      # Set missing CRPC and met dates to last visit date:
      crpc_date    = case_when(
        # Special case: Histologies that are castration resistant by definition
        # (e.g., neuroendocrine). Set CRPC date as the date of diagnosis.
        str_detect(string = is_crpc, pattern = "istology") |
          str_detect(string = is_crpc, pattern = "euroendo") ~ dxdate,
        is.na(crpc_date) ~ lastvisit,
        TRUE             ~ crpc_date),
      met_date = case_when(
        mstage == "M1" & abs(lubridate::interval(dxdate, met_date) / months(1)) <= 1 ~
          dxdate,  # Metastatic at diagnosis. Change only if <1 mo discrepancy between met_date and dxdate
        is.na(met_date) & is_met == "No"   ~ lastvisit,  # censor at last visit
        is.na(met_date) & is_met == "Yes"  ~ NA_real_,   # this case should not exist
        !is.na(met_date) & is_met == "Yes" ~ met_date),
      is_met   = factor(case_when(
        mstage == "M1" ~ "Yes",
        TRUE           ~ as.character(is_met))),
      # CRPC event indicator for survival analyses set to missing for variant histologies:
      crpc_event   = case_when(is_crpc == "Yes" ~ 1,
                               is_crpc == "No"  ~ 0),
      met_event    = case_when(is_met == "Yes"  ~ 1,
                               is_met == "No"   ~ 0),
      death_event  = if_else(is_dead == "Dead", 1, 0),
      # Time intervals
      agedx        = lubridate::interval(dob,      dxdate)    / lubridate::years(1),
      dx_bx_mos    = lubridate::interval(dxdate,   bxdate)    / months(1),
      dx_adt_mos   = lubridate::interval(dxdate,   adtstart)  / months(1),
      adt_crpc_mos = lubridate::interval(adtstart, crpc_date) / months(1),
      dx_crpc_mos  = lubridate::interval(dxdate,   crpc_date) / months(1),
      dx_met_mos   = lubridate::interval(dxdate,   met_date)  / months(1),
      dx_os_mos    = lubridate::interval(dxdate,   lastfu)    / months(1),
      adt_os_mos   = lubridate::interval(adtstart, lastfu)    / months(1),
      crpc_os_mos  = if_else(is_crpc == "Yes",
                             true = lubridate::interval(crpc_date, lastfu) / months(1),
                             false = NA_real_),
      # Combine race categories and restrict to 3 race categories
      race4        = fct_lump(f = race, n = 3),
      race3        = fct_recode(fct_recode(race4, NULL = "Other"),
                                "Black" = "Black or African American"),
      race3        = fct_relevel(race3, "Asian", "White", "Black"),
      lnpsa        = log(if_else(psadx == 0, true = 0.01, false = psadx)),
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
    tidyr::separate(col = "dmpid", into = paste0("dmpid_part", 1:4), sep = "-") %>%
    mutate(dmpid_part4 = if_else(dmpid_part4 %in% c("IM06", "iM6", "Im6"),
                                 true = "IM6", false = dmpid_part4)) %>%
    transmute(
      ptid      = ptid,
      complete_smp = complete_smp,
      dmpid     = paste(dmpid_part1, dmpid_part2, dmpid_part3, dmpid_part4, sep = "-"),
      smpdate   = guessdate(`Date of Collection`),
      hist_smp  = factor(unknowns(`Histology for Sample`)),
      hist_cmt  = `Histology (Other)`,
      dzextent  = factor(unknowns(`Extent of Disease at Collection`)),
      ext_pros  = factor(if_else(`Sites of Disease (choice=Prostate/Prostate Bed)` == "Checked",
                                 TRUE, FALSE)),
      ext_lndis = factor(if_else(`Sites of Disease (choice=LN (distant))`          == "Checked",
                                 TRUE, FALSE)),
      ext_bone  = factor(if_else(`Sites of Disease (choice=Bone)`                  == "Checked",
                                 TRUE, FALSE)),
      ext_vis   = factor(if_else(`Sites of Disease (choice=Liver)`                 == "Checked" |
                                   `Sites of Disease (choice=Lung)`                == "Checked" |
                                   `Sites of Disease (choice=Other Soft Tissue)`   == "Checked",
                                 TRUE, FALSE)),
      bonevol   = factor(unknowns(`Volume of Bone Metastases at Time of Collection`)),
      cntadt    = factor(unknowns(`Continuous ADT`)),
      tissue    = factor(unknowns(`Sample Type`)),
      smp_pros  = fct_other(tissue, keep = "Prostate", other_level = "Non-prostate"),
      smp_type  = fct_other(tissue, keep = c("Prostate", "Bone", "Lymph Node", "Liver", "Lung"),
                            other_level = "Other soft tissue"),
      smp_type  = fct_recode(tissue, `Lymph node` = "Lymph Node"),
      smp_type  = fct_collapse(smp_type, Visceral = c("Liver", "Lung")),
      smp_type  = fct_relevel(smp_type, "Prostate", "Lymph node", "Bone", "Visceral", "Other soft tissue"),
      # Unsure how these two variables are supposed to be used:
      pur_rev   = factor(`Reviewed for Tumor Purity`),
      pur_remov = factor(`Removed for Low Tumor Purity`)) %>%
    # Derive variables for which data from "pts" are needed:
    left_join(pts %>% select(ptid, stage, agedx, dxdate, met_date,
                             crpc_date, is_crpc, lastfu, adtstart),
              by = "ptid") %>%
    mutate(
      dzextent       = factor(case_when(
        dzextent %in% c("Localized", "Regional nodes") &
          is_crpc == "Yes" & crpc_date <= smpdate                   ~ "Non-metastatic castration-resistant",
        dzextent %in% c("Localized", "Regional nodes") &
          (is_crpc != "Yes" | crpc_date > smpdate)                  ~ as.character(dzextent),
        dzextent == "Metastatic" &
          grepl("N/A-", is_crpc, fixed = TRUE)                      ~ "Metastatic, variant histology",
        dzextent=="Metastatic" &
          (is_crpc == "No" | (is_crpc == "Yes" & crpc_date > smpdate))  ~ "Metastatic hormone-sensitive",
        dzextent == "Metastatic" &
          (is_crpc != "No" & crpc_date <= smpdate)                  ~ "Metastatic castration-resistant")),
      dzextent = fct_relevel(dzextent, "Localized", "Regional nodes", "Metastatic hormone-sensitive",
                             "Non-metastatic castration-resistant", "Metastatic castration-resistant",
                             "Metastatic, variant histology"),
      primmet = case_when(
        dzextent %in% c("Localized", "Regional nodes") ~ "Primary",
        dzextent %in% c("Metastatic castration-resistant",
                        "Metastatic hormone-sensitive",
                        "Metastatic, variant histology") ~ "Metastatic"),
      agesmp        = agedx + lubridate::interval(dxdate, smpdate) / lubridate::years(1),
      dx_smp_mos    = lubridate::interval(dxdate,   smpdate)   / months(1),
      adt_smp_mos   = lubridate::interval(adtstart, smpdate)   / months(1),
      smp_met_mos   = lubridate::interval(smpdate,  met_date)  / months(1),
      smp_crpc_mos  = lubridate::interval(smpdate,  crpc_date) / months(1),
      smp_os_mos    = lubridate::interval(smpdate,  lastfu)    / months(1),
      dzvol         = factor(case_when(
        str_detect(string = as.character(dzextent), pattern = "Metastatic") &
          (bonevol != "High-Volume Bone Metastases" | is.na(bonevol)) &
          ext_vis == FALSE                                             ~ "Low-volume disease",
        str_detect(string = as.character(dzextent), pattern = "Metastatic") &
          (bonevol == "High-Volume Bone Metastases" | ext_vis == TRUE) ~ "High-volume disease"),
        levels = c("Low-volume disease", "High-volume disease")),
      denovom      = factor(case_when(
        str_detect(string = as.character(dzextent), pattern = "Metastatic") & stage != "M1" ~
          "Metastatic recurrence",
        str_detect(string = as.character(dzextent), pattern = "Metastatic") & stage == "M1" ~
          "De-novo metastatic")),
      stage_for_qc = stage) %>%
    # remove variables from "pts":
    select(-stage, -agedx, -dxdate, -met_date, -crpc_date, -is_crpc, -lastfu, -adtstart)

  pts <- pts %>% labelled::set_variable_labels(
    ptid           = "Patient ID",
    complete_pts   = "Record completely entered",
    agedx          = "Age at diagnosis (years)",
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
    psadx          = "PSA at diagnosis (ng/ml)",
    psadxcat       = "PSA at diagnosis (ng/ml)",
    lnpsa          = "PSA at diagnosis (ng/ml, log)",
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
    is_dead        = "Death",
    death_event    = "Death",
    dx_bx_mos      = "Diagnosis to biopsy (months)",
    dx_adt_mos     = "Diagnosis to ADT (months)",
    adt_crpc_mos   = "ADT to castration resistance (months)",
    dx_crpc_mos    = "Diagnosis to castration resistance (months)",
    dx_met_mos     = "Diagnosis to metastasis (months)",
    dx_os_mos      = "Overall survival from diagnosis (months)",
    adt_os_mos     = "Overall survival from ADT initiation (months)",
    crpc_os_mos    = "Overall survival from castration resistance (months)") %>%
    # reorder variables
    select(ptid, complete_pts, agedx, race, race4, race3, ethnicity,
           smoking, smoke01, bx_gl_sum, bx_gl, bx_gl34, bx_gl_maj, bx_gl_min,
           psadx, psadxcat, lnpsa, clin_t, clin_n, clin_m, stage_detailed,
           stage, clin_tstage, clin_nstage, mstage, rxprim, rxprim_oth, rxprim_rp,
           rxprim_adt, rxprim_chemo, rxprim_xrt, rxprim_other, rp_gl_sum, rp_gl34,
           rp_gl_maj, rp_gl_min, path_t, path_n, is_crpc, crpc_event, is_met,
           met_event, is_dead, death_event, dx_bx_mos, dx_adt_mos, adt_crpc_mos,
           dx_crpc_mos, dx_met_mos, dx_os_mos, adt_os_mos, crpc_os_mos, everything())

  smp <- smp %>% labelled::set_variable_labels(
    ptid         = "Patient ID",
    complete_smp = "Record completely entered",
    dmpid        = "Sample/molecular pathology ID",
    hist_smp     = "Histology of sample",
    hist_cmt     = "Histology of sample (freetext)",
    dzextent     = "Extent of disease",
    primmet      = "Extent of disease",
    ext_pros     = "Extent of disease: Prostate",
    ext_lndis    = "Extent of disease: Distant lymph node",
    ext_bone     = "Extent of disease: Bone",
    ext_vis      = "Extent of disease: Visceral/soft tissue",
    bonevol      = "Volume of bone metastases",
    cntadt       = "Sample on continuous ADT",
    tissue       = "Sample tissue",
    smp_type     = "Sample tissue",
    smp_pros     = "Prostate sample",
    pur_rev      = "Reviewed for purity",
    pur_remov    = "Removed for purity",
    agesmp       = "Age at sample (years)",
    dx_smp_mos   = "Diagnosis to sample (months)",
    adt_smp_mos  = "ADT to sample (months)",
    smp_met_mos  = "Sample to metastases (months)",
    smp_crpc_mos = "Sample to castration resistance (months)",
    smp_os_mos   = "Overall survival from sample (months)",
    dzvol        = "Volume of disease",
    denovom      = "Timing of metastases")

  # Default: deidentify the dataset
  if(deidentify == TRUE) {
    pts <- pts %>%
      select(-dob, -dxdate, -met_date, -crpc_date, -lastfu, -adtstart, -bxdate,
             -freezedate, -lastvisit) %>%
      mutate(agedx = as.numeric(round(agedx, digits = 1)),
             psadx = as.numeric(round(psadx, digits = 1)),
             lnpsa = as.numeric(round(lnpsa, digits = 2)),
             ptid2 = row_number())
    smp <- smp %>%
      left_join(pts %>% select(ptid, ptid2), by = "ptid") %>%
      select(-ptid, -smpdate) %>%
      mutate(agesmp = as.numeric(round(agesmp, digits = 1))) %>%
      select(ptid = ptid2, everything()) %>%
      labelled::set_variable_labels(ptid   = "Patient ID",
                                    agesmp = "Age at sample (years)")
    pts <- pts %>%
      select(-ptid) %>%
      select(ptid = ptid2, everything()) %>%
      labelled::set_variable_labels(ptid  = "Patient ID",
                                    agedx = "Age at diagnosis (years)",
                                    psadx = "PSA at diagnosis (ng/ml)",
                                    lnpsa = "PSA at diagnosis (ng/ml, log)")
  }

  lst(pts, smp)
}
