# prostateredcap 0.3.2

* `load_prostate_redcap()`: {dplyr} 1.1.0 input check compatibility.
* Return types of visceral sites even if `recommended_only = TRUE`.
* Define time/event variables for metastasis-free survival.
* Add metastasis as time origin, in addition to sequencing, diagnosis, and ADT. 
  Document all available time origins for survival-type analyses in vignette. 
  Add argument `time_origin` to `load_prostate_redcap()` for easier selection. 


# prostateredcap 0.3.1

* `load_prostate_redcap()`: Fix name repair issues due to tidyverse changes
* Update manuscript reference.


# prostateredcap 0.3.0

* Add loading and deidentification of treatment data (`trt`).
* `deidentify_prostate_redcap()` Allow removal of ID from additional datasets
  for merging (`ext`).


# prostateredcap 0.2.1

* Add introductory vignette and (fake) example dataset.
* Export additional sample (`smp`) data elements; fix tissue type variable, and 
  add sample to `keep_also` list.


# prostateredcap 0.2.0

* Make use of newly available sequencing date and define entry time for survival
  analyses on that date; define time-varying variables separately at sample and
  at sequencing.
* Move `deidentify_prostate_redcap()` into separate function.
* Expand documentation, add pkgdown site.


# prostateredcap 0.1.2

* Fix sample QC filter 7: `primmet == "Primary" & smp_met_mos < 0.5` only in 
  samples with `is_met == "Yes"`.
* Separate the patient QC filter for missing censor date.


# prostateredcap 0.1.1

* Fix definitions of `tissue` (sample type), `dzextent` (variant histologies).
* Fix sample QC filter 8, `smp_met_mos > 0.5` in samples that are metastatic 
  according `dzextent`, to keep nmCRPC samples.
* Reorder `smp_type`, `dzextent`.


# prostateredcap 0.1.0

* First release.
