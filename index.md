# prostateredcap: Reproducible clinical-genomic research on prostate cancer

<!-- badges: start -->
<!-- badges: end -->

## Modules for reproducible clinical research on prostate cancer

1. **Data dictionary for REDCap:** Definitions for key data elements in clinical-genomic prostate cancer research. The field definitions can be [downloaded as a template](https://github.com/stopsack/prostateredcap/blob/main/inst/extdata/prostateredcap_datadictionary.csv) that serves for creating a [REDCap database](https://projectredcap.org/about/faq/) for data input. 

2. **Assessment of reproducibility:** Empirical assessment of the validity of key data fields for reproducibility in the specific clinical setting (as described in Keegan *et al.*).

3. **prostateredcap R package:** Load and reshape the raw data from REDCap, perform automated quality control checks, and deidentify, using an R package available via Github. The resulting sharable dataset is ready for statistical/bioinformatic analyses.


## Workflow

1. **Design and implement the REDCap database**  -- see Keegan *et al.*

   + Review rationale, create data definitions, assemble team.
   + Create a REDCap database using the [database template](https://github.com/stopsack/prostateredcap/blob/main/inst/extdata/prostateredcap_datadictionary.csv) and define access rights.
   + Pilot annotations, review pilot, and expand data extraction into REDCap. 
   
   
2. **Assessment of reproducibility** -- see Keegan *et al.*

   + Conducting a formal quality control study as described for implementation at Memorial Sloan Kettering Cancer Center is recommended.
   

3. **Data processing for analyses using the prostateredcap package** -- described in detail in the [Get Started](articles/prostateredcap.html) vignette.

   +  Export REDCap dataset in "Labeled CSV" format. (Because data collection involves protected health information, raw datasets cannot be shared. To test the process without having an actual dataset at hand yet, the [Get Started](articles/prostateredcap.html) vignette uses an [example dataset](https://github.com/stopsack/prostateredcap/blob/main/inst/extdata/SampleGUPIMPACTDatab_DATA_LABELS_2021-05-26.csv).)
   + In R, [install the prostateredcap package](articles/prostateredcap.html#install-the-prostateredcap-package-1) using `remotes::install_github("stopsack/prostateredcap")`.
   + Import the labelled CSV exported from REDCap using `load_prostate_redcap()`, which loads, merges, reformats, corrects, and labels the dataset (see [Details](reference/load_prostate_redcap.html#details)). By default, this function will also pass the dataset through `deidentify_prostate_redcap()` to remove protected health information.

       ``` r
       library(prostateredcap)
       datasets <- load_prostate_redcap(labeled_csv = "file_from_redcap.csv")
       ```
   
   + Run automated quality control checks and exclusions of records that fail these checks using `check_prostate_redcap()`. By using `recommended_only = TRUE`, only data elements recommended for analyses will be returned.

       ``` r
       datasets <- check_prostate_redcap(datasets, recommended_only = TRUE)
     
       # View quality control results and exclusion criteria:
       datasets$qc_pts
       datasets$qc_smp
       ```
   
  + Use the datasets for statistical and bioinformatic analyses: `datasets$pts` is the set of patient-level data; `datasets$smp` is the set of sample-level data; `datasets$trt` is the set of treatment data per line of treatment. See [data dictionary of elements recommended for analyses](articles/dataelements.html).


### References

Keegan NM, Vasselman SE, Barnett ES, Nweji B, Carbone EA, Blum A, Morris MJ, Rathkopf DE, Slovin SF, Danila DC, Autio KA, Scher HI, Kantoff PW, Abida W,\* Stopsack KH.\* Clinical annotations for prostate cancer research: Defining data elements, creating a reproducible analytical pipeline, and assessing data quality. *The Prostate*. 2022. doi: 10.1002/pros.24363. [Article](https://doi.org/10.1002/pros.24363) | [Preprint](https://doi.org/10.1101/2021.09.20.21263842)
