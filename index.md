# prostateredcap: Reproducible clinical-genomic research on prostate cancer

<!-- badges: start -->
<!-- badges: end -->

Three modules for reproducible clinical-genomic research on prostate cancer are available:

1. **Data dictionary for REDCap:** Definitions for key data elements in clinical-genomic prostate cancer research (~~described in this manuscript~~). The field definitions can be ~~[downloaded as a template](https://pendingurl)~~ that serves for creating a [REDCap database](https://projectredcap.org/about/faq/) for data input. 

2. **Assessment of reproducibility:** Empirical assessment of the validity of key data fields for reproducibility in the specific clinical setting (~~as described in this manuscript~~).

3. **prostateredcap R package:** Load and reshape the raw data from REDCap, perform automated quality control checks, and deidentify. The resulting sharable dataset is ready for statistical/bioinformatic analyses.


## Workflow

1. Review rationale and data definitions in the ~~accompanying manuscript~~. Assemble team.

2. Create a REDCap database using the ~~[database CSV template](https://pendingurl)~~ and define access rights.

3. Pilot annotations, review pilot, and expand data abstraction into REDCap. Conducting a formal quality control for reproducibility ~~as described for the MSK implementation~~ is recommended.

4. Export REDCap dataset in "Labeled CSV" format. (Because data collection involves protected health information, the raw dataset cannot be shared. To test the process without having an actual dataset at hand yet, use ~~this dummy dataset~~.)

5. In R, install (or update) the prostateredcap package from [GitHub](https://stopsack.github.io/prostateredcap). This step is only required once.

   ``` r
   install.packages("remotes")  # skip if 'remotes' package is already installed
   remotes::install_github("stopsack/prostateredcap")
   ```

6. Import the labelled CSV exported from REDCap using `load_prostate_redcap()`, which loads, merges, reformats, corrects, and labels the dataset (see [Details](reference/load_prostate_redcap.html#details)). By default, this function will also pass the dataset through `deidentify_prostate_redcap()` to remove protected health information.

   ``` r
   library(prostateredcap)
   datasets <- load_prostate_redcap(labeled_csv = "GUPIMPACTDatabaseFre_DATA_LABELS_2020-01-01_0001.csv")
   ```

7. Run automated quality control checks and exclusions of records that fail these checks. By using `recommended_only = TRUE`, only data elements recommended for analyses will be returned.

   ``` r
   datasets <- check_prostate_redcap(datasets, recommended_only = TRUE)
 
   # View quality control results and exclusion criteria:
   datasets$qc_pts
   datasets$qc_smp
   ```

8. Use the datasets for statistical and bioinformatic analyses: `datasets$pts` is the set of patient-level data; `datasets$smp` is the set of sample-level data. See [data dictionary of elements recommended for analyses](articles/dataelements.html).


