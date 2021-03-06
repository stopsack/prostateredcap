% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deidentify_prostate_redcap.R
\name{deidentify_prostate_redcap}
\alias{deidentify_prostate_redcap}
\title{Deidentify Prostate REDCap Data}
\usage{
deidentify_prostate_redcap(data)
}
\arguments{
\item{data}{List with elements \code{pts} and \code{smp}, as generated
by \code{\link{load_prostate_redcap}}.}
}
\value{
List:
\itemize{
\item \code{pts}: Deidentified patient-level data.
\item \code{smp}: Deidentified sample-level data.
}
}
\description{
Deidentifies the dataset for analysis and sharing:
\itemize{
\item Drops all dates, keeping time intervals only
\item Rounds age (in years, 1 decimal) and PSA (1 decimal)
\item Replaces patient ID by a sequential index number that still allows for
merging \code{pts} and \code{smp} datasets.
}

By default, this function is already being called automatically within
\code{\link{load_prostate_redcap}}.
}
