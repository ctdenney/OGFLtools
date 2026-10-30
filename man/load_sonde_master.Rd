\name{load_sonde_master}
\alias{load_sonde_master}
\title{Load the sonde master record with corrections applied}
\description{
Reads a sonde master CSV (the raw record as it came off the instruments),
converts \code{date} to Date and \code{time} to chron times, applies the
corrections file with \code{\link{apply_corrections}}, and adds a POSIXct
\code{datetime} column (UTC-labelled clock time) for joining to sonde times
written on field datasheets. Corrections run before \code{datetime} is built,
so retimed rows join on their corrected timestamp.
}
\usage{
load_sonde_master(master_path, corrections_path = NULL, verbose = TRUE)
}
\arguments{
  \item{master_path}{path to sonde_data_master_<date>.csv}
  \item{corrections_path}{path to data_corrections.csv; NULL applies no corrections}
  \item{verbose}{passed to \code{apply_corrections}}
}
\value{A data frame: the master record, corrected, with a \code{datetime} column.}
\seealso{\code{\link{apply_corrections}}}
