\name{apply_corrections}
\alias{apply_corrections}
\title{Apply declarative corrections to sonde data}
\description{
Applies each row of a corrections CSV to a wide sonde data frame: for rows of
\code{unit_id} whose timestamp falls in [start_datetime, end_datetime], set
\code{field} to \code{new_value} (optionally only where the cell currently
equals \code{original_value}). A \code{field} of "datetime" moves the row's
timestamp instead; those edits run first. The source data are never modified
on disk.
}
\usage{
apply_corrections(df, corrections_path = "corrections/data_corrections.csv",
                  date_col = "date", time_col = "time", unit_col = "unit_id",
                  origins = NULL, verbose = TRUE)
}
\arguments{
  \item{df}{wide sonde data frame with date (Date), time (chron times) and unit columns}
  \item{corrections_path}{path to the corrections CSV (columns unit_id, field,
    start_datetime, end_datetime, original_value, new_value, reason, author,
    date_logged, origin)}
  \item{date_col, time_col, unit_col}{names of the date, time and unit columns in \code{df}}
  \item{origins}{NULL to apply every row, or a vector of origins (e.g. "sentinel") to apply only those}
  \item{verbose}{print a line per manual correction applied}
}
\value{\code{df} with the corrections applied.}
\seealso{\code{\link{load_sonde_master}}}
