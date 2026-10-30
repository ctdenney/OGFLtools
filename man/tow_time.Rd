\name{parse_tow_time}
\alias{parse_tow_time}
\alias{format_tow_time}
\title{Parse and format datasheet tow times}
\description{
Field datasheet times are written as "HHMM" (e.g. "0915" or "915") or "HH:MM",
sometimes with stray separators or whitespace. \code{parse_tow_time} converts
them to minutes since midnight; \code{format_tow_time} converts them to a
canonical "HH:MM" string.
}
\usage{
parse_tow_time(x)
format_tow_time(x)
}
\arguments{
  \item{x}{a character (or numeric) vector of times}
}
\value{
\code{parse_tow_time}: an integer vector of minutes since midnight.
\code{format_tow_time}: a character vector of "HH:MM" strings.
Blank, malformed, or out-of-range (later than 23:59) values return NA.
}
\examples{
parse_tow_time(c("0915", "915", "09:15", "2460", NA))
format_tow_time(c("0915", "915", " 9:15", "abc"))
}
