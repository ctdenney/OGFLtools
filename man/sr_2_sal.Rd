\name{sr_2_sal}
\alias{sr_2_sal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert measured Strontium ratio to salinity

}
\description{
This function calculates salinity from measured Sr using a 2 end member mixing model. It requires the measured Sr ratio, Endmember 1 Strontium, Endmember 2 Strontium, concentrations for each, and salinity of the two endmembers. Default values only correct for Sacramento Delta
}
\usage{
sr_2_sal(sr, srfw = 0.705781, srmar = 0.70918, confw = 103.8, conmar = 7900, salfw = 0.1, salmar = 31.8)

}

\arguments{
  \item{sr}{Measured strontium. Should be giving as the ratio of Sr87/Sr86.}
  \item{srfw}{Strontium ratio of freshwater (or endmember 1).}
  \item{srmar}{Strontium ratio of the ocean (or endmember 2).}
  \item{confw}{Concentation (parts per million) of Strontium in freshwater (or endmember 1).}
  \item{conmar}{Concentration (parts per million) of Strontium in the ocean (or endmember 2).}
  \item{salfw}{Salinity of freshwater (or endmember 1).}
  \item{salmar}{Salinity of the ocean (or endmember 2).}
  \item{sallim_high}{Maximum allowed salinity value.}
  \item{high_fill}{what to replace calculated values that exceed the sallim_high value (e.g. "NA" or 7). Must be either "NA" (as a string) or numeric value.}
  \item{sallim_low}{Minimum allowed salinity value.}
  \item{lowfill}{what to replace calculated values that are lower than the sallim_low value (e.g. "NA" or 0.1). Must be either "NA" (as a string) or numeric value.}
  \item{suppress.warnings}{argument to suppress the warning that some values are outside the boundaries of the srfw and srmar values. Default is F and warnings will be displayed. Value of T will hide warnings. Note that, without the warnings, values outside the sallim values will STILL be replaced.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
    \item{}{Value returned will be the calculated salinity value unless that calculated value exceeds the bounds set by the user with sallim_high and sallim_low. In those cases, the calculated value will be replaced with the appropriate fill value.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}% use one of  RShowDoc("KEYWORDS")
\keyword{}% __ONLY ONE__ keyword per line
