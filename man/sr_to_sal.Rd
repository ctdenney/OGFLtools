\name{sr_to_sal}
\alias{sr_to_sal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert measured Strontium ratio to salinity

}
\description{
This function calculates salinity from measured Strontium isotope ratio (87Sr/86Sr) using a simple 2 end member mixing model. The function requires the measured sample Sr isotope ratio (87Sr/86Sr), and for each of the two endmembers the Sr isotope ratio (87Sr/86Sr), concentration, and salinity. Default values are based on measurements from the San Francisco Estuary. Values associated with Endmember 1 are sr1, con1, and sal1. Values associated with Endmember 2 are sr2, con2, and sal2. By convention (although not required), endmember 1 is the lower salinity endmember. The default endmember 1 values were collected from the Sacramento river at CDFW station 711 on 2012-03-08. The default endmember 2 values were collected at Muir Beach on 2007-04-23.
}
\usage{
sr_to_sal(sr, sr1 = 0.70576, sr2 = 0.70918, con1 = 95, con2 = 6819, sal1 = 0.11, sal2 = 31.8, sallim_high, sallim_low, fill_high, fill_low)

}

\arguments{
  \item{sr}{Measured strontium isotopic ratio. Must be given as the ratio of 87Sr/86Sr.}
  \item{sr1}{Strontium ratio of endmember 1 (by default, the Sacramento River value).}
  \item{sr2}{Strontium ratio of endmember 2 (by default, the Pacfic Ocean value).}
  \item{con1}{Concentation of Strontium in endmember 1. Default values are in parts per billion (ppb), but the function will accept and work with any units, as long as they are the same for con1 and con2}
  \item{con2}{Concentration of Strontium in the ocean (or endmember 2).Default values are in parts per billion (ppb), but the function will accept and work with any units, as long as they are the same for con1 and con2}
  \item{sal1}{Salinity of endmember 1 (by default, the Sacramento River value).}
  \item{sal2}{Salinity of the ocean (or endmember 2).}
  \item{sallim_high}{Maximum allowed salinity value. Default value is the higher salinity of the two endmembers}
  \item{fill_high}{what to replace calculated values that exceed the sallim_high value (e.g. "NA" or 7). Must be either "NA" (as a string) or numeric value. Default is "NA".}
  \item{sallim_low}{Minimum allowed salinity value.Default is the lower salinity of the two endmembers.}
  \item{fill_low}{what to replace calculated values that are lower than the sallim_low value (e.g. "NA" or 0.1). Must be either "NA" (as a string) or numeric value. Default is "NA".}
  \item{suppress.warnings}{argument to suppress the warning that some values are outside the boundaries of the sr1 and sr2 values. Default is F and warnings will be displayed. Value of T will hide warnings. Note that, without the warnings, values outside the sallim values will STILL be replaced.}
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
\item{}{Mixing model adapted from Phillis, Corey C., David J. Ostrach, B. Lynn Ingram, and Peter K. Weber. "Evaluating otolith Sr/Ca as a tool for reconstructing estuarine habitat use." Canadian Journal of Fisheries and Aquatic Sciences 68, no. 2 (2011): 360-373.}
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
