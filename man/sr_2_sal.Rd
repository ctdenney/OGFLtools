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
  \item{sallim}{Limit of ability to accurately distinguish different salinities.}
  \item{fill}{what to replace calculated values that exceed the sallim value (e.g. "NA" or 7).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
    \item{sallim}{The user designated salinity resolution limit. Calculated sainities above this value will be replaced by the fill argument}
    \item{fill }{Fill argument excepts either "NA" or "sallim". If NA, values above the sallim will be replaced with NA. If "sallim", valuese above the sallim with be replaced with the sallim. }
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
