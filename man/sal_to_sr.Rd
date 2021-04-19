\name{sal_to_sr}
\alias{sal_to_sr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert salinity to Sr isotope ratio
%%  ~~function to do ... ~~
}
\description{
Uses the Sr/Sal mixing model to take a measured salinity value of a sample and determine a Sr isotope ratio (87Sr/86Sr)
}
\usage{
sal_to_sr(sal, sr1 = 0.70576, sr2 = 0.70918, con1 = 95, con2 = 6819, sal1 = 0.11, sal2 = 31.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sal}{Measured salinity of the sample. Should be given as PSU.}
  \item{sr1}{Strontium isotope ratio (87Sr/86Sr) of freshwater (or endmember 1).}
  \item{sr2}{Strontium isotope ratio (87Sr/86Sr) of the ocean (or endmember 2).}
  \item{con1}{Concentation (parts per million) of Strontium in freshwater (or endmember 1).}
  \item{con2}{Concentration (parts per million) of Strontium in the ocean (or endmember 2).}
  \item{sal1}{Salinity of freshwater (or endmember 1).}
  \item{sal2}{Salinity of the ocean (or endmember 2).}
%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}% use one of  RShowDoc("KEYWORDS")
\keyword{}% __ONLY ONE__ keyword per line
