\name{sal_2_sr}
\alias{sal_2_sr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert salinity to Sr ratio
%%  ~~function to do ... ~~
}
\description{
Uses the Sr/Sal mixing model to take a measured salinity value and determine a Sr ratio
}
\usage{
sal_2_sr(sal, srfw = 0.705781, srmar = 0.70918, confw = 103.8, conmar = 7900, salfw = 0.1, salmar = 31.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sal}{Measured salinity. Should be giving as PSU.}
  \item{srfw}{Strontium ratio of freshwater (or endmember 1).}
  \item{srmar}{Strontium ratio of the ocean (or endmember 2).}
  \item{confw}{Concentation (parts per million) of Strontium in freshwater (or endmember 1).}
  \item{conmar}{Concentration (parts per million) of Strontium in the ocean (or endmember 2).}
  \item{salfw}{Salinity of freshwater (or endmember 1).}
  \item{salmar}{Salinity of the ocean (or endmember 2).}
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
