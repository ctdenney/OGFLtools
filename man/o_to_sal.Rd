\name{o_to_sal}
\alias{o_to_sal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert Oxygen isotope ratio to salinity
%%  ~~function to do ... ~~
}
\description{
This function uses published relationships of d18O and environmental salinity from either Ingram et al(1996) or McLaughlin et al(2006) to calculate salinity when given d18O. It is specific to the San Francisco Estuary.
}
\usage{
o_to_sal(d18O, source)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d18O}{is the measured d18O}
  \item{source}{is the desired source function. Currently only 'ingram' and 'mclg' are supported.}
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
Ingram BL, Conrad ME, Ingle JC. Stable isotope and salinity systematics
in estuarine waters and carbonates: San Francisco Bay. Geochim
Cosmochim Acta. 1996;60(3):455‐467. https://doi.org/10.1016/0016‐
7037(95)00398‐3
McLaughlin K, Kendall C, Silva SR, Young M, Paytan A. Phosphate
oxygen isotope ratios as a tracer for sources and cycling of
phosphate in North San Francisco Bay, California. J Geophys Res
Biogeo. 2006;111(3):1‐12. https://doi.org/10.1029/2005JG000079
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
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
