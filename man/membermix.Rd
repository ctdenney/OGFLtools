\name{membermix}
\alias{membermix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the Sr 87/86 ratio of a mixture.
%%  ~~function to do ... ~~
}
\description{
This function takes an arbitrary number of contributing endmembers, their respective Sr isotope ratios (87Sr/86Sr), their total Sr concentrations, their salinities, and their contribution to the mixutre (proportion) and calculates the ratio and concentration of the mixture
}
\usage{
membermix(sr, conc, sal, mix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sr}{a numeric vector of the Sr isotope ratios (87Sr/86Sr) of the endmembers}
  \item{conc}{a numeric vector of the Sr concentrations of the endmembers}
  \item{sal}{a numeric vector of the salinties of the endmembers}
  \item{mix}{a numeric vector of the relative contributions of the endmembers (needs to sum to 1)}
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
