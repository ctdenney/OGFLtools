\name{spc_2_ec}
\alias{spc_2_ec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to convert specific conductance to electrical conductivity.
%%  ~~function to do ... ~~
}
\description{
Specific conductance is a standardized measure as if conductivity had been taken at a water temperature of 25C. This function returns the actual conductivity value as meaured at the environtmental temperature. Formula taken from Standard Methods for the examination of water and wastewater.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spc_2_ec(temp, spc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{temp}{
  Temperature in degrees celsius of the sample
%%     ~~Describe \code{x} here~~
}
  \item{spc}{
  The specific conductance of the sample, in uS/cm
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns electrical conductivity in units of uS/cm
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Miner, Gary. "Standard methods for the examination of water and wastewater." American Water Works Association. Journal 98.1 (2006): 130.
%% ~put references to the literature/web site here ~
}
\author{
Christian Denney
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
