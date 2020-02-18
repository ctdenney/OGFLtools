\name{ec_2_spc}
\alias{ec_2_spc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to convert electrical conductivity to specific conductance
%%  ~~function to do ... ~~
}
\description{
Electrical conductivity is the measured conductivity of a sample. Specific conductance is a standardized measure representing what that conductivity *would have been* had it been taken at 25C. The function is described in Standard methods for examination of water and wastewater in section 2510 B under number 5 (19th edition)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ec_2_spc(temp, ec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{temp}{
  temperature at which the sample was taken, degrees Celsius
%%     ~~Describe \code{x} here~~
}
  \item{ec}{
  conductivity of the sample, should be in units of uS/cm (micro Siemens per centimeter)
%%     ~~Describe \code{x} here~~
}
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
