\name{ec_2_sal}
\alias{ec_2_sal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert electrical conductivity to salinity (PSU)
%%  ~~function to do ... ~~
}
\description{
This function converts electrical conductivity to salinity given temperature and the measured conductivity
}
\usage{
ec_2_sal(temp, cond)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{temp}{Temperature (degrees Celsius, ?C)}
  \item{cond}{Conductivitiy (microsiemens, uS)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{temp }{Temperature should be either the temperature at which the reading was taken in celsius, in the case of electrical conductivity, or 25C, in the case of specific conductance (which is electrical conductivity normalized to a a temperature readonf of 25C). It is extremely important that the user knows which their values are. The units are identical for both, but can result in wildly different calculated salinities if one is mistaken for the other. The SLS, 20mm and SKT surveys from CDFW have been specific conductance for the entirety of their run, EVEN WHEN THE OLD DATABASES REFFERRED TO ELECTRICAL CONDUCTIVITY. Bay study was electical conductivity from 1981-1990 (inclusive), but specific conductance after that period. Our own database is unfortunately a bit of both, but we have recorded calculated salinity directly, so this should always be used instead of conductance. This might change in the future.}
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

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
