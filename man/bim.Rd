\name{bim}
\alias{bim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculage growth trajectory using biological intercept model
%%  ~~function to do ... ~~
}
\description{
Uses a BIM to calculate the growth at age for the full life history of a fish. It requires the fork length of the fish at capture (fl), the length at hatch for the species (hl) and the otolith radius trajecotry (gt).
}
\usage{
bim(fl, hl, gt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fl}{Fork length of a fish as measured at time of capture.}
  \item{hl}{Length at hatch for the species of interest}
  \item{gt}{The full otolith radius trajectory for a given fish. At a minimum, the minimum radius, maximum radius, and radius of interest must be supplied (for example, total otolith radius, radius at age 0, and radius for a given age can be supplied if you only want to know the size at a specific age).}
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

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}% use one of  RShowDoc("KEYWORDS")
\keyword{}% __ONLY ONE__ keyword per line
