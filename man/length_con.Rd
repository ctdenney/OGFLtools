\name{length_con}
\alias{length_con}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to convert between various fish length measurements
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function uses various fish length measurement relationships to convert from an available measurement into a desired measurement. For example, if only standard length after ethanol preservation (sle) was recorded, but fresh fork length measurement (flf) is desired, then flf can be calculated using the relationship between sle and flf from fish which had both measurements taken.
}
\usage{
length_con(from, to, length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{
%%     ~~Describe \code{from} here~~
The measured length type.
}
  \item{to}{
%%     ~~Describe \code{to} here~~
The desired calculated length type.
}
  \item{length}{
The length or lengths to be converted.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
    \item{from }{Must be one of flf, fle, slf, sle, and tlf are supported (tle, total length ethanol, is not a supported input measurement). If not one of these, will return an error message.}
    \item{to }{Must be one of flf, slf, and tlf are supported. If not, will return an error message.}
    \item{length}{A single value or vector of numeric values. Conversion functions were created using fish >40mm. If a length is smaller than that, a warning will be produced.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Arthur Barros
Christian Denney
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
