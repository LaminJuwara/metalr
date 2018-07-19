\name{ici.or}
\alias{ici.or}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating Intrinsic confidence interval (ici) for odds ratio.
%%  ~~function to do ... ~~
}
\description{
Calculates traditional and intrinsic confidence intervals for odds ratio from an observational study.
}
\usage{
ici.or(idata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{idata}{ A dataframe of four columns: cases for treatment A, controls for treatment A, cases for treatment B and control for treatment B.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{OR}{The MLE estimate of the odds ratio}
\item{ll.ci}{lower 95 percent confidence limits}
\item{ul.ci}{upper 95 percent confidence limits}
\item{ll.ici}{lower 95 percent intrinsic confidence limits}
\item{ul.ici}{upper 95 percent intrinsic confidence limits}
}
\references{
%% ~put references to the literature/web site here ~
Dormuth, Colin R., Kristian B. Filion, and Robert W. Platt. "Likelihood ratio meta-analysis: New motivation and approach for an old method." Contemporary clinical trials 47 (2016): 259-265.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

# statin potency and acute kidney injury data
data("statindata")
ici.or(idata = statindata[1,2:5]) # odds ratio intrinsic ci for study 1

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Likelihood ratio }
\keyword{ Intrinsic confidence interval}
\keyword{ Meta-analytic methods}
