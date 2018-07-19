\name{ici.rr}
\alias{ici.rr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intrinsic confidence interval (ici) for rate ratio
%%  ~~function to do ... ~~
}
\description{
Calculates traditional and intrinsic confidence intervals for rate ratio from epidemiological studies.
}
\usage{
ici.rr(cases, patients, person_yrs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cases}{
  The number of outcome events among followed participants.
}
  \item{patients}{
The number of patients or participants followed.
}
  \item{person_yrs}{
%%     ~~Describe \code{person_yrs} here~~
The period for which the cases were followed.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{OR}{The MLE estimate of the rate ratio}
\item{ll.ci}{lower 95 percent confidence limits}
\item{ul.ci}{upper 95 percent confidence limits}
\item{ll.ici}{lower 95 percent intrinsic confidence limits}
\item{ul.ici}{upper 95 percent intrinsic confidence limits}
}
\references{
Dormuth, Colin R., Kristian B. Filion, and Robert W. Platt. "Likelihood ratio meta-analysis: New motivation and approach for an old method." Contemporary clinical trials 47 (2016): 259-265.
%% ~put references to the literature/web site here ~
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
## Dataset: eg. Clopidogrel vs Aspirin trial
cases<-c(939,1021)
person_yrs<-c(17636,17519)
patients<-c(9599,9586)

ici.rr(cases, patients, person_yrs)

}

\keyword{ Likelihood ratio }
\keyword{ Intrinsic confidence interval}
\keyword{ Meta-analytic methods}
