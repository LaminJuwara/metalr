\name{forest_metalr}
\alias{forest_metalr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forest plot for likeliihood ratio meta-analysis. Plots 95\% CIs and 95\% Intrinsic CIs.
%%  ~~function to do ... ~~
}
\description{
The function plots confidence limits of traditional 95\% CIs and ICIs for the studies included in the meta-analysis as well as confidence bars associated with the combined effect estimates.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
forest_metalr(metalr_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{metalr_obj}{
  An abject from the metalr functions metalr_or() or metalr_rr(). The metalr object is a list of results computed by the metalr functions which includes a dataframe of mle of the effect estimates and their corresponding 95\% CIs and ICIs. See the example below.
%%     ~~Describe \code{metalr_obj} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
Returns a forest plot of the 95\% CIs and 95\% ICIs.
%% ...
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

# the metalr object
metalr_obj<-metalr_or(idata=statindata[,2:5],refval=0,num_iter=3000,increm=0.001,method = "random")

# forest plot of the metalr object
forest_metalr(metalr_obj)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Likelihood ratio }
\keyword{ Intrinsic confidence interval}
\keyword{ Meta-analytic methods}
\keyword{ Forest plot}% __ONLY ONE__ keyword per line
