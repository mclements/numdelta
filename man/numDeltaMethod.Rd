\name{numDeltaMethod}
\alias{numDeltaMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility method to calculate variance using the numerical delta method.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
numDeltaMethod(object, fun, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{fun}{
%%     ~~Describe \code{fun} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, fun, ...) 
{
    coef <- coef(object)
    Sigma <- vcov(object)
    fit <- fun(coef, ...)
    gd <- grad(fun, coef, ...)
    se.fit <- as.vector(sqrt(diag(t(gd) \%*\% Sigma \%*\% gd)))
    names(se.fit) <- names(fit)
    if (all(se.fit == 0)) 
        warning("Zero variance estimated. Do you need to pass a newdata statement?")
    structure(list(fit = fit, se.fit = se.fit), class = "predictnl")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
