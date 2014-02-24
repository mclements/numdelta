\name{predictnl}
\alias{predictnl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of standard errors using the numerical delta method.
}
\description{
A simple, yet exceedingly useful, approach to estimate the variance of a
function using the numerical delta method. A number of packages provide
functions that analytically calculate the gradients; we use numerical
derivatives, which generalises to models that do not offer analytical
derivatives (e.g. ordinary differential equations, integration), or to
examples that are tedious or error-prone to calculate (e.g. sums of
predictions from GLMs).
}
\usage{
\method{predictnl}{default}(object, fun, newdata=NULL, ...)
\method{predictnl}{lm}(object, fun, newdata=NULL, ...)
\S4method{predictnl}{mle}(object, fun, ...)
\S4method{predictnl}{mle2}(object, fun, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object with \code{coef}, \code{vcov} and \code{`coef<-`}
methods (required). 
}
  \item{fun}{
A function that takes \code{object} as the first argument, possibly with
\code{newdata} and other arguments (required). See notes for why it is
often useful to include \code{newdata} as an argument to the function.
}
  \item{newdata}{
An optional argument that defines newdata to be passed to \code{fun}.
}
  \item{\dots}{
Other arguments that are passed to \code{fun}.
}
}
\details{
  
The signature for \code{fun}
is either \code{fun(object, ...)} or \code{fun(object, newdata=NULL,
  ...)}.

The different \code{predictnl} methods call the utility function
\code{numDeltaMethod}, which in turn calls the \code{grad} function for
numerical differentiation. The \code{numDeltaMethod} function calls the
standard \code{coef} and \code{vcov} methods, and the non-standard
\code{`coef<-`} method for changing the coefficients in a regression
object. This non-standard method has been provided for several
regression objects and essentially mirrors the \code{coef} method.

One potential issue is that some \code{predict} methods do not
re-calculate their predictions for the fitted dataset (i.e. when
\code{newdata=NULL}). As the \code{predictnl} function changes the
fitted coefficients, it is required that the predictions are
re-calculated. One solution is to pass \code{newdata} as an argument to
both \code{predictnl} and \code{fun}; alternatively, \code{newdata} can
be specified in \code{fun} or a specialised \code{numdelta:::predict}
method can be called. These approaches are described in the examples
below. The \code{numDeltaMethod} method called by \code{predictnl}
provides a warning when the variance estimates are zero, which may be
due to this cause.

For completeness, it is worth discussing why the example
\code{predictnl(fit,predict)} does not work for when \code{fit} is a
\code{glm} object. First, \code{predict.glm} does not update the
predictions for the fitted data. Second, the default \code{predict}
method has a signature \code{predict(object, ...)}, which does not
include a \code{newdata} argument. We could then either (i) require that
a \code{newdata} argument be passed to the \code{fun} function for all
examples, which would make this corner case work, or (ii) only pass the
\code{newdata} argument if it is non-null or in the formals for the
\code{fun} function, which would fail for this corner case. The current
API defaults to the latter case (ii).  To support this approach, the
\code{predictnl.lm} method replaces a null \code{newdata} with
\code{object$data}. We also provide a revised
\code{numdelta:::predict.lm} method that performs the same operation,
although its use is not encouraged due to its clumsiness.

}
\value{ Returns an object of class
\code{predictnl} with methods \code{print} and \code{confint}.
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
df <- data.frame(x=0:1, y=c(10, 20))
fit <- glm(y ~ x, df, family=poisson)

## Some of the prediction methods do not recalculate the predictions
## if the predictions are within the fitted data (that is, newdata=NULL). 
predictnl(fit,
          function(obj)
          diff(predict(obj,type="response"))) ## WRONG - with warning
## We can fix this in three ways.
## First, we can pass a newdata argument
predictnl(fit, function(obj,newdata)
          diff(predict(obj,newdata,type="response")))
## Second, we can specify the newdata argument in the prediction call
predictnl(fit,
          function(obj)
          diff(predict(obj,newdata=obj$data,type="response")))
## Third, we can use a specialised (unexported) predict method from numdelta
predictnl(fit,
          function(obj)
          diff(numdelta:::predict.lm(obj,type="response")))

## Simpler example
predict(fit,se.fit=TRUE)                # what we should get
predictnl(fit,predict)                  # WRONG - with warning 
predictnl(fit,predict.glm)              # okay 
predictnl(fit,predict,newdata=fit$data) # okay
predictnl(fit,numdelta:::predict.lm)    # okay

## A more complex example
\dontrun{
## log-hazard survival with natural splines
require(bbmle)
require(splines)
require(survival)
set.seed(12345)
y <- rweibull(100,2,10)
loghX <- ns(y,df=5,intercept=TRUE)
h <- function(coef,x) exp(drop(predict(loghX,x) \%*\% coef))
S <- function(coef,time)
     sapply(time, function(yi)
            exp(-integrate(function(x) h(coef,x),0,yi)$value))
negloglike <- function(coef) -sum(log(S(coef,y)) + log(h(coef,y)))
init <- structure(c(-1.4696, -2.1932, 1.2183, -5.1613, 0.607806),
            .Names = c("beta1", "beta2", "beta3", "beta4", "beta5"))
parnames(negloglike) <- names(init)
fit <- mle2(negloglike,init,vecpar=TRUE)

x <- seq(0,max(y),length=301)
pred1 <- predictnl(fit,function(obj) log(-log(S(coef(obj),x))))
plot(x,S(coef(fit),x),type="l")
lines(survfit(Surv(y,rep(TRUE,length(y)))~1),lty=2)
matlines(x,exp(-exp(confint(pred1))),col="blue",lty=1)
lines(x,pweibull(x,2,10,lower.tail=FALSE),col="red")

pred2 <- predictnl(fit,function(obj) log(h(coef(obj),x)))
plot(x,h(coef(fit),x),type="l")
matlines(x,exp(confint(pred2)),col="blue",lty=1)
lines(x,dweibull(x,2,10)/pweibull(x,2,10,lower.tail=FALSE),col="red")
}
 
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
