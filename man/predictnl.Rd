\name{predictnl}
\alias{predictnl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of standard errors using the numerical delta method.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predictnl(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
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
##
require(numdelta)
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
require(bbmle)
set.seed(12345)
y <- rweibull(100,2,10)
loghX <- ns(y,df=5,intercept=TRUE)
h <- function(coef,x) exp(drop(predict(loghX,x) %*% coef))
S <- function(coef,time) sapply(time, function(yi)
    exp(-integrate(function(x) h(coef,x),0,yi)$value))
negloglike <- function(coef) -sum(log(S(coef,y)) + log(h(coef,y)))
init <- structure(c(-1.4696, -2.1932, 1.2183, 
-5.1613, 0.607806), .Names = c("beta1", "beta2", 
"beta3", "beta4", "beta5"))
parnames(negloglike) <-names(init)
fit <- mle2(negloglike,init,vecpar=TRUE)

x <- seq(0,max(y),length=301)
pred1 <- predictnl(fit,function(obj) log(-log(S(coef(obj),x))))
plot(x,S(coef(fit),x),type="l")
lines(survfit(Surv(y,rep(T,length(y)))~1),lty=2)
matlines(x,exp(-exp(confint(pred1))),col="blue",lty=1)
lines(x,pweibull(x,2,10,lower.tail=FALSE),col="red")

pred2 <- predictnl(fit,function(obj) log(h(coef(obj),x)))
plot(x,h(coef(fit),x),type="l")
matlines(x,exp(confint(pred2)),col="blue",lty=1)
lines(x,dweibull(x,2,10)/pweibull(x,2,10,lower.tail=FALSE),col="red")



 
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
