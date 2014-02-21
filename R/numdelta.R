## package.skeleton(name="Rstpm2", path="c:/usr/src/R", force=T, namespace=T, code_files="pm2-3.R")
## Rtools.bat
## R CMD INSTALL "c:/usr/src/R/numdelta"
## R CMD build "c:/usr/src/R/numdelta"
## R CMD build --binary "c:/usr/src/R/numdelta"

## numerically calculate the gradient (func may return a vector)
grad <- function(func,x,...) # would shadow numDeriv::grad()
  {
    h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
    temp <- x+h
    h.hi <- temp-x
    temp <- x-h
    h.lo <- x-temp
    twoeps <- h.hi+h.lo
    nx <- length(x)
    ny <- length(func(x,...))
    if (ny==0L) stop("Length of function equals 0")
    df <- if(ny==1L) rep(NA, nx) else matrix(NA, nrow=nx,ncol=ny)
    for (i in 1L:nx) {
      hi <- lo <- x
      hi[i] <- x[i] + h.hi[i]
      lo[i] <- x[i] - h.lo[i]
      if (ny==1L)
        df[i] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
      else df[i,] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
      }
    return(df)
  }
## fun: takes coef as its first argument
## requires: coef() and vcov() on the object
numDeltaMethod <- function(object,fun,...) {
  coef <- coef(object)
  est <- fun(coef,...)
  Sigma <- vcov(object)
  gd <- grad(fun,coef,...)
  se.est <- as.vector(sqrt(diag(t(gd) %*% Sigma %*% gd)))
  data.frame(Estimate = est, SE = se.est)
}
predictnl <- function (object, ...) 
  UseMethod("predictnl")
predictnl.default <- function(object,fun,newdata=NULL,...)
  {
    ## link=c(I,log,sqrt),invlink=NULL
    ## link <- match.arg(link)
    ## if (is.null(invlink))
    ## invlink <- switch(deparse(substitute(link)),I=I,log=exp,sqrt=function(x) x^2)
    if (is.null(newdata) && !is.null(object$data))
      newdata <- object$data
    localf <- function(coef,...)
      {
        object$coefficients = coef
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,...)
  }
## setMethod("predictnl", "mle", function(object,fun,...)
##   {
##     localf <- function(coef,...)
##       {
##         object@fullcoef = coef # changed from predictnl.default()
##         fun(object,...)
##       }
##     numDeltaMethod(object,localf,...)
##   })
predictnl.glm <- function(object,fun,newdata=NULL,...)
  {
    ## stopifnot(require(stats))
    if (is.null(newdata) && !is.null(object$data))
      newdata <- object$data
    localf <- function(coef,...)
      {
        object$coefficients = coef
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,...)
  }

if (FALSE) {
  try(detach("package:numdelta",unload=TRUE),silent=TRUE)
}
