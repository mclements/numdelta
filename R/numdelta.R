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
  Sigma <- vcov(object)
  fit <- fun(coef,...)
  gd <- grad(fun,coef,...)
  se.fit <- as.vector(sqrt(colSums(gd* (Sigma %*% gd))))
  if (!is.null(names(fit)))
      names(se.fit) <- names(fit)
  if(all(se.fit==0)) warning("Zero variance estimated. Do you need to pass a newdata argument to fun()?")
  structure(list(fit = fit, se.fit = se.fit), # vcov=Sigma,
            class="predictnl")
}
confint.predictnl <- function(object,parm,level=0.95,...) {
    cf <- object$fit
    pnames <- names(cf)
    if (is.null(pnames))
        pnames <- 1:length(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- stats:::format.perc(a, 3)
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
        pct))
    ses <- object$se.fit[parm]
    ci[] <- as.vector(cf[parm]) + ses %o% fac
    ci
}
print.predictnl <- function(x, ...)
    print(structure(x,class=NULL),...)
predictnl <- function (object, ...) 
  UseMethod("predictnl")
`coef<-` <- function(obj,value) UseMethod("coef<-")
`coef<-.default` <- function(obj,value) {
    obj$coefficients <- value
    obj
}
`coef<-.mle` <- `coef<-.mle2` <- function(obj,value) {
    obj@fullcoef <- value
    obj
}
`coef<-.aov` <- function(obj,value) {
    obj$coefficients[!is.na(obj$coefficients)] <- value
    obj
}
`coef<-.Arima` <- function(obj,value) {
    obj$coef <- value
    obj
}


predictnl.default <- function(object,fun,newdata=NULL,...)
  {
      if (!is.null(newdata) || "newdata" %in% names(formals(fun))) {
          localf <- function(coef,newdata,...)
              {
                  coef(object) <- coef
                  fun(object,newdata=newdata,...)
              }
          numDeltaMethod(object,localf,newdata=newdata, ...)
      }
      else {
          localf <- function(coef,...)
              {
                  coef(object) <- coef
                  fun(object,...)
              }
          numDeltaMethod(object,localf,...)
      }
  }
setMethod("predictnl", "mle", function(object,fun,...)
          predictnl.default(object, fun, ...))
setMethod("predictnl", "mle2", function(object,fun,...)
          predictnl.default(object, fun, ...))
predictnl.lm <- function(object,fun,newdata=NULL,...)
  {
    ## Corner case: fun=predict with no newdata
    ##
    if (is.null(newdata) && "newdata" %in% names(formals(fun))) {
        stopifnot(!is.null(object$data))
        newdata <- object$data
    }
    predictnl.default(object,fun,newdata,...)
  }

predict.lm <- function(object,newdata=NULL,...) {
    if (is.null(newdata)) {
        stopifnot(!is.null(object$data))
        newdata <- object$data
    }
    stats::predict(object,newdata,...)
}
