# S3 methods from manuplating eDist objects. -----------------------------------------------------
#' @title S3 methods from manuplating eDist objects.
#' @description S3 methods from manuplating eDist objects
#' @rdname eDist
#' @name eDist
#' @aliases logLik.eDist AIC.eDist BIC.eDist vcov.eDist print.eDist

#' @param object,x a eDist object, which is the output of the parameter estimation functions.
#' @param ... other parameters
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param corr logic, if the vcov return correlation matrix (instead of variance-covariance matrix) 
#' or not.
#' @method logLik eDist
#' @method AIC eDist
#' @method BIC eDist
#' @method print eDist
#' @method vcov eDist

#' @author A. Jonathan R. Godfrey and Haizhen Wu

#' @examples \donttest{
#' X <- rnorm(20)
#' est.par <- eNormal(X, method ="numerical.MLE")
#' logLik(est.par)
#' AIC(est.par)
#' BIC(est.par)
#' vcov(est.par)
#' vcov(est.par,corr=TRUE)
#' print(est.par)
#' }

#' @rdname eDist
#' @export logLik.eDist
logLik.eDist <- function(object,...){
  lFoo <- get(paste0("l", attr(object, "distname")))
  ll <- lFoo(attr(object, "ob"), w=attr(object,"weights"), params=object)
  return(ll)
}

#' @rdname eDist
#' @export AIC.eDist
AIC.eDist <- function(object,..., k = 2){
  lFoo <- get(paste0("l", attr(object, "distname")))
  ll <- lFoo(attr(object, "ob"), w=attr(object,"weights"), params=object)
  AIC <- k*(-ll + length(object))
  return(AIC)
}

#' @rdname eDist
#' @export BIC.eDist
BIC.eDist <- function(object,...){
  lFoo <- get(paste0("l", attr(object, "distname")))
  ll <- lFoo(attr(object, "ob"), w=attr(object,"weights"), params=object)
  n <- length(attr(object, "ob"))
  BIC <- 2*(-ll) + length(object) *log(n) 
  return(BIC)
}

#' @rdname eDist
#' @export vcov.eDist
vcov.eDist <- function(object,..., corr=FALSE){
  vcov = solve(attr(object,"nll.hessian"))
  cor = cov2cor(vcov)
  if(corr){return(cor)} else {return(vcov)}
}

#' @rdname eDist
#' @export print.eDist
print.eDist <- function(x,...){
  cat("\nParameters for the", attr(x,"distname"), "distribution. \n(found using the ", attr(x,"method"), "method.)\n\n")
  print(data.frame(Parameter=attr(x,"par.name"), 
                   Type=attr(x,"par.type"), 
                   Estimate=attr(x,"par.vals"),
                   S.E. = attr(x,"par.s.e") ), 
        row.names=FALSE )
  cat("\n\n")
}

