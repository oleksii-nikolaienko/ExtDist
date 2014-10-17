#' S3 methods from manuplating eDist objects.
#' 
#' @title S3 methods from manuplating eDist objects.
#' @description S3 methods from manuplating eDist objects
#' @rdname eDist
#' @name eDist
#' @aliases logLik.eDist AIC.eDist AICc.eDist BIC.eDist MDL.eDist vcov.eDist print.eDist

#' @param object,x a eDist object, which is the output of the parameter estimation functions.
#' @param ... other parameters
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param corr logic, if the vcov return correlation matrix (instead of variance-covariance matrix) 
#' or not.
#' @method logLik eDist
#' @method AIC eDist
#' @method AICc eDist
#' @method BIC eDist
#' @method MDL eDist
#' @method print eDist
#' @method vcov eDist

#' @author A. Jonathan R. Godfrey and Haizhen Wu

#' @examples \donttest{
#' X <- rnorm(20)
#' est.par <- eNormal(X, method ="numerical.MLE")
#' logLik(est.par)
#' AIC(est.par)
#' AICc(est.par)
#' BIC(est.par)
#' MDL(est.par)
#' vcov(est.par)
#' vcov(est.par,corr=TRUE)
#' print(est.par)
#' }

#' @rdname eDist
#' @export logLik.eDist
logLik.eDist <- function(object,...){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)
  return(ll)
}

#' @rdname eDist
#' @export AIC.eDist
AIC.eDist <- function(object,..., k=2){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)
  p <- length(object)
  AIC <- k*(-ll + p)
  return(AIC)
}


#' @rdname eDist
#' @export AICc
AICc <- function(object) UseMethod("AICc")

#' @rdname eDist
#' @export AICc.eDist
AICc.eDist <- function(object,...){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)
  p <- length(object) #number of paramters
  n <- length(attributes(object)$ob)
  AICc <- 2*(-ll + p + p*(p+1)/(n-p-1))
  return(AICc)
}

#' @rdname eDist
#' @export BIC.eDist
BIC.eDist <- function(object,...){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)
  n <- length(attributes(object)$ob)
  BIC <- 2*(-ll) + length(object) *log(n) 
  return(BIC)
}

#' @rdname eDist
#' @export vcov.eDist
vcov.eDist <- function(object,..., corr=FALSE){
  vcov = solve(attributes(object)$nll.hessian)
  cor = cov2cor(vcov)
  if(corr){return(cor)} else {return(vcov)}
}

#' @references Myung, I. (2000). The Importance of Complexity in Model Selection. 
#' Journal of mathematical psychology, 44(1), 190-204.

#' @rdname eDist
#' @export MDL
MDL <- function(object) UseMethod("MDL")

#' @rdname eDist
#' @export MDL.eDist
MDL.eDist <- function(object,...){
  lFoo <- get(paste0("l", attributes(object)$distname))
  ll <- lFoo(attributes(object)$ob, w=attributes(object)$weights, params=object)  
  nll.hessian <- attributes(object)$nll.hessian
  MDL = ifelse(all(is.finite(nll.hessian)),
               -ll + 1/2*log(abs(det(attributes(object)$nll.hessian))),
               NA)
  return(MDL)
}


#' @rdname eDist
#' @export print.eDist
print.eDist <- function(x,...){
  cat("\nParameters for the", attributes(x)$distname, "distribution. \n(found using the ", attributes(x)$method, "method.)\n\n")
  if(any(is.na(attributes(x)$par.s.e))) {
    print(data.frame(Parameter=attributes(x)$par.name, 
                     Type=attributes(x)$par.type, 
                     Estimate=attributes(x)$par.vals), 
          row.names=FALSE )
  } else {
    print(data.frame(Parameter=attributes(x)$par.name, 
                     Type=attributes(x)$par.type, 
                     Estimate=attributes(x)$par.vals,
                     S.E. = attributes(x)$par.s.e ), 
          row.names=FALSE )
  }
  cat("\n\n")
}

