# Symmetric-Reflected Truncated Beta (SRTB)  Distribution -----------------------------------------------------
#' @title The Symmetric-Reflected Truncated Beta (SRTB)  Distribution.
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the Symmetric-Reflected Truncated Beta (SRTB). 
#' @rdname SRTB_ab
#' @name SRTB_ab
#' @aliases dSRTB_ab pSRTB_ab qSRTB_ab rSRTB_ab eSRTB_ab lSRTB_ab sSRTB_ab
#' @details See \href{../doc/Distributions-Four-Parameter-Beta.html}{Distributions-Four-Parameter-Beta}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param shape1,shape2 shape parameters.
#' @param a,b boundary parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lSRTB_ab gives log likelihood.
#' @param ... other parameters

#' @return dSRTB_ab gives the density; pSRTB_ab gives the distribution function;
#' qSRTB_ab gives the quantile function; rSRTB_ab generates random variables; 
#' eSRTB_ab estimate the parameters; sSRTB_ab gives observed scorn function 

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' a <- 1
#' b <- 2
#' shape1 <- 2
#' shape2 <- 10
#' X <- rSRTB_ab(n, shape1, shape2, a, b)
#' (est.par <- eSRTB_ab(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dSRTB_ab(den.x,params = est.par)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qSRTB_ab((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(a,b), ylim = c(a,b))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pSRTB_ab(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(shape1=2, shape2=10, a= 1, b=2)
#' X <- rSRTB_ab(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eSRTB_ab(X,w) # estimated parameters of weighted sample
#' eSRTB_ab(X) # estimated parameters of unweighted sample
#' 
#' # Extracting boundary and shape parameters
#' est.par[attributes(est.par)$par.type=="boundary"]
#' est.par[attributes(est.par)$par.type=="shape"]
#'  
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rSRTB_ab,edist=eSRTB_ab,n = 1000, rep.num = 1e3, 
#' params = list(shape1=2, shape2=10, a=0, b=1), method ="numerical.MLE")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rSRTB_ab(1000, shape1, shape2, a, b)
#' (est.par <- eSRTB_ab(X))
#' H <- attributes(eSRTB_ab(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix 
#' lSRTB_ab(X,param = est.par)
#' lSRTB_ab(X,param = est.par, logL=FALSE)
#' }

#' @rdname SRTB_ab
#' @export dSRTB_ab
dSRTB_ab <-
  function(x, shape1=2, shape2=3, a = 0, b=1, params = list(shape1, shape2, a, b)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out = dSSRTB((x-a)/(b-a),shape1,shape2)/(b-a)
    return(out)
  }

#' @rdname SRTB_ab
#' @export pSRTB_ab
pSRTB_ab <- 
  function(q, shape1=2, shape2=3, a = 0, b=1, params = list(shape1=2, shape2 = 5, a = 0, b = 1)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out = pSSRTB((q-a)/(b-a),shape1,shape2)
    return(out)
  }

#' @rdname SRTB_ab
#' @export qSRTB_ab
qSRTB_ab <- 
  function(p, shape1=2, shape2=3, a = 0, b=1, params = list(shape1=2, shape2 = 5, a = 0, b = 1)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    out = (b-a)*qSSRTB(p,shape1,shape2) + a
    return(out)
  }

#' @rdname SRTB_ab
#' @export rSRTB_ab
rSRTB_ab <- 
  function(n, shape1=2, shape2=3, a = 0, b = 1, params = list(shape1, shape2, a, b)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    x = rSSRTB(n,shape1,shape2)
    out = y = (b-a)*x + a
    return(out)
  }

#' @rdname SRTB_ab
#' @export eSRTB_ab
eSRTB_ab <-     
  function(X,w, method ="numerical.MLE"){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
{
  if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
  method = "numerical.MLE"  
  
  d <- max(X)-min(X)
  est.par <- wmle(X=X, w=w, distname = "SRTB_ab",
                  initial=list(shape1=3,shape2=3,a=min(X)-0.1*d,b=max(X)+0.1*d),
                  lower=list(shape1=1,shape2=1,a=-Inf,b=max(X)),
                  upper=list(shape1=Inf,shape2=Inf,a=min(X),b=Inf))
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
} 

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "SRTB_ab"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("shape1","shape2","a","b")
attributes(est.par)$par.type <- c("shape","shape","boundary","boundary")
attributes(est.par)$par.vals <- c(est.par$shape1, est.par$shape2, est.par$a, est.par$b)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname SRTB_ab
#' @export lSRTB_ab
## (weighted) (log) likelihood function
lSRTB_ab <- 
  function(X, w, shape1=2, shape2 =3, a = 0, b = 1,  params = list(shape1, shape2, a, b), logL = TRUE){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
      a <- params$a
      b <- params$b
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dSRTB_ab(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }