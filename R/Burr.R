#' Burr's Distribution
#' 
#' @title The Burr's Distribution.
#' @importFrom VGAM dparetoIV
#' @importFrom VGAM rparetoIV
#' @importFrom VGAM pparetoIV
#' @importFrom VGAM qparetoIV
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the Burr distribution 
#' @rdname Burr
#' @name Burr
#' @aliases dBurr pBurr qBurr rBurr eBurr lBurr
#' @details See \href{../doc/Distributions-Burr.html}{Distributions-Burr}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param b scale parameters.
#' @param g,s shape parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lBurr gives log likelihood.
#' @param ... other parameters

#' @return dBurr gives the density; pBurr gives the distribution function;
#' qBurr gives the quantile function; rBurr generates random variables; 
#' eBurr estimate the parameters
#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' b = 1; g = 2; s = 2
#' X <- rBurr(n, b = 1, g = 2, s = 2)
#' (est.par <- eBurr(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dBurr(den.x, b=est.par$b, g=est.par$g, s=est.par$s)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qBurr((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(0,5), ylim = c(0,5))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pBurr(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(b=1, g=2, s =2)
#' X <- rBurr(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eBurr(X,w) # estimated parameters of weighted sample
#' eBurr(X) # estimated parameters of unweighted sample
#' 
#' # Extracting shape or scale parameters
#' est.par[attributes(est.par)$par.type=="scale"]
#' est.par[attributes(est.par)$par.type=="shape"]
#' 
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rBurr,edist=eBurr,n = 1000, rep.num = 1e3, params = list(b=1, g=2, s =2))
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rBurr(1000, b = 1, g = 2, s = 2)
#' (est.par <- eBurr(X))
#' H <- attributes(eBurr(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix
#' lBurr(X,param = est.par)
#' lBurr(X,param = est.par, logL=FALSE)
#' }

#' @rdname Burr
#' @export dBurr
dBurr <-
  function(x, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2)){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    VGAM::dparetoIV(x, location = 0, scale = b, inequality = 1/g, shape = s)  }

#' @rdname Burr
#' @export pBurr
pBurr <- 
  function(q, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2)){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    q[q<=0] =0
    out = q
    c = (q/b)^g
    d = s*c
    out[q == Inf] =1  
    out[c>1e10 & q!=Inf] = 1- (q[c>1e10 & q!=Inf]/b)^(-g*s)
    out[c<1e-10 & q!=Inf] = 1- exp(-d[c<1e-10 & q!=Inf])
    out[c<=1e10 & c>=1e-10 & q!=Inf] = VGAM::pparetoIV(q[c<=1e10 & c>=1e-10 & q!=Inf], location = 0, scale = b, inequality = 1/g, shape = s)
    return(out)
  }

#' @rdname Burr
#' @export qBurr
qBurr <- 
  function(p, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2)){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    out = p
    c= (1/(1-p))^(1/s)
    d = -log(1-p)/s
    
    out[p== 0] =0
    
    c= (1/(1-p))^(1/s)
    out[c>1e10 & p!=0] = (1/(1-p[c>1e10 & p!=0]))^(1/(s*g))*b
    
    d = -log(1-p)/s
    out[d<1e-10 & p!=0] = b*(d[d<1e-10 & p!=0])^(1/g) 
    out[c<=1e10 & d>=(1e-10) & p!=0] = VGAM::qparetoIV(p[c<=1e10 & d>=(1e-10) & p!=0], location = 0, scale = b, inequality = 1/g, shape = s)
    return(out)
  }

#' @rdname Burr
#' @export rBurr
rBurr <- 
  function(n, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2)){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    VGAM::rparetoIV(n, location = 0, scale = b, inequality = 1/g, shape = s)
  }

#' @rdname Burr
#' @export eBurr
eBurr <-     
  function(X,w, method ="numerical.MLE"){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
{ if(method != "numerical.MLE") stop(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
  method = "numerical.MLE"  
  
  est.par <- wmle(X=X, w=w, distname = "Burr",
                  initial=list(b=2, g=2, s =2),
                  lower=list(b=0, g=0, s =0),
                  upper=list(b=Inf, g=Inf, s =Inf))
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
}

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "Burr"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("b","g","s")
attributes(est.par)$par.type <- c("scale","shape","shape")
attributes(est.par)$par.vals <- c(est.par$b, est.par$g, est.par$s)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname Burr
#' @export lBurr
## (weighted) (log) likelihood function
lBurr <- 
  function(X, w, b = 1, g = 2, s = 2, params = list(b = 1, g = 2, s = 2), logL = TRUE){
    if(!missing(params)){
      b = params$b; g = params$g; s = params$s
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dBurr(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }