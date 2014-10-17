#' Weibull Distribution
#' 
#' @title The Weibull Distribution.
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the Weibull distribution 
#' @rdname Weibull
#' @name Weibull

#' @aliases dWeibull pWeibull qWeibull rWeibull eWeibull lWeibull
#' @details See \href{../doc/Distributions-Weibull.html}{Distributions-Weibull}

#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param shape shape parameter.
#' @param scale scale parameter.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lWeibull gives log likelihood.
#' @param ... other parameters

#' @return dWeibull gives the density; pWeibull gives the distribution function;
#' qWeibull gives the quantile function; rWeibull generates random variables; 
#' eWeibull estimate the parameters
#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' shape <- 1.5
#' scale <- 0.5
#' X <- rWeibull(n, shape, scale)
#' (est.par <- eWeibull(X))
#' 

#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dWeibull(den.x,shape=est.par$shape,scale=est.par$scale)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qWeibull((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(0,5), ylim = c(0,5))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pWeibull(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(shape=1, scale=2)
#' X <- rWeibull(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eWeibull(X,w) # estimated parameters of weighted sample
#' eWeibull(X) # estimated parameters of unweighted sample
#' 
#' # Extracting shape or scale parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#' est.par[attributes(est.par)$par.type=="scale"]
#' 
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rWeibull,edist=eWeibull,n = 1000, rep.num = 1e3, 
#'    params = list(shape=1, scale=2))
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rWeibull(1000, shape, scale)
#' (est.par <- eWeibull(X))
#' H <- attributes(eWeibull(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix
#' lWeibull(X,param = est.par)
#' lWeibull(X,param = est.par, logL=FALSE)
#' }

#' @rdname Weibull
#' @export dWeibull

dWeibull <-
  function(x, shape = 2, scale = 2, params = list(shape = 2, scale = 2)){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::dgamma(x, shape, scale)
    return(out)
  }

#' @rdname Weibull
#' @export pWeibull

pWeibull <- 
  function(q, shape = 2, scale = 2, params = list(shape = 2, scale = 2)){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::pgamma(q,shape,scale)
    return(out)
}

#' @rdname Weibull
#' @export qWeibull

qWeibull <- 
  function(p, shape = 2, scale = 2, params = list(shape = 2, scale = 2)){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::qgamma(p,shape,scale)
    return(out)
}

#' @rdname Weibull
#' @export rWeibull
rWeibull <- 
  function(n, shape = 2, scale = 2, params = list(shape = 2, scale = 2)){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    out = stats::rgamma(n,shape,scale)
    return(out)
  }

#' @rdname Weibull
#' @export eWeibull
eWeibull <-     
  function(X,w, method ="numerical.MLE"){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
	{if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
      method = "numerical.MLE"  
      est.par <- wmle(X=X, w=w, distname = "Weibull",
                      initial=list(shape = 1, scale = 1),
                      lower=list(shape = 0, scale = 0),
                      upper=list(shape = Inf, scale = Inf))

      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Weibull"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("shape","scale")
    attributes(est.par)$par.type <- c("shape","scale")
    attributes(est.par)$par.vals <- c(est.par$shape, est.par$scale)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Weibull
#' @export lWeibull
## (weighted) (log) likelihood function
lWeibull <- 
  function(X, w, shape = 2, scale = 2, params = list(shape = 2, scale = 2), logL = TRUE){
    if(!missing(params)){
      shape <- params$shape
      scale <- params$scale
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dWeibull(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
