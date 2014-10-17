#' Logistic Distribution
#' 
#' @title The Logistic Distribution.

#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the Logistic distribution 
#' @rdname Logistic
#' @name Logistic

#' @aliases dLogistic pLogistic qLogistic rLogistic eLogistic lLogistic
#' @details See \href{../doc/Distributions-Logistic.html}{Distributions-Logistic}

#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lLogistic gives log likelihood.
#' @param ... other parameters

#' @return dLogistic gives the density; pLogistic gives the distribution function;
#' qLogistic gives the quantile function; rLogistic generates random variables; 
#' eLogistic estimate the parameters
#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' location <- 1.5
#' scale <- 0.5
#' X <- rLogistic(n, location, scale)
#' (est.par <- eLogistic(X))
#' 

#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dLogistic(den.x,location=est.par$location,scale=est.par$scale)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qLogistic((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(-5,5), ylim = c(-5,5))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pLogistic(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(location=1, scale=2)
#' X <- rLogistic(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eLogistic(X,w) # estimated parameters of weighted sample
#' eLogistic(X) # estimated parameters of unweighted sample
#' 
#' # Extracting location or scale parameters
#' est.par[attributes(est.par)$par.type=="location"]
#' est.par[attributes(est.par)$par.type=="scale"]
#' 
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rLogistic,edist=eLogistic,n = 1000, rep.num = 1e3, 
#'    params = list(location=1, scale=2))
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rLogistic(1000, location, scale)
#' (est.par <- eLogistic(X))
#' H <- attributes(eLogistic(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix
#' lLogistic(X,param = est.par)
#' lLogistic(X,param = est.par, logL=FALSE)
#' }

#' @rdname Logistic
#' @export dLogistic

dLogistic <-
  function(x, location = 0, scale = 1, params = list(location = 0, scale = 1)){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = dlogis(x, location, scale)
    return(out)
  }

#' @rdname Logistic
#' @export pLogistic

pLogistic <- 
  function(q, location = 0, scale = 1, params = list(location = 0, scale = 1)){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = plogis(q,location,scale)
    return(out)
}

#' @rdname Logistic
#' @export qLogistic

qLogistic <- 
  function(p, location = 0, scale = 1, params = list(location = 0, scale = 1)){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = qlogis(p,location,scale)
    return(out)
}

#' @rdname Logistic
#' @export rLogistic
rLogistic <- 
  function(n, location = 0, scale = 1, params = list(location = 0, scale = 1)){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    out = rlogis(n,location,scale)
    return(out)
  }

#' @rdname Logistic
#' @export eLogistic
eLogistic <-     
  function(X,w, method ="numerical.MLE"){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
	{if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
      method = "numerical.MLE"  
      est.par <- wmle(X=X, w=w, distname = "Logistic",
                      initial=list(location = 0, scale = 1),
                      lower=list(location = -Inf, scale = 0),
                      upper=list(location = Inf, scale = Inf))

      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Logistic"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("location","scale")
    attributes(est.par)$par.type <- c("location","scale")
    attributes(est.par)$par.vals <- c(est.par$location, est.par$scale)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Logistic
#' @export lLogistic
## (weighted) (log) likelihood function
lLogistic <- 
  function(X, w, location = 0, scale = 1, params = list(location = 0, scale = 1), logL = TRUE){
    if(!missing(params)){
      location <- params$location
      scale <- params$scale
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dLogistic(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
