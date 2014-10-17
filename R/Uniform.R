#' Uniform Distribution -----------------------------------------------------
#'
#' @title The Uniform Distribution.
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the Uniform distribution 
#' @rdname Uniform
#' @name Uniform
#' @aliases dUniform pUniform qUniform rUniform eUniform lUniform
#' @details See \href{../doc/Distributions-Uniform.html}{Distributions-Uniform}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param a,b boundary parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lUniform gives log likelihood.
#' @param ... other parameters

#' @return dUniform gives the density; pUniform gives the distribution function;
#' qUniform gives the quantile function; rUniform generates random variables; 
#' eUniform estimate the parameters; sUniform gives observed scorn function 

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' a <- 1
#' b <- 2
#' X <- rUniform(n, a, b)
#' (est.par <- eUniform(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dUniform(den.x,a=est.par$a,b=est.par$b)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qUniform((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", xlab="Theoretical Quantiles", 
#' ylab="Sample Quantiles", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pUniform(sort(X), params=est.par), main="P-P Plot", xlab="Theoretical Percentile", 
#' ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(a=1, b=2)
#' X <- rUniform(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eUniform(X,w) # estimated parameters of weighted sample
#' eUniform(X) # estimated parameters of unweighted sample
#' 
#' # Alternative parameter estimation methods
#' (est.par <- eUniform(X, method = "numerical.MLE"))
#' 
#' # Extracting boundary parameters
#' est.par[attributes(est.par)$par.type=="boundary"]
#'  
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rUniform,edist=eUniform,n = 1000, rep.num = 1e3, 
#' params = list(a=2, b=5), method ="numerical.MLE")
#' eval.estimation(rdist=rUniform,edist=eUniform,n = 1000, rep.num = 1e3, 
#' params = list(a=2, b=5), method ="MOM")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rUniform(1000, a, b)
#' (est.par <- eUniform(X))
#' H <- attributes(eUniform(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix 
#' lUniform(X,param = est.par)
#' lUniform(X,param = est.par, logL=FALSE)
#' sUniform(X,param = est.par)
#' iUniform(X,param = est.par)
#' }

#' @rdname Uniform
#' @export dUniform
dUniform <-
  function(x, a=2, b =3, params = list(a, b)){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    out <- dunif(x, min = a, max = b)
    return(out)
  }


#' @rdname Uniform
#' @export pUniform
pUniform <- 
  function(q, a=2, b =3, params = list(a, b)){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    out <- punif(q, min = a, max = b)
    return(out)
  }

#' @rdname Uniform
#' @export qUniform
qUniform <- 
  function(p, a=2, b =3, params = list(a, b)){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    out <- qunif(p, min = a, max = b)
    return(out)
  }

#' @rdname Uniform
#' @export rUniform
rUniform <- 
  function(n, a=2, b =3, params = list(a, b)){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    out <- runif(n, min = a, max = b)
    return(out)
  }

#' @rdname Uniform
#' @export eUniform
eUniform <-     
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
      
	  d <- (max(X)-min(X))
      est.par <- wmle(X=X, w=w, distname = "Uniform",
                      initial=list(a = min(X)-0.1*d, b = max(X)+0.1*d),
                      lower=list(a = -Inf, b = max(X)),
                      upper=list(a = min(X), b = Inf))
      
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Uniform"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("a","b")
    attributes(est.par)$par.type <- c("boundary","boundary")
    attributes(est.par)$par.vals <- c(est.par$a, est.par$b)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Uniform
#' @export lUniform
## (weighted) (log) likelihood function
lUniform <- 
  function(X, w, a=2, b =3, params = list(a, b), logL = TRUE){
    if(!missing(params)){
      a <- params$a
      b <- params$b
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dUniform(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }