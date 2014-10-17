#' Triangular Distribution
#' 
#' @title The Triangular Distribution.
#' @importFrom VGAM dtriangle
#' @importFrom VGAM rtriangle
#' @importFrom VGAM ptriangle
#' @importFrom VGAM qtriangle
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the Triangular distribution 
#' @rdname Triangular
#' @name Triangular
#' @aliases dTriangular pTriangular qTriangular rTriangular eTriangular lTriangular
#' @details See \href{../doc/Distributions-Four-Parameter-Beta.html}{Distributions-Four-Parameter-Beta}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param theta shape parameters.
#' @param a,b boundary parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lTriangular gives log likelihood.
#' @param ... other parameters

#' @return dTriangular gives the density; pTriangular gives the distribution function;
#' qTriangular gives the quantile function; rTriangular generates random variables; 
#' eTriangular estimate the parameters

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' a <- 1
#' b <- 2
#' theta <- 0.5
#' X <- rTriangular(n, a, b, theta)
#' (est.par <- eTriangular(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dTriangular(den.x,params = est.par)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qTriangular((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(a,b), ylim = c(a,b))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pTriangular(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(a=0, b=1, theta=0.5)
#' X <- rTriangular(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eTriangular(X,w) # estimated parameters of weighted sample
#' eTriangular(X) # estimated parameters of unweighted sample
#' 
#' # Extracting boundary and shape parameters
#' est.par[attributes(est.par)$par.type=="boundary"]
#' est.par[attributes(est.par)$par.type=="shape"]
#'  
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rTriangular,edist=eTriangular,n = 1000, rep.num = 1e3, 
#' params = list(a=0, b=1, theta=0.5), method ="numerical.MLE")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rTriangular(1000, a, b, theta)
#' (est.par <- eTriangular(X))
#' H <- attributes(eTriangular(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix 
#' lTriangular(X,param = est.par)
#' lTriangular(X,param = est.par, logL=FALSE)
#' sTriangular(X,param = est.par)
#' }

#' @rdname Triangular
#' @export dTriangular
dTriangular <-
  function(x, a=0, b=1, theta=0.5, params = list(a, b, theta)){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    out <- VGAM::dtriangle(x, theta = theta*(b-a)+a, lower = a, upper = b)
    return(out)
  }

#' @rdname Triangular
#' @export pTriangular
pTriangular <- 
  function(q, a=0, b=1, theta=0.5, params = list(a, b, theta)){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    out <- VGAM::ptriangle(q, theta = theta*(b-a)+a, lower = a, upper = b)
    return(out)
  }

#' @rdname Triangular
#' @export qTriangular
qTriangular <- 
  function(p, a=0, b=1, theta=0.5, params = list(a, b, theta)){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    out <- p
    out[p==0] <- a
    out[p!=0] <- VGAM::qtriangle(p[p!=0], theta = theta*(b-a)+a, lower = a, upper = b)
    return(out)
  }

#' @rdname Triangular
#' @export rTriangular
rTriangular <- 
  function(n, a=0, b=1, theta=0.5, params = list(a, b, theta)){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    out <- VGAM::rtriangle(n, theta = theta*(b-a)+a, lower = a, upper = b)
    return(out)
  }

#' @rdname Triangular
#' @export eTriangular
eTriangular <-     
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
  est.par <- wmle(X=X, w=w, distname = "Triangular",
            initial=list(a=min(X)-0.1*d,b=max(X)+0.1*d, theta= 0.5),
            lower=list( a=-Inf,b=max(X), theta= 1e-10),
            upper=list(a=min(X),b=Inf, theta= 1- 1e-10))
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
} 

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "Triangular"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("a","b","theta")
attributes(est.par)$par.type <- c("boundary","boundary","shape")
attributes(est.par)$par.vals <- c(est.par$a, est.par$b, est.par$theta)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname Triangular
#' @export lTriangular
## (weighted) (log) likelihood function
lTriangular <- 
  function(X, w, a=0, b=1, theta=0.5,  params = list(a, b, theta), logL = TRUE){
    if(!missing(params)){
      a <- params$a; b <- params$b; theta = params$theta  
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dTriangular(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }