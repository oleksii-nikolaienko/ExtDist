#' Laplace Distribution
#' 
#' @title The Laplace Distribution.
#' @description Density (d), distribution (p), and quantile (q), random number generation (r), 
#' and parameter estimation (e) functions for the Laplace distribution. Parameter estimation can 
#' only be based on an unweighted i.i.d. sample.

#' @rdname Laplace
#' @name Laplace
#' @aliases dLaplace pLaplace qLaplace rLaplace eLaplace lLaplace sLaplace iLaplace

#' @details no details yet.

#' @param params a list that includes all named parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param mu location parameter.
#' @param b scale parameter.
#' @param method parameter estimation method.
#' @param logL logical, it is assumed that the log likelihood is desired. Set to FALSE if the 
#' likelihood is wanted.
#' @param ... other parameters

#' @return dLaplace gives the density; pLaplace gives the distribution function, qLaplace gives the 
#' quantile function, rLaplace generates random variables and eLaplace estimates the parameters. 
#' lLaplace will provide the log-likelihood.

#' @note The estimation of the population mean is done using the median of the sample. Unweighted 
#' samples are not yet catered for in the eLaplace() function.

#' @author A. Jonathan R. Godfrey and Haizhen Wu

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' mu <- 1
#' b <- 2
#' X <- rLaplace(n, mu, b)
#' (est.par <- eLaplace(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dLaplace(den.x,mu=est.par$mu,b=est.par$b)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim=c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qLaplace((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim=c(-5,5), ylim=c(-5,5))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pLaplace(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim=c(0,1), ylim=c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(mu=1, b=2)
#' X <- rLaplace(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eLaplace(X,w) # estimated parameters of weighted sample
#' eLaplace(X) # estimated parameters of unweighted sample
#' 
#' # Alternative parameter estimation methods
#' eLaplace(X, method="numerical.MLE")
#' 
#' # Extracting location or scale parameters
#' est.par[attributes(est.par)$par.type=="location"]
#' est.par[attributes(est.par)$par.type=="scale"]
#' 
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rLaplace,edist=eLaplace,n=1000, rep.num=1e3, params=list(mu=1, b=2))
#' eval.estimation(rdist=rLaplace,edist=eLaplace,n=1000, rep.num=1e3, params=list(mu=1, b=2), 
#' method ="analytical.MLE")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rLaplace(1000, mu, b)
#' (est.par <- eLaplace(X))
#' H <- attributes(eLaplace(X, method="numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix
#' lLaplace(X,param=est.par)
#' lLaplace(X,param=est.par, logL=FALSE)
#' }

#' @rdname Laplace
#' @export dLaplace
dLaplace <- function(x, mu=0, b=1, params=list(mu, b)){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b
  }
  
  d <- exp(-abs(x-mu)/b) / (2*b)
}

#' @rdname Laplace
#' @export pLaplace
pLaplace <- function(q, mu=0, b=1,params=list(mu, b)){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b
  }
  
  x <- q-mu
  0.5 + 0.5*sign(x)*(1-exp(-abs(x)/b))
}

#' @rdname Laplace
#' @export qLaplace
qLaplace <- function(p, mu=0, b=1,params=list(mu, b)){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b
  }
  
  x <- p-0.5
  mu-b*sign(x)*log(1-2*abs(x))
}

#' @rdname Laplace
#' @export rLaplace
rLaplace <- function(n, mu=0, b=1,params=list(mu, b)){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b}
  
  u<-runif(n)-0.5
  x<-mu-b*sign(u)*log(1-2*abs(u))
}

#' @rdname Laplace
#' @export eLaplace
eLaplace <-     
  function(X,w, method ="numerical.MLE"){
    n<-length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    if(method == "analytical.MLE") {
      mu <- median(X, na.rm=T)
      b <- mean(abs(X-mu))
      est.par <- list(mu=mu, b=b) 
      est.par.se <- rep(NA, length(est.par))
    } else{
      if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
      method = "numerical.MLE"  
      
      est.par <- wmle(X=X, w=w, distname="Laplace",
                      initial=list(mu=0, b=1),
                      lower=list(mu=-Inf, b=0),
                      upper=list(mu=Inf, b=Inf))
      
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    } 
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Laplace"
    
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("mu","b")
    attributes(est.par)$par.type <- c("location", "scale")
    attributes(est.par)$par.vals <- c(est.par$mu, est.par$b)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Laplace
#' @export lLaplace
## (weighted) (log) likelihood function
lLaplace <- function(x, w=1, mu=0, b=1, params=list(mu, b), logL=TRUE){
  if(!missing(params)){
    mu <- params$mu
    b <- params$b}
  
  ll <- sum((-abs(x-mu)/b) - log(2*b))
  if(logL){return(ll)} else{return(exp(ll))}
}
