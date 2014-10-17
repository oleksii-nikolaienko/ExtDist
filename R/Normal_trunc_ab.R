#' truncated normal Distribution
#' 
#' @title The truncated normal distribution.
#' @importFrom truncdist dtrunc
#' @importFrom truncdist rtrunc
#' @importFrom truncdist ptrunc
#' @importFrom truncdist qtrunc
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the truncated normal distribution 
#' @rdname Normal_trunc_ab
#' @name Normal_trunc_ab
#' @aliases dNormal_trunc_ab pNormal_trunc_ab qNormal_trunc_ab rNormal_trunc_ab eNormal_trunc_ab lNormal_trunc_ab
#' @details See \href{../doc/Distributions-Truncted-Normal.html}{Distributions-Truncted-Normal}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param a,b boundary parameters.
#' @param mu,sigma shape parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lNormal_trunc_ab gives log likelihood.
#' @param ... other parameters

#' @return dNormal_trunc_ab gives the density; pNormal_trunc_ab gives the distribution function;
#' qNormal_trunc_ab gives the quantile function; rNormal_trunc_ab generates random variables; 
#' eNormal_trunc_ab estimate the parameters 

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' a <- 1
#' b <- 2
#' mu <- 2
#' sigma <- 5
#' X <- rNormal_trunc_ab(n, mu, sigma, a, b)
#' (est.par <- eNormal_trunc_ab(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dNormal_trunc_ab(den.x,params = est.par)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qNormal_trunc_ab((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(a,b), ylim = c(a,b))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pNormal_trunc_ab(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(mu=2, sigma=5, a= 1, b=2)
#' X <- rNormal_trunc_ab(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eNormal_trunc_ab(X,w) # estimated parameters of weighted sample
#' eNormal_trunc_ab(X) # estimated parameters of unweighted sample
#' 
#' # Extracting boundary and shape parameters
#' est.par[attributes(est.par)$par.type=="boundary"]
#' est.par[attributes(est.par)$par.type=="shape"]
#'  
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rNormal_trunc_ab,edist=eNormal_trunc_ab,n = 1000, rep.num = 1e3, 
#' params = list(mu=2, sigma=5, a=0, b=1), method ="numerical.MLE")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rNormal_trunc_ab(1000, mu, sigma, a, b)
#' (est.par <- eNormal_trunc_ab(X))
#' H <- attributes(eNormal_trunc_ab(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix 
#' lNormal_trunc_ab(X,param = est.par)
#' lNormal_trunc_ab(X,param = est.par, logL=FALSE)
#' }

#' @rdname Normal_trunc_ab
#' @export dNormal_trunc_ab
dNormal_trunc_ab <-
  function(x, mu=2, sigma=3, a = 0, b=1, params = list(mu, sigma, a, b)){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::dtrunc( x, spec="norm", a=a, b=b, mean = mu, sd = sigma)
    return(out)
  }

#' @rdname Normal_trunc_ab
#' @export pNormal_trunc_ab
pNormal_trunc_ab <- 
  function(q, mu=2, sigma=3, a = 0, b=1, params = list(mu=2, sigma = 5, a = 0, b = 1)){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::ptrunc( q, spec="norm", a=a, b=b, mean = mu, sd = sigma)
    return(out)
  }

#' @rdname Normal_trunc_ab
#' @export qNormal_trunc_ab
qNormal_trunc_ab <- 
  function(p, mu=2, sigma=3, a = 0, b=1, params = list(mu=2, sigma = 5, a = 0, b = 1)){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::qtrunc( p, spec="norm", a=a, b=b, mean = mu, sd = sigma)
    return(out)
  }

#' @rdname Normal_trunc_ab
#' @export rNormal_trunc_ab
rNormal_trunc_ab <- 
  function(n, mu=2, sigma=3, a = 0, b = 1, params = list(mu, sigma, a, b)){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    out <- truncdist::rtrunc(n, spec="norm", a=a, b=b, mean = mu, sd = sigma)
    return(out)
  }

#' @rdname Normal_trunc_ab
#' @export eNormal_trunc_ab
eNormal_trunc_ab <-     
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
  est.par <- wmle(X=X, w=w, distname = "Normal_trunc_ab",
                  initial=list(mu=mean(min(X),max(X)),sigma=1,a=min(X)-0.1*d,b=max(X)+0.1*d),
                  lower=list(mu=-Inf,sigma=0,a=-Inf,b=max(X)),
                  upper=list(mu=Inf,sigma=Inf,a=min(X),b=Inf)
				  )
  
  est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
  if(class(est.par.se) == "try-error") {
    est.par.se <- rep(NA, length(est.par))
  } 
} 

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "Normal_trunc_ab"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("mu","sigma","a","b")
attributes(est.par)$par.type <- c("shape","shape","boundary","boundary")
attributes(est.par)$par.vals <- c(est.par$mu, est.par$sigma, est.par$a, est.par$b)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname Normal_trunc_ab
#' @export lNormal_trunc_ab
## (weighted) (log) likelihood function
lNormal_trunc_ab <- 
  function(X, w, mu=2, sigma =3, a = 0, b = 1,  params = list(mu, sigma, a, b), logL = TRUE){
    if(!missing(params)){
      mu <- params$mu; sigma <- params$sigma; a <- params$a; b <- params$b
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dNormal_trunc_ab(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }