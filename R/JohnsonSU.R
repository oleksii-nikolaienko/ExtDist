#' Johnson SU Distribution
#' 
#' @title The Johnson SU Distribution.
#' @importFrom SuppDists dJohnson
#' @importFrom SuppDists rJohnson
#' @importFrom SuppDists pJohnson
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the 4 Parameter beta distribution 
#' @rdname JohnsonSU
#' @name JohnsonSU
#' @aliases dJohnsonSU pJohnsonSU qJohnsonSU rJohnsonSU eJohnsonSU lJohnsonSU
#' @details See \href{../doc/Distributions-Johnson-SU.html}{Distributions-Johnson-SU}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param xi,lambda location-scale parameters.
#' @param gamma,delta shape parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lJohnsonSU gives log likelihood.
#' @param ... other parameters

#' @return dJohnsonSU gives the density; pJohnsonSU gives the distribution function;
#' qJohnsonSU gives the quantile function; rJohnsonSU generates random variables; 
#' eJohnsonSU estimate the parameters
#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' gamma = -0.5; delta = 2; xi = -0.5; lambda = 2
#' X <- rJohnsonSU(n, gamma, delta, xi, lambda)
#' (est.par <- eJohnsonSU(X))
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dJohnsonSU(den.x,params = est.par)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qJohnsonSU((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(-5,5), ylim = c(-5,5))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pJohnsonSU(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)
#' X <- rJohnsonSU(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eJohnsonSU(X,w) # estimated parameters of weighted sample
#' eJohnsonSU(X) # estimated parameters of unweighted sample
#' 
#' # Extracting shape parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#'  
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rJohnsonSU,edist=eJohnsonSU,n = 1000, rep.num = 1e3, 
#' params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2), method ="numerical.MLE")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rJohnsonSU(1000, gamma, delta, xi, lambda)
#' (est.par <- eJohnsonSU(X))
#' H <- attributes(eJohnsonSU(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix 
#' lJohnsonSU(X,param = est.par)
#' lJohnsonSU(X,param = est.par, logL=FALSE)
#' }

#' @rdname JohnsonSU
#' @export dJohnsonSU
dJohnsonSU <-
  function(x, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }

    u = (x - xi)/lambda
    out = x
    out[abs(u)<1e-3] = dnorm(gamma+delta*u[abs(u)<1e-3])*delta/lambda  # good approximation when abs(u) is small. log(u+sqrt(1+u^2)) ~ u.
    out[u < (-1e10)] = dnorm(gamma+delta*(-Inf)) *(delta/lambda/sqrt(u[u<(-1e10)]^2+1))  # for u<-1e10, log(u+sqrt(1+u^2)) ~ -Inf . but if u = -1e200, R gives Inf
    out[abs(u)>=1e-3 & u>(-1e10)] = dnorm(gamma+delta*log(u[abs(u)>=1e-3 & u>(-1e10)]+sqrt(1+u[abs(u)>=1e-3 & u>(-1e10)]^2))) *(delta/lambda/sqrt(u[abs(u)>=1e-3 & u>(-1e10)]^2+1))    
    return(out)
  }
#' @rdname JohnsonSU
#' @export pJohnsonSU
pJohnsonSU <- 
  function(q, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }

    u = (q - xi)/lambda
    out = pnorm(delta*asinh(u)+gamma)
    return(out)
  }

#' @rdname JohnsonSU
#' @export qJohnsonSU
qJohnsonSU <- 
  function(p, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    u = (qnorm(p)-gamma)/delta
    out = lambda*sinh((qnorm(p)-gamma)/delta)+xi
    return(out)
  }

#' @rdname JohnsonSU
#' @export rJohnsonSU
rJohnsonSU <- 
  function(n, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out <- SuppDists::rJohnson(n, parms = list(gamma=gamma, delta=delta, xi=xi, lambda=lambda, type= "SU"))
    return(out)
  }

#' @rdname JohnsonSU
#' @export eJohnsonSU
eJohnsonSU <-     
  function(X,w, method ="numerical.MLE"){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
{if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
 method = "numerical.MLE"  
	
 est.par <- wmle(X=X, w=w, distname = "JohnsonSU",
          initial=list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2),
          lower=list(gamma = -Inf, delta = 1, xi = -Inf, lambda = 0),
          upper=list(gamma = Inf, delta = Inf, xi = Inf, lambda = Inf)
 )

 est.par.se <- try(suppressWarnings(sqrt(diag(solve(attributes(est.par)$nll.hessian )))),silent=TRUE)
 if(class(est.par.se) == "try-error" | any(is.nan(est.par.se))) {
   est.par.se <- rep(NA, length(est.par))
 } 
}

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "JohnsonSU"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("gamma","delta","xi","lambda")
attributes(est.par)$par.type <- c("shape","shape","boundary","boundary")
attributes(est.par)$par.vals <- c(est.par$gamma, est.par$delta, est.par$xi, est.par$lambda)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname JohnsonSU
#' @export lJohnsonSU
## (weighted) (log) likelihood function
lJohnsonSU <- 
  function(X, w, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2,  params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2), logL = TRUE){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dJohnsonSU(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
