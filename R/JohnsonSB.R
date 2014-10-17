#' Johnson SB Distribution
#' 
#' @title The Johnson SB Distribution.
#' @importFrom SuppDists dJohnson
#' @importFrom SuppDists rJohnson
#' @importFrom SuppDists pJohnson
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the 4 Parameter beta distribution 
#' @rdname JohnsonSB
#' @name JohnsonSB
#' @aliases dJohnsonSB pJohnsonSB qJohnsonSB rJohnsonSB eJohnsonSB lJohnsonSB
#' @details Johnson SB is bounded in [xi,xi+lambda].
#' See \href{../doc/Distributions-Johnson-SB.html}{Distributions-Johnson-SB}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param xi,lambda location-scale parameters.
#' @param gamma,delta shape parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lJohnsonSB gives log likelihood.
#' @param ... other parameters

#' @return dJohnsonSB gives the density; pJohnsonSB gives the distribution function;
#' qJohnsonSB gives the quantile function; rJohnsonSB generates random variables; 
#' eJohnsonSB estimate the parameters
#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' gamma = -0.5; delta = 2; xi = -0.5; lambda = 2
#' X <- rJohnsonSB(n, gamma, delta, xi, lambda)
#' (est.par <- eJohnsonSB(X))
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dJohnsonSB(den.x,params = est.par)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qJohnsonSB((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(a,b), ylim = c(a,b))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pJohnsonSB(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)
#' X <- rJohnsonSB(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eJohnsonSB(X,w) # estimated parameters of weighted sample
#' eJohnsonSB(X) # estimated parameters of unweighted sample
#' 
#' # Extracting location, scale and shape parameters
#' est.par[attributes(est.par)$par.type=="location"]
#' est.par[attributes(est.par)$par.type=="scale"]
#' est.par[attributes(est.par)$par.type=="shape"]
#'  
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rJohnsonSB,edist=eJohnsonSB,n = 1000, rep.num = 1e3, 
#' params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2), method ="numerical.MLE")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rJohnsonSB(1000, gamma, delta, xi, lambda)
#' (est.par <- eJohnsonSB(X))
#' H <- attributes(eJohnsonSB(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix 
#' lJohnsonSB(X,param = est.par)
#' lJohnsonSB(X,param = est.par, logL=FALSE)
#' }

#' @rdname JohnsonSB
#' @export dJohnsonSB
dJohnsonSB <-
  function(x, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out=SuppDists::dJohnson(x, parms = list(gamma, delta, xi, lambda, type= "SB"))
    out[is.nan(out)]=0
    return(out)
  }
dJohnsonSB_ab <-
  function(x, gamma = -0.5, delta = 2, a = -0.5, b = 1.5, params = list(gamma = -0.5, delta = 2, a = -0.5, b = 1.5)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      a <- params$a
      b <- params$b
    }
    xi <- a
    lambda <- (b-a)
    
    out = dJohnsonSB(x, params=list(gamma=gamma, delta=delta, xi=xi, lambda=lambda))
    return(out)
  }
#' @rdname JohnsonSB
#' @export pJohnsonSB
pJohnsonSB <- 
  function(q, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out = NULL
    for( i in 1:length(q)){
      out[i] = try(SuppDists::pJohnson(q[i], parms = list(gamma, delta, xi, lambda, type= "SB")),silent=T)
      if((!is.numeric(out[i]))&(q[i]<=xi)) out[i]= 0
      if((!is.numeric(out[i]))&(q[i]>=xi+lambda)) out[i]= 1
    }
    return(as.numeric(out))
  }

#' @rdname JohnsonSB
#' @export qJohnsonSB
qJohnsonSB <- 
  function(p, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out = lambda/(exp(-(qnorm(p)-gamma)/delta)+1)+xi
    return(out)
  }

#' @rdname JohnsonSB
#' @export rJohnsonSB
rJohnsonSB <- 
  function(n, gamma = -0.5, delta = 2, xi = -0.5, lambda = 2, params = list(gamma = -0.5, delta = 2, xi = -0.5, lambda = 2)){
    if(!missing(params)){
      gamma <- params$gamma
      delta <- params$delta
      xi <- params$xi
      lambda <- params$lambda
    }
    out <- SuppDists::rJohnson(n, parms = list(gamma=gamma, delta=delta, xi=xi, lambda=lambda, type= "SB"))
    return(out)
  }

#' @rdname JohnsonSB
#' @export eJohnsonSB
eJohnsonSB <-     
  function(X,w, method ="numerical.MLE"){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    # xi and lambda cannot be fitted directly.
    
{if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
 method = "numerical.MLE"  
 
#     l = min(X)
#     u = max(X)
#     Y=(X-l)/(u-l)
# 	
#  out <- wmle(X=Y, w=w, distname = "JohnsonSB_ab",
#              initial=list(gamma = -0.5, delta = 2, a = -1, b = 2),
#              lower=list(gamma = -Inf, delta = 0, a = -Inf, b = 1),
#              upper=list(gamma = Inf, delta = Inf, a = 0, b = Inf)
#  )
#     out$a = l+out$a*(u-l)
#     out$b = l+out$b*(u-l)
# 	
# 	deriv = c(1,1,1/(u-l),1/(u-l))
#     if(!is.null(attr(out,"nll.hessian"))) {    
#       attr(out,"nll.hessian") = attr(out,"nll.hessian") * deriv%*%t(deriv)
#     }
	
 out <- wmle(X=X, w=w, distname = "JohnsonSB_ab",
             initial=list(gamma = -0.5, delta = 2, a = min(X)-1, b = max(X)+1),
             lower=list(gamma = -Inf, delta = 0, a = -Inf, b = max(X)),
             upper=list(gamma = Inf, delta = Inf, a = min(X), b = Inf)
 )

 est.par = list(gamma = out$gamma , delta = out$delta, xi = out$a, lambda = out$b-out$a)
 trans.matrix = diag(rep(1,4))
 trans.matrix[4,3]=1 
 attributes(est.par)$nll.hessian = trans.matrix %*% attributes(out)$nll.hessian %*% t(trans.matrix)
 
 est.par.se <- try(suppressWarnings(sqrt(diag(solve(attributes(est.par)$nll.hessian )))),silent=TRUE)
 if(class(est.par.se) == "try-error" | any(is.nan(est.par.se))) {
   est.par.se <- rep(NA, length(est.par))
 } 
}

attributes(est.par)$ob <- X
attributes(est.par)$weights <- w
attributes(est.par)$distname <- "JohnsonSB"
attributes(est.par)$method <- method
attributes(est.par)$par.name <- c("gamma","delta","xi","lambda")
attributes(est.par)$par.type <- c("shape","shape","location","scale")
attributes(est.par)$par.vals <- c(est.par$gamma, est.par$delta, est.par$xi, est.par$lambda)
attributes(est.par)$par.s.e <-  est.par.se  

class(est.par) <- "eDist"

return(est.par)
  }

#' @rdname JohnsonSB
#' @export lJohnsonSB
## (weighted) (log) likelihood function
lJohnsonSB <- 
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
    
    ll <- sum(w*log(dJohnsonSB(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }
