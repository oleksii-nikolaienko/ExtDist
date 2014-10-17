#' Standard Symmetric-Reflected Truncated Beta (SRTB) Distribution 
#'
#' @title The standard Symmetric-Reflected Truncated Beta (SRTB) Distribution.
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the SSRTB distribution 
#' @rdname SSRTB
#' @name SSRTB
#' @aliases dSSRTB pSSRTB qSSRTB rSSRTB eSSRTB lSSRTB sSSRTB iSSRTB
#' @details See \href{../doc/Distributions-SSRTB.html}{Distributions-SSRTB}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param shape1,shape2 shape parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lSSRTB gives log likelihood.
#' @param ... other parameters

#' @return dSSRTB gives the density; pSSRTB gives the distribution function;
#' qSSRTB gives the quantile function; rSSRTB generates random variables; 
#' eSSRTB estimate the parameters; sSSRTB gives observed scorn function 

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' shape1 <- 2
#' shape2 <- 10
#' X <- rSSRTB(n, shape1, shape2)
#' (est.par <- eSSRTB(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dSSRTB(den.x,shape1=est.par$shape1,shape2=est.par$shape2)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qSSRTB((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", xlab="Theoretical Quantiles", 
#' ylab="Sample Quantiles", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pSSRTB(sort(X), params=est.par), main="P-P Plot", xlab="Theoretical Percentile", 
#' ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(shape1=1, shape2=2)
#' X <- rSSRTB(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eSSRTB(X,w) # estimated parameters of weighted sample
#' eSSRTB(X) # estimated parameters of unweighted sample
#' 
#' # Alternative parameter estimation methods
#' (est.par <- eSSRTB(X, method = "numerical.MLE"))
#' 
#' # Extracting shape parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#'  
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rSSRTB,edist=eSSRTB,n = 1000, rep.num = 1e3, 
#' params = list(shape1=2, shape2=5), method ="numerical.MLE")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rSSRTB(1000, shape1, shape2)
#' (est.par <- eSSRTB(X))
#' H <- attributes(eSSRTB(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix 
#' lSSRTB(X,param = est.par)
#' lSSRTB(X,param = est.par, logL=FALSE)
#' }

#' @rdname SSRTB
#' @export dSSRTB
dSSRTB <-
  function(x, shape1=2, shape2 =3, params = list(shape1, shape2)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    m = (1-1/shape1)/(1+(shape2-2)/shape1)
    out = m*dbeta(2*m*(1/2-abs(x-1/2)),shape1,shape2)/pbeta(m,shape1,shape2)

    return(out)
  }


#' @rdname SSRTB
#' @export pSSRTB
pSSRTB <- 
  function(q, shape1=2, shape2 =3, params = list(shape1, shape2)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    m = (1-1/shape1)/(1+(shape2-2)/shape1)
    out= ifelse(q<1/2, 
                pbeta(2*m*q,shape1,shape2)/2/pbeta(m,shape1,shape2), 
                1-pbeta(2*m*(1-q),shape1,shape2)/2/pbeta(m,shape1,shape2))

    return(out)
  }

#' @rdname SSRTB
#' @export qSSRTB
qSSRTB <- 
  function(p, shape1=2, shape2 =3, params = list(shape1, shape2)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    m = (1-1/shape1)/(1+(shape2-2)/shape1)
    out = ifelse(p<1/2, qbeta((p<1/2)*p*2*pbeta(m,shape1,shape2),shape1,shape2)/2/m,
                 1-qbeta((p>=1/2)*(1-p)*2*pbeta(m,shape1,shape2),shape1,shape2)/2/m)
    return(out)
  }

#' @rdname SSRTB
#' @export rSSRTB
rSSRTB <- 
  function(n, shape1=2, shape2 =3, params = list(shape1, shape2)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    qSSRTB(runif(n),shape1,shape2)
  }

#' @rdname SSRTB
#' @export eSSRTB
eSSRTB <-     
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
      
      est.par <- wmle(X=X, w=w, distname = "SSRTB",
                      initial=list(shape1 = 5, shape2 = 5),
                      lower=list(shape1 = 1.1, shape2 = 1),
                      upper=list(shape1 = Inf, shape2 = Inf))
      
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "SSRTB"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("shape1","shape2")
    attributes(est.par)$par.type <- c("shape","shape")
    attributes(est.par)$par.vals <- c(est.par$shape1, est.par$shape2)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname SSRTB
#' @export lSSRTB
## (weighted) (log) likelihood function
lSSRTB <- 
  function(X, w, shape1=2, shape2 =3, params = list(shape1, shape2), logL = TRUE){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*log(dSSRTB(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }