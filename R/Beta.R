#' standard Beta Distribution
#'
#' @title The standard Beta Distribution.
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the Beta distribution 
#' @rdname Beta
#' @name Beta
#' @aliases dBeta pBeta qBeta rBeta eBeta lBeta sBeta iBeta
#' @details See \href{../doc/Distributions-Beta.html}{Distributions-Beta}
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param shape1,shape2 shape parameters.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lBeta gives log likelihood.
#' @param ... other parameters

#' @return dBeta gives the density; pBeta gives the distribution function;
#' qBeta gives the quantile function; rBeta generates random variables; 
#' eBeta estimate the parameters; sBeta gives observed scorn function 

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' shape1 <- 1
#' shape2 <- 2
#' X <- rBeta(n, shape1, shape2)
#' (est.par <- eBeta(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dBeta(den.x,shape1=est.par$shape1,shape2=est.par$shape2)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qBeta((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", xlab="Theoretical Quantiles", 
#' ylab="Sample Quantiles", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pBeta(sort(X), params=est.par), main="P-P Plot", xlab="Theoretical Percentile", 
#' ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(shape1=1, shape2=2)
#' X <- rBeta(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eBeta(X,w) # estimated parameters of weighted sample
#' eBeta(X) # estimated parameters of unweighted sample
#' 
#' # Alternative parameter estimation methods
#' (est.par <- eBeta(X, method = "numerical.MLE"))
#' 
#' # Extracting shape parameters
#' est.par[attributes(est.par)$par.type=="shape"]
#'  
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rBeta,edist=eBeta,n = 1000, rep.num = 1e3, 
#' params = list(shape1=2, shape2=5), method ="numerical.MLE")
#' eval.estimation(rdist=rBeta,edist=eBeta,n = 1000, rep.num = 1e3, 
#' params = list(shape1=2, shape2=5), method ="MOM")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rBeta(1000, shape1, shape2)
#' (est.par <- eBeta(X))
#' H <- attributes(eBeta(X, method = "numerical.MLE"))$nll.hessian
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix 
#' lBeta(X,param = est.par)
#' lBeta(X,param = est.par, logL=FALSE)
#' sBeta(X,param = est.par)
#' iBeta(X,param = est.par)
#' }

#' @rdname Beta
#' @export dBeta
dBeta <-
  function(x, shape1=2, shape2 =3, params = list(shape1, shape2)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    out <- dbeta(x, shape1, shape2)
    return(out)
  }


#' @rdname Beta
#' @export pBeta
pBeta <- 
  function(q, shape1=2, shape2 =3, params = list(shape1, shape2)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    out <- pbeta(q, shape1, shape2)
    return(out)
  }

#' @rdname Beta
#' @export qBeta
qBeta <- 
  function(p, shape1=2, shape2 =3, params = list(shape1, shape2)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    qbeta(p, shape1 = shape1, shape2 = shape2)
  }

#' @rdname Beta
#' @export rBeta
rBeta <- 
  function(n, shape1=2, shape2 =3, params = list(shape1, shape2)){
    if(!missing(params)){
      shape1 <- params$shape1
      shape2 <- params$shape2
    }
    rbeta(n, shape1 = shape1, shape2 = shape2)
  }

#' @rdname Beta
#' @export eBeta
eBeta <-     
  function(X,w, method ="numerical.MLE"){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    if(method == "MOM") {
      sample.mean <- mean(w*X)
      sample.var <- (mean(w*X^2)-sample.mean^2)*n/(n-1)
      
      v <- sample.mean*(1-sample.mean)
      if(sample.var<v){
        shape1 <- sample.mean*(v/sample.var-1)
        shape2 <- (1-sample.mean)*(v/sample.var-1)
      } else {
        shape2 <- sample.mean*(v/sample.var-1)
        shape1 <- (1-sample.mean)*(v/sample.var-1)      
      } 
      
      est.par <- list(shape1 = shape1, shape2 = shape2)
      est.par.se <- rep(NA, length(est.par))
    } else {
      if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
      method = "numerical.MLE"  
      
      est.par <- wmle(X=X, w=w, distname = "Beta",
                      initial=list(shape1 = 3, shape2 = 3),
                      lower=list(shape1 = 0, shape2 = 0),
                      upper=list(shape1 = Inf, shape2 = Inf))
      
      est.par.se <- try(sqrt(diag(solve(attributes(est.par)$nll.hessian))),silent=TRUE)
      if(class(est.par.se) == "try-error") {
        est.par.se <- rep(NA, length(est.par))
      } 
    }
    
    attributes(est.par)$ob <- X
    attributes(est.par)$weights <- w
    attributes(est.par)$distname <- "Beta"
    attributes(est.par)$method <- method
    attributes(est.par)$par.name <- c("shape1","shape2")
    attributes(est.par)$par.type <- c("shape","shape")
    attributes(est.par)$par.vals <- c(est.par$shape1, est.par$shape2)
    attributes(est.par)$par.s.e <-  est.par.se  
    
    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Beta
#' @export lBeta
## (weighted) (log) likelihood function
lBeta <- 
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
    
    ll <- sum(w*((shape1-1)*log(X)+(shape2-1)*log(1-X)-log(beta(shape1,shape2))))
    #     ll <- sum(w*log(dBeta(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }

#' @rdname Beta
#' @export sBeta
## (weighted) score vectors
sBeta <- 
  function(X, w, shape1=2, shape2 =3, params = list(shape1, shape2)){
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
    
    score1 <- sum(w*(digamma(shape1+shape2)-digamma(shape1)+log(X)))
    score2 <- sum(w*(digamma(shape1+shape2)-digamma(shape2)+log(1-X)))
    
    score <- c(score1,score2)
    names(score) <- c("shape1","shape2")
    return(score)
  }

#' @rdname Beta
#' @export iBeta
## (weighted) (observed) information matrix
iBeta <- 
  function(X, w, shape1=2, shape2 =3, params = list(shape1, shape2)){
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
    
    info11 <- -sum(w*(trigamma(shape1+shape2)-trigamma(shape1)))
    info12 <- -sum(w*trigamma(shape1+shape2))
    info21 <- info12
    info22 <- -sum(w*(trigamma(shape1+shape2)-trigamma(shape2)))
    info <- matrix(c(info11,info12,info21,info22), nrow=2,ncol=2)
    rownames(info) <- colnames(info) <- c("shape1","shape2")
    return(info)
  }
