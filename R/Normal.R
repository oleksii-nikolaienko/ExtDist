# Normal Distribution -----------------------------------------------------
#' @title The Normal Distribution.
#' @description Density, distribution function, quantile function, random 
#' generation function and parameter estimation function (based on weighted or 
#' unweighted i.i.d. sample) for the Normal distribution 
#' @rdname Normal
#' @name Normal
#' @aliases dNormal pNormal qNormal rNormal eNormal lNormal sNormal iNormal enorm lnorm snorm inorm
#' @details 
#' \itemize{
#'  \item{Density function:}
#'  \deqn{f(x) = \frac 1 {\sqrt{2\pi}\sigma} e^{-((x - \mu)^2/(2 \sigma^2))}}
#'  with \eqn{\mu} the mean of the distribution and \eqn{\sigma} the standard deviation.
#'  \item{Log-likelihood function:}
#'  \deqn{L(\mu,\sigma;X_i)=-\frac 1 2 \ln(2\pi)-\ln(\sigma)-\frac{1}{2\sigma^2}(X_i-\mu)^2}
#'  \item{Score function vector:}
#'  \deqn{V(\mu,\sigma;X_i)
#'  =\left( \begin{array}{c}
#'  \frac{\partial L}{\partial \mu}  \\
#'  \frac{\partial L}{\partial \sigma}
#'  \end{array} \right)
#'  =\left( \begin{array}{c}
#'  \frac {X_i-\mu}{\sigma^2} \\
#'  \frac {(X_i-\mu)^2-\sigma^2}{\sigma^3} 
#'  \end{array} \right)
#'  }
#'  \item{Observed information matrix:}
#'  \deqn{\mathcal J (\mu,\sigma;X_i)=
#'  -\left( \begin{array}{cc}
#'   \frac{\partial^2 L}{\partial \mu^2} & \frac{\partial^2 L}{\partial \mu \partial \sigma} \\
#'   \frac{\partial^2 L}{\partial \sigma \partial \mu} & \frac{\partial^2 L}{\partial \sigma^2} \end{array} \right)
#'  =
#'  \left( \begin{array}{cc}
#'  \frac{1}{\sigma^2} & \frac{2(X_i-\mu)}{\sigma^3} \\
#'  \frac{2(X_i-\mu)}{\sigma^3} & \frac{3(X_i-\mu)^2-\sigma^2}{\sigma^4} \end{array} \right)
#'  }
#' }
#' @param params a list includes all parameters
#' @param x,q vector of quantiles.
#' @param w weights of sample.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param X sample observations.
#' @param mean location parameter.
#' @param sd scale parameter.
#' @param method parameter estimation method.
#' @param logL logical; if TRUE, lNormal gives log likelihood.
#' @param ... other parameters

#' @return dNormal gives the density; pNormal gives the distribution function;
#' qNormal gives the quantile function; rNormal generates random variables; 
#' eNormal estimate the parameters; sNormal gives observed scorn function 

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' # Parameter estimation
#' n <- 500
#' mean <- 1
#' sd <- 2
#' X <- rNormal(n, mean, sd)
#' (est.par <- eNormal(X))
#' 
#' # Histogram and fitted density
#' den.x <- seq(min(X),max(X),length=100)
#' den.y <- dNormal(den.x,mean=est.par$mean,sd=est.par$sd)
#' hist(X, breaks=10, col="red", probability=TRUE, ylim = c(0,1.1*max(den.y)))
#' lines(den.x, den.y, col="blue", lwd=2)
#' 
#' # Q-Q plot and P-P plot
#' plot(qNormal((1:n-0.5)/n, params=est.par), sort(X), main="Q-Q Plot", 
#' xlab="Theoretical Quantiles", ylab="Sample Quantiles", xlim = c(-5,5), ylim = c(-5,5))
#' abline(0,1)
#' 
#' plot((1:n-0.5)/n, pNormal(sort(X), params=est.par), main="P-P Plot", 
#' xlab="Theoretical Percentile", ylab="Sample Percentile", xlim = c(0,1), ylim = c(0,1))
#' abline(0,1)
#' 
#' # A weighted parameter estimation example
#' n <- 10
#' par <- list(mean=1, sd=2)
#' X <- rNormal(n, params=par)
#' w <- c(0.13, 0.06, 0.16, 0.07, 0.2, 0.01, 0.06, 0.09, 0.1, 0.12)
#' eNormal(X,w) # estimated parameters of weighted sample
#' eNormal(X) # estimated parameters of unweighted sample
#' 
#' # Alternative parameter estimation methods
#' eNormal(X, method = "numerical.MLE")
#' eNormal(X, method = "bias.adjusted.MLE")
#' mean(X); sd(X); sd(X)*sqrt((n-1)/n)
#' 
#' # Extracting location or scale parameters
#' est.par[attr(est.par,"par.type")=="location"]
#' est.par[attr(est.par,"par.type")=="scale"]
#' 
#' # evaluate the performance of the parameter estimation function by simulation
#' eval.estimation(rdist=rNormal,edist=eNormal,n = 1000, rep.num = 1e3, params = list(mean=1, sd=2))
#' eval.estimation(rdist=rNormal,edist=eNormal,n = 1000, rep.num = 1e3, params = list(mean=1, sd=2), 
#' method ="analytical.MLE")
#' 
#' # evaluate the precision of estimation by Hessian matrix
#' X <- rNormal(1000, mean, sd)
#' (est.par <- eNormal(X))
#' H <- attr(eNormal(X, method = "numerical.MLE"),"nll.hessian")
#' fisher_info <- solve(H)
#' sqrt(diag(fisher_info))
#' 
#' # log-likelihood, score vector and observed information matrix
#' lNormal(X,param = est.par)
#' lNormal(X,param = est.par, logL=FALSE)
#' sNormal(X,param = est.par)
#' iNormal(X,param = est.par)
#' }

#' @rdname Normal
#' @export dNormal
dNormal <-
  function(x, mean=0, sd =1, params = list(mean, sd)){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    dnorm(x, mean = mean, sd = sd)
  }

#' @rdname Normal
#' @export pNormal
pNormal <- 
  function(q, mean=0, sd =1, params = list(mean, sd)){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    pnorm(q, mean = mean, sd = sd)
  }

#' @rdname Normal
#' @export qNormal
qNormal <- 
  function(p, mean=0, sd =1, params = list(mean, sd)){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    qnorm(p, mean = mean, sd = sd)
  }

#' @rdname Normal
#' @export rNormal
rNormal <- 
  function(n, mean=0, sd =1, params = list(mean, sd)){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    rnorm(n, mean = mean, sd = sd)
  }

#' @rdname Normal
#' @export eNormal
eNormal <-     
  function(X,w, method ="analytical.MLE"){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    if(method == "analytical.MLE") {
      mu <- mean(w*X)
      sigma <- sqrt(mean(w*X^2)-mu^2)
      est.par <- list(mean = mu, sd = sigma)
      est.par.se <- rep(NA, length(est.par))
    } else if(method == "bias.adjusted.MLE"){
      mu <- mean(w*X)
      sigma <- sqrt((mean(w*X^2)-mu^2)*n/(n-1))
      est.par <- list(mean = mu, sd = sigma)
      est.par.se <- rep(NA, length(est.par))
    } else {
      if(method != "numerical.MLE") warning(paste("method ", method, " is not avaliable, use numerial.MLE instead."))  
      method = "numerical.MLE"  
      
      est.par <- wmle(X=X, w=w, distname = "Normal",
                      initial=list(mean = 0, sd = 1),
                      lower=list(mean = -Inf, sd = 0),
                      upper=list(mean = Inf, sd = Inf))
      est.par.se <- sqrt(diag(solve(attr(est.par,"nll.hessian"))))      
    }
    
    attr(est.par,"ob") <- X
    attr(est.par,"weights") <- w
    attr(est.par,"distname") <- "Normal"
    attr(est.par,"method") <- method
    attr(est.par,"par.name") <- c("mean","sd")
    attr(est.par,"par.type") <- c("location","scale")
    attr(est.par,"par.vals") <- c(est.par$mean, est.par$sd)
    attr(est.par,"par.s.e") <-  est.par.se  

    class(est.par) <- "eDist"
    
    return(est.par)
  }

#' @rdname Normal
#' @export enorm
enorm <- eNormal

#' @rdname Normal
#' @export lNormal
## (weighted) (log) likelihood function
lNormal <- 
  function(X, w, mean=0, sd =1, params = list(mean, sd), logL = TRUE){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    ll <- sum(w*(-log(2*pi)/2-log(sd)-(X-mean)^2/(2*sd^2)))
    #     ll <- sum(w*log(dNormal(x=X,params = params)))
    l <- exp(ll)
    
    if(logL) {return(ll)} else{return(l)}
  }

#' @rdname Normal
#' @export sNormal
## (weighted) score vectors
sNormal <- 
  function(X, w, mean=0, sd =1, params = list(mean, sd)){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    score1 <- sum(w*(X-mean)/sd^2)
    score2 <- sum(w*((X-mean)^2-sd^2)/sd^3)
    
    score <- c(score1,score2)
    names(score) <- c("mean","sd")
    return(score)
  }

#' @rdname Normal
#' @export iNormal
## (weighted) (observed) information matrix
iNormal <- 
  function(X, w, mean=0, sd =1, params = list(mean, sd)){
    if(!missing(params)){
      mean <- params$mean
      sd <- params$sd
    }
    
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    } else {
      w <- n*w/sum(w)
    }
    
    info11 <- sum(w*rep(1/sd^2,n))
    info12 <- sum(w*2*(X-mean)/sd^3)
    info21 <- info12
    info22 <- sum(w*(3*(X-mean)^2-sd^2)/sd^4)
    info <- matrix(c(info11,info12,info21,info22), nrow=2,ncol=2)
    rownames(info) <- colnames(info) <- c("mean","sd")
    return(info)
  }

#' @rdname Normal
#' @export snorm
lnorm <- sNormal

#' @rdname Normal
#' @export snorm
snorm <- sNormal

#' @rdname Normal
#' @export inorm
inorm <- iNormal
