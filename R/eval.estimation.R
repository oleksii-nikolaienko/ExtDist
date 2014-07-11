#' Parameter Estimation Evaluation
#'
#' @title Parameter Estimation Evaluation.
#' @description A function to evaluate the parameter estimation function.
#' @rdname eval.estimation
#' @name eval.estimation
#' @details Details
#' @param rdist random variable generating function.
#' @param edist parameter estimation function.
#' @param n sample size
#' @param rep.num number of replacates 
#' @param params true parameters of the distribution.
#' @param method esimation method

#' @return a list containing the mean, sd of the estimated parameters. na.cont 
#' is the number of "na"s appeared in the parameter estimation.

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' eval.estimation(rdist=rNormal,edist=eNormal,n = 10, rep.num = 1e3, 
#' params = list(mean = 1,sd = 5))
#' eval.estimation(rdist=rNormal,edist=eNormal,n = 1000, rep.num = 1e3, 
#' params = list(mean = 1,sd = 5))
#' }

#' @export eval.estimation
eval.estimation <- 
  function(rdist, edist, n = 20, rep.num = 1e3, params, method = "numerical.MLE"){
    k <- length(params)
    
    est.par <- array(NA, dim = c(rep.num,k))
    start.time <- proc.time()
    for(i in 1:rep.num){
      X <- rdist(n=n,params=params)
      est.par[i,] <- as.numeric(unlist(edist(X,method=method)))
      #       print(paste("i=",i))
    }
    return(list(method = method,
                est.mean = apply(est.par,2, mean, na.rm =T), 
                est.sd = apply(est.par,2, sd, na.rm =T), 
                time = proc.time() - start.time,
                na.cont = sum(is.na(est.par[,1])))
    )
  }
