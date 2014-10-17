# Best distribution for (weighted) sample ---------------------------------
#' @title Best distribution for (weighted) sample.
#' @description A function to choose the best fitted distibution based on 
#' specified criteria.
#' @rdname bestDist
#' @name bestDist
#' @details Details

#' @param X obersevations.
#' @param w weights of sample.
#' @param candDist a vector of names of candidate distributions.
#' @param criterion the criterion based on which the best fitted distribution is
#'  choosen.

#' @return the name of best fitted distribution and the parameter estimates.
#' @export bestDist

#' @examples \donttest{
#' X <- rBeta_ab(30, a = 0, b = 1, shape1 = 2, shape2 = 10 )
#' bestDist(X, candDist = c("Beta_ab","Laplace","Normal"), criterion = "logLik")
#' bestDist(X, candDist = c("Beta_ab","Laplace","Normal"), criterion = "AIC")
#' bestDist(X, candDist = c("Beta_ab","Laplace","Normal"), criterion = "AICc")
#' bestDist(X, candDist = c("Beta_ab","Laplace","Normal"), criterion = "BIC")
#' bestDist(X, candDist = c("Beta_ab","Laplace","Normal"), criterion = "MDL")
#' 
#' w <- c(0.32, 1.77, 1.22, 0.64, 0.38, 0.93, 1.63, 1.34, 0.57, 1.73, 1.67, 0.67, 
#' 0.09, 1, 1.55, 0.53, 0.76, 1.06, 1.13, 1.31, 1.18, 1.64, 0.07, 1.41, 1.18, 
#' 0.69, 0.28, 1.27, 0.9, 1.08)
#' bestDist(X, w, candDist = c("Beta_ab","Laplace","Normal"), criterion = "logLik")
#' bestDist(X, w, candDist = c("Beta_ab","Laplace","Normal"), criterion = "AIC")
#' bestDist(X, w, candDist = c("Beta_ab","Laplace","Normal"), criterion = "AICc")
#' bestDist(X, w, candDist = c("Beta_ab","Laplace","Normal"), criterion = "BIC")
#' bestDist(X, w, candDist = c("Beta_ab","Laplace","Normal"), criterion = "MDL")
#' 
#' # parameter for best distribution
#' best_dist <- bestDist(X, candDist = c("Beta_ab","Laplace","Normal"), criterion = "logLik")
#' attributes(best_dist)$best.dist.par
#' }

bestDist <- 
  function(X,
           w = rep(1,length(X))/length(X),
           candDist = c("Beta_ab","Laplace","Normal"),
           criterion = "AICc"
  ){
    if(!(criterion %in% c("logLik","AIC","AICc","BIC","MDL"))) {return("criterion unknown")}
    
    w <- w/sum(w)*length(X)
    
    est.pars <- lapply(paste0("e",candDist), do.call, args=list(X,w)) 
    
    criterion.value <- unlist(lapply(est.pars, getS3method(criterion, class = "eDist")))
    names(criterion.value) <- candDist
    
    if(criterion %in% c("logLik")){
      new.order <- order(criterion.value,decreasing=T)
      best <- candDist[new.order][1]
      est.par <- est.pars[new.order][[1]]
    } else {
      new.order <- order(criterion.value,decreasing=F)
      best <- candDist[new.order][1]
      est.par <- est.pars[new.order][[1]]
    }
    
    attributes(best)$best.dist.par <- est.par
    attributes(best)$criterion.value <- criterion.value
    
    return(best)
  }

