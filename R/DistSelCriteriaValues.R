# Distribution Selection Criteria Values ----------------------------------
#' @title Distribution Selection Criteria Values.
#' @description A function to calculate the distribution selection criteria 
#' values for a list of candidate distributions
#' @rdname DistSelCriteriaValues
#' @name DistSelCriteriaValues
#' @aliases DistSelCriteriaValues
#' @details Details

#' @param X obersevations.
#' @param w weights of sample.
#' @param candDist a vector of names of candidate distributions.
#' @param criteria a vector of criteria to be calculated

#' @return a table contains mle parameter estimates and corresponding 
#' distribution selection criteria values.
#' @export DistSelCriteriaValues

#' @examples \donttest{
#' X <- rBeta_ab(30, a = 0, b = 1, shape1 = 2, shape2 = 10 )
#' DistSelCriteriaValues(X, candDist = c("Beta_ab","Laplace","Normal"), 
#'                          criteria = c("logLik","AIC","AICc","BIC","MDL"))
#' 
#' w <- c(0.32, 1.77, 1.22, 0.64, 0.38, 0.93, 1.63, 1.34, 0.57, 1.73, 1.67, 0.67, 
#' 0.09, 1, 1.55, 0.53, 0.76, 1.06, 1.13, 1.31, 1.18, 1.64, 0.07, 1.41, 1.18, 
#' 0.69, 0.28, 1.27, 0.9, 1.08)
#' DistSelCriteriaValues(X, w, candDist = c("Beta_ab","Laplace","Normal"), 
#'                              criteria = c("logLik","AIC","AICc","BIC","MDL"))
#' }

#' @name DistSelCriteriaValues
#' @export DistSelCriteriaValues
DistSelCriteriaValues <- 
  function(X,
           w = rep(1,length(X))/length(X),
           candDist = c("Beta_ab","Laplace","Normal"),
           criteria = c("logLik","AIC","AICc","BIC","MDL")
  ){
    if(!(all(criteria %in% c("logLik","AIC","AICc","BIC","MDL")))) {return("some criteria unknown")}

    
    w <- w/sum(w)*length(X)   
    est.pars <- lapply(paste0("e",candDist), do.call, args=list(X,w)) 

    
    function.list <- lapply(criteria, getS3method, class = "eDist")
    fn <- function (est.par) (lapply(X=function.list, FUN=do.call, args=list(est.par)))
    criteria.values <- sapply(est.pars, fn )
    colnames(criteria.values) <- candDist
    rownames(criteria.values) <- criteria
    
    return(criteria.values)
  }

