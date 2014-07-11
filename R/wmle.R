#' Weighted Maximum Likelihood Estimation 
#' @title Weighted Maximum Likelihood Estimation.
#' @importFrom numDeriv grad
#' @importFrom optimx optimx
#' @importFrom matrixcalc is.positive.definite
#' @description A general weighted maximum likelihood estimation function.
#' @rdname wmle
#' @name wmle
#' @param X observation.
#' @param w frequency (or weights) of observation.
#' @param distname name of distribution to be estimated.
#' @param initial initial value of the parameters.
#' @param lower the lower bound of the parameters.
#' @param upper the upper bound of the parameters.
#' @param loglik.fn function to compute (weighted) log likelihood
#' @param score.fn  function to compute (weighted) score
#' @param obs.info.fn function to compute observed information matrix

#' @return weighted mle estimates.

#' @author Haizhen Wu and A. Jonathan R. Godfrey

#' @examples \donttest{
#' 
#' #if only density function available
#' n <- 200
#' X <- rnorm(n)
#' dFoo <- function(x, params = list(mean=0, sd=1)){
#'   mean <- params$mean
#'   sd <- params$sd
#'   return(dnorm(x, mean = mean, sd = sd))
#' }
#' wmle(X, w=rep(1,n), distname = "Foo",
#'      initial=list(mean = 0, sd = 1),
#'      lower=list(mean = -Inf, sd = 0),
#'      upper=list(mean = Inf, sd = Inf))
#' 
#' #if log likehood function available
#' lFoo2 <- function(X,w, params = list(mean=0, sd=1)){
#'   mean <- params$mean
#'   sd <- params$sd
#'   return(sum(w*log(dnorm(X, mean = mean, sd = sd))))
#' }
#' wmle(X, w=rep(1,n), loglik.fn = lFoo2,
#'      initial=list(mean = 0, sd = 1),
#'      lower=list(mean = -Inf, sd = 0),
#'      upper=list(mean = Inf, sd = Inf))
#' 
#' #if score function available
#' sFoo <- function(X,w, params = list(mean=0, sd=1)){
#'   mean <- params$mean
#'   sd <- params$sd
#'   score1 <- sum(w*(X-mean)/sd^2)
#'   score2 <- sum(w*((X-mean)^2-sd^2)/sd^3)  
#'   score <- c(score1,score2)
#'   return(score)
#' }
#' lFoo <- lFoo2
#' wmle(X, w=rep(1,n), loglik.fn = lFoo, score.fn= sFoo,
#'      initial=list(mean = 0, sd = 1),
#'      lower=list(mean = -Inf, sd = 0),
#'      upper=list(mean = Inf, sd = Inf))
#' 
#' # if score function & observed information matrix available
#' iFoo <- function(X,w, params = list(mean=0, sd=1)){
#'   mean <- params$mean
#'   sd <- params$sd
#'   
#'   n <- length(X)
#'   if(missing(w)){
#'     w <- rep(1,n)
#'   } else {
#'     w <- n*w/sum(w)
#'   }
#'   
#'   info11 <- sum(w*rep(1/sd^2,n))
#'   info12 <- sum(w*2*(X-mean)/sd^3)
#'   info21 <- info12
#'   info22 <- sum(w*(3*(X-mean)^2-sd^2)/sd^4)
#'   info <- matrix(c(info11,info12,info21,info22), nrow=2,ncol=2)
#'   rownames(info) <- colnames(info) <- c("mean","sd")
#'   return(info)
#' }
#' wmle(X, w=rep(1,n), loglik.fn = lFoo, score.fn= sFoo, obs.info.fn=iFoo,
#'      initial=list(mean = 0, sd = 1),
#'      lower=list(mean = -Inf, sd = 0),
#'      upper=list(mean = Inf, sd = Inf))
#' 
#' # Speed comparison 
#' N <- 1000
#' if(exists("sFoo2")) rm(sFoo2)
#' if(exists("iFoo2")) rm(iFoo2)
#' pt <- proc.time()
#' for(i in 1:N){
#'   wmle(X, w=rep(1,n), loglik.fn = lFoo2,
#'        initial=list(mean = 0, sd = 1),
#'        lower=list(mean = -Inf, sd = 0),
#'        upper=list(mean = Inf, sd = Inf))  
#' }
#' proc.time() - pt
#' 
#' sFoo2 <- sFoo
#' if(exists("iFoo2")) rm(iFoo2)
#' pt <- proc.time()
#' for(i in 1:N){
#'   wmle(X, w=rep(1,n), loglik.fn = lFoo2, score.fn= sFoo2,
#'        initial=list(mean = 0, sd = 1),
#'        lower=list(mean = -Inf, sd = 0),
#'        upper=list(mean = Inf, sd = Inf))  
#' }
#' proc.time() - pt
#' 
#' sFoo2 <- sFoo; iFoo2 <- iFoo
#' pt <- proc.time()
#' for(i in 1:N){
#'   wmle(X, w=rep(1,n), loglik.fn = lFoo2, score.fn= sFoo2, obs.info.fn=iFoo2,
#'        initial=list(mean = 0, sd = 1),
#'        lower=list(mean = -Inf, sd = 0),
#'        upper=list(mean = Inf, sd = Inf))  
#' }
#' proc.time() - pt
#' }

#' @export wmle
wmle <- 
  function(X, w, distname, initial, lower, upper, loglik.fn, score.fn, obs.info.fn
  ){
    n <- length(X)
    if(missing(w)){
      w <- rep(1,n)
    }
    w <- n*w/sum(w)
    
    if(!missing(loglik.fn)){
      ll <- function(params) {loglik.fn(X,w,params=as.list(params))}      
    } else {
      if(missing(distname)) {
        stop("Both distribution name and loglikelihood are missing!")
      } else{
        if(exists(paste("l",distname,sep=""), mode = "function")) {
          ldist <- get(paste("l",distname,sep=""))
          ll <- function(params) {ldist(X,w,params=as.list(params))}
        } else {
          ddist <- get(paste("d",distname,sep=""))
          ll <- function(params) {sum(w*log(ddist(x=X,params=as.list(params))))}
        }
      }
    }
    nll <- function(par) {-ll(as.list(par))}
    
    if(!missing(score.fn)){
      gr <- TRUE
      ll.gr <- function(params) {score.fn(X,w,params=params)}
    } else {
      if(missing(distname)) {gr <- FALSE} else{
        if(exists(paste("s",distname,sep=""), mode = "function")) {
          gr <- TRUE
          sdist <- get(paste("s",distname,sep=""))
          ll.gr <- function(params) {sdist(X,w,params=params)}
        } else {
          gr <- FALSE
        }
      }
    } 
    
    if(!missing(obs.info.fn)){
      hess <- TRUE
      ll.hess <- function(params) {-obs.info.fn(X,w,params=params)}
    } else {
      if(missing(distname)) {hess <- FALSE} else{
        if(exists(paste("i",distname,sep=""), mode = "function")) {
          hess <- TRUE
          idist <- get(paste("i",distname,sep=""))
          ll.hess <- function(params) {-idist(X,w,params=params)}
        } else {
          hess <- FALSE
        }
      }
    } 
    
    par.name <- names(initial)
    
    # set rules to transform all bounded intervals to whole real line
    num.par <- max(length(lower),length(upper))
    trans.fn <- vector("list",num.par)
    trans.fn.deriv <- vector("list",num.par)
    trans.fn.dd <- vector("list",num.par)    
    inv.trans.fn <- vector("list",num.par)
    
    for(k in 1:num.par) {
      if(lower[[k]] == -Inf & upper[[k]] == Inf) {
        trans.fn[[k]] <- function(y){y}
        trans.fn.deriv[[k]] <- function(y) {1}
        trans.fn.dd[[k]] <- function(y) {0}
        inv.trans.fn[[k]] <- function(y){y}
      } else if(lower[[k]] != -Inf & upper[[k]] == Inf) {
        trans.fn[[k]] <- function(y){log(y-lower[[k]])}
        body(trans.fn[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn[[k]])[2],fixed=T))[[1]]
        
        trans.fn.deriv[[k]] <- function(y){1/(y-lower[[k]])}
        body(trans.fn.deriv[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn.deriv[[k]])[2],fixed=T))[[1]]
        
        trans.fn.dd[[k]] <- function(y){-1/(lower[[k]]-y)^2}
        body(trans.fn.dd[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn.dd[[k]])[2],fixed=T))[[1]]
        
        inv.trans.fn[[k]] <- function(y){lower[[k]]+exp(y)}
        body(inv.trans.fn[[k]])[[2]] <- parse(text=gsub("k",k,body(inv.trans.fn[[k]])[2],fixed=T))[[1]]
        
      } else if(lower[[k]] == -Inf & upper[[k]] != Inf) {
        trans.fn[[k]] <- function(y){log(upper[[k]]-y)}
        body(trans.fn[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn[[k]])[2],fixed=T))[[1]]
        
        trans.fn.deriv[[k]] <- function(y){-1/(upper[[k]]-y)}
        body(trans.fn.deriv[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn.deriv[[k]])[2],fixed=T))[[1]]
        
        trans.fn.dd[[k]] <- function(y){-1/(upper[[k]]-y)^2}
        body(trans.fn.dd[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn.dd[[k]])[2],fixed=T))[[1]]
        
        inv.trans.fn[[k]] <- function(y){upper[[k]]-exp(y)}
        body(inv.trans.fn[[k]])[[2]] <- parse(text=gsub("k",k,body(inv.trans.fn[[k]])[2],fixed=T))[[1]]
        
      } else if(lower[[k]] != -Inf & upper[[k]] != Inf) {
        trans.fn[[k]] <- function(y){logit((upper[[k]]-y)/(upper[[k]]-lower[[k]]))}
        body(trans.fn[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn[[k]])[2],fixed=T))[[1]]
        
        trans.fn.deriv[[k]] <- function(y){ ifelse(y==upper[[k]] | y==lower[[k]], -Inf, 1/(lower[[k]]-y)+1/(y-upper[[k]]))}
        body(trans.fn.deriv[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn.deriv[[k]])[2],fixed=T))[[1]]
        
        trans.fn.dd[[k]] <- function(y){ ifelse(y==upper[[k]], -Inf, ifelse(y==lower[[k]], Inf, 1/(lower[[k]]-y)^2-1/(y-upper[[k]])^2))}
        body(trans.fn.dd[[k]])[[2]] <- parse(text=gsub("k",k,body(trans.fn.dd[[k]])[2],fixed=T))[[1]]
        
        inv.trans.fn[[k]] <- function(y){upper[[k]]-invlogit(y)*(upper[[k]]-lower[[k]])}
        body(inv.trans.fn[[k]])[[2]] <- parse(text=gsub("k",k,body(inv.trans.fn[[k]])[2],fixed=T))[[1]]
        
      } else {
        stop("the boundary parameters are not set appropriately!")
      }
    }
    
    # trans.initial - transformed initial point
    trans.initial <- initial
    for(k in 1:num.par) {
      trans.initial[[k]] <- trans.fn[[k]](as.numeric(initial[[k]]))
      if(!is.finite(trans.initial[[k]])) trans.initial[[k]] <- 0
    }
    
    # log-likelihood function for transformed parameters    
    trans.ll <- function(trans.arg) { 
      inv.trans.arg <- vector("list",length = num.par)
      for(k in 1:num.par){
        inv.trans.arg[[k]] <- as.numeric(inv.trans.fn[[k]](trans.arg[[k]]))
      }
      names(inv.trans.arg) <- names(initial)
      return(ll(params = as.list(inv.trans.arg)))
    }
    trans.nll <- function(trans.arg) {-trans.ll(trans.arg)}
    
    
    # gradient of ll function for transformed parameters 
    if(gr==TRUE){
      trans.ll.gr <- function(trans.arg) { 
        inv.trans.arg <- vector("list",length = num.par)
        deriv <- NULL
        for(k in 1:num.par){
          inv.trans.arg[[k]] <- as.numeric(inv.trans.fn[[k]](trans.arg[[k]]))
          deriv[k] <- as.numeric(trans.fn.deriv[[k]](inv.trans.arg[[k]]))
        }
        names(inv.trans.arg) <- names(initial)
        return(ll.gr(params <-  as.list(inv.trans.arg))/deriv)
      }
      trans.nll.gr <- function(trans.arg) {-trans.ll.gr(trans.arg)}
      
      # hessian of ll function for transformed parameters 
      if(hess==TRUE){
        trans.ll.hess <- function(trans.arg) { 
          inv.trans.arg <- vector("list",length = num.par)
          deriv <- NULL
          dd <- NULL
          for(k in 1:num.par){
            inv.trans.arg[[k]] <- as.numeric(inv.trans.fn[[k]](trans.arg[[k]]))
            deriv[k] <- as.numeric(trans.fn.deriv[[k]](inv.trans.arg[[k]]))
            dd[k] <- as.numeric(trans.fn.dd[[k]](inv.trans.arg[[k]]))
          }
          names(inv.trans.arg) <- names(initial)
          return( (ll.hess(params <- as.list(inv.trans.arg)) - diag(ll.gr(params <-  as.list(inv.trans.arg))/deriv*dd)) / 
                    deriv%*%t(deriv) )
        }
        trans.nll.hess <- function(trans.arg) {-trans.ll.hess(trans.arg)}
      }
      
    }
    
    approx.est.par <- NULL
    ## default - use nlm
    if(gr == TRUE){
      if(hess == TRUE) {
        f <- function(trans.arg) {
          res <- trans.nll(trans.arg)
          attr(res, "gradient") <- trans.nll.gr(trans.arg)
          attr(res, "hessian") <- trans.nll.hess(trans.arg)
          return(res)
        } 
        #         rst <- nlm(f, p=unlist(trans.initial), hessian=TRUE)
        rst <- try(suppressWarnings(nlm(f, p=unlist(trans.initial), hessian=TRUE)))
      } else {
        f <- function(trans.arg) {
          res <- trans.nll(trans.arg)
          attr(res, "gradient") <- trans.nll.gr(trans.arg)
          return(res)
        } 
        rst <- try(suppressWarnings(nlm(f, p=unlist(trans.initial), hessian=TRUE)))
      }
    } else {
      f <- function(trans.arg) {
        res <- trans.nll(trans.arg)
        return(res)
      } 
      rst <- try(suppressWarnings(nlm(f, p=unlist(trans.initial), hessian=TRUE)))
    }
    
    if(any(class(rst)!="try-error")){      
      if(any(rst$code == c(1,2,3))) {
        trans.est.par <- rst$estimate
        est.par <- NULL
        deriv <- NULL
        dd <- NULL
        for(k in 1:num.par) {
          est.par[k] <- inv.trans.fn[[k]](trans.est.par[[k]])
          deriv[k] <- trans.fn.deriv[[k]](est.par[k])
          dd[k] <- trans.fn.dd[[k]](est.par[k])
        }
        names(est.par) <- par.name
        est.par <- as.list(est.par)
        
        # rst$gradient == trans.nll.gr(trans.est.par)
        # rst$hessian == trans.nll.hess(trans.est.par)
        
        #         hessian <- -rst$hessian * deriv%*%t(deriv)
        nll.gr <- rst$gradient*deriv
        nll.hessian <- diag(rst$gradient*dd) + rst$hessian * (deriv%*%t(deriv))        
        
        colnames(nll.hessian) <- rownames(nll.hessian) <- par.name
        
        attributes(est.par)$nll.hessian <- nll.hessian
        attributes(est.par)$nll.gr <- nll.gr
        attributes(est.par)$nll <- rst$minimum
        
        attributes(est.par)$optim.fn <- "nlm"
        
        if(any(rst$code == c(1,2)) & matrixcalc::is.positive.definite(nll.hessian)) {
          attributes(est.par)$convergence <- 0; return(est.par)
        } else {
          approx.est.par[[1]] <- est.par
        }
      }
    }
    
    ## optim
    hess <- FALSE
    if(gr == TRUE){
      rst <- try(suppressWarnings(optim( par=unlist(trans.initial), fn=trans.nll, gr=trans.nll.gr, hessian = TRUE)))
    } else {
      rst <- try(suppressWarnings(optim( par=unlist(trans.initial), fn=trans.nll, hessian = TRUE)))
    }
    
    if(any(class(rst)!="try-error")){
      if( any(rst$convergence == c(0,1)) ) {
        if(!any(is.na(rst$par),is.na(rst$hessian))) {
          trans.est.par <- rst$par
          est.par <- NULL
          deriv <- NULL
          dd <- NULL
          for(k in 1:num.par) {
            est.par[k] <- inv.trans.fn[[k]](trans.est.par[[k]])
            deriv[k] <- trans.fn.deriv[[k]](est.par[k])
            dd[k] <- trans.fn.dd[[k]](est.par[k])
          }
          names(est.par) <- par.name
          est.par <- as.list(est.par)
          
          if(exists("trans.nll.gr", mode = "function")) {
            rst <- c(rst, list(gradient = trans.nll.gr(trans.est.par)))
          } else {
            rst <- c(rst, list(gradient = c(numDeriv::grad(trans.nll, trans.est.par))))
          }
          
          nll.gr <- rst$gradient*deriv
          names(nll.gr) <- par.name
          
          nll.hessian <- diag(rst$gradient*dd) + rst$hessian * (deriv%*%t(deriv))
          colnames(nll.hessian) <- rownames(nll.hessian) <- par.name
          
          attributes(est.par)$nll.hessian <- nll.hessian
          attributes(est.par)$nll.gr <- nll.gr
          attributes(est.par)$nll <- rst$value
          
          attributes(est.par)$optim.fn <- "optim"
          
          if( (rst$convergence == 0) & matrixcalc::is.positive.definite(nll.hessian)) {
            attributes(est.par)$convergence <- 0; return(est.par)
          } else {
            approx.est.par[[2]] <- est.par
          }      
        }
      }
    }
    
    ## optimx only use ll
    rst <- try(suppressWarnings(optimx::optimx( par=unlist(trans.initial), fn=trans.nll, 
                                                control=list(all.methods=TRUE), hessian = TRUE)))
    
    if(any(class(rst)!="try-error")) {
      rst <- rst[!is.infinite(rst$value),] 
      rst <- rst[rst$convcode==0,]
      
      good.rst.ind <- (rst$kkt1 & rst$kkt2)
      good.rst.ind[is.na(good.rst.ind)] <- FALSE
      
      acceptable.rst.ind <- (rst$kkt1)
      acceptable.rst.ind[is.na(acceptable.rst.ind)] <- FALSE
      
      if(any(sum(good.rst.ind)>0,sum(acceptable.rst.ind)>0)) {
        if(sum(good.rst.ind)>0) {
          good.rst <- rst[good.rst.ind,]
          rst <- good.rst[min(which(good.rst$value == min(good.rst$value))),] # chosse the first with minimum value
        } else if(sum(acceptable.rst.ind)>0) {
          acceptable.rst <- rst[acceptable.rst.ind,]
          rst <- acceptable.rst[min(which(acceptable.rst$value == min(acceptable.rst$value))),] # chosse the first with minimum value
        }
        
        trans.est.par <- coef(rst)
        est.par <- NULL
        deriv <- NULL
        dd <- NULL
        for(k in 1:num.par) {
          est.par[k] <- inv.trans.fn[[k]](trans.est.par[[k]])
          deriv[k] <- trans.fn.deriv[[k]](est.par[k])
          dd[k] <- trans.fn.dd[[k]](est.par[k])
        }
        names(est.par) <- par.name
        est.par <- as.list(est.par)
        
        trans.nll.gr <- attributes(rst)$details[[2]]
        nll.gr <- trans.nll.gr*deriv
        names(nll.gr) <- par.name
        
        trans.nll.hess <- attributes(rst)$details[[3]]
        nll.hessian <- diag(trans.nll.gr*dd) + trans.nll.hess * (deriv%*%t(deriv))
        colnames(nll.hessian) <- rownames(nll.hessian) <- par.name
        
        attributes(est.par)$nll.hessian <- nll.hessian
        attributes(est.par)$nll.gr <- nll.gr
        attributes(est.par)$nll <- rst$value
        
        attributes(est.par)$optim.fn <- paste0("optimx-",attributes(rst)$details[[1]])
        
        if(all(is.finite(nll.hessian))) {
          if(matrixcalc::is.positive.definite(nll.hessian)) {
            attributes(est.par)$convergence <- 0; return(est.par)
          } 
        } else {
          approx.est.par[[3]] <- est.par
        }
      }
    }
    
    nll.values <- NULL
    finite.hess <- NULL
    for(i in 1:length(approx.est.par)){
      if(is.null(approx.est.par[[i]])) {
        nll.values[i] <- NA
        finite.hess[i] <- NA
      } else {
        nll.values[i] = attributes(approx.est.par[[i]])$nll
        finite.hess[i] = all(is.finite(attributes(approx.est.par[[i]])$nll.hessian))
      }
    }
    
    finite.hess[is.na(finite.hess)] <- FALSE
    
    ind <- NULL
    if(sum(finite.hess)>1){
      if(any(!is.na(nll.values[finite.hess]))) {
        ind <- which(finite.hess & (nll.values %in% min(nll.values[finite.hess])))
      }
    }
    
    if(is.integer(ind)) { 
      est.par <- approx.est.par[[ind]]
      attributes(est.par)$convergence <- 1
    } else {
      est.par <- initial; for(k in 1:length(initial)) est.par[[k]] <- NA; 
      attributes(est.par)$convergence <- 20
    }
    
    return(est.par)  
  }      

###############################################################################
## Auxiliary Functions 
###############################################################################
## logit: a mathematical function transforming [0,1] to [-Inf, Inf]
###############################################################################
logit <- 
  function(x){
    log(x/(1-x))
  }
###############################################################################
## invlogit: Inverse of logit function; transform [-Inf, Inf] to [0,1]
###############################################################################
invlogit <- 
  function(x){
    1/(1+exp(-x))
  }
