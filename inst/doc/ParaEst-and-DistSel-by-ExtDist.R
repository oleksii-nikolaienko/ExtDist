## ----, message=FALSE-----------------------------------------------------
require(ExtDist)

## ----, echo=FALSE, message=FALSE-----------------------------------------
set.seed(1234)
X <- rWeibull(50, shape = 2, scale = 3)

## ------------------------------------------------------------------------
head(X)

## ------------------------------------------------------------------------
(est.par <- eWeibull(X))

## ------------------------------------------------------------------------
class(est.par)

## ------------------------------------------------------------------------
dWeibull(seq(0,2,0.4), params = est.par)
pWeibull(seq(0,2,0.4), params = est.par)
qWeibull(seq(0,1,0.2), params = est.par)
rWeibull(10, params = est.par)

## ----, results='hold'----------------------------------------------------
dWeibull(seq(0,2,0.4), shape = est.par$shape, scale = est.par$scale)
pWeibull(seq(0,2,0.4), shape = est.par$shape, scale = est.par$scale)
qWeibull(seq(0,1,0.2), shape = est.par$shape, scale = est.par$scale)
rWeibull(10, shape = est.par$shape, scale = est.par$scale)

## ------------------------------------------------------------------------
fit_Dist <- function(X, Dist){
  l <- min(X); u <- max(X); d <- u-l; n <- length(X)
  
  est.par <- get(paste0("e",Dist))(X)
  dDist <- function(X) get(paste0("d",Dist))(X,param = est.par)
  pDist <- function(X) get(paste0("p",Dist))(X,param = est.par)
  qDist <- function(X) get(paste0("q",Dist))(X,param = est.par)

  op <- par(mfrow=c(2,2)) 
  PerformanceAnalytics::textplot(capture.output(print(est.par)), valign = "top")

  hist(X, col="red", probability=TRUE, xlim=c(l-0.1*d,u+0.1*d))
  curve(dDist, add=TRUE, col="blue", lwd=2)
  
  plot(qDist((1:n-0.5)/n), sort(X), main="Q-Q Plot", xlim = c(l,u), ylim = c(l,u), 
       xlab="Theoretical Quantiles", ylab="Sample Quantiles")
  abline(0,1)

  plot((1:n-0.5)/n, pDist(sort(X)), main="P-P Plot", xlim = c(0,1), ylim = c(0,1),
       xlab="Theoretical Percentile", ylab="Sample Percentile")
  abline(0,1)
  
  par(op)
}

## ------------------------------------------------------------------------
X <- rBeta(100,2,5)
fit_Dist(X, "Beta")

## ------------------------------------------------------------------------
logLik(est.par) # log likihood
AIC(est.par) # Akaike information criterion
AICc(est.par) # corrected Akaike information criterion
BIC(est.par) # Bayesian Information Criterion. 
MDL(est.par) # minimum description length 
vcov(est.par) # variance-covariance matrix of the parameters of the fitted distribution

## ------------------------------------------------------------------------
set.seed(1234)
X <- rBeta(50, shape1 = 2, shape2 = 10 )
bestDist(X, candDist = c("Beta_ab","Laplace","Normal"), criterion = "AIC")

## ------------------------------------------------------------------------
set.seed(1234)
X <- rBeta(50, shape1 = 2, shape2 = 10 )
DistSelCriteriaValues(X, candDist = c("Beta_ab","Laplace","Normal"),
                         criteria = c("logLik","AIC","AICc","BIC","MDL"))

## ------------------------------------------------------------------------
Y <- c(0.1703, 0.4307, 0.6085, 0.0503, 0.4625, 0.479, 0.2695, 0.2744, 0.2713, 0.2177, 
       0.2865, 0.2009, 0.2359, 0.3877, 0.5799, 0.3537, 0.2805, 0.2144, 0.2261, 0.4016)
w <- c(0.85, 1.11, 0.88, 1.34, 1.01, 0.96, 0.86, 1.34, 0.87, 1.34, 0.84, 0.84, 0.83, 1.09, 
       0.95, 0.77, 0.96, 1.24, 0.78, 1.12)

## ------------------------------------------------------------------------
eBeta(Y,w)

bestDist(Y, w, candDist = c("Beta_ab","Laplace","Normal"), criterion = "AIC")

DistSelCriteriaValues(Y, w, candDist = c("Beta_ab","Laplace","Normal"),
                         criteria = c("logLik","AIC","AICc","BIC","MDL"))

