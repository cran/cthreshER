#'
#' test the existence of change point in the continuous threshold expectile regression
#'
#' This function for calculating the test statistics and p-value by wild bootstrap.
#'
#' @rdname cterTest
#' @param y A vector of response
#' @param x A scalar covariate with threshold
#' @param z A vector of covariates
#' @param tau the expectile level, 0.5 for default
#' @param NB resampling times, 1000 for default
#'
#' @return A list with the elements
#' \item{Tn}{The statistic based on original data.}
#' \item{Tn.NB}{The statistics by wild bootstrap.}
#' \item{p.value}{The p-value by wild bootstrap.}
#'
#' @author Feipeng Zhang and Qunhua Li
#' @keywords cterTest
#'
#' @importFrom stats rnorm lsfit sd
#' @importFrom Matrix nearPD
#'
#' @export
#'
#' @examples
#'
#' ## simulated data
#' ptm <- proc.time()
#' set.seed(1)
#' n <- 200
#' t0 <- 1.5
#' bet0 <- c(1, 3, 0, 1)
#' tau <- 0.3
#' modtype <- 1
#' errtype <- 1
#' dat <- cterSimData(n, bet0, t0, tau, modtype, errtype)
#' y <- dat[, 1]
#' x <- dat[, 2]
#' z <- dat[, 3]
#' fit.test <- cterTest(y, x, z, tau, NB = 30)
#' fit.test$p.value
#' 
#' ## The example of Baseball pitcher salary
#' data(data_bbsalaries)
#' y <- data_bbsalaries$y
#' x <- data_bbsalaries$x
#' z <- NULL
#' tau <- 0.5
#' fit.test <- cterTest(y, x, z, tau, NB = 30)
#' fit.test$p.value
#' proc.time() - ptm



cterTest <- function(y, x, z, tau = 0.5, NB = 1000){

  ## global variable
  max.iter = 100
  tol = 1e-5

  expTLlaws <- function(X, y, tau,  max.iter, tol)
  {
    # X is the model matrix
    # y is the response vector of observed proportion
    # maxIter is the maximum number of iterations
    # tol is a convergence criterion
    #X <- cbind(1, X) # add constant
    b <- bLast <- rep(0, ncol(X)) # initialize
    it <- 1 # iteration index
    while (it <= max.iter){

      ypred <- c(X %*% b)
      w <- as.vector(tau *(y>= ypred) + (1-tau)* (y<ypred))
      b <- lsfit(X, y, w, intercept=FALSE)$coef
      if (max(abs(b - bLast)/(abs(bLast) + 0.01*tol)) < tol) break
      bLast <- b
      it <- it + 1 # increment index
    }
    if (it > max.iter) warning('maximum iterations exceeded')

    ## the loss function
    loss <- sum(w*(y-c(X%*%b))^2)

    ## the variance
    Am <- t(X) %*% diag(w) %*% X
    if(min(eigen(Am)$values)<1e-8){Am=as.matrix(nearPD(Am)$mat)}

    A <- solve(Am)
    H <- w * diag(X %*% A %*% t(X))
    B <- t(X) %*% diag(w^2*(y-ypred)^2/(1-H)) %*% X
    Vb <- A %*% B %*% A  # variance function
    list(coefficients=b, variance= Vb, loss=loss, iterations=it)
  }


  wfun <- function(u, tau){
    abs(tau- (u <= 0))
  }

  ### test statistic based on origin data
  testFun <- function(y, x, z, tau, tt){

    ## under H0
    ## by asymmetric weight least square
    X <- cbind(1, x, z)
    fit <- expTLlaws(X, y, tau, max.iter, tol)
    bet <- fit$coefficients
    res <- as.vector(y - X %*% bet)
    wt <- wfun(res, tau)

    n <- length(y)

    Rn <-  rep(0, length(tt))
    for (kk in 1:length(tt)){

      Rn[kk] <- 1/sqrt(n)*sum(
        wt*res*(x-tt[kk])*ifelse(x<=tt[kk], 1, 0)
      )

    }

    Tn <- max(abs(Rn))

    return(Tn)
  }

  ## perturbed method to calculate the p-value
  testFun.resample <- function(y, x, z, tau, tt){

    #########################
    ## permutation random errors
    n<- length(y)
    u <- rnorm(n, 0, 1)

    #########################

    X <- cbind(1, x, z)
    fit <- expTLlaws(X, y, tau, max.iter, tol)
    bet <- fit$coefficients
    res <- as.vector(y - X %*% bet)
    wt <- wfun(res, tau)

    ## Sn
    xz <-  cbind(1, x, z)
    Sn <- (t(xz) %*% diag(wt) %*% xz)/n


    ## under H0
    Rn <-  rep(0, length(tt))
    for (kk in 1:length(tt)){

      Sn.t <- apply(xz*(wt*(x-tt[kk])*ifelse(x <= tt[kk], 1, 0)), 2, mean)
      Rn[kk] <- 1/sqrt(n)*sum(
        u * wt *res *((x-tt[kk])*ifelse(x<= tt[kk],1,0) -
                        xz %*% solve(Sn) %*% Sn.t)
      )


    }
    Tn <- max(abs(Rn))

    return(Tn)
  }


  #######################################################
  ###  calculate the p-value by wild bootstrap

  tt <- seq(min(x), max(x), length = 100)
  Tn <-  testFun(y, x, z, tau, tt)
  Tn.NB <- replicate(NB, testFun.resample(y, x, z, tau, tt))

  pv <- mean(Tn.NB > Tn,  na.rm = TRUE)

  return(list(Tn = Tn, Tn.NB = Tn.NB, p.value = pv))


}

