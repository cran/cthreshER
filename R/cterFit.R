#' Fit the continuous threshold expectile regression
#'
#' The grid search algorithm for the continuous threshold expectile regression
#'
#' @rdname cterFit
#' @param y A vector of response
#' @param x A scalar covariate with threshold
#' @param z A vector of covariates
#' @param tau the expectile level, 0.5 for default
#' @param tol  tolerance value, 1e-4 for default
#' @param max.iter the maximum iteration steps, 100 for default
#'
#' @return A list with the elements
#' \item{coef.est}{The estimated regression coefficients with intercept.}
#' \item{threshold.est}{The estimated threshold.}
#' \item{coef.se}{The estimated standard error of the regression coefficients.}
#' \item{threshold.se}{The estimated standard error of the threshold.}
#' \item{iter}{The iteration steps.}
#'
#' @author Feipeng Zhang and Qunhua Li
#' @keywords cterFit

#' @importFrom stats approxfun density lsfit sd
#' @importFrom Matrix nearPD
#'
#' @export
#'
#' @examples
#'
#'
#' ## simulated data
#' ptm <- proc.time()
#' n <- 200
#' t0 <- 1.5
#' bet0 <- c(1, 3, -2, 1)
#' tau <- 0.3
#' modtype <- 1
#' errtype <- 1
#' dat <- cterSimData(n, bet0, t0, tau, modtype, errtype)
#' y <- dat[, 1]
#' x <- dat[, 2]
#' z <- dat[, 3]
#' fit <- cterFit(y, x, z, tau)
#'
#' ## The example of Baseball pitcher salary
#' data(data_bbsalaries)
#' y <- data_bbsalaries$y
#' x <- data_bbsalaries$x
#' z <- NULL
#' tau <- 0.5
#' fit <- cterFit(y, x, z, tau)
#' proc.time() - ptm





cterFit <- function(y, x, z,  tau = 0.5, max.iter = 100, tol = 1E-4){

  ##############################################
  ##  some useful function
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
    list(coefficients = b, variance = Vb, loss = loss, it = it)
  }

  ### grid search method for estimate the coefficients
  expTLgrid <- function(y, x, z, tau, max.iter, tol){

    tt <- seq(min(x) + 0.01, max(x) - 0.01, length = 100)
    sst <- it <- NULL
    p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
    bet <- matrix(0, length(tt), p+3)

    for (k in 1:length(tt)){
      
      X <- cbind(1, x, pmax(x-tt[k],0), z)  # with intercept

      fit <- expTLlaws(X, y, tau, max.iter, tol)
      bet[k, ] <- fit$coefficients
      sst[k] <- fit$loss
      it[k] <- fit$it
    }

    bp.grid <- min(tt[sst== min(sst)])
    bet.grid <- bet[tt==bp.grid, ]
    it.grid <- it[tt==bp.grid]

    return(list(bet = bet.grid, bp = bp.grid, it = it.grid))
  }


  ### the estimated standard deviation
  wfun <- function(u, tau){
    abs(tau- (u <= 0))
  }

  # Epanechnikov kernel
  kerfun <- function(u){
    3/4*(1 - u^2)*(abs(u) <= 1)
  }

  expTLstd <- function(y, x, z, bet, t, tau){

    n <- length(y)
    Vt <- cbind(1, x, pmax(x-t, 0), z)  # the covariates
    res <- as.vector(y - Vt %*% bet)
    wt <- wfun(res, tau)
    bet2 <- bet[3]


    ## estimation for density value at t by Silverman's rule of thumb
    h <- 1.06* n^(-1/5)* sd(x)
    pdf<- approxfun(density(x, kernel= "epanechnikov", bw = h))
    ft <- pdf(t)



    Hn11 <- 2/n * t(Vt) %*% diag(wt) %*% Vt
    zmat0 <- 0*z
    Hn12 <- 2/n * apply(-Vt*(wt* bet2*(x>t)) +
                          cbind(0, 0, (x>t), zmat0)*(wt*res), 2, sum)

    Hn1 <- cbind(Hn11, Hn12)
    Hn21 <- t(Hn12)
    Hn22 <- 2/n * sum( wt* (bet2*(x>t))^2)  -
      2/n * sum( wt*res*bet2 * ft)
    Hn2 <- cbind(Hn21, Hn22)
    Hn <- rbind(Hn1, Hn2)

    Gn <- cbind(-2* Vt * res * wt,
                2*wt*res*bet2*(x>t))
    Sig <- 1/n * t(Gn) %*%  Gn
    p <- dim(Sig)[1]
    A <- solve(Hn + diag(1e-8, p)) %*% Sig %*% solve(Hn + diag(1e-8, p))
    thet.sd <- sqrt(diag(A)/n)
    return(thet.sd)
  }

  ##############################################
  ## the main fit
  fit <- expTLgrid(y, x, z, tau, max.iter, tol)
  coef.est <- fit$bet
  threshold.est <- fit$bp
  it <- fit$it
  ese <- expTLstd(y, x, z, fit$bet,  fit$bp, tau)
  coef.se <- ese[-length(ese)]
  threhold.se <- ese[length(ese)]

  return(list(coef.est = coef.est, threshold.est = threshold.est,
              coef.se = coef.se, threhold.se = threhold.se,
              iter = it))


}# end function



