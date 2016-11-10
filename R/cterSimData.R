#' Simulated data from the continuous threshold expectile regression
#'
#' The function for  simulating data from the continuous threshold expectile regression
#'
#' @rdname cterSimData
#' @param n sample size.
#' @param bet0 the vecotr of true regression coefficients.
#' @param t0 the true location of threshold.
#' @param tau the expectile level, 0.5 for default.
#' @param modtype type of model, 1 = IID for default, 2 = Heteroscedasticity,
#'  modtype = 1, \eqn{Y = beta_0 + beta_1 X + beta_2 (X-t)_+ gamma Z + e,}
#'  modtype = 1, \eqn{Y = beta_0 + beta_1 X + beta_2 (X-t)_+ gamma Z + (1+0.2Z)e,}
#   where X from Uniform(-2, 4), Z from N(1, 0.5^2).
#' @param errtype type of error, 1 for default,
#'  errtype = 1 for N(0, 1),
#'  errtype = 2 for t_4,
#'  errtype = 3 for 0.9 N(0, 1) + 0.1 t_4.
#'
#' @return A matrix with the elements
#' \item{y}{The response variable.}
#' \item{x}{The scalar covariate with threshold.}
#' \item{z}{A vector of covariates.}
#'
#' @author Feipeng Zhang and Qunhua Li
#' @keywords cterSimData

#' @importFrom stats integrate
#' @importFrom stats runif rnorm pnorm dnorm rt pt dt rbinom
#'
#' @export
#'
#' @examples
#' 
#' ## simulated data
#' ptm <- proc.time()
#' n <- 200
#' t0 <- 1.5
#' bet0 <- c(1, 3, -2, 1)
#' tau <- 0.5
#' modtype <- 1
#' errtype <- 1
#' dat <- cterSimData(n, bet0, t0, tau, modtype, errtype)
#' head(dat)
#' proc.time() - ptm

### generate the data
cterSimData <- function(n, bet0, t0, tau = 0.5, modtype = 1, errtype = 1)
  {

  ### calculate the true expectile
  pefun <- function(t, errtype){

    if (errtype==1){
      ## example1: normal distribution
      F <- pnorm(t, 0, 1)

      integrand <- function(x) {x*dnorm(x, 0, 1)}

    } else if (errtype==2){
      ## example2: t4 distribution
      F <- pt(t, 4)

      integrand <- function(x) {x*dt(x, 4)}

    } else if (errtype==3) {
      ## example3: mixtrue distribution of normal distributions
      prop <- 0.1  # mixture proportion
      F <- prop * pt(t, 4) + (1-prop)*pnorm(t, 0, 1)

      integrand <- function(x){
        x*(prop * dt(x, 4) +  (1-prop)*dnorm(x, 0, 1))
      }
    }

    G <- integrate(integrand, lower = -(1e+3), upper = t)$value
    gmean <-  integrate(integrand, lower = -(1e+3), upper = 1e+3)$value
    u <- G -t * F
    asy <- u/(2*u + t-gmean)

    return(asy)
  }

  efun <- function (tau, errtype){
    tau[tau > 1 | tau < 0] = NA
    zz = 0 * tau
    lower = rep(-10, length(tau))
    upper = rep(10, length(tau))
    diff = 1
    index = 1
    while (diff > 1e-10 && index < 1000) {
      root = pefun(zz, errtype) - tau
      root[is.na(root)] = 0
      lower[root < 0] = zz[root < 0]
      upper[root > 0] = zz[root > 0]
      zz = (upper + lower)/2
      diff = max(abs(root), na.rm = T)
      index = index + 1
    }
    zz[is.na(tau)] = NA
    return(zz)
  }


  ### parameter settings
  x <- runif(n, -2, 4)  # threshold scale variable
  z <- rnorm(n, 1, 0.5)  # covariates
  xz <- cbind(1, x, pmax(x-t0, 0), z)


  if (errtype==1){
    ## example1: normal distribution
    err <- rnorm(n, 0, 1)    # normal distribution
  } else if (errtype==2){
    ## example2: t4 distribution
    err <- rt(n, 4)
  } else if (errtype==3) {
    ## example3: mixtrue distribution of normal distributions
    prop <- 0.1  # mixture proportion of t4 distribution
    B <- rbinom(n, 1, prop)
    err <-  (B==1)* rt(n, 4) +  (B==0)*rnorm(n, 0, 1)
  }
  err0 <- err - efun(tau, errtype)   # the tau-th expectile of err is zero

  if (modtype == 1){
    y <- xz %*% bet0 + err0
  } else if (modtype == 2){
    y <- xz %*% bet0 + (1 + 0.2*z)*err0
  }

  dat <- cbind(y, x, z)

  return(dat)
}

