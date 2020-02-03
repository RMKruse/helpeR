#' Optimize weights for model averaging.
#'
#' Function to constructed an optimal vector of weights for model averaging of
#' Linear Mixed Models based on the proposal of Zhang et al. (2014) of using Stein's Formular
#' to derive a suitable criterion based on the conditional Akaike Information Criterion as
#' proposed by Greven and Kneib. The underlying optimization used is a customized version
#' of the Augmented Lagrangian Method.
#'
#' @param models An list object containing all considered candidate models fitted by
#' \code{\link[lme4]{lmer}} of the lme4-package or of class
#' \code{\link[nlme]{lme}}.
#' @return An updated object containing a vector of weights for the underlying candidate models, value
#' of the object given said weights as well as the time needed.
#' @author Rene-Marcel Kruse
#' @keywords weights, model averaging, augmented lagrangian method
#' @rdname MC_getWeights
#' @export MC_getWeights
#' @importFrom lme4 getME
#' @import parallel
#' @importFrom utils capture.output
#'
MC_getWeights <- function(models)
{
  m             <- models
  .envi         <- environment()
  # Creation of the variables required for the optimization of weights
  # TODO: Suppress anocAIC's automatic output
  modelcAIC     <- MC_anocAIC(m)
  df            <- modelcAIC[[2]]
  tempm         <- m[[which.max(modelcAIC$df)]]
  seDF          <- getME(tempm, "sigma")
  varDF         <- seDF * seDF
  y             <- getME(m[[1]], "y")
  mu            <- list()
  # TODO: that needs to be made more effective ...
  for(i in 1:length(m)){
    mu[[i]]     <- fitted(m[[i]])
  }
  mu            <- t(matrix(unlist(mu), nrow = length(m), byrow = TRUE))
  weights       <- rep(1/length(m), times = length(m))
  fun           <- find_weights <- function(w){
    (norm(y - matrix(mu %*% w)))^(2) + 2 * varDF * (w %*% df)}
  eqfun         <- equal <-function(w){sum(w)}
  equB           <- 1
  lowb          <- rep(0, times = length(m))
  uppb          <- rep(1, times = length(m))
  nw            <- length(weights)
  funv 	        <- find_weights(weights)
  eqv 	        <- (sum(weights)-equB)
  rho           <- 0
  maxit         <- 400
  minit         <- 800
  delta         <- (1.0e-7)
  tol           <- (1.0e-8)
  # Start of optimization:
  j             <- jh <- funv
  lambda        <- c(0)
  constraint    <- eqv
  p             <- c(weights)
  hess          <- diag(nw)
  mue           <- nw
  .iters        <- 0
  targets            <- c(funv, eqv)
  tic           <- Sys.time()
  while( .iters < maxit ){
    .iters <- .iters + 1
    scaler <- c( targets[ 1 ], rep(1, 1) * max( abs(targets[ 2:(1  + 1) ]) ) )
    scaler <- c(scaler, rep( 1, length.out = length(p) ) )
    scaler <- apply( matrix(scaler, ncol = 1), 1,
                     FUN = function( x ) min( max( abs(x), tol ), 1/tol ) )
    res    <- .MC_weightOptim(weights = p, lm = lambda, targets = targets,
                           hess = hess, lambda = mue, scaler = scaler, .envi = .envi)
    p      <- res$p
    lambda <- res$y
    hess   <- res$hess
    mue    <- res$lambda
    temp   <- p
    funv 	 <- find_weights(temp)
    eqv 	 <- (sum(temp)-equB)
    targets     <- c(funv, eqv)
    # Change of the target function through optimization
    tt     <- (j - targets[ 1 ]) / max(targets[ 1 ], 1)
    j      <- targets[ 1 ]
    constraint <- targets[ 2 ]
    if( abs(constraint) < 10 * tol ) {
      rho  <- 0
      mue  <- min(mue, tol)
    }
    if( c( tol + tt ) <= 0 ) {
      lambda <- 0
      hess  <- diag( diag ( hess ) )
    }
    if( sqrt(sum( (c(tt, eqv))^2 )) <= tol ) {
      maxit <- .iters
    }
  }
  toc <- Sys.time() - tic
  ans <- list("weights" = p, "functionvalue" = j,
              "duration" = toc)
  return( ans )
}
