#' Multi-Core Version von Model Averaging Funktionen aus cAIC4-Package
#'
#' Unterschied zur originalen Version, ist dass hier Multi-Core Support besteht.
#'
#'
#' @param models liste an kandidaten modellen fürs model averaging
#' @param opt Parameter ob die Augmented Lagrangian Optimierung genutzt werden soll. TRUE = JA.
#'
#' @author Rene-Marcel Kruse
#' @seealso \code{\link[cAIC4]{cAIC4-package}}
#' @references Saefken, B., Ruegamer, D., Kneib, T. and Greven, S. (2018):
#' Conditional Model Selection in Mixed-Effects Models with cAIC4.
#' \url{https://arxiv.org/abs/1803.05664}
#' @references Greven, S. and Kneib T. (2010) On the behaviour of marginal and
#' conditional AIC in linear mixed models. Biometrika 97(4), 773-789.
#' @keywords multicore, cAIC, model averaging
#' @rdname MC_modelAvg
#' @export MC_modelAvg
#' @importFrom lme4 getME ranef
#'
MC_modelAvg <- function (models, opt = TRUE){
  call <- match.call()
  if (opt == TRUE) {
    tempres <- MC_getWeights(models)
  }
  else {
    invisible(capture.output(tempres <- MC_anocAIC(models)))
    tempres$delta <- tempres$cAIC - min(tempres$cAIC)
    tempres$weights <- exp(-tempres$delta/2)/sum(exp(-tempres$delta/2))
  }
  betas <- list()
  for (i in 1:length(models)) {
    betas[[i]] <- getME(models[[i]], "fixef")
  }
  avg.betas <- list()
  for (i in 1:length(models)) {
    avg.betas[[i]] <- betas[[i]] * tempres$weights[i]
  }
  sum.avg.betas <- tapply((unlist(avg.betas)), names(unlist(avg.betas)),
                          FUN = sum)
  rand <- list()
  for (i in 1:length(models)) {
    rand[[i]] <- ranef(models[[i]])
  }
  avg.rand <- list()
  for (i in 1:length(models)) {
    dummy <- unlist(rand[[i]])
    avg.rand[[i]] <- dummy * tempres$weight[i]
  }
  sum.avg.rand <- tapply((unlist(avg.rand)), names(unlist(avg.rand)),
                         FUN = sum)
  res <- list(call = call, fixeff = sum.avg.betas, raneff = sum.avg.rand,
              optimresults = tempres, candidatmodels = models)
}
