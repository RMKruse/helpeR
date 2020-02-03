#' Multi-Core Version von Comparison of several Models via cAIC
#'
#' Unterschied zur originalen Version, ist dass hier Multi-Core Support besteht.
#' Desweiteren ist es möglich mit dieser Version auch lm-Objekte in der Liste an
#' Kandidaten Modellen zu nutzen.
#'
#' @param objects list() an Kandidaten Modelle
#' @param digits Anzahl der auszugebene Nachkommastelle
#'
#' @author Rene-Marcel Kruse
#' @seealso \code{\link[cAIC4]{cAIC4-package}}, \code{\link[lme4]{lmer}}
#' @references Saefken, B., Ruegamer, D., Kneib, T. and Greven, S. (2018):
#' Conditional Model Selection in Mixed-Effects Models with cAIC4.
#' \url{https://arxiv.org/abs/1803.05664}
#' @references Greven, S. and Kneib T. (2010) On the behaviour of marginal and
#' conditional AIC in linear mixed models. Biometrika 97(4), 773-789.
#' @keywords multicore, cAIC
#' @rdname MC_anocAIC
#' @export MC_anocAIC
#' @import parallel
#' @importFrom stats formula
#' @importFrom cAIC4 cAIC
#'
  MC_anocAIC <- function(objects, digits = 2){
    #Multicore Backend
    numberCores <- detectCores()
    tasks <- length(objects)
    if (tasks < numberCores) {
      threads <- tasks
    } else {
      threads <- numberCores -1
    }

    objs <- objects
    cAICs <- mclapply(objs, cAIC, mc.cores = threads)
    frms <- sapply(objs, function(x) Reduce(paste, deparse(formula(x))))
    refit <- sapply(cAICs, "[[", "new")
    if (any(refit))
      frms[which(refit)] <- sapply(cAICs[which(refit)],
                                   function(x) Reduce(paste,
                                                      deparse(formula(x$reducedModel))))
    ret <- as.data.frame(do.call("rbind",
                                 lapply(cAICs, function(x)
                                   round(unlist(x[c("loglikelihood", "df", "caic", "new")]),
                                         digits = digits))))
    ret[, 4] <- as.logical(ret[, 4])
    rownames(ret) <- make.unique(frms, sep = " % duplicate #")
    colnames(ret) <- c("cll", "df", "cAIC", "Refit")
    return(ret)
  }
