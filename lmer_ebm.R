# Dependency
library(lme4)

.get_lambda <- function(object, fpc2 = NULL, pop_clus_size = NULL) {
  # Computing reliability (lambda)
  nj <- ave(object@resp$y, object@flist[[1]], FUN = length)
  theta2x <- object@theta^2
  if (is.null(fpc2) & !is.null(pop_clus_size)) {
    fpc2 <- (pop_clus_size - nj) / (pop_clus_size - 1)
  } else if (is.null(fpc2) & is.null(pop_clus_size)) {
    fpc2 <- 1
  }
  theta2x / (theta2x + 1 / nj * fpc2)
}

compute_ebmc <- function(data, formulax,
                         fpc2 = NULL, pop_clus_size = NULL,
                         boundary_avoiding = FALSE) {
  fun_lmer_x <- if (boundary_avoiding) blme::blmer else lme4::lmer
  sample_mlmxcov <- fun_lmer_x(formulax,
      data = data,
      control = lmerControl(calc.derivs = FALSE)
    )
  lambdax_cm <- .get_lambda(sample_mlmxcov,
    fpc2 = fpc2,
    pop_clus_size = pop_clus_size
  )
  name_x <- paste(formulax)[2]
  data[[paste0(name_x, "_cm")]] <-
    ave(sample_mlmxcov@resp$y, sample_mlmxcov@flist[[1]], FUN = mean)
  fixed <- as.vector(sample_mlmxcov@pp$X %*% sample_mlmxcov@beta)
  data[[paste0(name_x, "_ebm")]] <-
    ave(fixed + lambdax_cm * (sample_mlmxcov@resp$y - fixed),
      sample_mlmxcov@flist[[1]],
      FUN = mean
    )
  data[[paste0(name_x, "_ebmc")]] <-
    data[[name_x]] - data[[paste0(name_x, "_ebm")]]
  list(data = data, relx_cm = lambdax_cm, 
       name_x = name_x, 
       taux = sample_mlmxcov@theta * sigma(sample_mlmxcov))
}

#' Perform empirical Bayes mean centering (EBM)
#'
#' This function performs 2-step EBM by first obtaining the empirical Bayes
#'   estimates of cluster means, and then perform centering using those means
#'   instead of the observed cluster means.
#'
#' @param formula A two-sided linear formula for the final model. Cluster mean
#'   variables should be denoted with a "_ebm" extension, whereas cluster-mean
#'   centered variables should be denoted with a "_ebmc" extension.
#' @param data A data frame containing the variables named in formula, to be
#'   passed to [lme4::lmer()].
#' @param formulax A two-sided linear formula for obtaining the empirical
#'   Bayes means for the specified predictor. All level-2 covariates in the
#'   final model should be included here, but no random slopes should be
#'   included.
#' @param fpc2 A numeric vector specifying the finite-population correction
#'   factors due to finite cluster sizes for each observation, with length
#'   equal to the length of the outcome variable (after dropping missing data).
#'   If NULL, fpc2 will be computed based on [pop_clus_size].
#' @param pop_clus_size A numeric vector specifying the population cluster sizes
#'   for each observation. Only used when [fpc2] is NULL. If both [fpc2] and
#'   [pop_clus_size] are NULL, [fpc2] will be assumed 1 (i.e., no finite
#'   population correction).
#' @param boundary_avoiding Logical indicating whether boundary avoiding priors
#'   recommended by Chung et al. (2013, doi: 10.1007/ S 11336-013-9328-2)
#'   should be used when estimating the cluster mean reliability. Useful when
#'   the estimated between-cluster variance of the predictor is zero with
#'   [lme4::lmer()].
#' @param ... Additional arguments passed to [lme4::lmer()].
#'
#' @return An object of class \linkS4class{merMod}.
lmer_ebm <- function(formula, data, formulax,
                     fpc2 = NULL, pop_clus_size = NULL,
                     boundary_avoiding = FALSE, ...) {
  if (boundary_avoiding && !require(blme)) {
    stop("Please install the `blme` package for boundary-avoiding priors.")
  }
  ebmc <- compute_ebmc(data,
    formulax = formulax,
    fpc2 = fpc2, pop_clus_size = pop_clus_size,
    boundary_avoiding = boundary_avoiding
  )
  out <- lmer(formula, data = ebmc$data, ...)
  est <- fixef(out)
  tau0sq_bias <- (1 - mean(ebmc$relx_cm)) *
    (est[paste0(ebmc$name_x, "_ebm")] - est[paste0(ebmc$name_x, "_ebmc")])^2 *
    ebmc$taux^2
  list(fit = out, tau0sq_bias)
}