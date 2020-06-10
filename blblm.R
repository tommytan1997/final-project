#' @import purrr
#' @import furrr
#' @import future
#' @import stats
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' linear regression with mini bag bootstrap
#'
#' @param formula an object of class "formula"
#' @param data a data frame
#' @param m times to split data
#' @param B number of bootstrap performed
#' @param core number of cores applied
#'
#' @return blblm object
#'
#' @export
#'
blblm <- function(formula, data, m = 10, B = 5000, core = 1) {
  if (is.data.frame(data)) {
    data_list <- split_data(data, m)
    if (core == 1) {
      estimates <- map(data_list,
                       ~ lm_each_subsample(
                         formula = formula,
                         data = .,
                         n = nrow(data),
                         B = B
                       ))
    } else{
      suppressWarnings(plan(multiprocess, workers = core))
      estimates <- future_map(data_list,
                              ~ lm_each_subsample(
                                formula = formula,
                                data = .,
                                n = nrow(data),
                                B = B
                              ))
    }
  } else{
    if (core == 1) {
      sample_size = data %>% map( ~ {
        df <- read.csv(., )
        nrow(df)
      }) %>% reduce(`+`)
      estimates = data %>% map( ~ {
        df <- read.csv(., )
        lm_each_subsample(
          formula = formula,
          data = df,
          n = sample_size,
          B = B
        )
      })
    } else{
      suppressWarnings(plan(multiprocess, workers = core))
      sample_size = data %>% map(~ {
        df <- read.csv(.,)
        nrow(df)
      }) %>% reduce(`+`)
      estimates = data %>% future_map(~ {
        df <- read.csv(.,)
        lm_each_subsample(
          formula = formula,
          data = df,
          n = sample_size,
          B = B
        )
      })
    }
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' split data to m parts
#'
#' @param data a data frame
#' @param m times to split data
#'
#' @export
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' call lm_each_boot
#'
#' @param formula formula class object
#' @param data a data frame
#' @param n sample size
#' @param B number of bootstrap
#' @export
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

#' linear regression for the dataset
#'
#' @param formula formula
#' @param data a data frame
#' @param n total sample size
#'
#' @export
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}

#' fit regression model
#'
#' @param formula formula
#' @param data a data frame
#' @param freqs times each rows are repeated
#'
#' @export
lm1 <- function(formula, data, freqs) {
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' extract coefficients
#'
#' @param fit fitted model
#'
#' @export
blbcoef <- function(fit) {
  coef(fit)
}

#' extract sigma
#'
#' @param fit fitted model
#'
#' @export
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' print method on blblm
#'
#' @param x x object
#' @param ... more parameters
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

##################

#' print method for blblm
#'
#' @param object blblm class
#' @param confidence T/F for confidence interval
#' @param level level of confidence intervals
#' @param ... more parameters
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' extract coefficients
#'
#' @param object blblm
#' @param ... more parameters
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' get confidence interval
#'
#' @param object blblm
#' @param parm parameters
#' @param level confidence level
#' @param ... more parameters
#'
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' prediciton
#'
#' @param object blblm
#' @param new_data new data frame
#' @param confidence show confidence interval
#' @param level confidence level
#' @param ... more parameters
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}

mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
